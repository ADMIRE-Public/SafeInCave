# Copyright (c) 2026, The SafeInCave Developers
#
# SPDX-License-Identifier: BSD-3-Clause

"""Simulation report logger for element-level monitoring.

Provides centralized logging infrastructure to track element-level variables during
nonlinear iterations. Uses a flexible variable registry system that allows easy
extension without modifying the core logger.

Supports two data-sourcing backends behind one public API:

- "push" (torch/Newton path): the simulator hands full per-element ``stress``/
  ``strain`` tensors and yield kwargs to ``log_step`` every call.
- "pull" (JAX/dolfinx external-operator path, selected by passing ``problem=``):
  the logger reaches directly into a ``PlasticityProblem``-like object's
  dolfinx quadrature state every call.

Both backends resolve each target point to the nearest mesh cell in an
MPI-partition-aware way (local nearest cell per rank, then a ``MPI.MINLOC``
reduction to find the globally-nearest cell's owning rank), and both support
tracking an arbitrary number of target points in one logger instance. All
points are written to a single CSV file, one table per point (separated by
a blank line, each preceded by a "Point: ..." line giving its coordinates).
See ``safeincave.PostProcessing.extract_variable_series`` to pull one
variable's time series across points back out of that file.
"""

import os
import csv
import sys
from typing import Any, Callable, Dict, List, Optional

import numpy as np
import torch as to
from mpi4py import MPI


# ============================================================================
# Variable Registry System
# ============================================================================

# Global registry for variable definitions
VARIABLE_REGISTRY: Dict[str, tuple] = {}


def register_variable(name: str, header: str, extractor: Callable) -> None:
    """
    Register a new variable in the global registry.

    This allows users to add custom variables without modifying the logger code.

    Parameters
    ----------
    name : str
        Variable name (used in variables_to_track list).
    header : str
        Column header for CSV output (e.g., 'Stress XX (Pa)').
    extractor : Callable
        Function that extracts the variable value. Signature:
        def extractor(elem_id: int, **kwargs) -> float:
            ...

    Examples
    --------
    >>> def extract_temperature(elem_id, **kwargs):
    ...     T = kwargs.get('temperature')
    ...     return T[elem_id].item() if T is not None else 0.0
    >>> register_variable('T', 'Temperature (K)', extract_temperature)
    """
    VARIABLE_REGISTRY[name] = (header, extractor)


def get_variable(name: str) -> Optional[tuple]:
    """
    Get variable definition from registry.

    Parameters
    ----------
    name : str
        Variable name.

    Returns
    -------
    tuple or None
        (header, extractor) tuple if found, None otherwise.
    """
    return VARIABLE_REGISTRY.get(name)


def list_registered_variables() -> List[str]:
    """
    List all registered variable names.

    Returns
    -------
    list of str
        Names of all registered variables.
    """
    return list(VARIABLE_REGISTRY.keys())


# ============================================================================
# Built-in Variable Extractors
# ============================================================================

def _extract_stress_component(component: str):
    """Create an extractor for stress tensor components."""
    component_map = {
        'xx': (0, 0), 'yy': (1, 1), 'zz': (2, 2),
        'xy': (0, 1), 'yx': (0, 1), 'yz': (1, 2), 'zy': (1, 2),
        'xz': (0, 2), 'zx': (0, 2)
    }

    if component not in component_map:
        raise ValueError(f"Unknown stress component: {component}")

    i, j = component_map[component]

    def extractor(elem_id: int, stress, **kwargs) -> float:
        """Extract stress component for given element."""
        return stress[elem_id, i, j].item()

    return extractor


def _extract_strain_component(component: str):
    """Create an extractor for strain tensor components."""
    component_map = {
        'xx': (0, 0), 'yy': (1, 1), 'zz': (2, 2),
        'xy': (0, 1), 'yx': (0, 1), 'yz': (1, 2), 'zy': (1, 2),
        'xz': (0, 2), 'zx': (0, 2)
    }

    if component not in component_map:
        raise ValueError(f"Unknown strain component: {component}")

    i, j = component_map[component]

    def extractor(elem_id: int, strain=None, **kwargs) -> float:
        """Extract strain component for given element."""
        if strain is None:
            return 0.0
        return strain[elem_id, i, j].item()

    return extractor


def _extract_mean_stress():
    """Create extractor for mean stress (hydrostatic)."""
    def extractor(elem_id: int, stress, **kwargs) -> float:
        """Extract mean stress for given element."""
        return (stress[elem_id, 0, 0] + stress[elem_id, 1, 1] + stress[elem_id, 2, 2]).item() / 3.0
    return extractor


def _extract_von_mises_stress():
    """Create extractor for von Mises stress."""
    def extractor(elem_id: int, stress, **kwargs) -> float:
        """Extract von Mises stress for given element."""
        s = stress[elem_id]
        s_xx, s_yy, s_zz = s[0, 0], s[1, 1], s[2, 2]
        s_xy, s_yz, s_xz = s[0, 1], s[1, 2], s[0, 2]

        # Compute mean stress (hydrostatic part)
        s_m = (s_xx + s_yy + s_zz) / 3.0

        # Compute deviatoric stress components
        s_xx_dev = s_xx - s_m
        s_yy_dev = s_yy - s_m
        s_zz_dev = s_zz - s_m

        # Compute second deviatoric invariant J2 using deviatoric components
        j2 = ((s_xx_dev - s_yy_dev)**2 + (s_yy_dev - s_zz_dev)**2 +
              (s_zz_dev - s_xx_dev)**2 + 6.0 * (s_xy**2 + s_yz**2 + s_xz**2)) / 6.0

        # Von Mises stress: σ_vm = √(3·J2)
        von_mises = (3.0 * j2) ** 0.5
        return float(von_mises)
    return extractor


def _extract_principal_stress(index: int):
    """Create extractor for principal stress at given index (0=s1, 1=s2, 2=s3).
    Index 0 gives the largest eigenvalue (least compressive / most tensile)."""
    def extractor(elem_id: int, stress, **kwargs) -> float:
        """Extract principal stress for given element."""
        s = stress[elem_id]
        s_np = s.numpy() if hasattr(s, "numpy") else np.asarray(s)
        eigvals = np.linalg.eigvalsh(s_np)
        if eigvals.ndim > 1:
            eigvals = eigvals[0]
        eigvals_desc = eigvals[::-1]
        return float(eigvals_desc[index])
    return extractor


def _extract_principal_strain(index: int):
    """Create extractor for principal strain at given index (0=e1, 1=e2, 2=e3).
    Index 0 gives the largest eigenvalue (greatest extension / least compression)."""
    def extractor(elem_id: int, strain=None, **kwargs) -> float:
        """Extract principal strain for given element."""
        if strain is None:
            return 0.0
        eps = strain[elem_id]
        eps_np = eps.numpy() if hasattr(eps, "numpy") else np.asarray(eps)
        eigvals = np.linalg.eigvalsh(eps_np)
        if eigvals.ndim > 1:
            eigvals = eigvals[0]
        eigvals_desc = eigvals[::-1]
        return float(eigvals_desc[index])
    return extractor


def _extract_yield_function():
    """Create extractor for yield function value."""
    def extractor(elem_id: int, **kwargs) -> float:
        """Extract yield function value for given element."""
        F = kwargs.get('yield_function')
        if F is None:
            return 0.0
        return F[elem_id].item()
    return extractor


def _extract_plastic_multiplier():
    """Create extractor for plastic multiplier."""
    def extractor(elem_id: int, **kwargs) -> float:
        """Extract plastic multiplier for given element."""
        delta_lambda = kwargs.get('plastic_multiplier')
        if delta_lambda is None:
            return 0.0
        return delta_lambda[elem_id].item()
    return extractor


def _extract_yield_indicator():
    """Create extractor for yield indicator (boolean)."""
    def extractor(elem_id: int, **kwargs) -> int:
        """Extract yield indicator for given element."""
        yield_indicator = kwargs.get('yield_indicator')
        if yield_indicator is None:
            return 0
        return int(yield_indicator[elem_id].item())
    return extractor


def _extract_yield_indicator_h():
    """Create extractor for the tension-cutoff (Rankine) yield indicator."""
    def extractor(elem_id: int, **kwargs) -> int:
        """Extract Rankine yield indicator for given element."""
        yield_indicator_h = kwargs.get('yield_indicator_h')
        if yield_indicator_h is None:
            return 0
        return int(yield_indicator_h[elem_id].item())
    return extractor


def _extract_generic_yield():
    """Create extractor for a generic yield indicator: 1 if any active yield
    surface (Drucker-Prager, Mohr-Coulomb, Rankine tension cutoff, ...) is
    yielding, regardless of which one."""
    def extractor(elem_id: int, **kwargs) -> int:
        yielding = 0
        yield_indicator = kwargs.get('yield_indicator')
        if yield_indicator is not None and int(yield_indicator[elem_id].item()):
            yielding = 1
        yield_indicator_h = kwargs.get('yield_indicator_h')
        if yield_indicator_h is not None and int(yield_indicator_h[elem_id].item()):
            yielding = 1
        return yielding
    return extractor


# ============================================================================
# Register Default Variables
# ============================================================================

def _register_default_variables():
    """Register the default set of variables."""
    register_variable('sxx', 'sxx (Pa)', _extract_stress_component('xx'))
    register_variable('syy', 'syy (Pa)', _extract_stress_component('yy'))
    register_variable('szz', 'szz (Pa)', _extract_stress_component('zz'))
    register_variable('sm', 'sm (Pa)', _extract_mean_stress())
    register_variable('svM', 'svM (Pa)', _extract_von_mises_stress())
    register_variable('s1', 's1 (Pa)', _extract_principal_stress(0))
    register_variable('s2', 's2 (Pa)', _extract_principal_stress(1))
    register_variable('s3', 's3 (Pa)', _extract_principal_stress(2))
    register_variable('exx', 'exx (-)', _extract_strain_component('xx'))
    register_variable('eyy', 'eyy (-)', _extract_strain_component('yy'))
    register_variable('ezz', 'ezz (-)', _extract_strain_component('zz'))
    register_variable('e1', 'e1 (-)', _extract_principal_strain(0))
    register_variable('e2', 'e2 (-)', _extract_principal_strain(1))
    register_variable('e3', 'e3 (-)', _extract_principal_strain(2))
    register_variable('F', 'F (Pa)', _extract_yield_function())
    register_variable('dl', 'dl (-)', _extract_plastic_multiplier())
    register_variable('YieldR', 'YieldR', _extract_yield_indicator_h())
    register_variable('YieldDP', 'YieldDP', _extract_yield_indicator())
    register_variable('Yield', 'Yield', _extract_generic_yield())


# Register default variables when module is loaded
_register_default_variables()

#: Variables tracked when `variables_to_track` is omitted. Deliberately
#: leaner than "every registered variable": `sm`, `svM`, `F`, `dl`, `YieldR`,
#: `YieldDP` remain available on request but are no longer defaults, in
#: favor of the generic `Yield` indicator (works identically on both the
#: push and pull backends).
DEFAULT_VARIABLES: List[str] = [
    "sxx", "syy", "szz", "exx", "eyy", "ezz", "Yield",
]

#: Variables the pull (JAX/dolfinx) backend can never derive: the
#: external-operator `PlasticityProblem` only persists a stress field and a
#: single generic committed-yield field, no per-mechanism yield function or
#: plastic multiplier. Requesting these in pull mode logs 0.0 (falls out of
#: the per-variable try/except in `log_step`) after a one-time notice.
_PULL_UNSUPPORTED_VARIABLES = frozenset({"F", "dl", "YieldR"})


def _resolve_owner(local_min_dist: float, local_cell_id: int) -> tuple:
    """
    Reduce a per-rank (nearest local cell, its distance) pair across all MPI
    ranks to determine which rank owns the globally-nearest cell.

    Returns
    -------
    (owner_rank, local_cell_id_or_minus_one)
        `local_cell_id_or_minus_one` is `local_cell_id` on the owning rank,
        -1 on every other rank.
    """
    comm = MPI.COMM_WORLD
    _, owner = comm.allreduce((local_min_dist, comm.rank), op=MPI.MINLOC)
    return owner, (local_cell_id if comm.rank == owner else -1)


# ============================================================================
# Simulation Logging Class
# ============================================================================

class SimulationLogging:
    """
    Central simulation report logger for element-level monitoring at one or
    more target points, all written to a single CSV file (one table per
    point).

    Tracks user-specified variables at the mesh cell(s) nearest to the given
    target point(s), writing results to CSV with full solver context
    (iteration count, nonlinear error, time step). Uses a flexible variable
    registry system for easy extension.

    Two data-sourcing backends share this one class:

    - push (default, `problem=None`): the torch/Newton simulators hand full
      per-element `stress`/`strain` tensors and yield kwargs to `log_step`
      every call.
    - pull (`problem=<PlasticityProblem-like>`): the logger reads the
      dolfinx quadrature state directly off `problem` every call; the
      `stress`/`strain`/`**kwargs` arguments to `log_step` are ignored.

    Parameters
    ----------
    target_point : list or np.ndarray, optional
        A single target point [x, y, z]. Exactly one of `target_point` /
        `target_points` must be given.
    target_points : sequence of (list or np.ndarray), optional
        Multiple target points, each tracked independently and written as
        its own table in the log file, each preceded by a "Point: (x=..,
        y=.., z=..)" line giving its coordinates.
    variables_to_track : list of str, optional
        Names of variables to track. Variables must be registered in the
        VARIABLE_REGISTRY. Defaults to `DEFAULT_VARIABLES`.
    time_conversion : float, optional
        Conversion factor for time units (e.g., 3600 for hours to seconds).
        Default 1.0 (no conversion).
    problem : Any, optional
        A dolfinx-based `PlasticityProblem`-like object. If given, selects
        the pull backend (JAX/external-operator path). Default None (push
        backend, torch/Newton path).

    Examples
    --------
    >>> # Single point, defaults
    >>> logger = sf.SimulationLogging(target_point=[0.0, 0.0, 1.0])
    >>> # Several points along a wall, explicit variables
    >>> logger = sf.SimulationLogging(
    ...     target_points=[[18.0, 0.0, 0.3], [0.0, 18.0, 0.3]],
    ...     variables_to_track=['sxx', 'syy', 'szz', 'exx', 'eyy', 'ezz', 'Yield'],
    ... )
    """

    def __init__(
        self,
        target_point: list | np.ndarray | None = None,
        target_points: list | None = None,
        variables_to_track: Optional[List[str]] = None,
        time_conversion: float = 1.0,
        problem: Any | None = None,
    ):
        if (target_point is None) == (target_points is None):
            raise ValueError(
                "SimulationLogging requires exactly one of target_point or "
                "target_points."
            )
        raw_points = [target_point] if target_point is not None else list(target_points)
        self._points: List[np.ndarray] = []
        for p in raw_points:
            pt = np.zeros(3)
            arr = np.asarray(p, dtype=float)
            pt[: len(arr)] = arr
            self._points.append(pt)

        self.problem = problem
        self._backend = "pull" if problem is not None else "push"
        self.time_conversion = time_conversion

        # Internal state - set by configure()
        self._owner: List[int] | None = None
        self._cell: List[int] | None = None
        self._log_file: str | None = None
        self._rows: List[list] | None = None
        self._sig_view = None
        self._eps_fn = None
        self._eps_expr = None

        if variables_to_track is None:
            self.variables_to_track = list(DEFAULT_VARIABLES)
        else:
            self.variables_to_track = list(variables_to_track)

        self.variables = {}
        for var_name in self.variables_to_track:
            var_def = get_variable(var_name)
            if var_def is not None:
                self.variables[var_name] = var_def
            else:
                print(f"WARNING: Variable '{var_name}' not found in registry. Skipping.")

        if self._backend == "pull" and MPI.COMM_WORLD.rank == 0:
            unsupported = sorted(set(self.variables) & _PULL_UNSUPPORTED_VARIABLES)
            if unsupported:
                print(
                    f"NOTE: variable(s) {unsupported} are not derivable from "
                    "the JAX/dolfinx pull-mode problem state; column(s) will "
                    "read 0.0."
                )

    # ------------------------------------------------------------- setup

    def configure(
        self,
        log_file: str | None = None,
        *,
        centroids: np.ndarray | None = None,
        time_conversion: float | None = None,
    ) -> None:
        """
        Resolve each target point to its nearest cell and open a single CSV
        file at `log_file` containing one table per point (separated by a
        blank line, each preceded by a "Point: ..." line giving its
        coordinates).

        Parameters
        ----------
        log_file : str, optional
            Path to the CSV file, e.g. ".../log.csv". If None, logging is
            disabled.
        centroids : np.ndarray, optional
            Local-rank element centroids, shape (N, 3). Required for the
            push backend; ignored for the pull backend (which resolves
            nearest cells directly from `problem.domain`).
        time_conversion : float, optional
            Time conversion factor. If provided, overrides the value from
            __init__.
        """
        if log_file is None:
            self._owner = None
            self._cell = None
            self._log_file = None
            self._rows = None
            return

        if time_conversion is not None:
            self.time_conversion = time_conversion

        if self._backend == "push":
            if centroids is None:
                raise ValueError(
                    "centroids is required to configure a push-backend "
                    "SimulationLogging (problem=None)."
                )
            owners, cells, locations = self._resolve_push(centroids)
        else:
            owners, cells, locations = self._resolve_pull()

        # Replace the requested points with the actual resolved (and
        # deduplicated) locations data is extracted from -- both what gets
        # logged in "Point: ..." lines and what downstream extraction uses.
        self._owner, self._cell, self._points = self._deduplicate(owners, cells, locations)

        self._log_file = log_file
        self._rows = [[] for _ in self._points]
        if MPI.COMM_WORLD.rank == 0:
            self._write_file()
            print(self._log_file)

    def _resolve_push(self, centroids: np.ndarray) -> tuple:
        return self._resolve_against(centroids)

    def _resolve_pull(self) -> tuple:
        from dolfinx.mesh import compute_midpoints

        domain = self.problem.domain
        tdim = domain.topology.dim
        n_local = domain.topology.index_map(tdim).size_local
        mids = (
            compute_midpoints(domain, tdim, np.arange(n_local, dtype=np.int32))
            if n_local > 0
            else np.zeros((0, 3))
        )
        return self._resolve_against(mids)

    def _resolve_against(self, candidates: np.ndarray) -> tuple:
        """
        Resolve each requested point to its nearest candidate location
        (element centroid or cell midpoint), MPI-owner-aware. Returns
        (owners, cells, locations) with every entry identical on every
        rank -- `cells`/`locations` are broadcast from the owning rank so
        deduplication (by the caller) is consistent across ranks.
        """
        comm = MPI.COMM_WORLD
        owners, cells, locations = [], [], []
        for pt in self._points:
            if candidates.shape[0] > 0:
                d = np.linalg.norm(candidates - pt[None, : candidates.shape[1]], axis=1)
                local_min, local_cell = float(np.min(d)), int(np.argmin(d))
            else:
                local_min, local_cell = np.inf, -1
            owner, cell = _resolve_owner(local_min, local_cell)
            loc = candidates[cell].copy() if (comm.rank == owner and cell >= 0) else None
            cell = comm.bcast(cell, root=owner)
            loc = comm.bcast(loc, root=owner)
            if loc is None:
                # Degenerate case (e.g. empty mesh on every rank): fall
                # back to the requested point rather than crash.
                loc = pt.copy()
            padded = np.zeros(3)
            padded[: len(loc)] = loc
            owners.append(owner)
            cells.append(cell)
            locations.append(padded)
        return owners, cells, locations

    @staticmethod
    def _deduplicate(owners: list, cells: list, locations: list) -> tuple:
        """
        Merge requested points that resolved to the same mesh cell (same
        owner + local cell id -- a valid unique identifier for a specific
        cell even without a global numbering), keeping first-seen order.
        Runs identically on every rank since `owners`/`cells` are already
        broadcast/consistent everywhere.
        """
        seen = {}
        u_owners, u_cells, u_locations = [], [], []
        n_duplicates = 0
        for owner, cell, loc in zip(owners, cells, locations):
            key = (owner, cell)
            if key in seen:
                n_duplicates += 1
                continue
            seen[key] = True
            u_owners.append(owner)
            u_cells.append(cell)
            u_locations.append(loc)
        if n_duplicates and MPI.COMM_WORLD.rank == 0:
            print(
                f"NOTE: {n_duplicates} requested point(s) resolved to a mesh "
                f"cell already claimed by another point; merged, "
                f"{len(u_owners)} unique point(s) tracked."
            )
        return u_owners, u_cells, u_locations

    def _write_file(self) -> None:
        """
        Rewrite the whole log file (MPI rank 0 only): one table per point,
        each preceded by a "Point: ..." coordinate line and followed by a
        blank line. Rewritten from the in-memory row buffer on every
        `log_step` call -- simple and crash-safe (the file is always a
        complete, valid snapshot), and cheap at the row/point counts these
        loggers see in practice.
        """
        try:
            folder = os.path.dirname(self._log_file)
            if folder:
                os.makedirs(folder, exist_ok=True)

            header = ['Step', 't (h)', 'dt (h)', 't/t_final', 'Iters', 'NL_Error']
            for var_name, (var_header, _) in self.variables.items():
                header.append(var_header)

            with open(self._log_file, 'w', newline='') as f:
                writer = csv.writer(f)
                for i, point in enumerate(self._points):
                    writer.writerow([
                        f"Point: (x={point[0]:.6f}, y={point[1]:.6f}, "
                        f"z={point[2]:.6f})"
                    ])
                    writer.writerow(header)
                    for row in self._rows[i]:
                        writer.writerow(row)
                    writer.writerow([])
        except Exception as e:
            print(f"ERROR writing log file: {e}", file=sys.stderr)

    # ----------------------------------------------------------- extract

    def _pull_context(self, cell_id: int) -> tuple:
        """Build single-cell (elem_id=0) `stress`/`strain`/`**kwargs`
        equivalents to what a push-mode caller would pass, sourced live from
        `self.problem` at `cell_id`. Reuses the same registry extractors as
        the push backend."""
        from safeincave.Output.field_utils import _QuadratureView
        from safeincave.Output.mandel import mandel_to_tensor3x3

        problem = self.problem
        if self._sig_view is None:
            self._sig_view = _QuadratureView(problem.sigma_n, problem.stress_dim)

        sig_mandel = self._sig_view.cell_values()[cell_id].mean(axis=0)
        stress = mandel_to_tensor3x3(sig_mandel)  # (1, 3, 3)
        strain = self._pull_strain(cell_id)  # (1, 3, 3)

        kwargs = {}
        if getattr(problem, "yield_n", None) is not None:
            qp_rows = self._sig_view.dm[cell_id]
            F = float(np.max(problem.yield_n.reshape(-1)[qp_rows]))
            kwargs["yield_indicator"] = np.array([1 if F > 0.0 else 0])

        return 0, stress, strain, kwargs

    def _pull_strain(self, cell_id: int) -> np.ndarray:
        from dolfinx import fem
        import ufl

        problem = self.problem
        g = problem.gdim
        if self._eps_fn is None:
            W = fem.functionspace(problem.domain, ("DG", 0, (g, g)))
            self._eps_fn = fem.Function(W)
            self._eps_expr = fem.Expression(
                ufl.sym(ufl.grad(problem.u)), W.element.interpolation_points()
            )
        self._eps_fn.interpolate(self._eps_expr)
        dof = self._eps_fn.function_space.dofmap.list[cell_id, 0]
        eps_local = self._eps_fn.x.array.reshape(-1, g, g)[dof]

        full = np.zeros((1, 3, 3))
        full[0, :g, :g] = eps_local
        return full

    def _extract_point_values(
        self, point_index: int, stress=None, strain=None, **kwargs
    ) -> Dict[str, float]:
        cell_id = self._cell[point_index]
        if self._backend == "pull":
            elem_id, stress, strain, kwargs = self._pull_context(cell_id)
        else:
            elem_id = cell_id

        vals = {}
        for var_name, (_, extractor) in self.variables.items():
            try:
                vals[var_name] = extractor(elem_id, stress=stress, strain=strain, **kwargs)
            except Exception:
                vals[var_name] = 0.0
        return vals

    # -------------------------------------------------------------- log

    def log_initial_state(
        self,
        time: float = 0.0,
        time_final: float = 1.0,
        stress: Optional[to.Tensor] = None,
        strain: Optional[to.Tensor] = None,
        **kwargs
    ) -> None:
        """
        Log the initial state (step 0) with zero dt and zero iterations.

        Parameters
        ----------
        time : float, optional
            Initial time (raw, SI seconds). Default 0.0.
        time_final : float, optional
            Final time endpoint. Default 1.0.
        stress : torch.Tensor, optional
            Full stress tensor, shape (N, 3, 3). Push backend only.
        strain : torch.Tensor, optional
            Total strain tensor, shape (N, 3, 3). Push backend only.
        **kwargs : dict
            Additional variables for extractors.
        """
        self.log_step(
            step=0,
            time=time,
            time_step=0.0,
            time_final=time_final,
            iteration=-1,          # displays as Iters=0
            nonlinear_error=0.0,
            stress=stress,
            strain=strain,
            **kwargs
        )

    def log_step(
        self,
        step: int,
        time: float = 0.0,
        time_step: float = 0.0,
        time_final: float = 1.0,
        iteration: int = 0,
        nonlinear_error: float = 0.0,
        stress: Optional[to.Tensor] = None,
        strain: Optional[to.Tensor] = None,
        **kwargs
    ) -> None:
        """
        Log simulation state at each monitored point.

        Parameters
        ----------
        step : int
            Time step number (raw).
        time : float, optional
            Current absolute time (raw, SI seconds). Default 0.0.
        time_step : float, optional
            Time step size (raw, SI seconds). Default 0.0.
        time_final : float, optional
            Final time endpoint. Default 1.0.
        iteration : int, optional
            Newton iteration number (raw 0-indexed). Default 0.
        nonlinear_error : float, optional
            Nonlinear error estimate. Default 0.0.
        stress : torch.Tensor, optional
            Full stress tensor, shape (N, 3, 3). Push backend only; ignored
            by the pull backend, which sources state from `self.problem`.
        strain : torch.Tensor, optional
            Total strain tensor, shape (N, 3, 3). Push backend only.
        **kwargs : dict
            Additional variables for extractors (e.g., yield_function,
            plastic_multiplier, yield_indicator, temperature, etc.). Push
            backend only.
        """
        if self._log_file is None:
            return

        comm = MPI.COMM_WORLD
        time_converted = time / self.time_conversion if self.time_conversion != 1.0 else time
        time_step_converted = time_step / self.time_conversion if self.time_conversion != 1.0 else time_step
        iteration_display = iteration + 1
        time_ratio = 0.0 if time_final == 0.0 else time / time_final
        base_row = [step, time_converted, time_step_converted, time_ratio, iteration_display, nonlinear_error]

        for i in range(len(self._points)):
            owner = self._owner[i]
            vals = None
            if comm.rank == owner:
                vals = self._extract_point_values(i, stress=stress, strain=strain, **kwargs)

            if owner != 0:
                if comm.rank == owner:
                    comm.send(vals, dest=0, tag=71 + i)
                elif comm.rank == 0:
                    vals = comm.recv(source=owner, tag=71 + i)

            if comm.rank != 0:
                continue

            row = list(base_row)
            for var_name in self.variables:
                row.append(vals.get(var_name, 0.0) if vals is not None else 0.0)
            self._rows[i].append(row)

        if comm.rank == 0:
            self._write_file()

    @staticmethod
    def get_centroids(grid):
        """
        Compute element centroids from mesh geometry.

        Parameters
        ----------
        grid : GridHandlerGMSH
            Grid object with mesh and topology information.

        Returns
        -------
        numpy.ndarray
            Element centroids, shape (n_elems, 3), where each row is the
            arithmetic mean of the four vertices of the corresponding tetrahedral cell.
        """
        conn = grid.mesh.topology.connectivity(3, 0).array.reshape(-1, 4)
        coord = grid.mesh.geometry.x
        centroids = np.zeros((grid.n_elems, 3))
        for i in range(grid.n_elems):
            nodes = conn[i]
            centroids[i, :] = (
                coord[nodes[0], :] + coord[nodes[1], :] + coord[nodes[2], :] + coord[nodes[3], :]
            ) / 4.0
        return centroids

    @staticmethod
    def auto_configure_from_simulator(logger, simulator, outputs) -> None:
        """
        Auto-configure logger from simulator and outputs.

        This is a convenience method that extracts configuration from a simulator
        and configures the logger. Use this in simulator __init__ methods.

        Parameters
        ----------
        logger : SimulationLogging
            Logger to configure.
        simulator : Simulator
            Simulator instance (Simulator_M, Simulator_TM, etc.).
        outputs : list of SaveFields
            Output handlers to extract output folder from.
        """
        if logger is None:
            return

        # Get centroids from grid
        try:
            grid = getattr(simulator, 'eq_mom', None)
            if grid is None:
                grid = getattr(simulator, 'eq_heat', None)
            if grid is None:
                return

            centroids = SimulationLogging.get_centroids(grid.grid)
        except Exception as e:
            print(f"WARNING: Could not extract centroids for logger: {e}")
            return

        # Auto-detect log file path from outputs
        log_file = None
        for output in outputs:
            if hasattr(output, 'output_folder') and output.output_folder:
                import os
                log_file = os.path.join(output.output_folder, "log.csv")
                break

        # Get time_conversion from time controller
        time_conversion = getattr(simulator, 't_control', None)
        time_conversion = getattr(time_conversion, 'time_conversion', 1.0) if time_conversion else 1.0

        # Configure logger
        if log_file is not None:
            logger.configure(
                log_file=log_file,
                centroids=centroids,
                time_conversion=time_conversion,
            )

    @staticmethod
    def extract_yield_variables(material) -> dict:
        """
        Extract yield-related variables from material for logging.

        Parameters
        ----------
        material : Material
            Material object with non-elastic elements.

        Returns
        -------
        dict
            Dictionary with yield_function, plastic_multiplier, and yield_indicator
            if available, otherwise empty dict.
        """
        log_kwargs = {}

        if not hasattr(material, 'elems_ne'):
            return log_kwargs

        # Try to extract yield indicators from non-elastic elements
        for elem_ne in material.elems_ne:
            if hasattr(elem_ne, 'F') and hasattr(elem_ne, 'delta_lambda'):
                log_kwargs['yield_function'] = elem_ne.F
                log_kwargs['plastic_multiplier'] = elem_ne.delta_lambda
                if hasattr(elem_ne, 'YieldDP'):
                    log_kwargs['yield_indicator'] = elem_ne.YieldDP
                elif hasattr(elem_ne, 'is_plastic'):
                    log_kwargs['yield_indicator'] = elem_ne.is_plastic
                break

        # Tension-cutoff (Rankine) indicator lives on its own model, which may
        # not be the element that supplied F/delta_lambda above.
        for elem_ne in material.elems_ne:
            if hasattr(elem_ne, 'YieldR'):
                log_kwargs['yield_indicator_h'] = elem_ne.YieldR
                break

        return log_kwargs
