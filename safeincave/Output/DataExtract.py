# Copyright (c) 2026, The SafeInCave Developers
#
# SPDX-License-Identifier: BSD-3-Clause

"""Point-wise data recording: live simulation logging and after-the-run extraction.

This is the single home for both halves of the YAML `outputs:` entry types
`extract_fixed_point` (one point, several variables) and
`extract_fixed_variable` (one variable, several points), selected by
`while_simulating`:

- **live** (`while_simulating: true`): `SimulationLogging` tracks
  user-specified variables at the mesh cell(s) nearest one or more target
  points while the simulation runs, writing CSV with full solver context
  (iteration count, nonlinear error, time step) every step. Two
  data-sourcing backends share this one class:

  - "push" (torch/Newton path): the simulator hands full per-element
    ``stress``/``strain`` tensors and yield kwargs to ``log_step`` every call.
  - "pull" (JAX/dolfinx external-operator path, selected by passing
    ``problem=``): the logger reaches directly into a
    ``PlasticityProblem``-like object's dolfinx quadrature state every call.

  Both backends resolve each target point to the nearest mesh cell in an
  MPI-partition-aware way (local nearest cell per rank, then a
  ``MPI.MINLOC`` reduction to find the globally-nearest cell's owning rank),
  and both support tracking an arbitrary number of target points in one
  logger instance. A step with more than one live profile fans a
  simulator's single ``simulation_logger`` hook out through `CompositeLogger`.

- **post-run** (`while_simulating: false`, the default): `extract_point` and
  `extract_variable` read a completed simulation's actual result files (the
  XDMF field output written during the run -- stress tensors, principal
  stresses, etc.) and pull out time series at the mesh cell(s) nearest the
  requested points. Reading whatever fields were saved via
  `SaveFields`/`SaveFieldsNewton` (torch path) or `AutoDiffJAX.SaveFields`
  (JAX path) -- both write per-field files at the identical
  ``{output_folder}/{field_name}/{field_name}.xdmf`` convention -- needs no
  backend-specific branching at all. The flip side is that solver state
  never written as a field can only be recorded live; see each variable's
  `post_run` mapping in the registry (`VARIABLE_REGISTRY`).

Both halves write the same two CSV layouts, so a plotting script never has
to know which one produced a given file:

- **per point** (`layout="per_point"` / `extract_point`): one file per
  point, named ``{name}_{x}_{y}_{z}.csv`` after the point's *resolved* cell
  centroid, one column per variable.
- **per variable** (`layout="per_variable"` / `extract_variable`): one file
  per variable, named ``{name}_{variable}.csv``, one column per point, each
  headed by that point's resolved coordinates.
"""

from __future__ import annotations

import csv
import os
import re
import sys
from typing import Any, Callable, Dict, List, NamedTuple, Optional, Sequence

import numpy as np
import torch as to
from mpi4py import MPI

from ..PostProcessing.Readers import (
    find_closest_point,
    read_cell_tensor,
    read_cell_scalar,
    read_node_scalar,
    read_node_vector,
)
from ..PostProcessing.VectorTensorReaders import (
    read_cell_vector,
    read_node_tensor,
)
from ..PostProcessing.MergedFieldReaders import (
    read_cell_tensor as read_cell_tensor_merged,
    read_cell_scalar as read_cell_scalar_merged,
    read_node_scalar as read_node_scalar_merged,
    read_node_vector as read_node_vector_merged,
)


# ============================================================================
# Shared file-naming conventions
# ============================================================================

#: Header of the leading time column in every post-run extracted CSV.
TIME_HEADER = "t"


def point_label(point: Sequence[float]) -> str:
    """
    Render a point as the compact ``x_y_z`` string used both as a column
    header (per-variable layout) and as a file-name suffix (per-point
    layout), e.g. ``"0.22_18.09_0.27"``.

    The point passed in should be the *resolved* cell centroid rather than
    the coordinates the user requested, so that the label always identifies
    the data that was actually extracted.
    """
    return "_".join(f"{int(round(c * 100))}" for c in point[:3])


def point_file_name(name: str, point: Sequence[float]) -> str:
    """Per-point layout file name: ``{name}_{x}_{y}_{z}.csv``."""
    return f"{name}_{point_label(point)}.csv"


def variable_file_name(name: str, variable: str) -> str:
    """Per-variable layout file name: ``{name}.csv``."""
    return f"{name}.csv"


# ============================================================================
# Variable Registry System
# ============================================================================
#
# Single source of truth for every name a YAML `variables:`/`fields:` list
# can contain: how to pull it live during a run (`extractor`, used by
# `SimulationLogging`), and -- if it's ever written as a saved XDMF field --
# how to read it back after the run (`post_run`, used by `extract_point`/
# `extract_variable`). A variable with `post_run=None` is solver-only state
# that can only ever be recorded live.


class VariableSpec(NamedTuple):
    header: str
    extractor: Callable
    post_run: Optional[tuple] = None  # (field_name, shape_kind, component)


# Global registry for variable definitions.
VARIABLE_REGISTRY: Dict[str, VariableSpec] = {}


def _strip_unit_suffix(header: str) -> str:
    """Drop a trailing ' (unit)' annotation from a CSV column header.

    E.g. 'sxx (Pa)' -> 'sxx', 'exx (-)' -> 'exx'. Headers without one are
    returned unchanged.
    """
    return re.sub(r"\s*\([^)]*\)\s*$", "", header)


def register_variable(
    name: str, header: str, extractor: Callable, post_run: Optional[tuple] = None
) -> None:
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
    post_run : tuple, optional
        `(field_name, shape_kind, component)` if this variable is also
        readable from a saved XDMF field via `extract_point`/
        `extract_variable` (see `_VARIABLE_SPECS` usage below for the
        meaning of `shape_kind`/`component`). Omit for solver-only state
        that's never persisted as a field.

    Examples
    --------
    >>> def extract_temperature(elem_id, **kwargs):
    ...     T = kwargs.get('temperature')
    ...     return T[elem_id].item() if T is not None else 0.0
    >>> register_variable('T', 'Temperature (K)', extract_temperature)
    """
    VARIABLE_REGISTRY[name] = VariableSpec(header, extractor, post_run)


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
    spec = VARIABLE_REGISTRY.get(name)
    return None if spec is None else (spec.header, spec.extractor)


def list_registered_variables() -> List[str]:
    """
    List all registered variable names.

    Returns
    -------
    list of str
        Names of all registered variables.
    """
    return list(VARIABLE_REGISTRY.keys())


def _post_run_spec(name: str) -> Optional[tuple]:
    """`(field_name, shape_kind, component)` for `name` if it's readable
    from a saved XDMF field, else None."""
    spec = VARIABLE_REGISTRY.get(name)
    return None if spec is None else spec.post_run


def _post_run_names() -> List[str]:
    """Names of every registered variable that's readable post-run."""
    return sorted(n for n, spec in VARIABLE_REGISTRY.items() if spec.post_run is not None)


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


#: kwargs names `_extract_generic_yield` checks for an active yield surface.
#: Extend via `register_yield_surface` for a custom mechanism's indicator.
_YIELD_SURFACE_KWARGS = {"yield_indicator", "yield_indicator_h"}


def register_yield_surface(name: str) -> None:
    """Register an additional kwargs name for `_extract_generic_yield` (the
    `yield_flag` variable) to check for an active yield surface, alongside
    the built-in `yield_indicator`/`yield_indicator_h`."""
    _YIELD_SURFACE_KWARGS.add(name.lower())


def _extract_generic_yield():
    """Create extractor for a generic yield indicator: 1 if any active yield
    surface (Drucker-Prager, Mohr-Coulomb, Rankine tension cutoff, ...) is
    yielding, regardless of which one."""
    def extractor(elem_id: int, **kwargs) -> int:
        for surface_name in _YIELD_SURFACE_KWARGS:
            indicator = kwargs.get(surface_name)
            if indicator is not None and int(indicator[elem_id].item()):
                return 1
        return 0
    return extractor


# ============================================================================
# Register Default Variables
# ============================================================================

def _register_default_variables():
    """Register the default set of variables."""
    register_variable('sxx', 'sxx (Pa)', _extract_stress_component('xx'),
                       post_run=("sig", "tensor", (0, 0)))
    register_variable('syy', 'syy (Pa)', _extract_stress_component('yy'),
                       post_run=("sig", "tensor", (1, 1)))
    register_variable('szz', 'szz (Pa)', _extract_stress_component('zz'),
                       post_run=("sig", "tensor", (2, 2)))
    register_variable('sxy', 'sxy (Pa)', _extract_stress_component('xy'),
                       post_run=("sig", "tensor", (0, 1)))
    register_variable('syz', 'syz (Pa)', _extract_stress_component('yz'),
                       post_run=("sig", "tensor", (1, 2)))
    register_variable('sxz', 'sxz (Pa)', _extract_stress_component('xz'),
                       post_run=("sig", "tensor", (0, 2)))
    register_variable('sm', 'sm (Pa)', _extract_mean_stress(),
                       post_run=("p_elems", "scalar", None))
    register_variable('svM', 'svM (Pa)', _extract_von_mises_stress(),
                       post_run=("q_elems", "scalar", None))
    register_variable('s1', 's1 (Pa)', _extract_principal_stress(0),
                       post_run=("principal_stresses", "vector", 0))
    register_variable('s2', 's2 (Pa)', _extract_principal_stress(1),
                       post_run=("principal_stresses", "vector", 1))
    register_variable('s3', 's3 (Pa)', _extract_principal_stress(2),
                       post_run=("principal_stresses", "vector", 2))
    register_variable('exx', 'exx (-)', _extract_strain_component('xx'),
                       post_run=("eps_tot", "tensor", (0, 0)))
    register_variable('eyy', 'eyy (-)', _extract_strain_component('yy'),
                       post_run=("eps_tot", "tensor", (1, 1)))
    register_variable('ezz', 'ezz (-)', _extract_strain_component('zz'),
                       post_run=("eps_tot", "tensor", (2, 2)))
    register_variable('exy', 'exy (-)', _extract_strain_component('xy'),
                       post_run=("eps_tot", "tensor", (0, 1)))
    register_variable('eyz', 'eyz (-)', _extract_strain_component('yz'),
                       post_run=("eps_tot", "tensor", (1, 2)))
    register_variable('exz', 'exz (-)', _extract_strain_component('xz'),
                       post_run=("eps_tot", "tensor", (0, 2)))
    register_variable('e1', 'e1 (-)', _extract_principal_strain(0),
                       post_run=("principal_strains", "vector", 0))
    register_variable('e2', 'e2 (-)', _extract_principal_strain(1),
                       post_run=("principal_strains", "vector", 1))
    register_variable('e3', 'e3 (-)', _extract_principal_strain(2),
                       post_run=("principal_strains", "vector", 2))
    register_variable('yield_flag', 'yield_flag', _extract_generic_yield(),
                       post_run=("yield_flag", "scalar", None))


# Register default variables when module is loaded
_register_default_variables()

#: Variables tracked when `variables_to_track` is omitted. Deliberately
#: leaner than "every registered variable": `sm`, `svM` remain available on
#: request but are no longer defaults, in favor of the generic `yield_flag`
#: indicator (works identically on both the push and pull backends; the
#: mechanism-specific `YieldDP`/`YieldR` indicators it superseded have been
#: removed). Name matches the saved XDMF field (`fields: {yield_flag: ...}`)
#: one-for-one, so there's a single spelling across both namespaces.
DEFAULT_VARIABLES: List[str] = [
    "sxx", "syy", "szz", "exx", "eyy", "ezz", "yield_flag",
]


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
# Simulation Logging Class (live path)
# ============================================================================

class SimulationLogging:
    """
    Central simulation report logger for element-level monitoring at one or
    more target points, written as plain CSV tables while the simulation
    runs.

    Tracks user-specified variables at the mesh cell(s) nearest to the given
    target point(s), writing results to CSV with full solver context
    (iteration count, nonlinear error, time step). Uses a flexible variable
    registry system for easy extension.

    Two file layouts are available, matching the two `outputs:` entry types
    (`extract_fixed_point` / `extract_fixed_variable`) and the two
    module-level functions below that produce the same shapes after a run:

    - ``layout="per_point"`` (default, the ``extract_fixed_point`` case --
      one point, many variables): one file per tracked point, named
      ``{name}_{x}_{y}_{z}.csv`` after the point's *resolved* cell centroid,
      with one column per variable.
    - ``layout="per_variable"`` (the ``extract_fixed_variable`` case -- one
      variable, many points): one file per tracked variable, named
      ``{name}_{variable}.csv``, with one column per point.

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
        Multiple target points, each tracked independently.
    variables_to_track : list of str, optional
        Names of variables to track. Variables must be registered in the
        VARIABLE_REGISTRY. Defaults to `DEFAULT_VARIABLES`.
    layout : {"per_point", "per_variable"}, optional
        How the tracked data is laid out across files; see above. Default
        "per_point".
    name : str, optional
        File-name stem, so several logging profiles can write into one
        output folder without colliding. Default "extract".
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
    >>> # Several points along a wall, one variable, one file
    >>> logger = sf.SimulationLogging(
    ...     target_points=[[18.0, 0.0, 0.3], [0.0, 18.0, 0.3]],
    ...     variables_to_track=['sxx'],
    ...     layout='per_variable', name='wall',
    ... )
    """

    #: The layouts `_write_file` knows how to emit.
    LAYOUTS = ("per_point", "per_variable")

    def __init__(
        self,
        target_point: list | np.ndarray | None = None,
        target_points: list | None = None,
        variables_to_track: Optional[List[str]] = None,
        layout: str = "per_point",
        name: str = "extract",
        time_conversion: float = 1.0,
        problem: Any | None = None,
    ):
        if (target_point is None) == (target_points is None):
            raise ValueError(
                "SimulationLogging requires exactly one of target_point or "
                "target_points."
            )
        if layout not in self.LAYOUTS:
            raise ValueError(
                f"SimulationLogging: unknown layout {layout!r}; "
                f"expected one of {self.LAYOUTS}."
            )
        self.layout = layout
        self.name = name

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
        self._output_folder: str | None = None
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

    # ------------------------------------------------------------- setup

    def configure(
        self,
        output_folder: str | None = None,
        *,
        centroids: np.ndarray | None = None,
        time_conversion: float | None = None,
    ) -> None:
        """
        Resolve each target point to its nearest cell and open the CSV
        file(s) this logger's `layout` calls for, inside `output_folder`.

        Parameters
        ----------
        output_folder : str, optional
            Folder the CSV file(s) are written into. If None, logging is
            disabled.
        centroids : np.ndarray, optional
            Local-rank element centroids, shape (N, 3). Required for the
            push backend; ignored for the pull backend (which resolves
            nearest cells directly from `problem.domain`).
        time_conversion : float, optional
            Time conversion factor. If provided, overrides the value from
            __init__.
        """
        if output_folder is None:
            self._owner = None
            self._cell = None
            self._output_folder = None
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
        # deduplicated) locations data is extracted from -- these are what
        # the file names and column headers report.
        self._owner, self._cell, self._points = self._deduplicate(owners, cells, locations)

        self._output_folder = output_folder
        self._rows = [[] for _ in self._points]
        if MPI.COMM_WORLD.rank == 0:
            for path in self._write_file():
                print(path)

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

    #: Solver-context columns every live table carries ahead of its data
    #: columns. 't'/'dt' are in whatever unit `time_conversion` produces.
    CONTEXT_HEADER = ['Step', 't', 'dt', 't/t_final', 'Iters', 'NL_Error']

    #: Number of leading `CONTEXT_HEADER` entries in a buffered row.
    _N_CONTEXT = len(CONTEXT_HEADER)

    def _write_file(self) -> List[str]:
        """
        Rewrite the whole set of log files (MPI rank 0 only), one plain CSV
        table each per `self.layout`. Rewritten from the in-memory row
        buffer on every `log_step` call -- simple and crash-safe (each file
        is always a complete, valid snapshot), and cheap at the row/point
        counts these loggers see in practice.

        Returns the paths written.
        """
        try:
            os.makedirs(self._output_folder, exist_ok=True)
            tables = (
                self._tables_per_point()
                if self.layout == "per_point"
                else self._tables_per_variable()
            )
            paths = []
            for file_name, header, rows in tables:
                path = os.path.join(self._output_folder, file_name)
                with open(path, 'w', newline='') as f:
                    writer = csv.writer(f)
                    writer.writerow(header)
                    writer.writerows(rows)
                paths.append(path)
            return paths
        except Exception as e:
            print(f"ERROR writing log file: {e}", file=sys.stderr)
            return []

    def _tables_per_point(self) -> List[tuple]:
        """One table per point: solver context, then one column per variable."""
        header = self.CONTEXT_HEADER + [
            _strip_unit_suffix(var_header)
            for var_header, _ in self.variables.values()
        ]
        return [
            (point_file_name(self.name, point), header, self._rows[i])
            for i, point in enumerate(self._points)
        ]

    def _tables_per_variable(self) -> List[tuple]:
        """
        One table per variable: solver context, then one column per point.

        The buffered rows are per point, so this transposes them; the
        context columns are identical across points (they come from the
        same `log_step` call), so point 0's are used.
        """
        labels = [point_label(p) for p in self._points]
        tables = []
        for v, var_name in enumerate(self.variables):
            rows = []
            for r, context_row in enumerate(self._rows[0]):
                rows.append(
                    list(context_row[: self._N_CONTEXT])
                    + [self._rows[i][r][self._N_CONTEXT + v] for i in range(len(self._points))]
                )
            tables.append(
                (variable_file_name(self.name, var_name), self.CONTEXT_HEADER + labels, rows)
            )
        return tables

    # ----------------------------------------------------------- extract

    def _pull_context(self, cell_id: int) -> tuple:
        """Build single-cell (elem_id=0) `stress`/`strain`/`**kwargs`
        equivalents to what a push-mode caller would pass, sourced live from
        `self.problem` at `cell_id`. Reuses the same registry extractors as
        the push backend."""
        from ..Output.field_utils import _QuadratureView
        from ..Output.mandel import mandel_to_tensor3x3

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
            Additional variables for extractors (e.g., yield_indicator,
            yield_indicator_h, temperature, etc.). Push backend only.
        """
        if self._output_folder is None:
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

        # Auto-detect the output folder to log into
        output_folder = None
        for output in outputs:
            if hasattr(output, 'output_folder') and output.output_folder:
                output_folder = output.output_folder
                break

        # Get time_conversion from time controller
        time_conversion = getattr(simulator, 't_control', None)
        time_conversion = getattr(time_conversion, 'time_conversion', 1.0) if time_conversion else 1.0

        # Configure logger
        if output_folder is not None:
            logger.configure(
                output_folder,
                centroids=centroids,
                time_conversion=time_conversion,
            )

    @staticmethod
    def extract_yield_variables(material, cell_average=None) -> dict:
        """
        Extract yield-related variables from material for logging.

        Parameters
        ----------
        material : Material
            Material object with non-elastic elements.
        cell_average : callable, optional
            Reduces a per-state-point tensor to per-cell (e.g.
            `mom_eq.cell_average`), for solvers whose non-elastic elements
            are sized by quadrature/state point rather than by cell (e.g.
            the Newton P2 path) -- the logger's `elem_id` is a cell index,
            matching the already-cell-averaged `stress`/`strain` it also
            receives. Omit when the material's non-elastic elements are
            already sized per cell (true here for the DG0 path).

        Returns
        -------
        dict
            Dictionary with yield_indicator (and yield_indicator_h, for a
            separate tension-cutoff mechanism) if available, otherwise
            empty dict.
        """
        log_kwargs = {}

        if not hasattr(material, 'elems_ne'):
            return log_kwargs

        def _reduce(values, is_indicator):
            if cell_average is None:
                return values
            reduced = cell_average(values.to(to.float64))
            return reduced > 0 if is_indicator else reduced

        # Try to extract a generic yield indicator from non-elastic elements
        for elem_ne in material.elems_ne:
            if hasattr(elem_ne, 'YieldDP'):
                log_kwargs['yield_indicator'] = _reduce(elem_ne.YieldDP, True)
                break
            elif hasattr(elem_ne, 'is_plastic'):
                log_kwargs['yield_indicator'] = _reduce(elem_ne.is_plastic, True)
                break

        # Tension-cutoff (Rankine) indicator lives on its own model, which may
        # not be the element that supplied F/delta_lambda above.
        for elem_ne in material.elems_ne:
            if hasattr(elem_ne, 'YieldR'):
                log_kwargs['yield_indicator_h'] = _reduce(elem_ne.YieldR, True)
                break

        return log_kwargs


# ============================================================================
# Composite Logger
# ============================================================================

class CompositeLogger:
    """
    Fan a simulator's single ``simulation_logger`` hook out to several
    `SimulationLogging` instances.

    A step may declare more than one live recording profile (say, a
    per-point table of every variable at the cavern wall plus a per-variable
    table of ``sxx`` along a line), but the simulators accept exactly one
    logger object. This wrapper presents the same three-method interface
    they use -- `configure`, `log_initial_state`, `log_step` -- and forwards
    every call to each wrapped logger in turn.

    Parameters
    ----------
    loggers : sequence of SimulationLogging
        The loggers to drive. Each keeps its own points, variables, layout
        and file names, so their outputs never collide as long as their
        `name`s differ.

    Examples
    --------
    >>> logger = sf.CompositeLogger([
    ...     sf.SimulationLogging(target_point=[0., 18., .3], name='wall'),
    ...     sf.SimulationLogging(target_points=WALL_LINE, name='line',
    ...                          variables_to_track=['sxx'],
    ...                          layout='per_variable'),
    ... ])
    """

    def __init__(self, loggers):
        self.loggers = list(loggers)

    def configure(self, output_folder: str | None = None, **kwargs) -> None:
        for logger in self.loggers:
            logger.configure(output_folder, **kwargs)

    def log_initial_state(self, *args, **kwargs) -> None:
        for logger in self.loggers:
            logger.log_initial_state(*args, **kwargs)

    def log_step(self, *args, **kwargs) -> None:
        for logger in self.loggers:
            logger.log_step(*args, **kwargs)


# ============================================================================
# Post-run extraction (from saved fields)
# ============================================================================

# variable name -> (field_name, shape_kind, component) for post-run
# extraction now lives on each variable's `VariableSpec.post_run` in the
# registry above (`_post_run_spec`/`_post_run_names`); shape_kind is
# "tensor" (component = (i, j)), "vector" (component = index), or "scalar"
# (component unused). A variable with no `post_run` (e.g. previously F, dl)
# is solver-only state that can only be recorded with 'while_simulating: true'.

_CELL_READERS = {
    "tensor": read_cell_tensor,
    "scalar": read_cell_scalar,
    "vector": read_cell_vector,
}
_NODE_READERS = {
    "tensor": read_node_tensor,
    "scalar": read_node_scalar,
    "vector": read_node_vector,
}
# field_name-selecting counterparts, for reading a merged multi-field file.
# read_cell_vector/read_node_tensor are already field_name-aware, so they're
# shared between both variants unchanged.
_CELL_READERS_MERGED = {
    "tensor": read_cell_tensor_merged,
    "scalar": read_cell_scalar_merged,
    "vector": read_cell_vector,
}
_NODE_READERS_MERGED = {
    "tensor": read_node_tensor,
    "scalar": read_node_scalar_merged,
    "vector": read_node_vector_merged,
}


def _read_any(path: str, field_name: Optional[str], shape_kind: str):
    """
    Read `field_name` (or the sole field, if None) from `path`, trying
    cell-centered data first and falling back to node-centered. Some output
    configurations (e.g. `SaveFieldsNewton` with `smooth_output=True`) write
    fields node-centered rather than as the raw per-element DG0 cell data --
    the underlying quantity is otherwise identical, so callers shouldn't
    need to know or care which representation a given run used.

    The public `safeincave.PostProcessing.Readers` functions only support a
    single field per file, so a real `field_name` (merged-file read) is
    dispatched to the `_merged` reader variants instead.
    """
    if field_name is not None:
        try:
            return _CELL_READERS_MERGED[shape_kind](path, field_name=field_name)
        except (KeyError, IndexError):
            return _NODE_READERS_MERGED[shape_kind](path, field_name=field_name)
    try:
        return _CELL_READERS[shape_kind](path)
    except (KeyError, IndexError):
        return _NODE_READERS[shape_kind](path)


def _read_field(output_folder: str, field_name: str, shape_kind: str):
    """
    Read `field_name` from the merged solution file if available (falling
    back on any failure), else from the per-field result folder -- the
    latter is always present on both the torch and JAX output paths, so
    this is what makes the function backend-agnostic without any explicit
    branching.
    """
    merged_path = os.path.join(output_folder, "solution", "solution.xdmf")
    if os.path.isfile(merged_path):
        try:
            return _read_any(merged_path, field_name, shape_kind)
        except (KeyError, IndexError):
            # Expected in practice: the merged file's internal keys are the
            # human-readable labels from OUTPUT_FIELDS (e.g. "Stress (Pa)"),
            # not the short field_name ("sig") -- fall through to the
            # per-field file rather than guess at a fuzzy label match.
            pass

    per_field_path = os.path.join(output_folder, field_name, f"{field_name}.xdmf")
    if not os.path.isfile(per_field_path):
        raise FileNotFoundError(
            f"No output found for field '{field_name}' under '{output_folder}' "
            f"-- looked for '{merged_path}' and '{per_field_path}'."
        )
    return _read_any(per_field_path, None, shape_kind)


def _has_field_output(folder: str, field_name: str) -> bool:
    merged_path = os.path.join(folder, "solution", "solution.xdmf")
    per_field_path = os.path.join(folder, field_name, f"{field_name}.xdmf")
    return os.path.isfile(merged_path) or os.path.isfile(per_field_path)


def _discover_folders(output_folder: str, field_name: str) -> List[str]:
    """
    Find every folder actually holding `field_name`'s output under
    `output_folder`: `output_folder` itself, if it directly has the field
    (a flat, single-run layout, e.g. `examples/mechanics/1_triaxial`), plus
    any of its immediate subdirectories that do (a multi-phase layout, e.g.
    this benchmark's `output/geostatic` and `output/loading`) -- mirroring
    how the live logger writes separately under each phase's own
    `output.output_folder`, just discovered rather than named explicitly.
    """
    folders = []
    if _has_field_output(output_folder, field_name):
        folders.append(output_folder)
    if os.path.isdir(output_folder):
        for entry in sorted(os.listdir(output_folder)):
            sub = os.path.join(output_folder, entry)
            if os.path.isdir(sub) and _has_field_output(sub, field_name):
                folders.append(sub)
    return folders


def _any_saved_field(folder: str) -> Optional[str]:
    """
    Name of any one field actually saved as a per-field file directly under
    `folder` (a subdirectory `{field}/{field}.xdmf`), or None if `folder`
    holds no per-field output at all. Used as a last resort to establish a
    time axis/point location when none of the requested variables have
    their own saved field, so a fully-NaN table can still be written
    instead of failing outright.
    """
    if not os.path.isdir(folder):
        return None
    for entry in sorted(os.listdir(folder)):
        if entry == "solution":
            continue
        if os.path.isfile(os.path.join(folder, entry, f"{entry}.xdmf")):
            return entry
    return None


def _has_any_output(folder: str) -> bool:
    return (
        os.path.isfile(os.path.join(folder, "solution", "solution.xdmf"))
        or _any_saved_field(folder) is not None
    )


def _discover_any_folders(output_folder: str) -> List[str]:
    """Like `_discover_folders`, but for "holds any saved field at all" rather
    than one specific field -- the fallback once no requested variable's
    field is found anywhere."""
    folders = []
    if _has_any_output(output_folder):
        folders.append(output_folder)
    if os.path.isdir(output_folder):
        for entry in sorted(os.listdir(output_folder)):
            sub = os.path.join(output_folder, entry)
            if os.path.isdir(sub) and _has_any_output(sub):
                folders.append(sub)
    return folders


def _reference_axis(folder: str) -> tuple:
    """
    Read a (centroids, time_list) pair from any one saved field in
    `folder`, shape-kind-agnostic -- enough to build a NaN-filled table's
    time axis and resolve a point's nearest cell, without caring which
    field it came from. Returns (None, None) if `folder` has no per-field
    output to read at all (a merged-solution-only folder can't be probed
    this way without knowing a real field key).
    """
    field_name = _any_saved_field(folder)
    if field_name is None:
        return None, None
    path = os.path.join(folder, field_name, f"{field_name}.xdmf")
    for shape_kind in ("tensor", "vector", "scalar"):
        try:
            centroids, time_list, _ = _read_any(path, None, shape_kind)
            return centroids, time_list
        except Exception:
            continue
    return None, None


def _resolve_indices(points, centroids) -> List[int]:
    """
    Resolve every requested point to its nearest actual location up front,
    deduplicating before any extraction happens (several requested points
    can land on the same cell/node if they're packed closer than the local
    mesh resolution).
    """
    seen = set()
    unique_indices = []
    for p in points:
        idx = find_closest_point(np.asarray(p, dtype=float), centroids)
        if idx not in seen:
            seen.add(idx)
            unique_indices.append(idx)
    return unique_indices


def _series(field_data: np.ndarray, idx: int, shape_kind: str, component):
    if shape_kind == "tensor":
        i, j = component
        return field_data[:, idx, i, j]
    if shape_kind == "vector":
        return field_data[:, idx, component]
    return field_data[:, idx]


def _write_csv(output_file: str, headers: Sequence[str], time_list, columns) -> None:
    """Write a wide CSV (a time column plus one column per series), rank 0 only."""
    try:
        rank = MPI.COMM_WORLD.rank
    except Exception:
        rank = 0
    if rank != 0:
        return
    with open(output_file, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([TIME_HEADER] + list(headers))
        for row_idx, t in enumerate(time_list):
            writer.writerow([t] + [c[row_idx] for c in columns])


def _require_known(variable: str, caller: str) -> tuple:
    spec = _post_run_spec(variable)
    if spec is None:
        raise ValueError(
            f"Unknown variable '{variable}' for {caller}; it is not written as "
            f"a field by SaveFields, so it can only be recorded while "
            f"simulating. Known: {_post_run_names()}"
        )
    return spec


def _extract_variable_in_folder(
    folder: str,
    points: Sequence[Sequence[float]],
    variable: str,
    field_name: str,
    shape_kind: str,
    component,
    name: str,
    output_file: Optional[str],
) -> str:
    centroids, time_list, field_data = _read_field(folder, field_name, shape_kind)
    indices = _resolve_indices(points, centroids)

    headers = [point_label(centroids[idx]) for idx in indices]
    columns = [_series(field_data, idx, shape_kind, component) for idx in indices]

    if output_file is None:
        output_file = os.path.join(folder, variable_file_name(name, variable))
    _write_csv(output_file, headers, time_list, columns)
    return output_file


def _extract_point_in_folder(
    folder: str,
    point: Sequence[float],
    variables: Sequence[str],
    name: str,
    output_file: Optional[str],
) -> str:
    """
    Every variable whose field was actually saved by this run is read as
    usual. A variable that either has no `post_run` mapping at all
    (solver-only state) or whose field just wasn't among this step's saved
    ``outputs`` gets a column of NaN
    instead of failing the whole profile, sized to match the time axis
    established by whichever variable is read first. If *no* requested
    variable is readable, the time axis/point location falls back to
    `_reference_axis` (any saved field, whatever it is) so the file still
    comes out -- entirely NaN -- rather than failing outright; only a
    folder with no saved output whatsoever raises.
    """
    time_list = None
    location = None
    headers = []
    columns = []
    for variable in variables:
        headers.append(variable)
        spec = _post_run_spec(variable)
        if spec is None:
            columns.append(None)  # filled with NaN once time_list is known
            continue
        field_name, shape_kind, component = spec
        try:
            centroids, times, field_data = _read_field(folder, field_name, shape_kind)
        except FileNotFoundError:
            columns.append(None)  # field not saved by this run; NaN column
            continue
        idx = find_closest_point(np.asarray(point, dtype=float), centroids)
        if time_list is None:
            time_list, location = times, centroids[idx]
        columns.append(_series(field_data, idx, shape_kind, component))

    if time_list is None:
        centroids, time_list = _reference_axis(folder)
        if time_list is None:
            raise ValueError(
                f"extract_point: '{folder}' has no saved output at all; "
                f"nothing to extract."
            )
        idx = find_closest_point(np.asarray(point, dtype=float), centroids)
        location = centroids[idx]
    columns = [c if c is not None else np.full(len(time_list), np.nan) for c in columns]

    if output_file is None:
        output_file = os.path.join(folder, point_file_name(name, location))
    _write_csv(output_file, headers, time_list, columns)
    return output_file


def extract_variable(
    output_folder: str,
    points: Sequence[Sequence[float]],
    variable: str,
    name: str = "extract",
    output_file: Optional[str] = None,
) -> List[str]:
    """
    Extract one variable's time series at the mesh cell(s) nearest each of
    `points`, from the XDMF field output already saved under
    `output_folder` -- the after-the-fact counterpart of a
    `SimulationLogging` profile with ``layout="per_variable"``.

    `output_folder` may be a single run's output folder, or a base folder
    containing several phase subfolders (e.g. "output" containing
    "output/geostatic" and "output/loading") -- either way, every folder
    actually holding the requested variable's field output is discovered
    and processed, writing one CSV into each (no need to know or loop over
    phase names yourself).

    Parameters
    ----------
    output_folder : str
        A run's output folder, or a base folder containing phase
        subfolders.
    points : sequence of (x, y, z)
        Target points; each is resolved to its nearest cell via
        `find_closest_point`.
    variable : str
        A name with a `post_run` mapping in the variable registry, e.g.
        "sxx".
    name : str, optional
        File-name stem, so several extraction profiles can write into the
        same folder without colliding. Default "extract".
    output_file : str, optional
        Explicit output CSV path, used only when exactly one folder is
        processed. Default (and always, when multiple phase folders are
        found): "{name}_{variable}.csv" in each processed folder.

    Returns
    -------
    list of str
        The output file path(s) written, one per discovered folder.

    Examples
    --------
    >>> sf.PostProcessing.extract_variable(output_folder, POINTS, "sxx", name="wall")
    ['output/geostatic/wall_sxx.csv', 'output/loading/wall_sxx.csv']
    """
    field_name, shape_kind, component = _require_known(variable, "extract_variable")

    folders = _discover_folders(output_folder, field_name)
    if not folders:
        raise FileNotFoundError(
            f"No output found for field '{field_name}' (variable '{variable}') "
            f"under '{output_folder}' or its immediate subfolders."
        )

    per_call_output_file = output_file if len(folders) == 1 else None
    return [
        _extract_variable_in_folder(
            folder, points, variable, field_name, shape_kind, component,
            name, per_call_output_file,
        )
        for folder in folders
    ]


def extract_point(
    output_folder: str,
    point: Sequence[float],
    variables: Sequence[str],
    name: str = "extract",
    output_file: Optional[str] = None,
) -> List[str]:
    """
    Extract several variables' time series at the mesh cell nearest
    `point`, from the XDMF field output already saved under
    `output_folder` -- the after-the-fact counterpart of a
    `SimulationLogging` profile with ``layout="per_point"``.

    Folder discovery prefers folders holding one of the requested
    variables' own saved field, falling back to any folder with any saved
    output at all when none of them are found anywhere (e.g. every
    requested variable is solver-only, or this run just didn't save the
    relevant fields). A variable with no data to read -- whether it has no
    `post_run` mapping at all (solver-only state) or its field simply
    wasn't saved by this run -- gets a column of NaN instead of failing
    the whole profile; only a folder with no saved output whatsoever
    raises.

    Parameters
    ----------
    output_folder : str
        A run's output folder, or a base folder containing phase
        subfolders.
    point : (x, y, z)
        Target point, resolved to its nearest cell via `find_closest_point`.
    variables : sequence of str
        Names with a `post_run` mapping in the variable registry, e.g.
        ["sxx", "syy"]. Each becomes one column.
    name : str, optional
        File-name stem. Default "extract".
    output_file : str, optional
        Explicit output CSV path, used only when exactly one folder is
        processed. Default: "{name}_{x}_{y}_{z}.csv" (the *resolved* cell
        centroid) in each processed folder.

    Returns
    -------
    list of str
        The output file path(s) written, one per discovered folder.

    Examples
    --------
    >>> sf.PostProcessing.extract_point(output_folder, [0., 18., .3],
    ...                                 ["sxx", "syy"], name="wall")
    ['output/loading/wall_0.22_18.09_0.27.csv']
    """
    variables = list(variables)
    if not variables:
        raise ValueError("extract_point requires at least one variable.")
    known = [v for v in variables if _post_run_spec(v) is not None]

    # A folder is a candidate if it holds *any* requested variable's field
    # -- not necessarily the same one for every phase subfolder, and not
    # necessarily every requested variable's (those missing from a given
    # folder become NaN columns in `_extract_point_in_folder`).
    fields_tried = []
    folders, seen = [], set()
    for variable in known:
        field_name = _post_run_spec(variable)[0]
        if field_name in fields_tried:
            continue
        fields_tried.append(field_name)
        for folder in _discover_folders(output_folder, field_name):
            if folder not in seen:
                seen.add(folder)
                folders.append(folder)
    if not folders:
        # None of the requested variables have their own saved field
        # anywhere -- fall back to any folder with any saved output at
        # all, so the profile still comes out (fully NaN) instead of
        # failing; only a genuinely empty output_folder raises.
        folders = _discover_any_folders(output_folder)
    if not folders:
        raise FileNotFoundError(
            f"No output found at all under '{output_folder}' or its "
            f"immediate subfolders."
        )

    per_call_output_file = output_file if len(folders) == 1 else None
    return [
        _extract_point_in_folder(folder, point, variables, name, per_call_output_file)
        for folder in folders
    ]
