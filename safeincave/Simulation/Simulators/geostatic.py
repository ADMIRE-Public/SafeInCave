# Copyright (c) 2026, The SafeInCave Developers
#
# SPDX-License-Identifier: BSD-3-Clause

"""
Conventional geostatic equilibrium step -- the SafeInCave equivalent of
Abaqus' ``*GEOSTATIC`` procedure.

A geostatic step answers one question: *is the specified initial stress in
equilibrium with the applied loads and boundary conditions?* It then leaves
the model in a state fit to be the reference configuration for the loading
stages that follow:

1. **Diagnose.** Assemble the out-of-balance force at the current state
   (``u = 0`` with ``sigma_0`` applied) *without* solving, and report
   ``||R||`` and ``||R|| / ||f||``. A self-equilibrated ``sigma_0`` with
   matching boundary tractions gives a ratio at round-off; anything larger is
   a real force imbalance that the step is about to resolve into
   deformation.
2. **Equilibrate.** Solve to equilibrium.
3. **Commit as reference.** Set ``sigma_0 <- sigma_converged`` and zero every
   deformation quantity, so subsequent stages measure deformation relative to
   the equilibrated in-situ state while the stress field is preserved and
   still self-equilibrated.

Not every geostatic step is meant to have a zero residual at step 1. Leaving
an excavated surface traction-free *is* the excavation load, so an imbalance
there is the physics, not a mistake -- hence the diagnostic warns and
proceeds rather than failing.

Works with any compatible momentum equation/simulator
--------------------------------------------------------------
This class is solver-agnostic: it drives ``eq_mom`` through an internal
``simulator_cls`` (default :class:`Simulator_M`, the staggered path), but any
compatible ``Simulator`` subclass can be passed instead -- e.g. an
incremental-Newton simulator, for a momentum equation that supports it. It
never hard-imports a Newton-specific class; it only dispatches on which
members ``eq_mom``/the report happen to expose. Two things adapt
automatically:

* **The stress/count fields read throughout** (``check_equilibrium``,
  ``solve``, ``commit_as_reference``) prefer ``sig_state``/``n_state_points``
  when present (a higher-resolution quadrature/state-point field some
  formulations carry alongside the plain cell-averaged ``sig``/``n_elems``),
  falling back to ``sig``/``n_elems`` otherwise.
* **The equilibrium diagnostic** prefers ``eq_mom.assemble_residual()`` when
  present (a formulation's own native, already-tested residual), falling
  back to the generic force-residual assembly in
  :mod:`..Convergence.residuals` otherwise.
* **The rate reseed on commit** prefers ``eq_mom.initialize_rates()`` when
  present, falling back to the ``compute_eps_ne_rate``/
  ``update_eps_ne_rate_old`` idiom :class:`Simulator_M` itself uses at
  ``t=0``.
* **The deformation fields zeroed on commit** always include ``u``, ``X``,
  ``eps_tot``, plus ``du_sol``/``eps_tot_state`` when present.

Why this always runs with ``compute_elastic_response=False``
--------------------------------------------------------------
Two defects in the public elastic pre-solve corrupt exactly the state a
geostatic step exists to produce:

* ``Equations/Momentum/standard.py``'s ``solve_elastic_response`` assembles
  the initial-stress term with a ``+`` sign where the consistent form (the
  one ``compute_stress`` and the transient ``solve`` both use) is a ``-``,
  so a self-equilibrated ``sigma_0`` yields roughly twice the wrong elastic
  displacement instead of ``u = 0``;
* ``compute_elastic_stress`` (``stress = C:eps_e``, no ``sigma_0`` term) then
  drops the initial stress from the committed stress state entirely.

Both are only reachable when ``compute_elastic_response=True``. This module
never sets it, so the bugs are made unreachable here, not patched -- they
should be fixed upstream separately.

Scope: mechanical only -- thermo-mechanical simulators
(``Simulator_TM``/-Newton) are not supported.
"""

from __future__ import annotations

import inspect
import warnings
from dataclasses import dataclass
from typing import Any, Sequence

import numpy as np
import torch as to
from mpi4py import MPI

from .base import Simulator
from .mechanical import Simulator_M
from ..Convergence.residuals import _compute_force_residual, _compute_vector_norm
from ...Cavern import CavernHandler
from ...Utils import numpy2torch

__all__ = ["GeostaticStep", "GeostaticReport", "check_equilibrium"]

#: Members a momentum equation must expose to be driven by a geostatic step.
_REQUIRED_MEMBERS = (
    "apply_initial_stress",
    "compute_eps_ne_rate",
    "update_eps_ne_rate_old",
    "n_elems",
    "get_uV",
)

#: Reference-scale floor below which a residual ratio carries no information.
_REF_FLOOR = 1e-15

#: Constitutive fields snapshotted/restored around a native residual probe
#: (``assemble_residual`` mutates some of these as a side effect).
_CONSTITUTIVE_FIELDS = ("sig_state", "sig", "Dep", "eps_tot", "eps_tot_state")

#: Deformation fields zeroed on commit; the last two are Newton/quadrature-only.
_DEFORMATION_FIELDS = ("u", "X", "eps_tot", "du_sol", "eps_tot_state")


def _require_momentum_interface(eq_mom: Any) -> None:
    """Raise unless `eq_mom` implements the momentum-equation interface."""
    missing = [name for name in _REQUIRED_MEMBERS if not hasattr(eq_mom, name)]
    if missing:
        raise TypeError(
            f"{type(eq_mom).__name__} does not implement the momentum-equation "
            f"interface required by GeostaticStep (missing: {', '.join(missing)}). "
            f"GeostaticStep needs a LinearMomentumBase-family equation, whether "
            f"driven by the staggered path or an incremental-Newton one."
        )


def _stress_attr(eq_mom: Any) -> str:
    """Preferred stress-field attribute name: state-point if present."""
    return "sig_state" if hasattr(eq_mom, "sig_state") else "sig"


def _count_attr(eq_mom: Any) -> str:
    """Preferred count attribute name: state-point if present."""
    return "n_state_points" if hasattr(eq_mom, "n_state_points") else "n_elems"


def _snapshot_constitutive(eq_mom: Any) -> dict:
    """Copy every field a native residual probe might mutate."""
    snap: dict = {"sigma_k": None, "functions": {}}
    for name in _CONSTITUTIVE_FIELDS:
        fn = getattr(eq_mom, name, None)
        if fn is not None:
            snap["functions"][name] = fn.x.array.copy()
    sigma_k = getattr(eq_mom, "_sigma_k", None)
    snap["sigma_k"] = None if sigma_k is None else sigma_k.clone()
    return snap


def _restore_constitutive(eq_mom: Any, snap: dict) -> None:
    """Undo :func:`_snapshot_constitutive`, ghosts included."""
    for name, values in snap["functions"].items():
        fn = getattr(eq_mom, name)
        fn.x.array[:] = values
        fn.x.scatter_forward()
    if hasattr(eq_mom, "_sigma_k"):
        eq_mom._sigma_k = snap["sigma_k"]


def _filter_kwargs(cls: type, kwargs: dict) -> dict:
    """Keep only the keys `cls.__init__` actually accepts."""
    params = inspect.signature(cls.__init__).parameters
    if any(p.kind is inspect.Parameter.VAR_KEYWORD for p in params.values()):
        return dict(kwargs)
    return {k: v for k, v in kwargs.items() if k in params}


def _u_max(eq_mom: Any) -> float:
    """``||u||_inf`` over locally owned DOFs, reduced across ranks."""
    V = eq_mom.get_uV()
    n_owned = V.dofmap.index_map.size_local * V.dofmap.index_map_bs
    owned = eq_mom.u.x.array[:n_owned]
    local = float(np.abs(owned).max()) if owned.size else 0.0
    return float(eq_mom.grid.mesh.comm.allreduce(local, op=MPI.MAX))


def _verdict(r_norm: float, ref_norm: float, tolerance: float) -> tuple:
    """``(ratio, in_equilibrium, degenerate)`` from a residual/reference pair."""
    degenerate = ref_norm <= _REF_FLOOR
    ratio = r_norm / ref_norm
    in_equilibrium = (r_norm <= tolerance) if degenerate else (ratio <= tolerance)
    return ratio, in_equilibrium, degenerate


@dataclass(frozen=True)
class GeostaticReport:
    """
    Result of a geostatic diagnostic, solve, or commit.

    Attributes
    ----------
    r_norm, ref_norm, ratio
        Out-of-balance force norm, reference force scale, and their ratio.
    in_equilibrium
        Whether ``ratio`` (or ``r_norm`` when the scale is degenerate) is
        within ``tolerance``.
    tolerance
        Threshold used for the verdict.
    stage
        ``"check"`` or ``"solve"``.
    iterations
        Nonlinear iterations used by the last solved time step (``solve``
        stage only).
    converged
        Whether the solve converged.
    u_max
        ``||u||_inf`` at the time this report was produced.
    initial_ratio
        The pre-solve ``ratio``, carried through from ``check()``.
    degenerate_scale
        Whether the reference force scale was at the numerical floor.
    """

    r_norm: float
    ref_norm: float
    ratio: float
    in_equilibrium: bool
    tolerance: float
    stage: str = "check"
    iterations: int = 0
    converged: bool = True
    u_max: float = 0.0
    initial_ratio: float | None = None
    degenerate_scale: bool = False

    def summary(self) -> str:
        """Human-readable one-block digest of this report."""
        lines = [f"Geostatic {self.stage} report:"]
        lines.append(
            f"  ||R|| = {self.r_norm:.6e}, ||f|| = {self.ref_norm:.6e}, "
            f"ratio = {self.ratio:.6e} (tolerance = {self.tolerance:.1e})"
        )
        if self.degenerate_scale:
            lines.append("  reference force scale is at the numerical floor")
        lines.append(f"  in_equilibrium = {self.in_equilibrium}")
        if self.stage == "solve":
            lines.append(
                f"  converged = {self.converged}, iterations = {self.iterations}, "
                f"u_max = {self.u_max:.6e}"
            )
            if self.initial_ratio is not None:
                lines.append(f"  initial ratio (pre-solve) = {self.initial_ratio:.6e}")
        return "\n".join(lines)


def check_equilibrium(
    eq_mom: Any,
    tolerance: float = 1e-8,
    caverns: Any | None = None,
    t: float = 0.0,
    warn: bool = True,
    restore: bool = True,
) -> GeostaticReport:
    """
    Abaqus-style geostatic diagnostic: report the force imbalance at the
    current state without solving.

    Parameters
    ----------
    eq_mom
        A ``LinearMomentum``-like momentum equation with its material, body
        force, boundary conditions and initial stress already set.
    tolerance
        Threshold on ``||R|| / ||f||`` for the `in_equilibrium` verdict.
    caverns
        Cavern handler whose loads participate in the external force.
        Defaults to an empty handler -- correct when cavern surfaces carry
        plain ``NeumannBC`` tractions, as a geostatic stage normally wants.
    t
        Time at which to evaluate the boundary conditions.
    warn
        Emit a ``RuntimeWarning`` when not in equilibrium.
    restore
        Leave the constitutive state untouched (a pure probe). The residual
        assembly writes the probed stress into the equation as a side
        effect; pass False only when the caller wants that write to persist.

    Returns
    -------
    GeostaticReport
        With ``stage="check"``.
    """
    eq_mom.bc.update_dirichlet(t)
    eq_mom.bc.update_neumann(t)
    eq_mom.bc.update_cavern_bcs(caverns if caverns is not None else CavernHandler())

    if hasattr(eq_mom, "assemble_residual"):
        # Native diagnostic: the formulation's own residual, already tested
        # against its own solve path. constitutive_update mutates several
        # fields as a side effect, so snapshot/restore around it.
        snap = _snapshot_constitutive(eq_mom) if restore else None
        try:
            eq_mom.constitutive_update(0.0)
            r_norm, ref_norm = eq_mom.assemble_residual()
        finally:
            if snap is not None:
                _restore_constitutive(eq_mom, snap)
    else:
        stress_attr, count_attr = _stress_attr(eq_mom), _count_attr(eq_mom)
        stress_field = getattr(eq_mom, stress_attr)
        n = getattr(eq_mom, count_attr)
        snap = stress_field.x.array.copy() if restore else None
        try:
            stress_to = numpy2torch(stress_field.x.array.reshape((n, 3, 3)))
            residual_vector, ref_norm = _compute_force_residual(eq_mom, stress_to)
            r_norm = _compute_vector_norm(residual_vector)
        finally:
            if snap is not None:
                stress_field.x.array[:] = snap
                stress_field.x.scatter_forward()

    ratio, in_equilibrium, degenerate = _verdict(r_norm, ref_norm, tolerance)

    if warn and not in_equilibrium:
        warnings.warn(
            f"Geostatic equilibrium check: initial stress is NOT in equilibrium "
            f"with the applied loads/BCs (||R||/||f|| = {ratio:.3e} > "
            f"{tolerance:.1e}). This is expected when a boundary is "
            f"deliberately unbalanced -- e.g. a traction-free excavated "
            f"surface, where the imbalance IS the excavation load -- and the "
            f"step will resolve it into deformation. It indicates a mistake "
            f"only if you intended a zero-deformation in-situ state.",
            RuntimeWarning,
            stacklevel=2,
        )

    return GeostaticReport(
        r_norm=r_norm,
        ref_norm=ref_norm,
        ratio=ratio,
        in_equilibrium=in_equilibrium,
        tolerance=tolerance,
        stage="check",
        degenerate_scale=degenerate,
    )


class GeostaticStep(Simulator):
    """
    Conventional geostatic equilibrium step. Works with any compatible
    momentum equation/simulator (see the module docstring for how dispatch
    works) -- the staggered path by default, or e.g. an incremental-Newton
    one via ``simulator_cls``.

    Examples
    --------
    >>> mom_eq.apply_initial_stress(stress_0)
    >>> geo = GeostaticStep(mom_eq, t_control, outputs)
    >>> report = geo.run()          # check -> solve -> commit_as_reference
    >>> print(report.summary())

    After :meth:`run` the model carries the equilibrated stress with zero
    displacement, strain and inelastic history, so the next stage should be
    driven with ``compute_elastic_response=False`` to continue from it.

    Constructor parameter names deliberately mirror :class:`Simulator_M` so
    that the YAML transpiler's auto-wiring (``eq_mom``, ``t_control``,
    ``outputs``, ``caverns``, ``simulation_logger``) and its
    ``_check_simulator_wiring`` validation work with no registry changes;
    this class is usable directly as a ``stages[].simulator:`` block.

    Parameters
    ----------
    eq_mom
        Momentum equation (e.g. ``LinearMomentum`` or a compatible
        incremental-Newton formulation) implementing ``apply_initial_stress``,
        ``compute_eps_ne_rate``, ``update_eps_ne_rate_old``, ``n_elems``,
        ``get_uV``.
    t_control
        Time controller for the internal equilibrium solve.
    outputs
        ``SaveFields``-like output handlers for the step.
    caverns
        Cavern handler shared by the diagnostic and the solve. Defaults to a
        fresh empty handler; keep it empty so a geostatic stage does not run
        cavern thermodynamics against a relaxing volume.
    tolerance
        Equilibrium threshold on ``||R|| / ||f||``.
    maxiter
        Maximum nonlinear iterations per time step.
    simulator_cls
        The ``Simulator`` subclass used internally to drive ``eq_mom`` to
        equilibrium (default :class:`Simulator_M`, the staggered path). Pass
        a compatible incremental-Newton simulator to drive a Newton-capable
        ``eq_mom`` through genuine Newton convergence instead. Only the
        constructor parameters ``simulator_cls`` actually declares are
        forwarded to it (extras below are silently dropped if unsupported).
    convergence_criterion
        Forwarded to the internal simulator when given; left unset (``None``,
        the default) so each ``simulator_cls`` picks its own natural default
        (e.g. ``Simulator_M``'s ``"strain_based"`` vs. ``Simulator_MNewton``'s
        ``NewtonResidualCriterion``) -- do not hardcode one here, it would
        silently override the other simulator's correct default.
    simulation_logger, merged_solutions, smooth_output,
    plastic_consistency_tolerance
        Forwarded to the internal simulator.
    warn_on_imbalance
        Whether :meth:`check` warns when the initial state is out of
        equilibrium. Set False for excavation-style stages where the
        imbalance is intentional and the warning is just noise.

    Notes
    -----
    ``solve``'s reported ``iterations`` is the last time step's nonlinear
    iteration count (neither ``Simulator_M`` nor most simulators expose a
    per-step trace hook to track a peak across a multi-step ramp); for the
    recommended single-step geostatic usage this is the same number.
    """

    def __init__(
        self,
        eq_mom: Any,
        t_control: Any,
        outputs: Sequence[Any] = (),
        caverns: Any | None = None,
        tolerance: float = 1e-8,
        maxiter: int = 40,
        simulator_cls: type = Simulator_M,
        convergence_criterion: Any | None = None,
        simulation_logger: Any | None = None,
        merged_solutions: bool = False,
        smooth_output: bool = False,
        plastic_consistency_tolerance: float = 1e-4,
        warn_on_imbalance: bool = True,
    ):
        _require_momentum_interface(eq_mom)
        self.eq_mom = eq_mom
        self.t_control = t_control
        self.outputs = list(outputs)
        self.caverns = caverns if caverns is not None else CavernHandler()
        self.tolerance = float(tolerance)
        self.maxiter = int(maxiter)
        self.simulator_cls = simulator_cls
        self.convergence_criterion = convergence_criterion
        self.simulation_logger = simulation_logger
        self.merged_solutions = merged_solutions
        self.smooth_output = smooth_output
        self.plastic_consistency_tolerance = plastic_consistency_tolerance
        self.warn_on_imbalance = warn_on_imbalance
        #: Reason the last solve failed, when it did.
        self.failure_reason: str | None = None

    # ------------------------------------------------------------------

    def check(self) -> GeostaticReport:
        """Equilibrium diagnostic at the current state; mutates nothing."""
        return check_equilibrium(
            self.eq_mom,
            tolerance=self.tolerance,
            caverns=self.caverns,
            t=self.t_control.t,
            warn=self.warn_on_imbalance,
        )

    def solve(self, _pre: GeostaticReport | None = None) -> GeostaticReport:
        """
        Solve to equilibrium through ``simulator_cls``.

        A self-equilibrated initial stress converges with ``iterations == 0``
        and ``u_max`` at round-off; a genuine imbalance (excavation) is a
        normal nonlinear solve.
        """
        eq = self.eq_mom
        if _pre is None:
            _pre = check_equilibrium(
                eq,
                tolerance=self.tolerance,
                caverns=self.caverns,
                t=self.t_control.t,
                warn=False,
            )

        # The whole point: this branch performs no elastic pre-solve, so the
        # buggy solve_elastic_response/compute_elastic_stress path is never
        # reached and the already-correct sigma_0 state is carried straight
        # in. Only the kwargs simulator_cls actually declares are forwarded,
        # so any compatible Simulator subclass (staggered or Newton) works
        # from this one call.
        options = dict(
            caverns=self.caverns,
            compute_elastic_response=False,
            maxiter=self.maxiter,
            simulation_logger=self.simulation_logger,
            merged_solutions=self.merged_solutions,
            smooth_output=self.smooth_output,
            plastic_consistency_tolerance=self.plastic_consistency_tolerance,
        )
        if self.convergence_criterion is not None:
            # Only forward when explicitly set -- otherwise let simulator_cls
            # apply its own correct default (see __init__'s docstring).
            options["convergence_criterion"] = self.convergence_criterion
        kwargs = _filter_kwargs(self.simulator_cls, options)
        sim = self.simulator_cls(eq, self.t_control, self.outputs, **kwargs)

        converged = True
        self.failure_reason = None
        try:
            sim.run()
        except RuntimeError as exc:
            converged = False
            self.failure_reason = str(exc)

        iterations = int(sim.convergence_handler.ite)
        post = check_equilibrium(
            eq,
            tolerance=self.tolerance,
            caverns=self.caverns,
            t=self.t_control.t,
            warn=False,
        )

        return GeostaticReport(
            r_norm=post.r_norm,
            ref_norm=post.ref_norm,
            ratio=post.ratio,
            in_equilibrium=post.in_equilibrium,
            tolerance=self.tolerance,
            stage="solve",
            iterations=iterations,
            converged=converged,
            u_max=_u_max(eq),
            initial_ratio=_pre.ratio,
            degenerate_scale=post.degenerate_scale,
        )

    def commit_as_reference(self) -> GeostaticReport:
        """
        Make the equilibrated state the reference configuration.

        Sets ``sigma_0 <- sigma_converged`` and zeros displacement and total
        strain, so the next stage starts from zero deformation while
        carrying the geostatic stress.

        Returns
        -------
        GeostaticReport
            Post-commit diagnostic (``stage="check"``).
        """
        eq = self.eq_mom

        # The equilibrated stress. Prefer the state-point field when present
        # (full resolution), not the cell-averaged mirror. Clone before
        # handing it back in: numpy2torch may alias the dolfinx array that
        # apply_initial_stress is about to overwrite.
        stress_field = getattr(eq, _stress_attr(eq))
        n = getattr(eq, _count_attr(eq))
        sigma = numpy2torch(stress_field.x.array.reshape((n, 3, 3))).clone()
        # Writes every mirror: eps_0, sig0 (torch and FE) and sig itself.
        eq.apply_initial_stress(sigma)

        # Zero the deformation side. u and X may be the same object
        # (split_solution aliases self.u = self.X on the staggered path);
        # zeroing both is harmless either way. du_sol/eps_tot_state only
        # exist on Newton/quadrature formulations.
        for name in _DEFORMATION_FIELDS:
            fn = getattr(eq, name, None)
            if fn is not None:
                fn.x.array[:] = 0.0
                fn.x.scatter_forward()

        # Zero inelastic history. A no-op for an elastic-only geostatic
        # stage, but keeps the reference state honest when plasticity/creep
        # was active during equilibration.
        for elem in getattr(eq.mat, "elems_ne", []):
            for attr in (
                "eps_ne_old",
                "eps_ne_k",
                "eps_ne_rate",
                "eps_ne_rate_old",
            ):
                tensor = getattr(elem, attr, None)
                if tensor is not None:
                    setattr(elem, attr, to.zeros_like(tensor))

        # Reseed creep rate carryover at the new sigma_0. Prefer a native
        # initialize_rates() when the formulation has one; otherwise fall
        # back to the two-call idiom Simulator_M.run() itself uses at t=0.
        initialize_rates = getattr(eq, "initialize_rates", None)
        if callable(initialize_rates):
            initialize_rates()
        else:
            eq.compute_eps_ne_rate(sigma, self.t_control.t)
            eq.update_eps_ne_rate_old()

        # Refresh derived output fields so anything written after the commit
        # reflects the reference state. Not every momentum equation exposes
        # every one of these (e.g. compute_yield_mode is plasticity-only).
        for name in (
            "compute_p_elems",
            "compute_q_elems",
            "compute_yield_mode",
            "compute_principal_stresses",
            "compute_principal_strains",
        ):
            fn = getattr(eq, name, None)
            if callable(fn):
                fn()

        return check_equilibrium(
            eq,
            tolerance=self.tolerance,
            caverns=self.caverns,
            t=self.t_control.t,
            warn=False,
        )

    def run(self, commit: bool = True) -> GeostaticReport:
        """
        Full geostatic step: diagnose, equilibrate, commit as reference.

        Returns the post-commit diagnostic, carrying the solve's iteration
        count, converged ``u_max`` and the pre-solve imbalance. On
        non-convergence the solve report is returned unchanged and nothing
        is committed -- inspect ``converged`` and ``failure_reason``.
        """
        pre = self.check()
        solved = self.solve(_pre=pre)
        if not solved.converged or not commit:
            report = solved
        else:
            post = self.commit_as_reference()
            report = GeostaticReport(
                r_norm=post.r_norm,
                ref_norm=post.ref_norm,
                ratio=post.ratio,
                in_equilibrium=post.in_equilibrium,
                tolerance=self.tolerance,
                stage="solve",
                iterations=solved.iterations,
                converged=True,
                # commit_as_reference() zeros u, so the pre-commit solved
                # value is the informative one to report here.
                u_max=solved.u_max,
                initial_ratio=pre.ratio,
                degenerate_scale=post.degenerate_scale,
            )

        if self.eq_mom.grid.mesh.comm.rank == 0:
            print(report.summary())
        return report
