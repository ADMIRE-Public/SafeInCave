# Copyright (c) 2026, The SafeInCave Developers
#
# SPDX-License-Identifier: BSD-3-Clause

"""Force-displacement composite convergence criterion."""

from typing import TYPE_CHECKING, Optional, List, Dict, Any
from .base import ConvergenceCriterion
from .force_residual import ForceResidualCriterion
from .displacement_increment import DisplacementIncrementCriterion

if TYPE_CHECKING:
    pass

class ForceDisplacementCriterion(ConvergenceCriterion):
    """
    Dual force/displacement criterion. Converged when BOTH hold:

    - ``‖R‖ / ref ≤ force_tolerance`` — free-DOF residual over the reference
      force scale (max of free-DOF external norm and internal-force norm
      including reactions), and
    - ``‖δu‖ / ‖Δu_step‖ ≤ displacement_tolerance`` — last Newton correction
      over the accumulated step displacement increment. The default (1e-2)
      matches Abaqus's displacement-correction check: the force residual is
      the physical equilibrium gate; the displacement check only guards
      against acceptance while corrections are still large. A much tighter
      value stalls steps where the active plastic set oscillates even
      though force equilibrium and yield consistency are satisfied.

    The combined scalar error is
    ``max(error_force / tol_force, error_displacement / tol_displacement)``,
    converged at ≤ 1.0, which is equivalent to satisfying both component
    tolerances.

    Where the two errors come from
    -----------------------------
    Two sources, and they must not be swapped for one another:

    - **Driver-supplied.** If the momentum equation exposes ``_newton_state`` —
      a dict carrying ``force_error``, ``disp_error`` and optionally ``ite`` —
      those values are used VERBATIM and nothing is recomputed. The incremental
      Newton driver already has both norms in hand from its residual assembly,
      so recomputing them here would cost an extra assembly and risk disagreeing
      with the Jacobian the step was solved with. That driver's ``disp_error``
      also includes the line-search factor ``alpha`` and counts owned DOFs only
      (MPI-correct).
    - **Self-computed.** Otherwise the two component criteria
      :class:`ForceResidualCriterion` and :class:`DisplacementIncrementCriterion`
      are composed, each assembling what it needs from ``momentum_eq``.

    Iteration 0 is gated separately, on ``initial_force_tolerance``
    ------------------------------------------------------------------
    ``force_tolerance`` is relative to a reference force scale that
    ``assemble_residual`` takes as ``max(‖f_ext‖, ‖f_int‖)`` over the whole
    body — a TOTAL, accumulated quantity. A step's *initial* residual, by
    contrast, is the out-of-balance force created by that step's load
    INCREMENT, so it scales with dt. The two are not comparable, and using
    ``force_tolerance`` at iteration 0 makes the solver silently skip steps:
    on the cavern2D ramp (in-situ 8.5 MPa fixing ‖f_int‖, 0.5 MPa of fresh
    cavern pressure per step) one increment gives a relative initial residual
    of 9e-7 — just under 1e-6 — so every other step was accepted with ZERO
    Newton iterations and the displacement lagged the load by one increment.
    Refining dt makes it worse, not better: halving dt halves the initial
    residual, so MORE steps fall under the tolerance and are skipped, and the
    solution stops converging under refinement altogether.

    The fix is Abaqus's own rule — an increment always takes at least one
    iteration — with one exception kept for the geostatic/held-load steps the
    Newton path relies on: a step whose load genuinely did not change starts at
    round-off (measured at ~4e-16 relative, nine orders below a real
    increment), and solving it would be pure waste. ``initial_force_tolerance``
    (1e-12 on the Newton path, i.e. round-off and nothing else) draws that line.
    Do not raise it towards ``force_tolerance``; that reinstates the skipping.

    The gate is OPT-IN: ``initial_force_tolerance=None`` (the default) disables
    it entirely, so the staggered path behaves exactly as it always has.

    The plastic-consistency gate of :class:`ConvergenceErrorHandler` applies on
    top of this criterion, unchanged.

    Parameters
    ----------
    force_tolerance : float, optional
        Tolerance for the force residual. Default: 1e-3. The incremental-Newton
        path uses 1e-6 (see ``Simulator_MNewton``).
    displacement_tolerance : float, optional
        Tolerance for the displacement correction. Default: 1e-2.
    initial_force_tolerance : float or None, optional
        Force tolerance applied at iteration 0 only. ``None`` (default) means
        no separate iteration-0 gate. The Newton path passes 1e-12.
    name : str, optional
        Name for diagnostics. Default: "force_displacement".
    """

    # Progress-table description, read by Simulator_MNewton. `compute_error`
    # returns max(force_err/force_tol, disp_err/disp_tol) -- already divided BY
    # the tolerances and converged at <= 1 -- so printing it under the stock
    # "Non-linear error" heading makes a converged step read as "94% error"
    # when the force residual is actually ~1e-10. Report the two physical
    # relative errors the criterion is built from instead.
    progress_columns = ("force err", "disp err")
    progress_formats = ("%.2e", "%.2e")

    def __init__(
        self,
        force_tolerance: float = 1e-3,
        displacement_tolerance: float = 1e-2,
        initial_force_tolerance: Optional[float] = None,
        name: Optional[str] = None,
    ):
        criterion_name = name or "force_displacement"
        # tolerance=1 means both criteria satisfy their own tolerances
        super().__init__(tolerance=1.0, name=criterion_name)
        self.force_tolerance = float(force_tolerance)
        self.displacement_tolerance = float(displacement_tolerance)
        self.initial_force_tolerance = (
            None if initial_force_tolerance is None else float(initial_force_tolerance)
        )
        self.force_criterion = ForceResidualCriterion(tolerance=force_tolerance)
        self.displacement_criterion = DisplacementIncrementCriterion(
            tolerance=displacement_tolerance
        )
        self.component_history: List[Dict[str, float]] = []
        # Last physical errors, for the progress table. Populated on both
        # source paths so the display does not care which one ran.
        self.last_force_error: float = 0.0
        self.last_disp_error: float = 0.0

    @staticmethod
    def _driver_state(momentum_eq: Any) -> Optional[Dict[str, float]]:
        """The driver-supplied error dict, or None if this equation has none."""
        state = getattr(momentum_eq, "_newton_state", None)
        if isinstance(state, dict) and "force_error" in state and "disp_error" in state:
            return state
        return None

    def _supplies_errors(self, momentum_eq: Any) -> bool:
        """Whether this equation is a driver that reports errors itself.

        Checked by attribute presence rather than by value, because at step
        start the driver has not assembled anything yet and ``_newton_state``
        is legitimately None.
        """
        return hasattr(momentum_eq, "_newton_state")

    def initialize(self, momentum_eq: Any) -> None:
        """
        Initialize at step start.

        Parameters
        ----------
        momentum_eq : LinearMomentumBase
            Momentum equation solver.

        Returns
        -------
        None
        """
        self.reset_history()
        self.component_history = []
        if self._supplies_errors(momentum_eq):
            # Driver-supplied path: clear last step's state so iteration 0 of
            # this step cannot read a stale error. The component criteria are
            # never consulted, and initializing them would copy the whole
            # displacement array into `momentum_eq.u_step_start` every step for
            # nothing (the Newton driver tracks its own step start locally).
            momentum_eq._newton_state = None
            return
        self.force_criterion.initialize(momentum_eq)
        self.displacement_criterion.initialize(momentum_eq)

    def compute_error(self, momentum_eq: Any) -> float:
        """
        Compute combined normalized error from both components.

        Parameters
        ----------
        momentum_eq : LinearMomentumBase
            Momentum equation solver.

        Returns
        -------
        float
            Combined normalized error, converged at <= 1.0.
        """
        state = self._driver_state(momentum_eq)
        if state is not None:
            force_error = float(state["force_error"])
            displacement_error = float(state["disp_error"])
            iteration = int(state.get("ite", 0))
        elif self._supplies_errors(momentum_eq):
            # Driver present but nothing assembled yet (step start).
            force_error = float("inf")
            displacement_error = 0.0
            iteration = 0
        else:
            force_error = self.force_criterion.compute_error(momentum_eq)
            displacement_error = self.displacement_criterion.compute_error(momentum_eq)
            # Before this call's append, len(history) is the iteration index.
            iteration = len(self.history)

        force_tolerance = self.force_tolerance
        if self.initial_force_tolerance is not None and iteration == 0:
            force_tolerance = self.initial_force_tolerance

        force_ratio = force_error / max(force_tolerance, 1e-16)
        displacement_ratio = displacement_error / max(
            self.displacement_tolerance, 1e-16
        )
        combined_error = max(force_ratio, displacement_ratio)

        self.last_force_error = force_error
        self.last_disp_error = displacement_error
        self.component_history.append(
            {
                "force_error": float(force_error),
                "displacement_error": float(displacement_error),
                "force_ratio": float(force_ratio),
                "displacement_ratio": float(displacement_ratio),
            }
        )
        self.history.append(float(combined_error))
        return float(combined_error)

    def is_converged(self, error: float) -> bool:
        """
        Check composite convergence based on logical combination rule.

        Parameters
        ----------
        error : float
            Combined normalized error.

        Returns
        -------
        bool
            True when error <= 1.0.
        """
        return error <= self.tolerance

    def progress_values(self, momentum_eq: Any) -> tuple:
        """The two physical relative errors, for the progress table."""
        return (self.last_force_error, self.last_disp_error)

    def get_convergence_info(self) -> Dict[str, Any]:
        """
        Return detailed convergence diagnostics.

        Returns
        -------
        dict
            Composite convergence information with per-criterion history.
        """
        return {
            "n_iterations": len(self.history),
            "criterion": self.name,
            "combined_tolerance": self.tolerance,
            "force_tolerance": self.force_tolerance,
            "displacement_tolerance": self.displacement_tolerance,
            "initial_force_tolerance": self.initial_force_tolerance,
            "combined_history": self.history,
            "component_history": self.component_history,
        }
