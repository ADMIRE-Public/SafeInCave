# Copyright (c) 2026, The SafeInCave Developers
#
# SPDX-License-Identifier: BSD-3-Clause

from __future__ import annotations
from typing import TYPE_CHECKING, Optional, Any
import math
import os
from mpi4py import MPI
from .base import ConvergenceCriterion
from .strain_based import StrainBasedCriterion
from .force_residual import ForceResidualCriterion
from .displacement_increment import DisplacementIncrementCriterion
from .force_displacement import ForceDisplacementCriterion
from .newton_residual import NewtonResidualCriterion

if TYPE_CHECKING:
    pass

"""Handler and factory for convergence criteria."""

def resolve_convergence_criterion(
    convergence_criterion: str | "ConvergenceCriterion" = "strain_based",
) -> "ConvergenceCriterion":
    """
    Resolve convergence criterion selector into concrete strategy.

    Supported names:
    - "strain_based"
    - "force_residual"
    - "displacement_increment"
    - "force_displacement" (combined force residual + displacement increment)

    Parameters
    ----------
    convergence_criterion : str
        Selector name.

    Returns
    -------
    ConvergenceCriterion
        Concrete criterion strategy instance.
    """
    if hasattr(convergence_criterion, "compute_error") and hasattr(
        convergence_criterion, "is_converged"
    ):
        return convergence_criterion

    criterion_key = str(convergence_criterion).strip().lower()

    if criterion_key == "strain_based":
        return StrainBasedCriterion()
    if criterion_key == "force_residual":
        return ForceResidualCriterion()
    if criterion_key == "displacement_increment":
        return DisplacementIncrementCriterion()
    if criterion_key == "force_displacement":
        return ForceDisplacementCriterion()
    if criterion_key == "newton_residual":
        return NewtonResidualCriterion()

    raise ValueError(
        "Unknown convergence_criterion. Supported values: "
        "'strain_based', 'force_residual', "
        "'displacement_increment', 'force_displacement', 'newton_residual'."
    )


def initialize_convergence_state(
    momentum_eq: Any,
    convergence_criterion: "ConvergenceCriterion",
) -> None:
    """
    Initialize criterion state at time-step start.

    Parameters
    ----------
    momentum_eq : LinearMomentumBase
        Momentum equation solver instance.
    convergence_criterion : Any
        Resolved criterion strategy.
    """
    convergence_criterion.initialize(momentum_eq)


def compute_error_from_criterion(
    momentum_eq: Any,
    convergence_criterion: "ConvergenceCriterion",
) -> float:
    """
    Compute convergence error from the selected criterion strategy.

    Parameters
    ----------
    momentum_eq : LinearMomentumBase
        Momentum equation solver instance.
    convergence_criterion : Any
        Resolved criterion strategy.

    Returns
    -------
    float
        Raw criterion error value.
    """
    scalar_error = float(convergence_criterion.compute_error(momentum_eq))
    return scalar_error


class ConvergenceErrorHandler:
    """
    Facade for convergence strategy selection and error evaluation.

    This class centralizes all convergence-related selection/dispatch logic so
    simulator drivers only need one import from this module.
    """

    def __init__(
        self,
        convergence_criterion: str | "ConvergenceCriterion" = "strain_based",
        tol: Optional[float] = None,
        plastic_consistency_tolerance: float = 1e-4,
    ):
        self.criterion = resolve_convergence_criterion(convergence_criterion)
        self.error: float = 0.0
        self.not_converged_error: bool = True
        self.below_max_iterations: bool = True
        self.last_raw_error: float = 0.0
        # If tol is None, convergence is evaluated with criterion.is_converged(error).
        # If tol is provided, it overrides criterion tolerance for loop control.
        self.tol: Optional[float] = tol
        self.ite: int = 0
        self.maxiter: Optional[int] = None
        # Independent equilibrium gate, applied regardless of the named
        # criterion: a plastic material's yield-consistency residual
        # (max(F)+/f_c across elems_ne exposing `consistency_error`) must
        # also be small before the step is accepted. Without this, a metric
        # like strain-based can be satisfied trivially on a small load
        # increment while the stress still violates the yield surface by a
        # non-negligible margin. Tolerance calibration on the cavern2D
        # benchmark: 1e-3 lets per-step residuals compound to a 10%+ bias
        # over a long load ramp (do not loosen); 1e-4 is validated stable
        # across the full ramp (~3.4% vs COMSOL, plateaued, ~20-25
        # iters/step); 1e-6 adds ~1.5-2x more iterations with no measurable
        # accuracy gain. Kept out of individual criteria so their reported
        # error stays a truthful, single-purpose metric -- this is a safety
        # net on top, not folded into any one criterion's number.
        self.plastic_consistency_tolerance = plastic_consistency_tolerance
        self.plastic_consistency_error: float = 0.0
        # Optional per-iteration CSV trace (SAFEINCAVE_CONV_TRACE=<path>):
        # one row per nonlinear iteration with the criterion error, the
        # consistency residual and per-model plastic diagnostics. Meant for
        # plastic-model debugging; costs nothing when the variable is unset.
        self._trace_path: Optional[str] = os.environ.get("SAFEINCAVE_CONV_TRACE")
        self._trace_step: int = 0
        self._trace_header_written: bool = False

    def initialize_step(
        self,
        momentum_eq: Any,
        maxiter: Optional[int] = None,
    ) -> None:
        """Initialize criterion state at the start of a time step."""
        self.error = 0.0
        self.not_converged_error = True
        self.ite = 0
        self.maxiter = maxiter
        self.below_max_iterations = True if maxiter is None else (self.ite < maxiter)
        self._trace_step += 1
        initialize_convergence_state(
            momentum_eq,
            self.criterion,
        )

    def compute_error(self, momentum_eq: Any) -> float:
        """Compute raw criterion error and store for diagnostics."""
        raw_error = compute_error_from_criterion(momentum_eq, self.criterion)
        self.error = float(raw_error)
        self.last_raw_error = self.error

        local_consistency = 0.0
        for elem_ne in getattr(momentum_eq.mat, "elems_ne", []):
            local_consistency = max(
                local_consistency, float(getattr(elem_ne, "consistency_error", 0.0))
            )
        self.plastic_consistency_error = float(
            momentum_eq.grid.mesh.comm.allreduce(local_consistency, op=MPI.MAX)
        )

        if os.environ.get("SAFEINCAVE_DEBUG_CONV"):
            print(
                f"[conv] ite={self.ite:3d} error={self.error:.3e} "
                f"consistency={self.plastic_consistency_error:.3e}",
                flush=True,
            )

        if self._trace_path and momentum_eq.grid.mesh.comm.rank == 0:
            self._write_trace_row(momentum_eq)

        return self.error

    def _write_trace_row(self, momentum_eq: Any) -> None:
        """Append one per-iteration diagnostics row to the trace CSV."""
        columns = ["step", "ite", "error", "plastic_consistency_error"]
        values = [
            self._trace_step,
            self.ite,
            self.error,
            self.plastic_consistency_error,
        ]
        for elem_ne in getattr(momentum_eq.mat, "elems_ne", []):
            name = getattr(elem_ne, "name", type(elem_ne).__name__)
            consistency = getattr(elem_ne, "consistency_error", None)
            if consistency is not None:
                columns.append(f"{name}_consistency")
                values.append(float(consistency))
            is_plastic = getattr(elem_ne, "is_plastic", None)
            if is_plastic is not None:
                columns.append(f"{name}_n_plastic")
                values.append(int(is_plastic.sum()))
            corner = getattr(elem_ne, "corner", None)
            if corner is not None:
                columns.append(f"{name}_n_corner")
                values.append(int(corner.sum()))
        with open(self._trace_path, "a") as f:
            if not self._trace_header_written:
                if f.tell() == 0:
                    f.write(",".join(columns) + "\n")
                self._trace_header_written = True
            f.write(",".join(str(v) for v in values) + "\n")

    def evaluate(
        self,
        momentum_eq: Any,
        ite: Optional[int] = None,
        maxiter: Optional[int] = None,
    ) -> bool:
        """Compute error and refresh the convergence-state boolean."""
        self.compute_error(momentum_eq)
        if ite is not None or maxiter is not None:
            self.update_below_max_iterations(ite=ite, maxiter=maxiter)
        return self.update_not_converged_error()

    def update_not_converged_error(self) -> bool:
        """Update and return the current convergence-state boolean."""
        if not math.isfinite(self.error) or not math.isfinite(self.plastic_consistency_error):
            # NaN/Inf means the nonlinear iteration is unstable and must continue/retry.
            self.not_converged_error = True
        elif self.tol is None:
            self.not_converged_error = not self.criterion.is_converged(self.error)
        else:
            self.not_converged_error = self.error > self.tol

        if not self.not_converged_error:
            # Equilibrium gate: don't accept a step where the criterion's
            # metric (e.g. strain change) looks converged but a plastic
            # material's yield-consistency residual is still large.
            if self.plastic_consistency_error > self.plastic_consistency_tolerance:
                self.not_converged_error = True

        return self.not_converged_error

    def update_below_max_iterations(
        self,
        ite: Optional[int] = None,
        maxiter: Optional[int] = None,
    ) -> bool:
        """Update and return whether the current iteration is below max iterations."""
        if maxiter is not None:
            self.maxiter = maxiter
        if ite is not None:
            self.ite = ite

        if self.maxiter is None:
            self.below_max_iterations = True
        else:
            self.below_max_iterations = self.ite < self.maxiter

        return self.below_max_iterations

    def increment_iteration(self) -> int:
        """Increment iteration counter and refresh max-iteration state."""
        self.ite += 1
        self.update_below_max_iterations(ite=self.ite)
        return self.ite

    def not_converged(self) -> bool:
        """Return unified nonlinear-loop continuation condition."""
        return self.not_converged_error and self.below_max_iterations

    def get_tolerance(self) -> float:
        """Return active criterion tolerance used for convergence checks."""
        if self.tol is not None:
            return float(self.tol)
        return float(self.criterion.tolerance)

    def get_last_raw_error(self) -> float:
        """Return most recent raw criterion value (not zeroed by convergence logic)."""
        return self.last_raw_error

    def set_criterion(self, convergence_criterion: str) -> None:
        """Swap convergence strategy at runtime."""
        self.criterion = resolve_convergence_criterion(convergence_criterion)
