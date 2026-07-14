# Copyright (c) 2026, The SafeInCave Developers
#
# SPDX-License-Identifier: BSD-3-Clause

"""Displacement increment convergence criterion."""

from typing import Any,  TYPE_CHECKING, Optional
import torch as to
import numpy as np
from mpi4py import MPI
from .base import ConvergenceCriterion, _is_momentum_solver_instance

if TYPE_CHECKING:
    from ...MomentumEquation import LinearMomentumBase

def _initialize_step_displacement(momentum_eq: Any) -> to.Tensor:
    """
    Capture and store displacement field at time step start.

    Captures the displacement state at the beginning of a time step to enable
    computation of the cumulative displacement increment throughout the step.
    Used by displacement correction criterion for robust convergence checking.

    Parameters
    ----------
    momentum_eq : LinearMomentumBase
        Momentum equation solver instance.

    Returns
    -------
    torch.Tensor
        Copy of displacement array at time step start.

    Notes
    -----
    - Must be called ONCE per time step, before nonlinear iteration loop
    - Stores displacement as a FEniCS vector array copy (numpy format)
    - Used internally by DisplacementIncrementCriterion
    - MPI-safe: each rank stores its own DOFs plus ghost elements
    """
    # Capture displacement state at step start for total increment tracking
    displacement_step_start = momentum_eq.u.x.array.copy()
    momentum_eq.u_step_start = displacement_step_start
    return displacement_step_start

class DisplacementIncrementCriterion(ConvergenceCriterion):
    """
    Displacement increment convergence criterion (Newton step measure).

    Error metric (Newton step relative to total increment):
    $$\\text{error}_\\text{displacement} = \\frac{\\|\\Delta u_\\text{correction}\\|}{\\|\\Delta u_\\text{total}\\|}$$

    where:
    - Δu_correction = ||u^(k+1) - u^(k)|| (Newton iteration step in current iteration)
    - Δu_total = ||u^(k+1) - u_step_start|| (cumulative displacement since step start)

    **Purpose**:
    - Measures iteration progress relative to step size
    - Scale-independent (correction vs. total, not vs. old)
    - Detects convergence stalling (very small corrections)

    **Advantages**:
    - Scale-robust (independent of problem size)
    - Captures iteration difficulty relative to step magnitude
    - Robust for large deformations

    **Disadvantages**:
    - Requires tracking u_step_start (requires initialize_step_displacements())
    - Can be strict for small steps (ratio may remain large)
    - Should be combined with residual criterion

    **Typical Tolerance**: 1e-2 (1% of total step increment)

    Parameters
    ----------
    tolerance : float, optional
        Displacement correction tolerance. Default: 1e-2.
    """

    def __init__(self, tolerance: float = 1e-2):
        super().__init__(tolerance=tolerance, name="DisplacementIncrement")
        self.u_step_start: Optional[to.Tensor] = None
        self.u_previous: Optional[to.Tensor] = None

    def initialize(self, momentum_eq: Any) -> None:
        """
        Initialize displacement criterion and capture starting displacement.

        Called once per time step to store the initial displacement state
        for tracking the cumulative increment throughout the step.

        Parameters
        ----------
        momentum_eq : LinearMomentumBase
            Momentum equation solver instance.

        Returns
        -------
        None

        Side Effects
        -----------
        Stores initial displacement via _initialize_step_displacement().
        """
        self.u_step_start = _initialize_step_displacement(momentum_eq)
        self.u_previous = to.as_tensor(momentum_eq.u.x.array.copy())
        self.reset_history()

    def compute_error(self, *args, **kwargs) -> float:
        """
        Compute displacement correction error.

        Parameters
        ----------
        Supports two call styles:
        1) compute_error(momentum_eq, u_new, u_old)
        2) compute_error(momentum_eq)  # u_old tracked internally

        Returns
        -------
        float
            Displacement correction error. Dimensionless, in [0, ∞).
        """
        if len(args) == 3 and _is_momentum_solver_instance(args[0]):
            momentum_eq = args[0]
            u_new = args[1]
            u_old = args[2]
            self.u_previous = u_new.clone()
        elif len(args) == 1 and _is_momentum_solver_instance(args[0]):
            momentum_eq = args[0]
            u_new = to.as_tensor(momentum_eq.u.x.array.copy())
            if self.u_previous is None:
                u_old = u_new
            else:
                u_old = self.u_previous
            self.u_previous = u_new.clone()
        else:
            raise ValueError(
                "DisplacementIncrementCriterion.compute_error expects "
                "(momentum_eq) or (momentum_eq, u_new, u_old)."
            )

        # Compute Newton iteration correction (single step)
        displacement_correction = u_new - u_old
        correction_norm = float(to.linalg.norm(displacement_correction))

        # Compute cumulative displacement increment from step start
        displacement_total = u_new - self.u_step_start
        total_increment_norm = float(to.linalg.norm(displacement_total))

        # Degenerate step: displacement did not move at all relative to the
        # overall solution scale (e.g. a constant-load step already in
        # equilibrium). Both norms are ~zero; flooring both to 1e-16 would
        # yield a spurious ratio of 1.0 that never converges. Declare
        # converged instead.
        solution_scale = max(float(to.linalg.norm(u_new)), 1e-16)
        if total_increment_norm <= 1e-12 * solution_scale or total_increment_norm <= 1e-30:
            error_value = 0.0
        else:
            # Ratio: single step vs. cumulative (scale-independent)
            error_value = correction_norm / total_increment_norm

        self.history.append(error_value)
        return float(error_value)

    def is_converged(self, error: float) -> bool:
        """Check if displacement correction error is below tolerance."""
        return error <= self.tolerance
