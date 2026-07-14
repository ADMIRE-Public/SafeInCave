# Copyright (c) 2026, The SafeInCave Developers
#
# SPDX-License-Identifier: BSD-3-Clause

"""Force residual convergence criterion."""

from typing import Any,  TYPE_CHECKING, Optional
import torch as to
import numpy as np
from mpi4py import MPI
from .base import ConvergenceCriterion, _is_momentum_solver_instance
from .residuals import (
    _compute_force_residual,
    _compute_vector_norm,
)

if TYPE_CHECKING:
    from ...MomentumEquation import LinearMomentumBase

class ForceResidualCriterion(ConvergenceCriterion):
    """
    Force residual convergence criterion (equilibrium check).

    Error metric (normalized out-of-balance forces):
    $$\\text{error}_\\text{residual} = \\frac{\\|R\\|}{\\|P_\\text{ext}\\|}$$

    where R = P_ext - q_int (out-of-balance forces).

    **Purpose**:
    - Checks mechanical equilibrium
    - Detects stalled iterations (large residual despite small displacement)
    - Robust for problems with varying load magnitudes

    **Advantages**:
    - Direct equilibrium measure
    - Scale-robust (normalized by load)
    - Detects residual stalling

    **Disadvantages**:
    - Requires residual assembly (expensive)
    - Must be combined with displacement criterion for full convergence check

    **Typical Tolerance**: 1e-3 to 1e-4 (0.1% to 0.01% of load magnitude)

    Parameters
    ----------
    tolerance : float, optional
        Residual tolerance (as fraction of external load). Default: 1e-3.
    """

    def __init__(self, tolerance: float = 1e-3):
        super().__init__(tolerance=tolerance, name="ForceResidual")

    def initialize(self, momentum_eq: Any) -> None:
        """Initialize residual criterion."""
        self.reset_history()

    def compute_error(self, *args, **kwargs) -> float:
        """
        Compute residual-based error: ||R|| / ||P_ext||.

        Parameters
        ----------
        Supports two call styles:
        1) compute_error(momentum_eq, stress)
        2) compute_error(momentum_eq)  # stress inferred from momentum_eq.sig

        Returns
        -------
        float
            Residual error. Dimensionless, typically in [0, ∞).
        """
        if len(args) == 2 and _is_momentum_solver_instance(args[0]):
            momentum_eq = args[0]
            stress = args[1]
        elif len(args) == 1 and _is_momentum_solver_instance(args[0]):
            momentum_eq = args[0]
            stress = to.as_tensor(
                momentum_eq.sig.x.array.reshape((momentum_eq.n_elems, 3, 3))
            )
        else:
            raise ValueError(
                "ForceResidualCriterion.compute_error expects (momentum_eq) or "
                "(momentum_eq, stress)."
            )

        # Compute residual R = P_ext - q_int and the reference force scale
        # (max of free-DOF external load and full internal force incl.
        # reactions -- robust for displacement-controlled and fully
        # confined load cases where the free-DOF external load is ~zero).
        residual_vector, reference_norm = _compute_force_residual(momentum_eq, stress)
        residual_norm = _compute_vector_norm(residual_vector)

        error_value = residual_norm / max(reference_norm, 1e-16)
        self.history.append(error_value)
        return float(error_value)

    def is_converged(self, error: float) -> bool:
        """Check if residual error is below tolerance."""
        return error <= self.tolerance
