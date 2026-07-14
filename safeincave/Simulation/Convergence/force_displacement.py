# Copyright (c) 2026, The SafeInCave Developers
#
# SPDX-License-Identifier: BSD-3-Clause

"""Force-displacement composite convergence criterion."""

from typing import TYPE_CHECKING, Optional, List, Dict, Any
import torch as to
from .base import ConvergenceCriterion
from .force_residual import ForceResidualCriterion
from .displacement_increment import DisplacementIncrementCriterion

if TYPE_CHECKING:
    from ...MomentumEquation import LinearMomentumBase

class ForceDisplacementCriterion(ConvergenceCriterion):
    """
    Explicit combined criterion: force residual + displacement increment.

    The combined scalar error is defined as:
    max(error_force / tol_force, error_displacement / tol_displacement)

    Convergence is achieved when the combined error is <= 1.0, which is
    equivalent to satisfying both component tolerances.

    Parameters
    ----------
    force_tolerance : float, optional
        Tolerance for force residual criterion. Default: 1e-3.
    displacement_tolerance : float, optional
        Tolerance for displacement increment criterion. Default: 1e-2.
    name : str, optional
        Name for diagnostics. Default: "force_displacement".
    """

    def __init__(
        self,
        force_tolerance: float = 1e-3,
        displacement_tolerance: float = 1e-2,
        name: Optional[str] = None,
    ):
        criterion_name = name or "force_displacement"
        # tolerance=1 means both criteria satisfy their own tolerances
        super().__init__(tolerance=1.0, name=criterion_name)
        self.force_criterion = ForceResidualCriterion(tolerance=force_tolerance)
        self.displacement_criterion = DisplacementIncrementCriterion(
            tolerance=displacement_tolerance
        )
        self.component_history: List[Dict[str, float]] = []

    def initialize(self, momentum_eq: Any) -> None:
        """
        Initialize component criteria.

        Parameters
        ----------
        momentum_eq : LinearMomentumBase
            Momentum equation solver.

        Returns
        -------
        None
        """
        self.force_criterion.initialize(momentum_eq)
        self.displacement_criterion.initialize(momentum_eq)
        self.reset_history()
        self.component_history = []

    def compute_error(self, momentum_eq: Any) -> float:
        """
        Compute combined normalized error from both component criteria.

        Parameters
        ----------
        momentum_eq : LinearMomentumBase
            Momentum equation solver.

        Returns
        -------
        float
            Combined normalized error.
        """
        force_error = self.force_criterion.compute_error(momentum_eq)
        displacement_error = self.displacement_criterion.compute_error(momentum_eq)

        force_ratio = force_error / max(self.force_criterion.tolerance, 1e-16)
        displacement_ratio = displacement_error / max(
            self.displacement_criterion.tolerance, 1e-16
        )
        combined_error = max(force_ratio, displacement_ratio)

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
            "force_tolerance": self.force_criterion.tolerance,
            "displacement_tolerance": self.displacement_criterion.tolerance,
            "combined_history": self.history,
            "component_history": self.component_history,
        }

