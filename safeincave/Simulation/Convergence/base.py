# Copyright (c) 2026, The SafeInCave Developers
#
# SPDX-License-Identifier: BSD-3-Clause

"""Abstract base class for convergence criteria."""

from typing import TYPE_CHECKING, Optional, Any
from abc import ABC, abstractmethod

if TYPE_CHECKING:
    pass


def _is_momentum_solver_instance(obj: Any) -> bool:
    """Duck-type check: is this a momentum solver instance?"""
    return hasattr(obj, 'compute_total_strain') and hasattr(obj, 'u') and hasattr(obj, 'mat')


class ConvergenceCriterion(ABC):
    """
    Abstract base class for convergence criteria.

    Each criterion computes a single error metric that can be checked against
    a tolerance to create flexible convergence strategies.

    Parameters
    ----------
    tolerance : float, optional
        Convergence tolerance for this criterion. Default: 1e-7.
    name : str, optional
        Descriptive name for diagnostics/logging. Default: class name.

    Attributes
    ----------
    tolerance : float
        Convergence tolerance.
    name : str
        Criterion name.
    history : list
        Record of error values per iteration.
    """

    def __init__(self, tolerance: float = 1e-7, name: Optional[str] = None):
        self.tolerance = tolerance
        self.name = name or self.__class__.__name__
        self.history: list = []

    @abstractmethod
    def initialize(self, momentum_eq: Any) -> None:
        """
        Initialize criterion at step start.

        Called once per time step, before entering the nonlinear iteration loop.

        Parameters
        ----------
        momentum_eq : LinearMomentumBase
            Momentum equation solver instance.

        Returns
        -------
        None
        """
        pass

    @abstractmethod
    def compute_error(self, *args, **kwargs) -> float:
        """
        Compute error metric for current iteration.

        Returns
        -------
        float
            Dimensionless error value. Returns 0 if criterion not applicable
            (e.g., elastic-only problem for strain criterion).
        """
        pass

    @abstractmethod
    def is_converged(self, error: float) -> bool:
        """
        Check if error is below tolerance.

        Parameters
        ----------
        error : float
            Error from compute_error().

        Returns
        -------
        bool
            True if error ≤ tolerance.
        """
        pass

    def get_history(self) -> list:
        """Return iteration history."""
        return self.history

    def reset_history(self) -> None:
        """Reset iteration history."""
        self.history = []

    def __repr__(self) -> str:
        return f"{self.name}(tol={self.tolerance:.2e})"
