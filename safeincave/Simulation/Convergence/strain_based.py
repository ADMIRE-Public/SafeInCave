# Copyright (c) 2026, The SafeInCave Developers
#
# SPDX-License-Identifier: BSD-3-Clause

"""Strain-based convergence criterion."""

from typing import TYPE_CHECKING, Optional, Any
import torch as to
import numpy as np
from mpi4py import MPI
from .base import ConvergenceCriterion, _is_momentum_solver_instance

if TYPE_CHECKING:
    pass

class StrainBasedCriterion(ConvergenceCriterion):
    """
    Strain-based convergence criterion.

    Error metric:
    $$\\text{error} = \\frac{\\|\\varepsilon_\\text{new} - \\varepsilon_\\text{old}\\|}{\\|\\varepsilon_\\text{new}\\|}$$

    This implementation preserves prior behavior,
    including MPI SUM reduction and momentum-equation-tracked previous strain.

    Parameters
    ----------
    tolerance : float, optional
        Strain error tolerance. Default: 1e-7.
    """

    def __init__(
        self,
        tolerance: float = 1e-5,
        name: Optional[str] = None,
    ):
        super().__init__(tolerance=tolerance, name=name or "strain_based")

    def initialize(self, momentum_eq: Any) -> None:
        """Initialize tracked strain state at step start."""
        momentum_eq._strain_previous = momentum_eq.compute_total_strain().clone()
        self.reset_history()

    def compute_error(self, *args, **kwargs) -> float:
        """
        Compute strain-based error: ||Δε|| / ||ε_current||.

        Parameters
        ----------
        Supports call styles:
        1) compute_error(momentum_eq)
        2) compute_error(momentum_eq, strain_previous, strain_current)

        Returns
        -------
        float
            Relative strain change error. Returns 0 if not applicable.
        """
        if len(args) == 1 and _is_momentum_solver_instance(args[0]):
            momentum_eq = args[0]
            strain_current = momentum_eq.compute_total_strain()
            if hasattr(momentum_eq, "_strain_previous"):
                strain_previous = momentum_eq._strain_previous
            else:
                strain_previous = strain_current
        elif len(args) == 3 and _is_momentum_solver_instance(args[0]):
            momentum_eq = args[0]
            strain_previous = args[1]
            strain_current = args[2]
        else:
            raise ValueError(
                "StrainBasedCriterion.compute_error expects (momentum_eq) or "
                "(momentum_eq, strain_previous, strain_current)."
            )

        # Early exit for explicit time integration
        if momentum_eq.theta == 1.0:
            error_value = 0.0
            self.history.append(error_value)
            return error_value

        # Early exit for purely elastic material
        if len(momentum_eq.mat.elems_ne) == 0:
            error_value = 0.0
            self.history.append(error_value)
            return error_value

        strain_previous_flat = to.flatten(strain_previous)
        strain_current_flat = to.flatten(strain_current)

        denom = max(np.linalg.norm(strain_current_flat), 1e-16)
        local_error = np.linalg.norm(strain_previous_flat - strain_current_flat) / denom
        error = momentum_eq.grid.mesh.comm.allreduce(local_error, op=MPI.SUM)

        # Update tracked state for next call.
        momentum_eq._strain_previous = strain_current.clone()

        error_value = float(error)

        self.history.append(error_value)
        return float(error_value)

    def is_converged(self, error: float) -> bool:
        """Check if strain error is below tolerance."""
        return error <= self.tolerance
