# Copyright (c) 2026, The SafeInCave Developers
#
# SPDX-License-Identifier: BSD-3-Clause

"""Convergence criterion for the incremental Newton solver path."""

from __future__ import annotations
from typing import Any
from .base import ConvergenceCriterion


class NewtonResidualCriterion(ConvergenceCriterion):
    """
    Dual force/displacement criterion for the Newton path.

    Converged when BOTH hold:

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

    The Newton driver stores the two normalized errors on the momentum
    equation (``momentum_eq._newton_state``) after each residual assembly;
    this criterion only reads them. The reported error is
    ``max(force_err/force_tol, disp_err/disp_tol)``, converged at ≤ 1. The
    plastic-consistency gate of ConvergenceErrorHandler applies on top,
    unchanged.
    """

    def __init__(
        self,
        force_tolerance: float = 1e-6,
        displacement_tolerance: float = 1e-2,
        name: str | None = None,
    ):
        super().__init__(tolerance=1.0, name=name or "NewtonResidual")
        self.force_tolerance = float(force_tolerance)
        self.displacement_tolerance = float(displacement_tolerance)

    def initialize(self, momentum_eq: Any) -> None:
        momentum_eq._newton_state = None

    def compute_error(self, momentum_eq: Any) -> float:
        state = getattr(momentum_eq, "_newton_state", None)
        if state is None:
            return float("inf")
        force_error = float(state.get("force_error", float("inf")))
        # No correction solved yet: the force residual alone decides (a step
        # with an unchanged load may converge with zero solves).
        disp_error = float(state.get("disp_error", 0.0))
        return max(
            force_error / self.force_tolerance,
            disp_error / self.displacement_tolerance,
        )

    def is_converged(self, error: float) -> bool:
        return error <= 1.0
