# Copyright (c) 2026, The SafeInCave Developers
#
# SPDX-License-Identifier: BSD-3-Clause

"""Mechanical simulator output variant (deprecated)."""

from __future__ import annotations
from typing import Any
import warnings
from .base import Simulator
from .mechanical import Simulator_M
from ...Equations.Momentum import LinearMomentum

class Simulator_Mout(Simulator):
    """
    Mechanical-only simulator (linear momentum).

    .. deprecated:: 3.0.3
        Use :class:`Simulator_M` instead. This class is maintained for
        backward compatibility but will be removed in a future version.

    Solves the momentum equation with possible non-elastic behavior using a
    θ-method loop and fixed-point iterations per step. No thermal coupling.

    Parameters
    ----------
    eq_mom : LinearMomentum
        Configured momentum equation (materials, BCs, solver set).
    t_control : TimeControllerBase
        Time controller providing `t`, `dt`, and loop control.
    outputs : list of SaveFields
        Output writers to initialize and use at each saved time.
    compute_elastic_response : bool, default=True
        If True, starts with a purely elastic solve to initialize fields.

    Attributes
    ----------
    eq_mom : LinearMomentum
    t_control : TimeControllerBase
    outputs : list[SaveFields]
    compute_elastic_response : bool
    """

    def __init__(
        self,
        eq_mom: LinearMomentum,
        t_control: Any,
        outputs: Any,
        compute_elastic_response: bool = True,
        convergence_criterion: str = "strain_based",
        merged_solutions: bool = False,
        smooth_output: bool = False,
        simulation_logger: Any | None = None,
        plastic_consistency_tolerance: float = 1e-4,
    ):
        warnings.warn(
            "Simulator_Mout is deprecated and will be removed in a future version. "
            "Use Simulator_M instead, which supports the same functionality.",
            DeprecationWarning,
            stacklevel=2
        )
        
        # Delegate to Simulator_M with caverns=None (no cavern support in Mout)
        self._simulator_m = Simulator_M(
            eq_mom=eq_mom,
            t_control=t_control,
            outputs=outputs,
            caverns=None,  # Simulator_Mout doesn't support caverns
            compute_elastic_response=compute_elastic_response,
            convergence_criterion=convergence_criterion,
            maxiter=40,  # Default from Simulator_M
            simulation_logger=simulation_logger,
            merged_solutions=merged_solutions,
            smooth_output=smooth_output,
            plastic_consistency_tolerance=plastic_consistency_tolerance,
        )

    def run(self) -> None:
        """
        Run the mechanical simulation.

        .. deprecated:: 3.0.3
            Use :meth:`Simulator_M.run` instead.
        """
        # Delegate to the wrapped Simulator_M instance
        self._simulator_m.run()
