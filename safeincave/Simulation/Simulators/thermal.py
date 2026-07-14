# Copyright (c) 2026, The SafeInCave Developers
#
# SPDX-License-Identifier: BSD-3-Clause

"""Thermal simulator."""

from __future__ import annotations
from typing import Any
from .base import Simulator
from ...Output.SimLogging import SimulationLogging
from ...Output.Screen import ScreenPrinter
from ...Cavern import CavernHandler

class Simulator_T(Simulator):
    """
    Thermal-only simulator (heat diffusion).

    Advances the heat equation with fully-implicit time loop and writes fields.

    Parameters
    ----------
    eq_heat : HeatDiffusion
        Configured heat equation (materials, BCs, solver set).
    t_control : TimeControllerBase
        Time controller providing `t`, `dt`, and loop control.
    outputs : list of SaveFields
        Output writers to initialize and use at each saved time.
    compute_elastic_response : bool, default=True
        Unused placeholder kept for interface parity.

    Attributes
    ----------
    eq_heat : HeatDiffusion
    t_control : TimeControllerBase
    outputs : list[SaveFields]
    """

    def __init__(
        self,
        eq_heat: Any,
        t_control: Any,
        outputs: Any,
        caverns: Any | None = None,
        merged_solutions: bool = False,
        smooth_output: bool = False,
        simulation_logger: Any | None = None,
    ):
        self.eq_heat = eq_heat
        self.t_control = t_control
        self.outputs = outputs
        self.caverns = caverns
        self.simulation_logger = simulation_logger

        if caverns is None:
            self.caverns = CavernHandler()

        # Apply merged_solutions and smooth_output flags to all output handlers
        for output in self.outputs:
            output.merged_solutions = merged_solutions
            output.smooth_output = smooth_output

        # Auto-configure simulation logger if provided
        if self.simulation_logger is not None:
            SimulationLogging.auto_configure_from_simulator(
                self.simulation_logger, self, self.outputs
            )

        ScreenPrinter.reset_instance()
        self.screen = ScreenPrinter(
            self.eq_heat.grid,
            self.eq_heat.solver,
            self.eq_heat.mat,
            self.outputs,
            t_control.time_unit,
        )

    def run(self) -> None:
        """
        Run the thermal simulation.

        Workflow
        --------
        1. Initialize outputs.
        2. (Optionally) solve an initial step.
        3. Time loop: update BCs, solve heat equation for `(t, dt)`, and save.

        Returns
        -------
        None

        Notes
        -----
        Printing of progress occurs on rank 0 only.
        """
        # Output field
        for output in self.outputs:
            output.initialize()

        # # Solve initial T field
        # self.eq_heat.solve(0, self.t_control.dt)

        # Save fields
        output.save_fields(0)

        # Log initial state (step 0) to simulation logger
        if self.simulation_logger is not None:
            self.simulation_logger.log_initial_state(
                time=self.t_control.t,
                time_final=self.t_control.t_final,
            )

        # Print initial state header row
        self.screen.print_initial_state(
            self.t_control.t,
            self.t_control.t_final,
            self.t_control.time_conversion,
        )

        # Time loop
        while self.t_control.keep_looping():

            # Advance time
            self.t_control.advance_time()

            t = self.t_control.t
            dt = self.t_control.dt

            # Update boundary conditions
            self.eq_heat.bc.update_bcs(t)
            self.eq_heat.bc.update_cavern_bcs(self.caverns)

            # Solve heat
            self.eq_heat.solve(t, dt)

            # Save fields
            for output in self.outputs:
                output.save_fields(t)

            # Log simulation state
            if self.simulation_logger is not None:
                self.simulation_logger.log_step(
                    step=self.t_control.step_counter,
                    iteration=0,
                    nonlinear_error=0.0,
                    time=t,
                    time_step=dt,
                    time_final=self.t_control.t_final,
                )

            # Print stuff
            current_time = "%.3f" % (t / self.t_control.time_conversion)
            screen_output_row = [
                self.t_control.step_counter,
                self.t_control.dt_used / self.t_control.time_conversion,
                f"{current_time} / {self.t_control.t_final / self.t_control.time_conversion}",
                0,
                0,
            ]
            self.screen.print_row(screen_output_row)

        self.screen.close()

        self._finalize_outputs()
