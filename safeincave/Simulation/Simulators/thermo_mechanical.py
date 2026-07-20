# Copyright (c) 2026, The SafeInCave Developers
#
# SPDX-License-Identifier: BSD-3-Clause

"""Thermo-mechanical simulator."""

from __future__ import annotations
from typing import Any
from .base import Simulator
from ...Utils import numpy2torch
from ...Output.SimLogging import SimulationLogging
from ...Output.Screen import ScreenPrinter
from ...Cavern import CavernHandler
from ..Convergence import ConvergenceErrorHandler

class Simulator_TM(Simulator):
    """
    Run the coupled thermo–mechanical simulation.

    Workflow
    --------
    1. Initialize outputs.
    2. Initialize momentum temperature from the heat solution and update BCs.
    3. Optionally solve a purely elastic response.
    4. Initialize non-elastic rates.
    5. For each time step:

       - Advance time and update boundary conditions for both equations.
       - Solve the heat equation for ``(t, dt)`` and set temperatures in momentum.
       - Iterate the momentum step (assemble/solve, update internal variables and rates).
       - Save requested fields.

    Returns
    -------
    None
    """

    def __init__(
        self,
        eq_mom: Any,
        eq_heat: Any,
        t_control: Any,
        outputs: Any,
        caverns: Any | None = CavernHandler(),
        compute_elastic_response: bool = True,
        merged_solutions: bool = False,
        smooth_output: bool = False,
        simulation_logger: Any | None = None,
        convergence_criterion: str = "strain_based",
    ):
        self.eq_mom = eq_mom
        self.eq_heat = eq_heat
        self.t_control = t_control
        self.outputs = outputs
        self.caverns = caverns
        self.compute_elastic_response = compute_elastic_response
        self.simulation_logger = simulation_logger
        self.convergence_handler = ConvergenceErrorHandler(convergence_criterion)

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
            self.eq_mom.grid,
            self.eq_mom.solver,
            self.eq_mom.mat,
            self.outputs,
            t_control.time_unit,
        )

    def run(self) -> None:
        """
        Run the coupled thermo–mechanical simulation.

        Workflow
        --------
        1. Initialize outputs.
        2. Initialize momentum temperature history from the heat solution.
        3. Update BCs and optionally solve a purely elastic step.
        4. Initialize non-elastic rates.
        5. Time loop:

           - Advance time, update BCs, solve heat step.
           - Fixed-point (or single-pass) iterate the momentum step with the
             current temperature, updating internal variables and rates.
           - Save requested fields.

        Convergence
        -----------
        Uses a relative change in total strain between iterations as error.
        If ``theta == 1.0`` (backward Euler) or there are no non-elastic
        elements, iteration terminates immediately.

        Returns
        -------
        None

        Notes
        -----
        - Calls ``SaveFields.initialize()`` once and ``save_fields(t)`` at each
          saved time, followed by ``save_mesh()`` after the loop.
        - Printing of progress occurs on rank 0 only.
        - The first ``output.save_fields(0)`` call targets the last ``output``
          from the preceding loop variable; ensure all outputs are saved by
          iterating over ``self.outputs``.
        """
        # Output field
        for output in self.outputs:
            output.initialize()

        # Update boundary conditions
        self.eq_mom.bc.update_dirichlet(self.t_control.t)
        self.eq_mom.bc.update_neumann(self.t_control.t)
        self.eq_mom.bc.update_cavern_bcs(self.caverns)

        if self.compute_elastic_response:
            # Solve elasticity
            self.eq_mom.solve_elastic_response()

            # Calculate total (elastic) strain
            eps_tot_to = self.eq_mom.compute_total_strain()

            # Compute stress
            stress_to = self.eq_mom.compute_elastic_stress(eps_tot_to)

        else:
            # Calculate total strain
            eps_tot_to = self.eq_mom.compute_total_strain()

            # Retrieve stress
            stress_to = numpy2torch(
                self.eq_mom.sig.x.array.reshape((self.eq_mom.n_elems, 3, 3))
            )

        # Calculate initial cavern volumes
        self.caverns.calculate_volumes(self.eq_mom.u)

        # # Calculate initial heat flux in caverns
        # self.caverns.calculate_total_heat(self.eq_heat.T, self.eq_heat.k)

        # Calculate initial cavern masses and volumes
        self.caverns.calculate_initial_conditions()

        # Record initial cavern data
        self.caverns.record_cavern_data(self.t_control.t)

        # Set initial temperature to momentum equation
        T_elems = self.eq_heat.get_T_elems()
        self.eq_mom.set_T(T_elems)
        self.eq_mom.set_T0(T_elems)

        # Calculate and eps_ie_rate_old
        self.eq_mom.compute_eps_ne_rate(stress_to, self.t_control.t)
        self.eq_mom.update_eps_ne_rate_old()

        # Save fields
        self.eq_mom.compute_p_elems()
        self.eq_mom.compute_q_elems()
        self.eq_mom.compute_principal_stresses()
        self.eq_mom.compute_principal_strains()

        # Save initial fields
        for output in self.outputs:
            output.save_fields(0)

        # Log initial state (step 0) to simulation logger
        if self.simulation_logger is not None:
            self.simulation_logger.log_initial_state(
                time=self.t_control.t,
                time_final=self.t_control.t_final,
                stress=stress_to,
                strain=eps_tot_to,
            )

        # Print initial state header row
        self.screen.print_initial_state(
            self.t_control.t,
            self.t_control.t_final,
            self.t_control.time_conversion,
        )

        # First step has no previous nonlinear iteration count.
        ite = 0

        # Time loop
        while self.t_control.keep_looping():
            max_bisections = int(getattr(self.t_control, "max_bisections", 5))
            n_bisections = 0
            step_completed = False

            while not step_completed:
                step_state = self._capture_step_state(include_heat=True, include_caverns=True)
                stress_to_step_start = stress_to.clone()

                # Advance time
                self.t_control.advance_time()

                t = self.t_control.t
                dt = self.t_control.dt

                # Update boundary conditions
                self.eq_mom.bc.update_dirichlet(t)
                self.eq_mom.bc.update_neumann(t)
                self.eq_heat.bc.update_bcs(t)

                # Initialize criterion at step start
                maxiter = int(getattr(self.t_control, "maxiter", 80) or 80)
                self.convergence_handler.initialize_step(
                    self.eq_mom,
                    maxiter=maxiter,
                )

                while self.convergence_handler.not_converged():
                    # Update cavern boundary conditions for heat diffusion equation
                    self.eq_heat.bc.update_cavern_bcs(self.caverns)

                    # Solve heat
                    self.eq_heat.solve(t, dt)

                    # Calculate total heat transfered through cavern walls
                    self.caverns.calculate_total_heat(dt, self.eq_heat.T, self.eq_heat.k)

                    # Update thermodynamic state of caverns
                    self.caverns.update_caverns(t, dt)

                    # Update cavern boundary conditions for momentum equation
                    self.eq_mom.bc.update_cavern_bcs(self.caverns)

                    # Set new temperature to momentum equation
                    T_elems = self.eq_heat.get_T_elems()
                    self.eq_mom.set_T(T_elems)

                    # Update stress
                    stress_k_to = stress_to.clone()

                    # Build bi-linear form
                    self.eq_mom.solve(stress_k_to, t, dt)

                    # Compute total strain
                    eps_tot_to = self.eq_mom.compute_total_strain()

                    # Compute stress
                    stress_to = self.eq_mom.compute_stress(eps_tot_to)

                    # Increment internal variables
                    self.eq_mom.increment_internal_variables(stress_to, stress_k_to, dt)

                    # Compute inelastic strain rates
                    self.eq_mom.compute_eps_ne_rate(stress_to, dt, eps_tot_to)

                    # Recalculate volumes of caverns
                    self.caverns.calculate_volumes(self.eq_mom.u)

                    # Compute error via active convergence criterion
                    self.convergence_handler.evaluate(self.eq_mom)

                    self.convergence_handler.increment_iteration()

                ite = self.convergence_handler.ite
                converged = not self.convergence_handler.not_converged_error
                retry_scale = 0.5

                if not converged:
                    if n_bisections >= max_bisections:
                        raise RuntimeError(
                            f"Failed to converge after {max_bisections} retries "
                            f"(ite={ite}, maxiter={self.convergence_handler.maxiter})."
                        )
                    self._restore_step_state(
                        step_state,
                        include_heat=True,
                        include_caverns=True,
                    )
                    stress_to = stress_to_step_start.clone()
                    dt_floor = float(getattr(self.t_control, "dt_min", 0.0))
                    self.t_control.dt = max(
                        step_state["time"]["dt"] * retry_scale, dt_floor
                    )
                    n_bisections += 1
                    continue

                # Closing solve: consume the final iteration's inelastic rate
                # so the committed internal state matches the converged rate.
                stress_to, stress_k_to, eps_tot_to = self._final_consuming_solve(
                    stress_to, t, dt
                )

                # Adaptive dt integration via relative convergence ratio.
                if hasattr(self.t_control, "get_next_dt"):
                    maxiter = max(int(self.convergence_handler.maxiter or 1), 1)
                    convergence_ratio = ite / maxiter
                    self.t_control.dt = self.t_control.get_next_dt(
                        convergence_ratio=convergence_ratio,
                        n_bisections=n_bisections,
                        converged=True,
                    )

                # Calculate next cavern masses and volumes
                self.caverns.calculate_initial_conditions()

                # Record thermodynamic data for caverns
                self.caverns.record_cavern_data(t)

                # Update internal variables
                self.eq_mom.update_internal_variables()

                # Update strain rates
                self.eq_mom.update_eps_ne_rate_old()

                # Update strain
                self.eq_mom.update_eps_ne_old(stress_to, stress_k_to, dt)

                # Update old temperature field
                self.eq_heat.update_T_old()

                # Save fields
                self.eq_mom.compute_p_elems()
                self.eq_mom.compute_q_elems()
                self.eq_mom.compute_principal_stresses()
                self.eq_mom.compute_principal_strains()
                for output in self.outputs:
                    output.save_fields(t)

                # Log simulation state
                if self.simulation_logger is not None:
                    self.simulation_logger.log_step(
                        step=self.t_control.step_counter,
                        stress=stress_to,
                        iteration=ite,
                        nonlinear_error=self.convergence_handler.error,
                        time=t,
                        time_step=dt,
                        time_final=self.t_control.t_final,
                        strain=eps_tot_to,
                    )

                # Print stuff
                current_time = "%.3f" % (t / self.t_control.time_conversion)
                screen_output_row = [
                    self.t_control.step_counter,
                    self.t_control.dt_used / self.t_control.time_conversion,
                    f"{current_time} / {self.t_control.t_final / self.t_control.time_conversion}",
                    ite,
                    self.convergence_handler.error,
                ]
                self.screen.print_row(screen_output_row)
                step_completed = True

        self.caverns.save_caverns_data()

        self.screen.close()

        self._finalize_outputs()
