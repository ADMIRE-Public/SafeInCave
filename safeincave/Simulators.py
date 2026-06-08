# Copyright (c) 2026, The SafeInCave Developers
#
# SPDX-License-Identifier: BSD-3-Clause

from abc import ABC, abstractmethod
import torch as to
import numpy as np
from mpi4py import MPI
from .Utils import numpy2torch
from .HeatEquation import HeatDiffusion
from .MomentumEquation import LinearMomentum, LinearMomentumBase
from .TimeHandler import TimeControllerBase
from .OutputHandler import SaveFields
from .ScreenOutput import ScreenPrinter
from .CavernBC import CavernHandler


class Simulator(ABC):
    """
    Abstract simulation driver interface.

    Subclasses implement a concrete `run()` method that advances one or more
    coupled PDE solvers in time, handles I/O, and updates material/internal
    variables as needed.
    """

    @abstractmethod
    def run(self):
        """
        Execute the simulation.

        Returns
        -------
        None
        """
        pass


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
        eq_mom: LinearMomentumBase,
        eq_heat: HeatDiffusion,
        t_control: TimeControllerBase,
        outputs: list[SaveFields],
        caverns: CavernHandler | None = CavernHandler(),
        compute_elastic_response: bool = True,
    ):
        self.eq_mom = eq_mom
        self.eq_heat = eq_heat
        self.t_control = t_control
        self.outputs = outputs
        self.caverns = caverns
        self.compute_elastic_response = compute_elastic_response

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
        self.eq_mom.compute_p_nodes()
        self.eq_mom.compute_q_nodes()

        # Save initial fields
        for output in self.outputs:
            output.save_fields(0)

        # Time loop
        while self.t_control.keep_looping():
            # Advance time
            self.t_control.advance_time()
            t = self.t_control.t
            dt = self.t_control.dt

            # Update boundary conditions
            self.eq_mom.bc.update_dirichlet(t)
            self.eq_mom.bc.update_neumann(t)
            self.eq_heat.bc.update_bcs(t)

            # Iterative loop settings
            tol = 1e-7
            error = 2 * tol
            ite = 0
            maxiter = 80

            while error > tol and ite < maxiter:
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

                # Update total strain of previous iteration (eps_tot_k <-- eps_tot)
                eps_tot_k_to = eps_tot_to.clone()

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
                self.eq_mom.compute_eps_ne_rate(stress_to, dt)

                # Recalculate volumes of caverns
                self.caverns.calculate_volumes(self.eq_mom.u)

                # Compute error
                if self.eq_mom.theta == 1.0:
                    error = 0.0
                elif len(self.eq_mom.mat.elems_ne) == 0:
                    error = 0.0
                else:
                    eps_tot_k_flat = to.flatten(eps_tot_k_to)
                    eps_tot_flat = to.flatten(eps_tot_to)
                    local_error = np.linalg.norm(
                        eps_tot_k_flat - eps_tot_flat
                    ) / np.linalg.norm(eps_tot_flat)
                    error = self.eq_mom.grid.mesh.comm.allreduce(
                        local_error, op=MPI.SUM
                    )

                ite += 1

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
            self.eq_mom.compute_p_nodes()
            self.eq_mom.compute_q_nodes()
            for output in self.outputs:
                output.save_fields(t)

            # Print stuff
            current_time = "%.3f" % (t / self.t_control.time_conversion)
            screen_output_row = [
                self.t_control.step_counter,
                self.t_control.dt / self.t_control.time_conversion,
                f"{current_time} / {self.t_control.t_final / self.t_control.time_conversion}",
                ite,
                error,
            ]
            self.screen.print_row(screen_output_row)

        self.caverns.save_caverns_data()

        self.screen.close()

        for output in self.outputs:
            output.save_mesh()


class Simulator_M(Simulator):
    """
    Mechanical-only simulator (linear momentum).

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
        eq_mom: LinearMomentumBase,
        t_control: TimeControllerBase,
        outputs: list[SaveFields],
        caverns: CavernHandler | None = CavernHandler(),
        compute_elastic_response: bool = True,
    ):
        self.eq_mom = eq_mom
        self.t_control = t_control
        self.outputs = outputs
        self.caverns = caverns
        self.compute_elastic_response = compute_elastic_response

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
        Run the mechanical simulation.

        Workflow
        --------
        1. Initialize outputs and boundary conditions.
        2. Optionally solve a purely elastic step.
        3. Initialize non-elastic rates.
        4. For each time step: assemble/solve, update internal variables and
           rates, compute relevant quantities, and save fields.

        Convergence
        -----------
        Uses a relative change in total strain between iterations as error.
        If `theta == 1.0` or no non-elastic elements exist, iteration ends
        immediately.

        Returns
        -------
        None

        Notes
        -----
        - Printing occurs on rank 0 only.
        - The first `output.save_fields(0)` call uses the last `output`
          from the preceding loop variable.
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
            eps_tot_to = numpy2torch(
                self.eq_mom.eps_tot.x.array.reshape((self.eq_mom.n_elems, 3, 3))
            )

            # Retrieve stress
            stress_to = numpy2torch(
                self.eq_mom.sig.x.array.reshape((self.eq_mom.n_elems, 3, 3))
            )

        # Calculate initial cavern volumes
        self.caverns.calculate_volumes(self.eq_mom.u)

        # Calculate initial cavern masses and volumes
        self.caverns.calculate_initial_conditions()

        # Record initial cavern data
        self.caverns.record_cavern_data(self.t_control.t)

        # Calculate and eps_ie_rate_old
        self.eq_mom.compute_eps_ne_rate(stress_to, self.t_control.t)
        self.eq_mom.update_eps_ne_rate_old()

        # Save fields
        self.eq_mom.compute_p_elems()
        self.eq_mom.compute_q_elems()
        self.eq_mom.compute_p_nodes()
        self.eq_mom.compute_q_nodes()

        # Save initial fields
        for output in self.outputs:
            output.save_fields(0)

        # Time loop
        while self.t_control.keep_looping():
            # Advance time
            self.t_control.advance_time()
            t = self.t_control.t
            dt = self.t_control.dt

            # Update boundary conditions
            self.eq_mom.bc.update_dirichlet(t)
            self.eq_mom.bc.update_neumann(t)

            # Iterative loop settings
            tol = 1e-8
            error = 2 * tol
            ite = 0
            maxiter = 40

            while error > tol and ite < maxiter:
                # Update thermodynamic state of caverns
                self.caverns.update_caverns(t, dt)

                # Update cavern boundary conditions
                self.eq_mom.bc.update_cavern_bcs(self.caverns)

                # Update total strain of previous iteration (eps_tot_k <-- eps_tot)
                eps_tot_k_to = eps_tot_to.clone()

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
                self.eq_mom.compute_eps_ne_rate(stress_to, dt)

                # Recalculate volumes of caverns
                self.caverns.calculate_volumes(self.eq_mom.u)

                # Compute error
                if self.eq_mom.theta == 1.0:
                    error = 0.0
                elif len(self.eq_mom.mat.elems_ne) == 0:
                    error = 0.0
                else:
                    eps_tot_k_flat = to.flatten(eps_tot_k_to)
                    eps_tot_flat = to.flatten(eps_tot_to)
                    local_error = np.linalg.norm(
                        eps_tot_k_flat - eps_tot_flat
                    ) / np.linalg.norm(eps_tot_flat)
                    error = self.eq_mom.grid.mesh.comm.allreduce(
                        local_error, op=MPI.SUM
                    )

                ite += 1

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

            # Save fields
            self.eq_mom.compute_p_nodes()
            self.eq_mom.compute_p_elems()
            self.eq_mom.compute_q_elems()
            self.eq_mom.compute_q_nodes()
            for output in self.outputs:
                output.save_fields(t)

            # Print stuff
            current_time = "%.3f" % (t / self.t_control.time_conversion)
            screen_output_row = [
                self.t_control.step_counter,
                self.t_control.dt / self.t_control.time_conversion,
                f"{current_time} / {self.t_control.t_final / self.t_control.time_conversion}",
                ite,
                error,
            ]
            self.screen.print_row(screen_output_row)

        self.caverns.save_caverns_data()

        self.screen.close()

        for output in self.outputs:
            output.save_mesh()


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
        eq_heat: HeatDiffusion,
        t_control: TimeControllerBase,
        outputs: list[SaveFields],
        caverns: CavernHandler | None = None,
    ):
        self.eq_heat = eq_heat
        self.t_control = t_control
        self.outputs = outputs
        self.caverns = caverns

        if caverns is None:
            self.caverns = CavernHandler()

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

            # Print stuff
            current_time = "%.3f" % (t / self.t_control.time_conversion)
            screen_output_row = [
                self.t_control.step_counter,
                self.t_control.dt / self.t_control.time_conversion,
                f"{current_time} / {self.t_control.t_final / self.t_control.time_conversion}",
                0,
                0,
            ]
            self.screen.print_row(screen_output_row)

        self.screen.close()

        for output in self.outputs:
            output.save_mesh()


class Simulator_Mout(Simulator):
    """
    Mechanical-only simulator (linear momentum).

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
        t_control: TimeControllerBase,
        outputs: list[SaveFields],
        compute_elastic_response: bool = True,
    ):
        self.eq_mom = eq_mom
        self.t_control = t_control
        self.outputs = outputs
        self.compute_elastic_response = compute_elastic_response

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
        Run the mechanical simulation.

        Workflow
        --------
        1. Initialize outputs and boundary conditions.
        2. Optionally solve a purely elastic step.
        3. Initialize non-elastic rates.
        4. For each time step: assemble/solve, update internal variables and
           rates, compute relevant quantities, and save fields.

        Convergence
        -----------
        Uses a relative change in total strain between iterations as error.
        If `theta == 1.0` or no non-elastic elements exist, iteration ends
        immediately.

        Returns
        -------
        None

        Notes
        -----
        - Printing occurs on rank 0 only.
        - The first `output.save_fields(0)` call uses the last `output`
          from the preceding loop variable.
        """
        # Output field
        for output in self.outputs:
            output.initialize()

        # Update boundary conditions
        self.eq_mom.bc.update_dirichlet(self.t_control.t)
        self.eq_mom.bc.update_neumann(self.t_control.t)

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

        # Calculate and eps_ie_rate_old
        self.eq_mom.compute_eps_ne_rate(stress_to, self.t_control.t)
        self.eq_mom.update_eps_ne_rate_old()

        # Save fields
        self.eq_mom.compute_p_elems()
        self.eq_mom.compute_q_elems()
        self.eq_mom.compute_p_nodes()
        self.eq_mom.compute_q_nodes()
        output.save_fields(0)

        # Time loop
        while self.t_control.keep_looping():
            # Advance time
            self.t_control.advance_time()
            t = self.t_control.t
            dt = self.t_control.dt

            # Update boundary conditions
            self.eq_mom.bc.update_dirichlet(t)
            self.eq_mom.bc.update_neumann(t)

            # Iterative loop settings
            tol = 1e-8
            error = 2 * tol
            ite = 0
            maxiter = 40

            while error > tol and ite < maxiter:
                # Update total strain of previous iteration (eps_tot_k <-- eps_tot)
                eps_tot_k_to = eps_tot_to.clone()

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
                self.eq_mom.compute_eps_ne_rate(stress_to, dt)

                # Compute error
                if self.eq_mom.theta == 1.0:
                    error = 0.0
                elif len(self.eq_mom.mat.elems_ne) == 0:
                    error = 0.0
                else:
                    eps_tot_k_flat = to.flatten(eps_tot_k_to)
                    eps_tot_flat = to.flatten(eps_tot_to)
                    local_error = np.linalg.norm(
                        eps_tot_k_flat - eps_tot_flat
                    ) / np.linalg.norm(eps_tot_flat)
                    error = self.eq_mom.grid.mesh.comm.allreduce(
                        local_error, op=MPI.SUM
                    )

                ite += 1

            # Update internal variables
            self.eq_mom.update_internal_variables()

            # Update strain rates
            self.eq_mom.update_eps_ne_rate_old()

            # Update strain
            self.eq_mom.update_eps_ne_old(stress_to, stress_k_to, dt)

            # Save fields
            self.eq_mom.compute_p_elems()
            self.eq_mom.compute_q_elems()
            self.eq_mom.compute_p_nodes()
            self.eq_mom.compute_q_nodes()
            for output in self.outputs:
                output.save_fields(t)

            # Print stuff
            screen_output_row = [
                self.t_control.step_counter,
                self.t_control.dt / self.t_control.time_conversion,
                f"{t / self.t_control.time_conversion} / {self.t_control.t_final / self.t_control.time_conversion}",
                ite,
                error,
            ]
            self.screen.print_row(screen_output_row)

            # if self.eq_mom.grid.mesh.comm.rank == 0:
            # 	print(t/self.t_control.time_unit, ite, error)
            # 	sys.stdout.flush()
            # 	try:
            # 		print(float(self.eq_mom.mat.elems_ne[-1].Fvp.max()))
            # 		sys.stdout.flush()
            # 	except:
            # 		pass

        self.screen.close()
        for output in self.outputs:
            output.save_mesh()
