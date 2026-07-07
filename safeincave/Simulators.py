# Copyright (c) 2026, The SafeInCave Developers
#
# SPDX-License-Identifier: BSD-3-Clause

from abc import ABC, abstractmethod
import copy
from .Utils import numpy2torch
from .HeatEquation import HeatDiffusion
from .MomentumEquation import LinearMomentum, LinearMomentumBase
from .TimeHandler import TimeControllerBase
from .OutputHandler import SaveFields
from .SimulationLogging import SimulationLogging
from .ScreenOutput import ScreenPrinter
from .CavernBC import CavernHandler
from .ConvergenceCriteria import ConvergenceErrorHandler


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

    def _compute_error(self) -> float:
        """
        Compute current raw convergence error.

        Criterion is responsible for fetching required state from momentum_eq.
        Enables any error metric (strain, residual, displacement, composite)
        without loop coupling to specific implementations.

        Returns
        -------
        float
            Raw criterion error value for convergence checks and reporting.
        """
        return self.convergence_handler.compute_error(self.eq_mom)

    def _initialize_convergence_criterion(self) -> None:
        """
        Initialize active convergence strategy at time-step start.

        Parameters
        ----------
        None
        """
        self.convergence_handler.initialize_step(self.eq_mom)

    def _finalize_outputs(self) -> None:
        """
        Close and finalize all output handlers.

        Closes all XDMF file handles and copies mesh files to the output
        directory for provenance. This should be called at the end of a
        simulation's `run()` method.

        Returns
        -------
        None
        """
        for output in self.outputs:
            output.close()
            output.save_mesh()

    @staticmethod
    def _is_array_like(value) -> bool:
        """Return True for array-like values that expose copy/shape/dtype."""
        return (
            hasattr(value, "copy")
            and hasattr(value, "shape")
            and hasattr(value, "dtype")
        )

    @staticmethod
    def _clone_value(value):
        """Best-effort cloning utility for tensors/arrays/containers/scalars."""
        if hasattr(value, "clone"):
            try:
                return value.clone()
            except Exception:
                pass
        if hasattr(value, "copy"):
            try:
                return value.copy()
            except Exception:
                pass
        try:
            return copy.deepcopy(value)
        except Exception:
            return value

    def _snapshot_function_arrays(self, obj) -> dict:
        """Capture all dolfinx-like Function arrays found in object attributes."""
        snapshot = {}
        for attr_name, attr_value in getattr(obj, "__dict__", {}).items():
            if hasattr(attr_value, "x") and hasattr(attr_value.x, "array"):
                try:
                    snapshot[attr_name] = attr_value.x.array.copy()
                except Exception:
                    continue
        return snapshot

    def _restore_function_arrays(self, obj, snapshot: dict) -> None:
        """Restore dolfinx-like Function arrays captured by _snapshot_function_arrays."""
        for attr_name, saved_array in snapshot.items():
            attr_value = getattr(obj, attr_name, None)
            if hasattr(attr_value, "x") and hasattr(attr_value.x, "array"):
                attr_value.x.array[:] = saved_array

    def _snapshot_object_state(self, obj) -> dict:
        """Capture mutable object attributes using best-effort cloning."""
        snapshot = {}
        for attr_name, attr_value in getattr(obj, "__dict__", {}).items():
            if callable(attr_value):
                continue
            if hasattr(attr_value, "x") and hasattr(attr_value.x, "array"):
                continue
            if isinstance(attr_value, (int, float, bool, str, type(None), tuple, list, dict)):
                snapshot[attr_name] = self._clone_value(attr_value)
            elif hasattr(attr_value, "clone") or self._is_array_like(attr_value):
                snapshot[attr_name] = self._clone_value(attr_value)
        return snapshot

    def _restore_object_state(self, obj, snapshot: dict) -> None:
        """Restore object attributes captured by _snapshot_object_state."""
        for attr_name, saved_value in snapshot.items():
            setattr(obj, attr_name, self._clone_value(saved_value))

    def _snapshot_material_internal_state(self, eq_mom) -> list:
        """Capture non-elastic internal variable state from material elements."""
        material = getattr(eq_mom, "mat", None)
        if material is None or not hasattr(material, "elems_ne"):
            return []
        return [self._snapshot_object_state(elem) for elem in material.elems_ne]

    def _restore_material_internal_state(self, eq_mom, snapshot: list) -> None:
        """Restore non-elastic internal variable state for material elements."""
        material = getattr(eq_mom, "mat", None)
        if material is None or not hasattr(material, "elems_ne"):
            return
        for elem, elem_snapshot in zip(material.elems_ne, snapshot):
            self._restore_object_state(elem, elem_snapshot)

    def _snapshot_caverns_state(self) -> dict:
        """Capture cavern model mutable states."""
        caverns = getattr(self, "caverns", None)
        if caverns is None:
            return {}
        snapshot = {}
        for cavern_list_name in ("caverns_T", "caverns_PT", "caverns_MFlux"):
            cavern_list = getattr(caverns, cavern_list_name, [])
            snapshot[cavern_list_name] = [
                self._snapshot_object_state(cavern) for cavern in cavern_list
            ]
        return snapshot

    def _restore_caverns_state(self, snapshot: dict) -> None:
        """Restore cavern model mutable states."""
        caverns = getattr(self, "caverns", None)
        if caverns is None:
            return
        for cavern_list_name in ("caverns_T", "caverns_PT", "caverns_MFlux"):
            cavern_list = getattr(caverns, cavern_list_name, [])
            saved_list = snapshot.get(cavern_list_name, [])
            for cavern, cavern_snapshot in zip(cavern_list, saved_list):
                self._restore_object_state(cavern, cavern_snapshot)

    def _capture_step_state(self, include_heat: bool = False, include_caverns: bool = False) -> dict:
        """Capture state needed to rollback a failed nonlinear step attempt."""
        state = {
            "time": {
                "t": self.t_control.t,
                "dt": self.t_control.dt,
                "step_counter": self.t_control.step_counter,
            },
            "eq_mom_functions": self._snapshot_function_arrays(self.eq_mom),
            "material_state": self._snapshot_material_internal_state(self.eq_mom),
        }
        if include_heat and hasattr(self, "eq_heat"):
            state["eq_heat_functions"] = self._snapshot_function_arrays(self.eq_heat)
        if include_caverns:
            state["caverns_state"] = self._snapshot_caverns_state()
        return state

    def _restore_step_state(self, state: dict, include_heat: bool = False, include_caverns: bool = False) -> None:
        """Restore state captured by _capture_step_state."""
        time_state = state["time"]
        self.t_control.t = time_state["t"]
        self.t_control.dt = time_state["dt"]
        self.t_control.step_counter = time_state["step_counter"]

        self._restore_function_arrays(self.eq_mom, state["eq_mom_functions"])
        self._restore_material_internal_state(self.eq_mom, state["material_state"])

        if include_heat and hasattr(self, "eq_heat") and "eq_heat_functions" in state:
            self._restore_function_arrays(self.eq_heat, state["eq_heat_functions"])
        if include_caverns and "caverns_state" in state:
            self._restore_caverns_state(state["caverns_state"])


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
        merged_solutions: bool = False,
        simulation_logger: SimulationLogging | None = None,
    ):
        self.eq_mom = eq_mom
        self.eq_heat = eq_heat
        self.t_control = t_control
        self.outputs = outputs
        self.caverns = caverns
        self.compute_elastic_response = compute_elastic_response
        self.simulation_logger = simulation_logger

        # Apply merged_solutions flag to all output handlers
        for output in self.outputs:
            output.merged_solutions = merged_solutions

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
        self.eq_mom.compute_p_nodes()
        self.eq_mom.compute_q_nodes()

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
                    self.eq_mom.compute_eps_ne_rate(stress_to, dt)

                    # Recalculate volumes of caverns
                    self.caverns.calculate_volumes(self.eq_mom.u)

                    # Compute error via active convergence criterion
                    self.convergence_handler.evaluate(self.eq_mom)

                    self.convergence_handler.increment_iteration()

                ite = self.convergence_handler.ite
                converged = not self.convergence_handler.not_converged_error

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
                    self.t_control.dt = max(step_state["time"]["dt"] * 0.5, dt_floor)
                    n_bisections += 1
                    continue

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
                self.eq_mom.compute_p_nodes()
                self.eq_mom.compute_q_nodes()
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
        convergence_criterion: str = "strain_based",
        maxiter: int = 40,
        simulation_logger: SimulationLogging | None = None,
        merged_solutions: bool = False,
    ):
        self.eq_mom = eq_mom
        self.t_control = t_control
        self.outputs = outputs
        self.caverns = caverns
        self.compute_elastic_response = compute_elastic_response
        self.convergence_handler = ConvergenceErrorHandler(convergence_criterion)
        self.maxiter = int(maxiter)
        self.simulation_logger = simulation_logger

        # Apply merged_solutions flag to all output handlers
        for output in self.outputs:
            output.merged_solutions = merged_solutions

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

        # Log initial state (step 0) to simulation logger
        if self.simulation_logger is not None:
            log_kwargs = SimulationLogging.extract_yield_variables(self.eq_mom.mat)
            self.simulation_logger.log_initial_state(
                time=self.t_control.t,
                time_final=self.t_control.t_final,
                stress=stress_to,
                strain=eps_tot_to,
                **log_kwargs
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
                step_state = self._capture_step_state(include_heat=False, include_caverns=True)
                stress_to_step_start = stress_to.clone()

                # Advance time
                self.t_control.advance_time()

                t = self.t_control.t
                dt = self.t_control.dt

                # Update boundary conditions
                self.eq_mom.bc.update_dirichlet(t)
                self.eq_mom.bc.update_neumann(t)

                # Initialize criterion at step start
                maxiter = int(
                    getattr(self.t_control, "maxiter", self.maxiter)
                    or self.maxiter
                )
                self.convergence_handler.initialize_step(
                    self.eq_mom,
                    maxiter=maxiter,
                )

                while self.convergence_handler.not_converged():
                    # Update thermodynamic state of caverns
                    self.caverns.update_caverns(t, dt)

                    # Update cavern boundary conditions
                    self.eq_mom.bc.update_cavern_bcs(self.caverns)

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

                    # Compute error via active convergence criterion
                    self.convergence_handler.evaluate(self.eq_mom)

                    self.convergence_handler.increment_iteration()

                ite = self.convergence_handler.ite
                converged = not self.convergence_handler.not_converged_error

                if not converged:
                    if n_bisections >= max_bisections:
                        raise RuntimeError(
                            f"Failed to converge after {max_bisections} retries "
                            f"(ite={ite}, maxiter={self.convergence_handler.maxiter})."
                        )
                    self._restore_step_state(
                        step_state,
                        include_heat=False,
                        include_caverns=True,
                    )
                    stress_to = stress_to_step_start.clone()
                    dt_floor = float(getattr(self.t_control, "dt_min", 0.0))
                    self.t_control.dt = max(step_state["time"]["dt"] * 0.5, dt_floor)
                    n_bisections += 1
                    continue

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

                # Save fields
                self.eq_mom.compute_p_nodes()
                self.eq_mom.compute_p_elems()
                self.eq_mom.compute_q_elems()
                self.eq_mom.compute_q_nodes()

                # Persist simulation log row at each converged time step.
                if self.simulation_logger is not None:
                    # Extract yield variables if material has non-elastic elements
                    log_kwargs = SimulationLogging.extract_yield_variables(self.eq_mom.mat)
                    
                    self.simulation_logger.log_step(
                        step=self.t_control.step_counter,
                        stress=stress_to,
                        iteration=ite,
                        nonlinear_error=self.convergence_handler.error,
                        time=t,
                        time_step=dt,
                        time_final=self.t_control.t_final,
                        strain=eps_tot_to,
                        **log_kwargs
                    )

                for output in self.outputs:
                    output.save_fields(t)

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
        merged_solutions: bool = False,
        simulation_logger: SimulationLogging | None = None,
    ):
        self.eq_heat = eq_heat
        self.t_control = t_control
        self.outputs = outputs
        self.caverns = caverns
        self.simulation_logger = simulation_logger

        if caverns is None:
            self.caverns = CavernHandler()

        # Apply merged_solutions flag to all output handlers
        for output in self.outputs:
            output.merged_solutions = merged_solutions

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

        # Thermal-only solver has no nonlinear mechanical iterations.
        ite = 0

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
        convergence_criterion: str = "strain_based",
        merged_solutions: bool = False,
        simulation_logger: SimulationLogging | None = None,
    ):
        self.eq_mom = eq_mom
        self.t_control = t_control
        self.outputs = outputs
        self.compute_elastic_response = compute_elastic_response
        self.convergence_handler = ConvergenceErrorHandler(convergence_criterion)
        self.simulation_logger = simulation_logger

        # Apply merged_solutions flag to all output handlers
        for output in self.outputs:
            output.merged_solutions = merged_solutions

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
                step_state = self._capture_step_state(include_heat=False, include_caverns=False)
                stress_to_step_start = stress_to.clone()

                # Advance time
                self.t_control.advance_time()

                t = self.t_control.t
                dt = self.t_control.dt

                # Update boundary conditions
                self.eq_mom.bc.update_dirichlet(t)
                self.eq_mom.bc.update_neumann(t)

                # Initialize criterion at step start
                maxiter = int(getattr(self.t_control, "maxiter", 40) or 40)
                self.convergence_handler.initialize_step(
                    self.eq_mom,
                    maxiter=maxiter,
                )

                while self.convergence_handler.not_converged():
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

                    # Compute error via active convergence criterion
                    self.convergence_handler.evaluate(self.eq_mom)

                    self.convergence_handler.increment_iteration()

                ite = self.convergence_handler.ite
                converged = not self.convergence_handler.not_converged_error

                if not converged:
                    if n_bisections >= max_bisections:
                        raise RuntimeError(
                            f"Failed to converge after {max_bisections} retries "
                            f"(ite={ite}, maxiter={self.convergence_handler.maxiter})."
                        )
                    self._restore_step_state(
                        step_state,
                        include_heat=False,
                        include_caverns=False,
                    )
                    stress_to = stress_to_step_start.clone()
                    dt_floor = float(getattr(self.t_control, "dt_min", 0.0))
                    self.t_control.dt = max(step_state["time"]["dt"] * 0.5, dt_floor)
                    n_bisections += 1
                    continue

                # Adaptive dt integration via relative convergence ratio.
                if hasattr(self.t_control, "get_next_dt"):
                    maxiter = max(int(self.convergence_handler.maxiter or 1), 1)
                    convergence_ratio = ite / maxiter
                    self.t_control.dt = self.t_control.get_next_dt(
                        convergence_ratio=convergence_ratio,
                        n_bisections=n_bisections,
                        converged=True,
                    )

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
                screen_output_row = [
                    self.t_control.step_counter,
                    self.t_control.dt_used / self.t_control.time_conversion,
                    f"{t / self.t_control.time_conversion} / {self.t_control.t_final / self.t_control.time_conversion}",
                    ite,
                    self.convergence_handler.error,
                ]
                self.screen.print_row(screen_output_row)
                step_completed = True

            # if self.eq_mom.grid.mesh.comm.rank == 0:
            # 	print(t/self.t_control.time_unit, ite, error)
            # 	sys.stdout.flush()
            # 	try:
            # 		print(float(self.eq_mom.mat.elems_ne[-1].Fvp.max()))
            # 		sys.stdout.flush()
            # 	except:
            # 		pass

        self.screen.close()
        self._finalize_outputs()
