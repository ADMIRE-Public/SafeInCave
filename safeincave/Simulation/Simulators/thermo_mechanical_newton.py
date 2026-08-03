# Copyright (c) 2026, The SafeInCave Developers
#
# SPDX-License-Identifier: BSD-3-Clause

"""Thermo-mechanical simulator with the Newton mechanical substep."""

from __future__ import annotations
from typing import Any

from .mechanical_newton import Simulator_MNewton
from ...Utils import numpy2torch
from ...Output.SimLogging import SimulationLogging


class Simulator_TMNewton(Simulator_MNewton):
    """
    Thermo-mechanical simulator on the incremental Newton path.

    Mirrors Simulator_TM's staggered thermal coupling with the mechanical
    substep replaced by the Newton iteration of Simulator_MNewton: per step
    attempt, the heat equation is solved ONCE for ``(t, dt)`` (the thermal
    problem does not depend on displacement), the resulting temperature is
    set on the momentum equation (expanded to state points in quadrature
    mode), cavern thermodynamics are updated, and the mechanical step then
    iterates Newton to equilibrium. Thermal strain enters the elastic trial
    through ``compute_eps_th`` as on the staggered path.

    Cutback restores both the mechanical and the thermal state.
    """

    def __init__(
        self,
        eq_mom: Any,
        eq_heat: Any,
        t_control: Any,
        outputs: Any,
        caverns: Any | None = None,
        compute_elastic_response: bool = True,
        convergence_criterion: Any = "newton_residual",
        maxiter: int = 15,
        simulation_logger: Any | None = None,
        merged_solutions: bool = False,
        smooth_output: bool = False,
        plastic_consistency_tolerance: float = 1e-4,
        line_search: Any = True,
    ):
        super().__init__(
            eq_mom,
            t_control,
            outputs,
            caverns=caverns,
            compute_elastic_response=compute_elastic_response,
            convergence_criterion=convergence_criterion,
            maxiter=maxiter,
            simulation_logger=simulation_logger,
            merged_solutions=merged_solutions,
            smooth_output=smooth_output,
            plastic_consistency_tolerance=plastic_consistency_tolerance,
            line_search=line_search,
        )
        self.eq_heat = eq_heat

    def _set_momentum_temperature(self, initial: bool = False) -> None:
        """Push the heat solution's per-cell temperature into the momentum
        equation, expanded to state points in quadrature mode."""
        T_elems = self.eq_heat.get_T_elems()
        expand = getattr(self.eq_mom, "expand_cell_field", None)
        if expand is not None:
            T_elems = expand(T_elems)
        self.eq_mom.set_T(T_elems)
        if initial:
            self.eq_mom.set_T0(T_elems)

    def run(self) -> None:
        """Run the coupled simulation (see class docstring)."""
        for output in self.outputs:
            output.initialize()

        self.eq_mom.bc.update_dirichlet(self.t_control.t)
        self.eq_mom.bc.update_neumann(self.t_control.t)
        self.eq_mom.bc.update_cavern_bcs(self.caverns)

        self._set_momentum_temperature(initial=True)

        if self.compute_elastic_response:
            self.eq_mom.solve_elastic_response()
            eps_tot_to = self.eq_mom.compute_total_strain()
            stress_to = self.eq_mom.compute_elastic_stress(eps_tot_to)
        else:
            eps_tot_to = self.eq_mom.compute_total_strain()
            stress_to = numpy2torch(
                self.eq_mom.sig.x.array.reshape((self.eq_mom.n_elems, 3, 3))
            )

        self.caverns.calculate_volumes(self.eq_mom.u)
        self.caverns.calculate_initial_conditions()
        self.caverns.record_cavern_data(self.t_control.t)

        self.eq_mom.initialize_rates()

        self.eq_mom.compute_p_elems()
        self.eq_mom.compute_q_elems()
        self.eq_mom.compute_yield_mode()
        self.eq_mom.compute_principal_stresses()
        self.eq_mom.compute_principal_strains()

        for output in self.outputs:
            output.save_fields(0)

        if self.simulation_logger is not None:
            log_kwargs = SimulationLogging.extract_yield_variables(self.eq_mom.mat)
            self.simulation_logger.log_initial_state(
                time=self.t_control.t,
                time_final=self.t_control.t_final,
                stress=stress_to,
                strain=eps_tot_to,
                **log_kwargs,
            )

        self.screen.print_initial_state(
            self.t_control.t,
            self.t_control.t_final,
            self.t_control.time_conversion,
        )

        while self.t_control.keep_looping():
            max_bisections = int(getattr(self.t_control, "max_bisections", 5))
            n_bisections = 0
            step_completed = False

            while not step_completed:
                step_state = self._capture_step_state(
                    include_heat=True, include_caverns=True
                )

                self.t_control.advance_time()
                t = self.t_control.t
                dt = self.t_control.dt

                self.eq_mom.bc.update_dirichlet(t)
                self.eq_mom.bc.update_neumann(t)
                self.eq_heat.bc.update_bcs(t)
                self.eq_heat.bc.update_cavern_bcs(self.caverns)

                # Thermal substep (independent of displacement): solve once,
                # then hold the temperature fixed through the Newton loop.
                self.eq_heat.solve(t, dt)
                self.caverns.calculate_total_heat(
                    dt, self.eq_heat.T, self.eq_heat.k
                )
                self.caverns.update_caverns(t, dt)
                self.eq_mom.bc.update_cavern_bcs(self.caverns)
                self._set_momentum_temperature()

                maxiter = int(
                    getattr(self.t_control, "maxiter", self.maxiter) or self.maxiter
                )
                self.convergence_handler.initialize_step(
                    self.eq_mom, maxiter=maxiter
                )

                converged, stress_to, eps_tot_to = self._newton_step(t, dt)
                ite = self.convergence_handler.ite

                if not converged:
                    if n_bisections >= max_bisections:
                        raise RuntimeError(
                            f"Newton failed to converge after {max_bisections} "
                            f"retries (ite={ite}, maxiter={maxiter})."
                        )
                    self._restore_step_state(
                        step_state, include_heat=True, include_caverns=True
                    )
                    dt_floor = float(getattr(self.t_control, "dt_min", 0.0))
                    cutback = float(
                        getattr(self.t_control, "cutback_factor", 0.25)
                    )
                    self.t_control.dt = max(
                        step_state["time"]["dt"] * cutback, dt_floor
                    )
                    n_bisections += 1
                    continue

                self.eq_mom.commit_step(dt)
                self.eq_heat.update_T_old()

                if hasattr(self.t_control, "get_next_dt"):
                    maxiter_ref = max(int(self.convergence_handler.maxiter or 1), 1)
                    dt_kwargs = dict(
                        convergence_ratio=ite / maxiter_ref,
                        n_bisections=n_bisections,
                        converged=True,
                    )
                    if getattr(self.t_control, "accepts_iterations", False):
                        dt_kwargs["n_iterations"] = ite
                    self.t_control.dt = self.t_control.get_next_dt(**dt_kwargs)

                self.caverns.calculate_volumes(self.eq_mom.u)
                self.caverns.calculate_initial_conditions()
                self.caverns.record_cavern_data(t)

                self.eq_mom.compute_p_elems()
                self.eq_mom.compute_q_elems()
                self.eq_mom.compute_yield_mode()
                self.eq_mom.compute_principal_stresses()
                self.eq_mom.compute_principal_strains()

                if self.simulation_logger is not None:
                    log_kwargs = SimulationLogging.extract_yield_variables(
                        self.eq_mom.mat
                    )
                    self.simulation_logger.log_step(
                        step=self.t_control.step_counter,
                        stress=stress_to,
                        iteration=ite,
                        nonlinear_error=self.convergence_handler.error,
                        time=t,
                        time_step=dt,
                        time_final=self.t_control.t_final,
                        strain=eps_tot_to,
                        **log_kwargs,
                    )

                for output in self.outputs:
                    output.save_fields(t)

                current_time = "%.3f" % (t / self.t_control.time_conversion)
                self.screen.print_row(
                    [
                        self.t_control.step_counter,
                        self.t_control.dt_used / self.t_control.time_conversion,
                        f"{current_time} / "
                        f"{self.t_control.t_final / self.t_control.time_conversion}",
                        ite,
                        self.convergence_handler.error,
                    ]
                )
                step_completed = True

        self.caverns.save_caverns_data()
        self.screen.close()
        self._finalize_outputs()
