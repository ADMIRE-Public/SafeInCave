# Copyright (c) 2026, The SafeInCave Developers
#
# SPDX-License-Identifier: BSD-3-Clause

"""Mechanical simulator driving the incremental Newton solve path."""

from __future__ import annotations
import os
from typing import Any

import numpy as np

from .base import Simulator
from ...Utils import numpy2torch
from ...Output.SimLogging import SimulationLogging
from ...Output.Screen import ScreenPrinter
from ...Cavern import CavernHandler
from ..Convergence import ConvergenceErrorHandler
from ..Convergence.line_search import BacktrackingLineSearch


class Simulator_MNewton(Simulator):
    """
    Mechanical-only simulator on the incremental Newton path.

    Per time step: snapshot state, advance time, update BC/cavern loads,
    then iterate constitutive update → residual assembly → Newton correction
    with the consistent tangent until the NewtonResidualCriterion (plus the
    plastic-consistency gate) is satisfied; commit the plastic increments
    once per accepted step. Non-convergence restores the snapshot and
    retries with dt cut back (×0.25, clamped to dt_min), reusing the same
    bisection scaffold as Simulator_M.

    Differences from Simulator_M (staggered):

    - Requires every non-elastic element to support the incremental-Newton
      constitutive interface (validated at construction).
    - No ``_final_consuming_solve``: the residual is evaluated at the state
      that gets committed, by construction.
    - Cavern thermodynamics/loads are updated once per step attempt (not per
      iteration); volumes are recomputed from the converged displacement.

    Set ``SAFEINCAVE_NEWTON_TRACE=<path>`` to append a per-iteration CSV
    (step, ite, force_error, disp_error, consistency).
    """

    def __init__(
        self,
        eq_mom: Any,
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
        self.eq_mom = eq_mom
        self.t_control = t_control
        self.outputs = outputs
        self.caverns = caverns if caverns is not None else CavernHandler()
        self.compute_elastic_response = compute_elastic_response
        self.convergence_handler = ConvergenceErrorHandler(
            convergence_criterion,
            plastic_consistency_tolerance=plastic_consistency_tolerance,
        )
        self.maxiter = int(maxiter)
        self.simulation_logger = simulation_logger
        # Backtracking line search on Newton corrections: True (default
        # instance), an instance, or False/None to disable.
        if line_search is True:
            self.line_search = BacktrackingLineSearch()
        elif line_search:
            self.line_search = line_search
        else:
            self.line_search = None

        # Fail fast on materials the Newton path cannot drive.
        self.eq_mom._newton_elements()

        for output in self.outputs:
            output.merged_solutions = merged_solutions
            output.smooth_output = smooth_output

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

        self._trace_path = os.environ.get("SAFEINCAVE_NEWTON_TRACE")
        self._trace_header_written = False

    # ------------------------------------------------------------------

    def _trace(self, step: int, ite: int, force_error: float, disp_error: float):
        if not self._trace_path or self.eq_mom.grid.mesh.comm.rank != 0:
            return
        with open(self._trace_path, "a") as f:
            if not self._trace_header_written:
                if f.tell() == 0:
                    f.write("step,ite,force_error,disp_error,consistency\n")
                self._trace_header_written = True
            f.write(
                f"{step},{ite},{force_error:.6e},{disp_error:.6e},"
                f"{self.convergence_handler.plastic_consistency_error:.6e}\n"
            )

    def _newton_step(self, t: float, dt: float) -> tuple:
        """
        Run the Newton iteration for one step attempt. The first correction
        of a step is always applied in full (Dirichlet-gap increment); later
        corrections go through the backtracking line search when enabled
        (solve-free trials — the accepted trial's constitutive/residual
        evaluation is reused as the next iteration's evaluation).

        Returns ``(converged, stress_to, eps_tot_to)``.
        """
        eq = self.eq_mom
        handler = self.convergence_handler
        u_start = eq.u.x.array.copy()
        du_norm = None
        alpha = 1.0

        stress_to, eps_tot_to = eq.constitutive_update(dt)
        r_norm, ref_norm = eq.assemble_residual()

        while True:
            force_error = r_norm / ref_norm
            if du_norm is None:
                disp_error = 0.0
            else:
                n_owned = (
                    eq.V.dofmap.index_map.size_local * eq.V.dofmap.index_map_bs
                )
                step_increment = (eq.u.x.array - u_start)[:n_owned]
                step_norm = eq._global_norm(np.ascontiguousarray(step_increment))
                disp_error = (alpha * du_norm) / max(step_norm, 1e-16)

            eq._newton_state = {
                "force_error": force_error,
                "disp_error": disp_error,
            }
            handler.evaluate(eq)
            self._trace(
                self.t_control.step_counter, handler.ite, force_error, disp_error
            )
            if not handler.not_converged_error:
                break
            if not handler.below_max_iterations:
                break

            use_line_search = self.line_search is not None and du_norm is not None
            try:
                du_norm = eq.solve_increment(apply=not use_line_search)
            except RuntimeError:
                # Failed direct solve: treat as non-convergence → cutback.
                handler.not_converged_error = True
                break

            if use_line_search:
                alpha, r_norm, ref_norm = self.line_search.search(eq, dt, r_norm)
                stress_to = numpy2torch(
                    eq.sig.x.array.reshape((eq.n_elems, 3, 3))
                )
                eps_tot_to = numpy2torch(
                    eq.eps_tot.x.array.reshape((eq.n_elems, 3, 3))
                )
            else:
                alpha = 1.0
                stress_to, eps_tot_to = eq.constitutive_update(dt)
                r_norm, ref_norm = eq.assemble_residual()
            handler.increment_iteration()

        converged = not handler.not_converged_error
        return converged, stress_to, eps_tot_to

    # ------------------------------------------------------------------

    def run(self) -> None:
        """Run the simulation (see class docstring for the step algorithm)."""
        for output in self.outputs:
            output.initialize()

        self.eq_mom.bc.update_dirichlet(self.t_control.t)
        self.eq_mom.bc.update_neumann(self.t_control.t)
        self.eq_mom.bc.update_cavern_bcs(self.caverns)

        if self.compute_elastic_response:
            self.eq_mom.solve_elastic_response()
            eps_tot_to = self.eq_mom.compute_total_strain()
            stress_to = self.eq_mom.compute_elastic_stress(eps_tot_to)
        else:
            eps_tot_to = numpy2torch(
                self.eq_mom.eps_tot.x.array.reshape((self.eq_mom.n_elems, 3, 3))
            )
            stress_to = numpy2torch(
                self.eq_mom.sig.x.array.reshape((self.eq_mom.n_elems, 3, 3))
            )

        self.caverns.calculate_volumes(self.eq_mom.u)
        self.caverns.calculate_initial_conditions()
        self.caverns.record_cavern_data(self.t_control.t)

        # Seed creep rate carryover at the initial stress (theta scheme).
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
                    include_heat=False, include_caverns=True
                )

                self.t_control.advance_time()
                t = self.t_control.t
                dt = self.t_control.dt

                self.eq_mom.bc.update_dirichlet(t)
                self.eq_mom.bc.update_neumann(t)

                # Cavern loads held fixed within the step (updated from the
                # state at the start of the attempt).
                self.caverns.update_caverns(t, dt)
                self.eq_mom.bc.update_cavern_bcs(self.caverns)

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
                        step_state, include_heat=False, include_caverns=True
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

                # Commit the accepted increments (creep elements first, then
                # each Newton-capable element's staged share).
                self.eq_mom.commit_step(dt)

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
