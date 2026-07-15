# Copyright (c) 2026, The SafeInCave Developers
#
# SPDX-License-Identifier: BSD-3-Clause

"""
Integration tests for the incremental Newton solver path
(LinearMomentumNewton + Simulator_MNewton + NewtonResidualCriterion).
"""

import os

import numpy as np
import pytest
import torch as to

import safeincave as sf
import safeincave.BC.Momentum as momBC

MPa = 1.0e6
HAS_PLASTIC = hasattr(sf, "PlasticDruckerPrager")

E_MOD = 20.0e9
NU = 0.3
COHESION = 3.0 * MPa
FRICTION = 30.0
DILATION = 10.0
# Uniaxial-stress plateau: fc = 2 c cos(phi) / (1 - sin(phi))
PLATEAU = 2.0 * COHESION * np.cos(np.radians(FRICTION)) / (1.0 - np.sin(np.radians(FRICTION)))

T_FINAL = 1.0
TOP_DISPLACEMENT = -2.0e-3  # well beyond first yield (~5.2e-4)


def build_equation(order=1):
    test_dir = os.path.dirname(__file__)
    grid = sf.GridHandlerGMSH("geom", os.path.join(test_dir, "files", "cube_coarse"))
    return grid, sf.LinearMomentumNewton(grid, theta=0.5, order=order)


def elastic_material(n_elems):
    mat = sf.Material(n_elems)
    mat.set_density(2000.0 * to.ones(n_elems, dtype=to.float64))
    ones = to.ones(n_elems)
    mat.add_to_elastic(sf.Spring(E_MOD * ones, NU * ones, "spring"))
    return mat


def dp_material(n_elems):
    mat = elastic_material(n_elems)
    ones = to.ones(n_elems, dtype=to.float64)
    plastic = sf.PlasticDruckerPrager(
        cohesion=COHESION * ones,
        friction_angle=FRICTION * ones,
        dilation_angle=DILATION * ones,
        youngs_modulus=E_MOD * ones,
        poissons_ratio=NU * ones,
        acceleration=None,
    )
    mat.add_to_non_elastic(plastic)
    return mat


def uniaxial_bcs(t_final=T_FINAL, top_displacement=TOP_DISPLACEMENT):
    """Rollers on WEST(x)/SOUTH(y)/BOTTOM(z), z-displacement ramp on TOP."""
    bc = momBC.BcHandler()
    for name, comp in (("WEST", 0), ("SOUTH", 1), ("BOTTOM", 2)):
        bc.add_boundary_condition(
            momBC.DirichletBC(
                boundary_name=name,
                component=comp,
                values=[0.0, 0.0],
                time_values=[0.0, t_final],
            )
        )
    bc.add_boundary_condition(
        momBC.DirichletBC(
            boundary_name="TOP",
            component=2,
            values=[0.0, top_displacement],
            time_values=[0.0, t_final],
        )
    )
    return bc


def setup(material_factory, order=1):
    grid, eq = build_equation(order)
    mat = material_factory(eq.n_state_points)
    eq.set_material(mat)
    eq.build_body_force([0.0, 0.0, 0.0])
    T0 = 293.0 * to.ones(eq.n_state_points)
    eq.set_T0(T0)
    eq.set_T(T0)
    eq.set_boundary_conditions(uniaxial_bcs())
    return grid, eq


class TestNewtonKernel:
    @pytest.mark.parametrize("order", [1, 2])
    def test_elastic_single_iteration(self, order):
        """The cheapest total-pipeline check: for a linear elastic problem
        the first Newton correction must close the residual to roundoff
        (validates Dirichlet lifting, Voigt conventions, f_ext and — for
        order=2 — the quadrature state space at once)."""
        _, eq = setup(elastic_material, order=order)
        t = T_FINAL / 2.0
        eq.bc.update_dirichlet(t)
        eq.bc.update_neumann(t)
        eq.bc.update_cavern_bcs(sf.Cavern.CavernHandler())

        eq.constitutive_update(dt=0.5)
        r0, ref = eq.assemble_residual()
        assert r0 / ref > 1e-3  # BC gap produces a real residual

        eq.solve_increment()
        eq.constitutive_update(dt=0.5)
        r1, ref = eq.assemble_residual()
        assert r1 / ref < 1e-10, f"elastic residual after 1 iteration: {r1 / ref:.3e}"

        # Dirichlet values are met exactly.
        u_z = eq.u.x.array.reshape((-1, 3))[:, 2]
        assert np.isclose(u_z.min(), TOP_DISPLACEMENT / 2.0, rtol=1e-10)

    @pytest.mark.skipif(not HAS_PLASTIC, reason="plastic extensions not installed")
    @pytest.mark.parametrize("order", [1, 2])
    def test_plastic_quadratic_convergence(self, order):
        """Beyond yield, Newton with the consistent tangent must converge in
        a handful of iterations with a superlinear residual tail."""
        _, eq = setup(dp_material, order=order)
        eq.bc.update_dirichlet(T_FINAL)  # full displacement, deep plastic
        eq.bc.update_neumann(T_FINAL)
        eq.bc.update_cavern_bcs(sf.Cavern.CavernHandler())

        history = []
        for _ in range(12):
            eq.constitutive_update(dt=1.0)
            r, ref = eq.assemble_residual()
            history.append(r / ref)
            if r / ref < 1e-10:
                break
            eq.solve_increment()

        assert history[-1] < 1e-10, f"residual history: {history}"
        assert len(history) <= 8, f"too many Newton iterations: {history}"


@pytest.mark.skipif(not HAS_PLASTIC, reason="plastic extensions not installed")
class TestSimulatorNewton:
    def test_uniaxial_plateau(self):
        _, eq = setup(dp_material)
        t_control = sf.TimeController(
            dt=0.2, initial_time=0.0, final_time=T_FINAL, time_unit="second"
        )
        sim = sf.Simulator_MNewton(
            eq,
            t_control,
            outputs=[],
            compute_elastic_response=False,
            maxiter=15,
        )
        sim.run()

        sig_zz = to.from_numpy(
            eq.sig.x.array.reshape((eq.n_elems, 3, 3))[:, 2, 2].copy()
        )
        assert float(sig_zz.mean()) == pytest.approx(-PLATEAU, rel=1e-6), (
            f"plateau: expected {-PLATEAU:.1f}, got {float(sig_zz.mean()):.1f}"
        )
        # Committed plastic strain is nonzero and consistent
        plastic = eq.mat.elems_ne[0]
        assert float(plastic.eps_ne_old.abs().max()) > 1e-4
        assert plastic.consistency_error < 1e-8

    def test_snapshot_restores_committed_state(self):
        """Cutback restore must round-trip eps_ne_old and staged increments
        (risk R6: incomplete state restore on bisection)."""
        _, eq = setup(dp_material)
        t_control = sf.TimeController(
            dt=0.5, initial_time=0.0, final_time=T_FINAL, time_unit="second"
        )
        sim = sf.Simulator_MNewton(
            eq, t_control, outputs=[], compute_elastic_response=False
        )
        plastic = eq.mat.elems_ne[0]
        plastic.eps_ne_old = to.rand_like(plastic.eps_ne_old)
        reference = plastic.eps_ne_old.clone()

        state = sim._capture_step_state(include_heat=False, include_caverns=True)
        plastic.eps_ne_old += 1.0
        plastic._staged_increment = to.ones_like(reference)
        sim._restore_step_state(state, include_heat=False, include_caverns=True)

        assert to.equal(plastic.eps_ne_old, reference)


class TestTimeControllerNewton:
    def make(self, **kw):
        return sf.TimeControllerNewton(
            initial_time=0.0, initial_dt=0.1, final_time=1.0,
            time_unit="second", dt_min=0.001, dt_max=0.5, **kw,
        )

    def test_growth_after_two_easy_steps(self):
        tc = self.make()
        dt0 = tc.dt
        # first easy step: wait; second easy step: grow x1.5
        assert tc.get_next_dt(converged=True, n_iterations=3) == pytest.approx(dt0)
        assert tc.get_next_dt(converged=True, n_iterations=3) == pytest.approx(
            1.5 * dt0
        )

    def test_cutback_and_streak_reset(self):
        tc = self.make()
        dt0 = tc.dt
        tc.get_next_dt(converged=True, n_iterations=3)  # easy 1
        assert tc.get_next_dt(converged=False, n_iterations=15) == pytest.approx(
            0.25 * dt0
        )
        # streak was reset: one easy step is not enough to grow again
        assert tc.get_next_dt(converged=True, n_iterations=3) == pytest.approx(dt0)

    def test_clamping(self):
        tc = self.make()
        tc.dt = 0.4
        tc._easy_streak = 5
        assert tc.get_next_dt(converged=True, n_iterations=1) == pytest.approx(0.5)
        tc.dt = 0.002
        assert tc.get_next_dt(converged=False, n_iterations=15) == pytest.approx(0.001)


@pytest.mark.skipif(not HAS_PLASTIC, reason="plastic extensions not installed")
class TestRobustness:
    def test_duvaut_lions_fd_consistency(self):
        """The Duvaut-Lions blend must stay tangent-consistent (FD gate)."""
        import sys as _sys
        from pathlib import Path as _Path
        _sys.path.insert(0, str(_Path(__file__).parent))
        from helpers.fd_tangent import compare_element_dep
        from safeincave.Utils import dotdot_torch

        ones = to.ones(2, dtype=to.float64)
        elem = sf.PlasticDruckerPrager(
            cohesion=COHESION * ones,
            friction_angle=FRICTION * ones,
            dilation_angle=DILATION * ones,
            youngs_modulus=E_MOD * ones,
            poissons_ratio=NU * ones,
            acceleration=None,
            viscous_regularization_tau=0.3,
        )
        stress = to.stack([
            to.diag(to.tensor([-40.0, -5.0, -5.0], dtype=to.float64)) * MPa,
            to.diag(to.tensor([30.0, 29.0, 28.0], dtype=to.float64)) * MPa,
        ])
        eps_trial = dotdot_torch(elem._C_inv_voigt, stress)
        err = compare_element_dep(elem, eps_trial, dt=1.0)
        assert float(err.max()) < 2e-6, f"rel error: {err.tolist()}"

    def test_forced_cutback_recovers(self):
        """An absurd first dt with a tiny iteation budget must exercise the
        cutback path and still finish the ramp."""
        _, eq = setup(dp_material)
        t_control = sf.TimeControllerNewton(
            initial_time=0.0, initial_dt=1.0, final_time=T_FINAL,
            time_unit="second", dt_min=0.01, dt_max=1.0, max_bisections=8,
        )
        sim = sf.Simulator_MNewton(
            eq, t_control, outputs=[], compute_elastic_response=False, maxiter=3,
        )
        sim.run()
        sig_zz = eq.sig.x.array.reshape((eq.n_elems, 3, 3))[:, 2, 2]
        assert float(np.mean(sig_zz)) == pytest.approx(-PLATEAU, rel=1e-5)


class TestCreepComposition:
    """Phase 4: rate-based creep composed staggered-within-Newton."""

    def _creep_material(self, n):
        mat = elastic_material(n)
        ones = to.ones(n, dtype=to.float64)
        creep = sf.DislocationCreep(1.9e-20 * ones, 51600.0 * ones, 3.0 * ones, "creep")
        mat.add_to_non_elastic(creep)
        return mat

    def _bcs_const_load(self, t_final):
        bc = momBC.BcHandler()
        for name, comp in (("WEST", 0), ("SOUTH", 1), ("BOTTOM", 2)):
            bc.add_boundary_condition(momBC.DirichletBC(
                boundary_name=name, component=comp,
                values=[0.0, 0.0], time_values=[0.0, t_final]))
        bc.add_boundary_condition(momBC.NeumannBC(
            boundary_name="TOP", direction=2, density=0.0, ref_pos=0.0,
            values=[12.0 * MPa, 12.0 * MPa], time_values=[0.0, t_final], g=0.0))
        return bc

    def _run(self, newton: bool):
        test_dir = os.path.dirname(__file__)
        grid = sf.GridHandlerGMSH("geom", os.path.join(test_dir, "files", "cube_coarse"))
        if newton:
            eq = sf.LinearMomentumNewton(grid, theta=0.5, order=1)
        else:
            eq = sf.LinearMomentum(grid, theta=0.5, solver_name="preonly", preconditioner="lu")
        eq.set_material(self._creep_material(eq.n_elems))
        eq.build_body_force([0.0, 0.0, 0.0])
        T0 = 293.0 * to.ones(eq.n_elems)
        eq.set_T0(T0)
        eq.set_T(T0)
        t_control = sf.TimeController(dt=0.25, initial_time=0.0, final_time=1.0, time_unit="hour")
        eq.set_boundary_conditions(self._bcs_const_load(t_control.t_final))
        cls = sf.Simulator_MNewton if newton else sf.Simulator_M
        sim = cls(eq, t_control, outputs=[], compute_elastic_response=True)
        sim.run()
        creep = eq.mat.elems_ne[0]
        return (
            to.from_numpy(eq.sig.x.array.reshape((eq.n_elems, 3, 3)).copy()),
            creep.eps_ne_old.clone(),
        )

    def test_creep_newton_matches_staggered(self):
        sig_stag, epsne_stag = self._run(newton=False)
        sig_newt, epsne_newt = self._run(newton=True)
        scale = float(sig_stag.abs().max())
        assert float((sig_newt - sig_stag).abs().max()) < 1e-4 * scale
        eps_scale = max(float(epsne_stag.abs().max()), 1e-30)
        assert float((epsne_newt - epsne_stag).abs().max()) < 1e-3 * eps_scale
