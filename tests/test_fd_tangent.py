# Copyright (c) 2026, The SafeInCave Developers
#
# SPDX-License-Identifier: BSD-3-Clause

"""
Finite-difference validation of the plastic models' algorithmic tangents,
plus strain-driven single-element sanity checks. Skips cleanly when the
plastic extension models are not installed.
"""

import sys
from pathlib import Path

import pytest
import torch as to

sys.path.insert(0, str(Path(__file__).parent))

import safeincave as sf  # noqa: E402
from safeincave.Utils import dotdot_torch  # noqa: E402
from helpers.fd_tangent import (  # noqa: E402
    compare_element_tangent,
    compare_element_dep,
)
from helpers.single_element_driver import SingleElementDriver  # noqa: E402

pytestmark = pytest.mark.skipif(
    not hasattr(sf, "PlasticDruckerPrager") or not hasattr(sf, "PlasticRankine"),
    reason="plastic extension models not installed",
)

MPa = 1.0e6
E_MOD = 20.0e9
NU = 0.3


def make_dp(n_elems: int):
    ones = to.ones(n_elems, dtype=to.float64)
    return sf.PlasticDruckerPrager(
        cohesion=5.0 * MPa * ones,
        friction_angle=30.0 * ones,
        dilation_angle=10.0 * ones,
        youngs_modulus=E_MOD * ones,
        poissons_ratio=NU * ones,
        acceleration=None,
    )


def make_rankine(n_elems: int):
    ones = to.ones(n_elems, dtype=to.float64)
    return sf.PlasticRankine(
        tensile_strength=1.0 * MPa * ones,
        youngs_modulus=E_MOD * ones,
        poissons_ratio=NU * ones,
        acceleration=None,
    )


def diag_stress(*states) -> to.Tensor:
    return to.stack(
        [to.diag(to.tensor(s, dtype=to.float64)) * MPa for s in states]
    )


def rotate(stress: to.Tensor, angle_deg: float) -> to.Tensor:
    a = to.tensor(angle_deg * to.pi / 180.0, dtype=to.float64)
    R = to.tensor(
        [
            [to.cos(a), -to.sin(a), 0.0],
            [to.sin(a), to.cos(a), 0.0],
            [0.0, 0.0, 1.0],
        ],
        dtype=to.float64,
    )
    return to.einsum("ij,njk,lk->nil", R, stress, R)


class TestDruckerPragerTangent:
    def test_fd_matches_analytic(self):
        # elastic interior / smooth-cone face / apex, plus a rotated face
        # state exercising the shear columns.
        states = diag_stress(
            (-5.0, -4.0, -3.0),  # elastic
            (-40.0, -5.0, -5.0),  # face
            (30.0, 29.0, 28.0),  # apex
        )
        rotated_face = rotate(diag_stress((-40.0, -5.0, -5.0)), 30.0)
        stress = to.cat([states, rotated_face])
        elem = make_dp(stress.shape[0])

        err = compare_element_tangent(elem, stress)
        expected_modes = to.tensor([0, 1, 2, 1])
        assert to.equal(elem.yield_mode, expected_modes)
        assert float(err.max()) < 2e-6, f"per-element rel error: {err.tolist()}"


class TestRankineTangent:
    def test_fd_matches_analytic(self):
        states = diag_stress(
            (-2.0, -3.0, -4.0),  # elastic
            (3.0, -1.0, -2.0),  # face return
            (4.0, 3.5, -3.0),  # edge return
            (5.0, 4.6, 4.2),  # apex return
        )
        rotated_face = rotate(diag_stress((3.0, -1.0, -2.0)), 30.0)
        stress = to.cat([states, rotated_face])
        elem = make_rankine(stress.shape[0])

        err = compare_element_tangent(elem, stress)
        expected_modes = to.tensor([0, 3, 4, 5, 3])
        assert to.equal(elem.yield_mode, expected_modes)
        assert float(err.max()) < 2e-6, f"per-element rel error: {err.tolist()}"


class TestNewtonInterfaceDP:
    def test_dep_fd_matches(self):
        states = diag_stress(
            (-5.0, -4.0, -3.0),  # elastic
            (-40.0, -5.0, -5.0),  # face
            (30.0, 29.0, 28.0),  # apex
        )
        stress = to.cat([states, rotate(diag_stress((-40.0, -5.0, -5.0)), 30.0)])
        elem = make_dp(stress.shape[0])
        eps_trial = dotdot_torch(elem._C_inv_voigt, stress)

        err = compare_element_dep(elem, eps_trial)
        assert to.equal(elem.yield_mode, to.tensor([0, 1, 2, 1]))
        assert float(err.max()) < 2e-6, f"per-element rel error: {err.tolist()}"

        # Structural checks: elastic D_ep == C, apex D_ep == 0.
        _, D_ep = elem.compute_stress_and_tangent(
            eps_trial, to.zeros(4, dtype=to.float64), 1.0
        )
        assert to.allclose(D_ep[0], elem._C_voigt[0])
        assert float(D_ep[2].abs().max()) < 1e-3 * float(elem._C_voigt[2].abs().max())

    def test_matches_staggered_fixed_point(self, tmp_path):
        # Strain-driven: the Newton update committed once per step must
        # reproduce the staggered exact-trial driver exactly (same return
        # mapping, same commit).
        eye = to.eye(3, dtype=to.float64).unsqueeze(0)
        eps_path = [f * 2.0e-3 * eye for f in (0.25, 0.5, 0.75, 1.0)]

        elem_a = make_dp(1)
        driver = SingleElementDriver([elem_a])
        rec = [driver.step(e) for e in eps_path]
        sigma_staggered = rec[-1]["sigma"]

        elem_b = make_dp(1)
        T = to.zeros(1, dtype=to.float64)
        for eps_tot in eps_path:
            eps_trial = eps_tot - elem_b.eps_ne_old
            sigma_newton, _ = elem_b.compute_stress_and_tangent(eps_trial, T, 1.0)
            elem_b.commit_increment()

        assert to.allclose(sigma_newton, sigma_staggered, rtol=1e-12, atol=1e-3)


class TestNewtonInterfaceRankine:
    def test_dep_fd_matches(self):
        states = diag_stress(
            (-2.0, -3.0, -4.0),  # elastic
            (3.0, -1.0, -2.0),  # face
            (4.0, 3.5, -3.0),  # edge
            (5.0, 4.6, 4.2),  # apex
        )
        elem = make_rankine(states.shape[0])
        eps_trial = dotdot_torch(elem._C_inv_voigt, states)

        err = compare_element_dep(elem, eps_trial)
        assert to.equal(elem.yield_mode, to.tensor([0, 3, 4, 5]))
        assert float(err.max()) < 2e-6, f"per-element rel error: {err.tolist()}"


class TestNewtonInterfacePair:
    def make_pair(self, n):
        rankine = make_rankine(n)
        ones = to.ones(n, dtype=to.float64)
        dp = sf.PlasticDruckerPrager(
            cohesion=5.0 * MPa * ones,
            friction_angle=30.0 * ones,
            dilation_angle=10.0 * ones,
            youngs_modulus=E_MOD * ones,
            poissons_ratio=NU * ones,
            acceleration=None,
            cede_to=rankine,
        )
        return dp, rankine

    def test_pair_classification_and_fd(self):
        # elastic / DP face / Rankine-governed (ceded) / two-surface corner
        states = diag_stress(
            (-5.0, -4.0, -3.0),
            (-40.0, -5.0, -5.0),
            (2.5, 0.5, -3.0),
            (3.0, -0.5, -18.0),
        )
        dp, rankine = self.make_pair(states.shape[0])
        eps_trial = dotdot_torch(dp._C_inv_voigt, states)

        err = compare_element_dep(dp, eps_trial)
        # classification: DP mode reports corner via 6; rankine covers ceded
        assert to.equal(dp.yield_mode, to.tensor([0, 1, 0, 6]))
        assert int(rankine.yield_mode[2]) == 3  # ceded element: Rankine face
        assert int(rankine.yield_mode[3]) == 6

        # smooth regions must be tight; the corner tangent is the derivative
        # of the linearized (frozen-coefficient) one-shot corner return, so
        # its FD gap reflects the neglected coefficient variation.
        assert float(err[:3].max()) < 2e-6, f"rel error: {err.tolist()}"
        assert float(err[3]) < 0.05, f"corner rel error: {float(err[3])}"

    def test_pair_returned_stress_admissible(self):
        states = diag_stress((2.5, 0.5, -3.0), (3.0, -0.5, -18.0))
        dp, rankine = self.make_pair(states.shape[0])
        eps_trial = dotdot_torch(dp._C_inv_voigt, states)
        sigma, _ = dp.compute_stress_and_tangent(
            eps_trial, to.zeros(2, dtype=to.float64), 1.0
        )
        # Rankine cutoff satisfied on both (linearization tolerance)
        sig1 = to.linalg.eigvalsh(sigma)[:, -1]
        assert float((sig1 - 1.0 * MPa).max()) < 5e-2 * MPa
        # consistency errors reported small
        assert dp.consistency_error < 5e-2
        assert rankine.consistency_error < 5e-2


class TestSingleElementDriver:
    def test_dp_hydrostatic_apex_cap(self, tmp_path):
        # A hydrostatic tension strain ramp must cap the stress at the DP
        # apex sigma_m = -f_c/beta once the cone apex is reached.
        elem = make_dp(1)
        eye = to.eye(3, dtype=to.float64).unsqueeze(0)
        eps_path = to.stack([f * 2.0e-3 * eye for f in (0.25, 0.5, 0.75, 1.0)])
        driver = SingleElementDriver([elem])
        records = driver.run(eps_path, csv_path=str(tmp_path / "dp_apex.csv"))

        sigma_m_apex = float(-elem.fc[0] / elem.beta[0])
        sigma_final = records[-1]["sigma"]
        sigma_m_final = float(to.diagonal(sigma_final[0]).mean())
        assert sigma_m_final == pytest.approx(sigma_m_apex, rel=1e-9)
        assert int(records[-1]["drucker_prager_yield_mode"][0]) == 2

    def test_rankine_tension_cap(self, tmp_path):
        # A uniaxial-strain tension ramp must cap the largest principal
        # stress at sigma_t.
        elem = make_rankine(1)
        eps = to.zeros((1, 3, 3), dtype=to.float64)
        eps[0, 0, 0] = 1.0
        eps_path = to.stack([f * 5.0e-4 * eps for f in (0.2, 0.4, 0.6, 0.8, 1.0)])
        driver = SingleElementDriver([elem])
        records = driver.run(eps_path, csv_path=str(tmp_path / "rk_cap.csv"))

        sigma_final = records[-1]["sigma"][0]
        sigma_1 = float(to.linalg.eigvalsh(sigma_final)[-1])
        assert sigma_1 == pytest.approx(float(elem.tensile_strength[0]), rel=1e-9)
        assert int(records[-1]["rankine_yield_mode"][0]) in (3, 4)
