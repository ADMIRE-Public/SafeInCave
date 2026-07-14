# Copyright (c) 2026, The SafeInCave Developers
#
# SPDX-License-Identifier: BSD-3-Clause

"""
Regression tests for the Constant-backed boundary-condition loads.

The failure mode these tests guard against: a compiled form that references
a Neumann/Dirichlet term with the load baked in as a Python float freezes
that load at its compile-time value — time advances, the load doesn't. With
loads in ``fem.Constant`` objects updated in place, a form compiled ONCE must
track the load ramp exactly.
"""

import os

import numpy as np
import pytest
import dolfinx as do
import dolfinx.fem.petsc as fem_petsc

from safeincave import GridHandlerGMSH, LinearMomentum
import safeincave.BC.Momentum as momBC

MPa = 1.0e6
T_FINAL = 10.0


@pytest.fixture(scope="module")
def eq():
    test_dir = os.path.dirname(__file__)
    grid = GridHandlerGMSH("geom", os.path.join(test_dir, "files", "cube_coarse"))
    eq = LinearMomentum(grid, theta=0.5)

    bc_handler = momBC.BcHandler()
    bc_handler.add_boundary_condition(
        momBC.DirichletBC(
            boundary_name="BOTTOM",
            component=2,
            values=[0.0, -1.0e-3],
            time_values=[0.0, T_FINAL],
        )
    )
    bc_handler.add_boundary_condition(
        momBC.NeumannBC(
            boundary_name="TOP",
            direction=2,
            density=0.0,
            ref_pos=0.0,
            values=[1.0 * MPa, 5.0 * MPa],
            time_values=[0.0, T_FINAL],
            g=0.0,
        )
    )
    eq.set_boundary_conditions(bc_handler)
    return eq


def assemble_neumann(eq, form=None):
    """Assemble the Neumann surface terms into a fresh vector."""
    if form is None:
        # eq.bc.normal already carries the test function; each term is a
        # complete linear form (same construction as LinearMomentum.solve).
        form = do.fem.form(sum(eq.bc.neumann_bcs))
    b = fem_petsc.assemble_vector(form)
    b.ghostUpdate()
    return form, b.array.copy()


def reference_neumann(eq, t):
    """Float-baked reference assembly (the pre-Constant construction)."""
    terms = []
    for bc in eq.bc.neumann_boundaries:
        p = -np.interp(t, bc.time_values, bc.values)
        value = p + bc.density * bc.gravity * (bc.ref_pos - eq.bc.x[bc.direction])
        terms.append(
            value
            * eq.bc.normal
            * eq.bc.ds(eq.bc.dolfin_tags[eq.bc.boundary_dim][bc.boundary_name])
        )
    form = do.fem.form(sum(terms))
    b = fem_petsc.assemble_vector(form)
    b.ghostUpdate()
    return b.array.copy()


class TestNeumannConstants:
    def test_cached_form_tracks_load_ramp(self, eq):
        eq.bc.update_neumann(0.0)
        terms_t0 = list(eq.bc.neumann_bcs)
        form, b0 = assemble_neumann(eq)
        assert np.array_equal(b0, reference_neumann(eq, 0.0))

        # advance time WITHOUT recompiling the form
        eq.bc.update_neumann(T_FINAL / 2.0)
        # UFL terms and list must be identity-stable
        assert all(a is b for a, b in zip(terms_t0, eq.bc.neumann_bcs))
        _, b1 = assemble_neumann(eq, form)
        assert np.array_equal(b1, reference_neumann(eq, T_FINAL / 2.0))
        # the load actually changed (guards against a trivially-passing test)
        assert not np.array_equal(b0, b1)

    def test_ramp_is_exact_interpolation(self, eq):
        eq.bc.update_neumann(T_FINAL)
        _, b_end = assemble_neumann(eq)
        eq.bc.update_neumann(0.0)
        _, b_start = assemble_neumann(eq)
        nz = np.abs(b_start) > 1e-30
        # 5 MPa at t_final vs 1 MPa at t0
        assert np.allclose(b_end[nz] / b_start[nz], 5.0, rtol=1e-12)


class TestDirichletConstants:
    def test_bc_objects_stable_and_values_track(self, eq):
        eq.bc.update_dirichlet(0.0)
        bcs_t0 = list(eq.bc.dirichlet_bcs)
        u = do.fem.Function(eq.get_uV())
        do.fem.set_bc(u.x.array, eq.bc.dirichlet_bcs)
        assert np.isclose(u.x.array.min(), 0.0)

        eq.bc.update_dirichlet(T_FINAL)
        assert all(a is b for a, b in zip(bcs_t0, eq.bc.dirichlet_bcs))
        do.fem.set_bc(u.x.array, eq.bc.dirichlet_bcs)
        assert np.isclose(u.x.array.min(), -1.0e-3)
