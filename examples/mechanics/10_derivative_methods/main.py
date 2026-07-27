"""
Comparison of three derivative evaluation backends: finite_difference, torch_ad, jax_ad.

This example runs the same mechanical problem (cube under uniaxial loading with
dislocation creep) three times, each with a different derivative_method for
computing the tangent operator E = d(eps_ne)/d(stress). Results are written to
separate output folders for direct comparison.

The three backends all produce correct results (validated in the test suite);
differences in this run demonstrate numerical precision and performance
characteristics across finite-difference, torch autodiff, and JAX autodiff
approaches.
"""

import os
import time

import dolfinx as do
import numpy as np
import torch as to

import safeincave as sf
import safeincave.BC.Momentum as momBC
import safeincave.Utils as ut


class LinearMomentumMod(sf.LinearMomentum):
    def __init__(self, grid, theta):
        super().__init__(grid, theta)
        self.eps_cr = do.fem.Function(self.DG0_3x3)

    def run_after_solve(self):
        self.eps_cr.x.array[:] = to.flatten(self.mat.elems_ne[0].eps_ne_k)


def run(derivative_method):
    """
    Run the mechanics simulation with the specified derivative evaluation method.

    Parameters
    ----------
    derivative_method : str
        One of "finite_difference" (default, central differences),
        "automatic_differentiation_torch" (exact via torch.func.jvp), or
        "automatic_differentiation_jax" (exact via jax.jvp).

    Returns
    -------
    tuple of float
        ``(elapsed, u_mag)`` — simulation solve time (seconds), and the
        final displacement magnitude (m) at node 0.
    """

    # Read grid
    grid_path = os.path.join("..", "..", "..", "grids", "cube")
    grid = sf.GridHandlerGMSH("geom", grid_path)

    # Define output folder
    output_folder = os.path.join("output", "case_0", f"{derivative_method}")

    # Time settings
    unit = "hour"
    t_0 = 0.0
    dt = 0.5
    t_final = 6  # Short run for quick comparison
    t_control = sf.TimeController(
        dt=dt, initial_time=t_0, final_time=t_final, time_unit=unit
    )

    # Define momentum equation
    theta = 0.5
    mom_eq = LinearMomentumMod(grid, theta=theta)

    # Define material properties
    mat = sf.Material(grid.n_elems)

    # Set material density
    rho = 2000.0 * to.ones(grid.n_elems, dtype=to.float64)
    mat.set_density(rho)

    # Elastic element (Spring)
    E_elastic = 102e9 * to.ones(grid.n_elems)
    nu_elastic = 0.3 * to.ones(grid.n_elems)
    spring = sf.Spring(E_elastic, nu_elastic, "spring")

    # Inelastic creep element with selectable derivative method
    A = 1.9e-20 * to.ones(grid.n_elems)
    Q = 51600 * to.ones(grid.n_elems)
    n = 3.0 * to.ones(grid.n_elems)
    creep = sf.DislocationCreep(A, Q, n, "creep", derivative_method=derivative_method)

    # Add elements to material
    mat.add_to_elastic(spring)
    mat.add_to_non_elastic(creep)

    # Set constitutive model
    mom_eq.set_material(mat)

    # Set body forces
    g_vec = [0.0, 0.0, 0.0]
    mom_eq.build_body_force(g_vec)

    # Set initial temperature field
    T0_field = 293 * to.ones(grid.n_elems)
    mom_eq.set_T0(T0_field)
    mom_eq.set_T(T0_field)

    # Boundary conditions (same as 1_triaxial)
    bc_west = momBC.DirichletBC(
        boundary_name="WEST",
        component=0,
        values=[0.0, 0.0],
        time_values=[0.0, t_control.t_final],
    )

    bc_bottom = momBC.DirichletBC(
        boundary_name="BOTTOM",
        component=2,
        values=[0.0, 0.0],
        time_values=[0.0, t_control.t_final],
    )

    bc_south = momBC.DirichletBC(
        boundary_name="SOUTH",
        component=1,
        values=[0.0, 0.0],
        time_values=[0.0, t_control.t_final],
    )

    bc_east = momBC.NeumannBC(
        boundary_name="EAST",
        direction=2,
        density=0.0,
        ref_pos=0.0,
        values=[4.0 * ut.MPa, 4.0 * ut.MPa],
        time_values=[0.0, t_control.t_final],
        g=g_vec[2],
    )

    bc_north = momBC.NeumannBC(
        boundary_name="NORTH",
        direction=2,
        density=0.0,
        ref_pos=0.0,
        values=[4.0 * ut.MPa, 4.0 * ut.MPa],
        time_values=[0.0, t_control.t_final],
        g=g_vec[2],
    )

    bc_top = momBC.NeumannBC(
        boundary_name="TOP",
        direction=2,
        density=0.0,
        ref_pos=0.0,
        values=[4.1 * ut.MPa, 16 * ut.MPa],
        time_values=[0 * ut.hour, t_control.t_final],
        g=g_vec[2],
    )

    bc_handler = momBC.BcHandler()
    bc_handler.add_boundary_condition(bc_west)
    bc_handler.add_boundary_condition(bc_bottom)
    bc_handler.add_boundary_condition(bc_south)
    bc_handler.add_boundary_condition(bc_east)
    bc_handler.add_boundary_condition(bc_north)
    bc_handler.add_boundary_condition(bc_top)

    # Set boundary conditions
    mom_eq.set_boundary_conditions(bc_handler)

    # Create output handlers
    output_mom = sf.SaveFields(mom_eq)
    output_mom.set_output_folder(output_folder)
    output_mom.add_output_field("u", "Displacement (m)")
    output_mom.add_output_field("eps_tot", "Total strain (-)")
    output_mom.add_output_field("eps_cr", "Creep strain (-)")
    outputs = [output_mom]

    # Define and run simulator
    sim = sf.Simulator_M(mom_eq, t_control, outputs, compute_elastic_response=True)
    t_start = time.perf_counter()
    sim.run()
    t_elapsed = time.perf_counter() - t_start

    # Record final displacement magnitude at a representative node
    u = mom_eq.u.x.array.reshape(-1, 3)
    u_mag = float(np.linalg.norm(u[0]))

    return t_elapsed, u_mag


def main():
    """Run the comparison with all three derivative methods and report timing."""
    methods = ["finite_difference", "automatic_differentiation_torch", "automatic_differentiation_jax"]
    results = {}

    for method in methods:
        results[method] = run(method)

    # Print comparison table
    print("\n" + "=" * 70)
    print("Derivative Method Comparison")
    print("=" * 70)
    print(f"{'Method':<32} {'Final_displacement (m)':>23} {'Time (s)':>10}")
    print("-" * 70)
    for method in methods:
        elapsed, u_mag = results[method]
        print(f"{method:<32} {u_mag:>23.6e} {elapsed:>10.3f}")
    print("=" * 70)


if __name__ == "__main__":
    main()
