import safeincave as sf
from safeincave.Utils import day, GPa, create_field_elems
import safeincave.MomentumBC as momBC
import safeincave.CavernBC as caveBC
from mpi4py import MPI
import torch as to
import os
import sys


def get_geometry_parameters(path_to_grid):
    f = open(os.path.join(path_to_grid, "geom.geo"), "r")
    data = f.readlines()
    ovb_thickness = float(data[10][len("ovb_thickness = ") : -2])
    hanging_wall = float(data[12][len("hanging_wall = ") : -2])
    return ovb_thickness, hanging_wall


def main():
    # Read grid
    grid_path = os.path.join("..", "..", "..", "grids", "cavern_regular")
    grid = sf.GridHandlerGMSH("geom", grid_path)

    # Define output folder
    output_folder = os.path.join("output", "case_0")

    # Define momentum equation
    theta = 0.0
    mom_eq = sf.LinearMomentumMixed(grid, theta=theta, stab_scaling=1.0)

    # Define material properties
    mat = sf.Material(mom_eq.n_elems)

    # Set material density
    salt_density = 2200
    rho = salt_density * to.ones(mom_eq.n_elems, dtype=to.float64)
    mat.set_density(rho)

    # Constitutive model
    E0 = 102 * GPa * to.ones(mom_eq.n_elems)
    nu0 = 0.3 * to.ones(mom_eq.n_elems)
    spring_0 = sf.Spring(E0, nu0, "spring")

    # Create creep
    A = 1.9e-20 * to.ones(mom_eq.n_elems)
    Q = 51600 * to.ones(mom_eq.n_elems)
    n = 3.0 * to.ones(mom_eq.n_elems)
    creep_0 = sf.DislocationCreep(A, Q, n, "creep")

    # Create constitutive model
    mat.add_to_elastic(spring_0)
    mat.add_to_non_elastic(creep_0)

    # Set constitutive model
    mom_eq.set_material(mat)

    # Set body forces
    g = -9.81
    g_vec = [0.0, 0.0, g]
    mom_eq.build_body_force(g_vec)

    # Set initial temperature field
    def T_field_fun(x, y, z):
        km = 1000
        dTdZ = 27 / km
        T_surface = 20 + 273
        return T_surface - dTdZ * z

    T0_field = create_field_elems(grid, T_field_fun)
    mom_eq.set_T0(T0_field)
    mom_eq.set_T(T0_field)

    # Time settings for equilibrium stage
    tc_eq = sf.TimeControllerParabolic(
        n_time_steps=100, initial_time=0.0, final_time=365, time_unit="day"
    )

    # Boundary conditions
    bc_equilibrium = momBC.BcHandler()

    # Apply Dirichlet boundary conditions
    boundaries = [("West", 0), ("East", 0), ("South", 1), ("North", 1), ("Bottom", 2)]
    for b_name, component in boundaries:
        bc = momBC.DirichletBC(
            boundary_name=b_name,
            component=component,
            values=[0.0, 0.0],
            time_values=[0.0, tc_eq.t_final],
        )
        bc_equilibrium.add_boundary_condition(bc)

    # Apply overburden
    overburden = 10 * sf.Utils.MPa
    bc_top = momBC.NeumannBC(
        boundary_name="Top",
        direction=2,
        density=0.0,
        ref_pos=0.0,
        values=[overburden, overburden],
        time_values=[0 * day, tc_eq.t_final],
        g=g_vec[2],
    )
    bc_equilibrium.add_boundary_condition(bc_top)

    # Calculate lithostatic pressure at the cavern's roof
    hanging_wall = 430
    p_roof = overburden + salt_density * abs(g) * hanging_wall

    # Define cavern conditions
    cavern_handler = caveBC.CavernHandler()
    cave_1 = caveBC.Cavern_PT(
        grid=grid,
        cavern_name="Cavern",
        fluid="Hydrogen",
        sym_scale=4,
        reference_point=[0.0, 0.0, hanging_wall],
        P_values=[0.8 * p_roof, 0.8 * p_roof],
        T_values=[303, 303],
        time_values=[0 * day, tc_eq.t_final],
        ref_pos=hanging_wall,
        direction=2,
        g=g_vec[2],
    )
    cavern_handler.add_cavern(cave_1)
    cavern_handler.set_output_folder(output_folder)

    # Set boundary conditions
    mom_eq.set_boundary_conditions(bc_equilibrium)

    # Print output folder
    if MPI.COMM_WORLD.rank == 0:
        print(output_folder)
        sys.stdout.flush()

    # Create output handlers
    output_mom = sf.SaveFields(mom_eq)
    output_mom.set_output_folder(output_folder)
    output_mom.add_output_field("u", "Displacement (m)")
    output_mom.add_output_field("sig", "Stress (Pa)")
    output_mom.add_output_field("p_elems", "Mean stress (Pa)")
    output_mom.add_output_field("q_elems", "Von Mises stress (Pa)")
    outputs = [output_mom]

    # Define simulator
    sim = sf.Simulator_M(mom_eq, tc_eq, outputs, cavern_handler, True)
    sim.run()


if __name__ == "__main__":
    main()
