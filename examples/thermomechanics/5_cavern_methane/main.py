import safeincave as sf
from safeincave.Utils import GPa, create_field_elems, create_field_nodes
import safeincave.MomentumBC as momBC
import safeincave.HeatBC as heatBC
import safeincave.CavernBC as caveBC
from petsc4py import PETSc
import torch as to
import os


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

    # Define solver
    mom_solver = PETSc.KSP().create(grid.mesh.comm)
    mom_solver.setType("gmres")
    mom_solver.getPC().setType("asm")
    mom_solver.setTolerances(rtol=1e-12, max_it=100)
    mom_eq.set_solver(mom_solver)

    # Define heat diffusion equation
    heat_eq = sf.HeatDiffusion(grid)

    # Define solver
    solver_heat = PETSc.KSP().create(grid.mesh.comm)
    solver_heat.setType("cg")
    solver_heat.getPC().setType("asm")
    solver_heat.setTolerances(rtol=1e-12, max_it=100)
    heat_eq.set_solver(solver_heat)

    # Define material properties
    mat = sf.Material(mom_eq.n_elems)

    # Set material density
    salt_density = 2200
    rho = salt_density * to.ones(mom_eq.n_elems, dtype=to.float64)
    mat.set_density(rho)

    # Constitutive model
    E0 = 102 * GPa * to.ones(mom_eq.n_elems, dtype=to.float64)
    nu0 = 0.3 * to.ones(mom_eq.n_elems, dtype=to.float64)
    spring_0 = sf.Spring(E0, nu0, "spring")

    # Create creep
    A = 1.9e-20 * to.ones(mom_eq.n_elems, dtype=to.float64)
    Q = 51600 * to.ones(mom_eq.n_elems, dtype=to.float64)
    n = 3.0 * to.ones(mom_eq.n_elems, dtype=to.float64)
    creep_ds = sf.DislocationCreep(A, Q, n, "ds_creep")

    # Create pressure solution creep
    A = 1.29e-29 * to.ones(mom_eq.n_elems, dtype=to.float64)
    Q = 13184 * to.ones(mom_eq.n_elems, dtype=to.float64)
    d = 0.01 * to.ones(mom_eq.n_elems, dtype=to.float64)
    creep_ps = sf.PressureSolutionCreep(A, d, Q, "ps_creep")

    # Thermo-elastic element
    # alpha = 44e-7*to.ones(mom_eq.n_elems)
    alpha = 120e-6 * to.ones(mom_eq.n_elems)
    thermo = sf.Thermoelastic(alpha, "thermo")

    # Set specific heat capacity
    cp = 850 * to.ones(heat_eq.n_elems, dtype=to.float64)
    mat.set_specific_heat_capacity(cp)

    # Set thermal conductivity
    k = 7 * to.ones(heat_eq.n_elems, dtype=to.float64)
    mat.set_thermal_conductivity(k)

    # Create constitutive model
    mat.add_to_elastic(spring_0)
    mat.add_to_thermoelastic(thermo)
    mat.add_to_non_elastic(creep_ds)
    mat.add_to_non_elastic(creep_ps)

    # Set material properties to governing equations
    mom_eq.set_material(mat)
    heat_eq.set_material(mat)

    # Set body forces
    g = -9.81
    g_vec = [0.0, 0.0, g]
    mom_eq.build_body_force(g_vec)

    # Set initial temperature field
    km = 1000
    dTdZ = 27 / km
    T_top = 273 + 20

    def T_field_fun(x, y, z):
        return T_top + dTdZ * (660 - z)

    T0_field = create_field_elems(grid, T_field_fun)
    mom_eq.set_T0(T0_field)
    mom_eq.set_T(T0_field)

    # Time settings
    time_ctrl = sf.TimeController(
        dt=1.0, initial_time=0.0, final_time=365, time_unit="day"
    )
    time_values = [time_ctrl.t_initial, time_ctrl.t_final]
    nt = len(time_values)

    # Boundary conditions for momentum balance equation
    bc_mom = momBC.BcHandler()

    # Apply Dirichlet boundary conditions
    boundaries = [("West", 0), ("East", 0), ("South", 1), ("North", 1), ("Bottom", 2)]
    for b_name, component in boundaries:
        bc = momBC.DirichletBC(
            boundary_name=b_name,
            component=component,
            values=nt * [0.0],
            time_values=time_values,
        )
        bc_mom.add_boundary_condition(bc)

    # Apply overburden
    overburden = 10 * sf.Utils.MPa
    bc_top = momBC.NeumannBC(
        boundary_name="Top",
        direction=2,
        density=0.0,
        ref_pos=0.0,
        values=nt * [overburden],
        time_values=time_values,
        g=g_vec[2],
    )
    bc_mom.add_boundary_condition(bc_top)

    # Set boundary conditions for momentum balance equation
    mom_eq.set_boundary_conditions(bc_mom)

    # Set initial temperature
    T0_field_nodes = create_field_nodes(grid, T_field_fun)
    heat_eq.set_initial_T(T0_field_nodes)

    # Boundary conditions for heat diffusion equation
    bc_heat = heatBC.BcHandler()
    bc_top = heatBC.DirichletBC("Top", nt * [T_top], time_values)
    bc_bottom = heatBC.NeumannBC("Bottom", nt * [dTdZ], time_values)
    bc_west = heatBC.NeumannBC("West", nt * [0.0], time_values)
    bc_east = heatBC.NeumannBC("East", nt * [0.0], time_values)
    bc_north = heatBC.NeumannBC("North", nt * [0.0], time_values)
    bc_south = heatBC.NeumannBC("South", nt * [0.0], time_values)
    bc_heat.add_boundary_condition(bc_top)
    bc_heat.add_boundary_condition(bc_bottom)
    bc_heat.add_boundary_condition(bc_west)
    bc_heat.add_boundary_condition(bc_east)
    bc_heat.add_boundary_condition(bc_north)
    bc_heat.add_boundary_condition(bc_south)
    heat_eq.set_boundary_conditions(bc_heat)

    # Calculate lithostatic pressure at the cavern's roof
    hanging_wall = 430
    p_roof = overburden + salt_density * abs(g) * hanging_wall
    print(0.5 * p_roof)

    # Read flow rate values
    data_cavern = sf.Utils.read_json("input_massflux.json")

    # Define cavern conditions
    cavern_handler = caveBC.CavernHandler()
    cave_0 = caveBC.Cavern_MassFlux(
        grid=grid,
        cavern_name="Cavern",
        fluid="Methane",
        sym_scale=4,
        reference_point=[0.0, 0.0, hanging_wall],
        P_init=0.5 * p_roof,
        T_init=T_top,
        T_in=T_top,
        Mflux_values=data_cavern["flow"],
        time_values=data_cavern["time"],
        direction=2,
        h_conv=5.0,
        g=g_vec[2],
    )
    cavern_handler.add_cavern(cave_0)
    cavern_handler.set_output_folder(output_folder)

    # Create output handlers
    output_mom = sf.SaveFields(mom_eq)
    output_mom.set_output_folder(output_folder)
    output_mom.add_output_field("u", "Displacement (m)")
    output_mom.add_output_field("sig", "Stress (Pa)")
    output_mom.add_output_field("p_elems", "Mean stress (Pa)")
    output_mom.add_output_field("q_elems", "Von Mises stress (Pa)")

    output_heat = sf.SaveFields(heat_eq)
    output_heat.set_output_folder(output_folder)
    output_heat.add_output_field("T", "Temperature (K)")

    outputs = [output_mom, output_heat]

    # Define simulator
    sim = sf.Simulator_Full(
        eq_mom=mom_eq,
        eq_heat=heat_eq,
        t_control=time_ctrl,
        outputs=outputs,
        caverns=cavern_handler,
        compute_elastic_response=True,
    )
    sim.run()


if __name__ == "__main__":
    main()
