import safeincave as sf
import safeincave.Utils as ut
from safeincave.Utils import GPa, MPa, day, hour, create_field_elems, create_field_nodes
import safeincave.HeatBC as heatBC
import safeincave.MomentumBC as momBC
from petsc4py import PETSc
from mpi4py import MPI
import os
import sys
import torch as to


def get_geometry_parameters(path_to_grid):
	f = open(os.path.join(path_to_grid, "geom.geo"), "r")
	data = f.readlines()
	ovb_thickness = float(data[10][len("ovb_thickness = "):-2])
	salt_thickness = float(data[11][len("salt_thickness = "):-2])
	hanging_wall = float(data[12][len("hanging_wall = "):-2])
	return ovb_thickness, salt_thickness, hanging_wall


def main():
	# Read grid
	grid_path = os.path.join("..", "..", "..", "grids", "cavern_overburden_coarse")
	grid = sf.GridHandlerGMSH("geom", grid_path)

	# Define output folder
	output_folder = os.path.join("output", "case_1")

	# Extract region indices
	ind_salt = grid.region_indices["Salt"]
	ind_ovb = grid.region_indices["Overburden"]

	# Define momentum equation
	mom_eq = sf.LinearMomentum(grid, theta=0.0)

	# Define solver
	mom_solver = PETSc.KSP().create(grid.mesh.comm)
	mom_solver.setType("bcgs")
	mom_solver.getPC().setType("asm")
	mom_solver.setTolerances(rtol=1e-12, max_it=100)
	mom_eq.set_solver(mom_solver)

	# Define material properties
	mat = sf.Material(mom_eq.n_elems)

	# Set material density
	gas_density = 0.082
	salt_density = 2200
	ovb_density = 2800
	rho = to.zeros(mom_eq.n_elems, dtype=to.float64)
	rho[ind_salt] = salt_density
	rho[ind_ovb] = ovb_density
	mat.set_density(rho)

	# Constitutive model
	E0 = to.zeros(mom_eq.n_elems)
	E0[ind_salt] = 102*GPa
	E0[ind_ovb] = 180*GPa
	nu0 = 0.3*to.ones(mom_eq.n_elems)
	spring_0 = sf.Spring(E0, nu0, "spring")

	# Create Kelvin-Voigt viscoelastic element
	# eta = 105e11*to.ones(mom_eq.n_elems)
	E1 = 10*GPa*to.ones(mom_eq.n_elems)
	nu1 = 0.32*to.ones(mom_eq.n_elems)

	eta = to.zeros(mom_eq.n_elems)
	eta[ind_salt] = 105e11
	eta[ind_ovb] = 105e21
	kelvin = sf.Viscoelastic(eta, E1, nu1, "kelvin")

	# Create dislocation creep
	A = to.zeros(mom_eq.n_elems)
	A[ind_salt] = 1.9e-20
	A[ind_ovb] = 0.0
	Q = 51600*to.ones(mom_eq.n_elems)
	n = 3.0*to.ones(mom_eq.n_elems)
	creep_ds = sf.DislocationCreep(A, Q, n, "ds_creep")

	# Create pressure solution creep
	A = to.zeros(mom_eq.n_elems)
	A[ind_salt] = 1.29e-19
	A[ind_ovb] = 0.0
	Q = 13184*to.ones(mom_eq.n_elems)
	d = 0.01*to.ones(mom_eq.n_elems)
	creep_ps = sf.PressureSolutionCreep(A, d, Q, "ps_creep")

	# Thermo-elastic element
	alpha = to.zeros(mom_eq.n_elems)
	alpha[ind_salt] = 44e-6
	alpha[ind_ovb] = 0.0
	thermo = sf.Thermoelastic(alpha, "thermo")

	# Create constitutive model
	mat.add_to_elastic(spring_0)
	mat.add_to_thermoelastic(thermo)
	mat.add_to_non_elastic(kelvin)
	mat.add_to_non_elastic(creep_ds)
	mat.add_to_non_elastic(creep_ps)

	# Set constitutive model
	mom_eq.set_material(mat)

	# Set body forces
	g = -9.81
	g_vec = [0.0, 0.0, g]
	mom_eq.build_body_force(g_vec)

	# Set initial temperature field
	km = 1000
	dTdZ = 27/km
	T_top = 273 + 20
	T_field_fun = lambda x,y,z: T_top + dTdZ*(660 - z)
	T0_field = create_field_elems(grid, T_field_fun)
	mom_eq.set_T0(T0_field)
	mom_eq.set_T(T0_field)

	# Time settings for equilibrium stage
	tc_eq = sf.TimeControllerParabolic(n_time_steps=20, initial_time=0.0, final_time=10, time_unit="day")
	# tc_eq = sf.TimeController(dt=0.1, final_time=5, initial_time=0.0, time_unit="day")

	# Boundary conditions
	time_values = [0*hour,  1*hour]
	nt = len(time_values)

	bc_west_salt = momBC.DirichletBC(boundary_name="West_salt", component=0, values=[0.0, 0.0], time_values=[0.0, tc_eq.t_final])
	bc_west_ovb = momBC.DirichletBC(boundary_name = "West_ovb", component=0, values=[0.0, 0.0], time_values=[0.0, tc_eq.t_final])

	bc_east_salt = momBC.DirichletBC(boundary_name="East_salt", component=0, values=[0.0, 0.0], time_values=[0.0, tc_eq.t_final])
	bc_east_ovb = momBC.DirichletBC(boundary_name = "East_ovb", component=0, values=[0.0, 0.0], time_values=[0.0, tc_eq.t_final])

	bc_bottom = momBC.DirichletBC(boundary_name="Bottom", component=2, values=[0.0, 0.0], time_values=[0.0, tc_eq.t_final])

	bc_south_salt = momBC.DirichletBC(boundary_name="South_salt", component=1, values=[0.0, 0.0], time_values=[0.0, tc_eq.t_final])
	bc_south_ovb = momBC.DirichletBC(boundary_name="South_ovb", component=1, values=[0.0, 0.0], time_values=[0.0, tc_eq.t_final])

	bc_north_salt = momBC.DirichletBC(boundary_name="North_salt", component=1, values=[0.0, 0.0], time_values=[0.0, tc_eq.t_final])
	bc_north_ovb = momBC.DirichletBC(boundary_name="North_ovb", component=1, values=[0.0, 0.0], time_values=[0.0, tc_eq.t_final])

	# Extract geometry dimensions
	Lx = grid.Lx
	Ly = grid.Ly
	Lz = grid.Lz
	z_surface = 0.0

	g = 9.81
	ovb_thickness, salt_thickness, hanging_wall = get_geometry_parameters(grid_path)
	cavern_roof = ovb_thickness + hanging_wall
	p_roof = 0 + salt_density*g*hanging_wall + ovb_density*g*ovb_thickness

	# Pressure at the top of the salt layer (bottom of overburden)
	p_top = ovb_density*g*ovb_thickness

	bc_top = momBC.NeumannBC(boundary_name = "Top",
						direction = 2,
						density = 0.0,
						ref_pos = z_surface,
						values = [0*MPa, 0*MPa],
						time_values = [0*day,  10*day],
						g = g_vec[2])

	bc_cavern = momBC.NeumannBC(boundary_name = "Cavern",
						direction = 2,
						density = gas_density,
						ref_pos = cavern_roof,
						values = [0.8*p_roof, 0.8*p_roof],
						time_values = [0*day,  tc_eq.t_final],
						g = g_vec[2])

	bc_equilibrium = momBC.BcHandler(mom_eq)
	bc_equilibrium.add_boundary_condition(bc_west_salt)
	bc_equilibrium.add_boundary_condition(bc_west_ovb)
	bc_equilibrium.add_boundary_condition(bc_east_salt)
	bc_equilibrium.add_boundary_condition(bc_east_ovb)
	bc_equilibrium.add_boundary_condition(bc_bottom)
	bc_equilibrium.add_boundary_condition(bc_south_salt)
	bc_equilibrium.add_boundary_condition(bc_south_ovb)
	bc_equilibrium.add_boundary_condition(bc_north_salt)
	bc_equilibrium.add_boundary_condition(bc_north_ovb)
	bc_equilibrium.add_boundary_condition(bc_top)
	bc_equilibrium.add_boundary_condition(bc_cavern)

	# Set boundary conditions
	mom_eq.set_boundary_conditions(bc_equilibrium)

	# Equilibrium output folder
	ouput_folder_equilibrium = os.path.join(output_folder, "equilibrium")

	# Print output folder
	if MPI.COMM_WORLD.rank == 0:
		print(ouput_folder_equilibrium)

	# Create output handlers
	output_mom = sf.SaveFields(mom_eq)
	output_mom.set_output_folder(ouput_folder_equilibrium)
	output_mom.add_output_field("u", "Displacement (m)")
	output_mom.add_output_field("eps_tot", "Total strain (-)")
	output_mom.add_output_field("sig", "Stress (Pa)")
	output_mom.add_output_field("p_elems", "Mean stress (Pa)")
	output_mom.add_output_field("q_elems", "Von Mises stress (Pa)")
	outputs = [output_mom]

	# Define simulator
	sim = sf.Simulator_M(mom_eq, tc_eq, outputs, True)
	sim.run()






	# Time settings for operation stage
	tc_op = sf.TimeController(dt=0.5, initial_time=0.0, final_time=240, time_unit="day")

	# Define heat diffusion equation
	heat_eq = sf.HeatDiffusion(grid)

	# Define solver
	solver_heat = PETSc.KSP().create(grid.mesh.comm)
	solver_heat.setType("cg")
	solver_heat.getPC().setType("asm")
	solver_heat.setTolerances(rtol=1e-12, max_it=100)
	heat_eq.set_solver(solver_heat)

	# Set specific heat capacity
	cp = 850*to.ones(heat_eq.n_elems, dtype=to.float64)
	mat.set_specific_heat_capacity(cp)

	# Set thermal conductivity
	k = 7*to.ones(heat_eq.n_elems, dtype=to.float64)
	mat.set_thermal_conductivity(k)

	# Set material properties to heat_equation
	heat_eq.set_material(mat)

	# Set initial temperature
	T0_field_nodes = create_field_nodes(grid, T_field_fun)
	heat_eq.set_initial_T(T0_field_nodes)

	# Define boundary conditions for heat diffusion
	time_values = [tc_op.t_initial, tc_op.t_final]
	nt = len(time_values)

	bc_handler = heatBC.BcHandler(heat_eq)

	bc_top = heatBC.DirichletBC("Top", nt*[T_top], time_values)
	bc_bottom = heatBC.NeumannBC("Bottom", nt*[dTdZ], time_values)
	bc_east_salt = heatBC.NeumannBC("East_salt", nt*[0.0], time_values)
	bc_east_ovb = heatBC.NeumannBC("East_ovb", nt*[0.0], time_values)
	bc_west_salt = heatBC.NeumannBC("West_salt", nt*[0.0], time_values)
	bc_west_ovb = heatBC.NeumannBC("West_ovb", nt*[0.0], time_values)
	bc_south_salt = heatBC.NeumannBC("South_salt", nt*[0.0], time_values)
	bc_south_ovb = heatBC.NeumannBC("South_ovb", nt*[0.0], time_values)
	bc_north_salt = heatBC.NeumannBC("North_salt", nt*[0.0], time_values)
	bc_north_ovb = heatBC.NeumannBC("North_ovb", nt*[0.0], time_values)

	bc_handler.add_boundary_condition(bc_top)
	bc_handler.add_boundary_condition(bc_bottom)
	bc_handler.add_boundary_condition(bc_east_salt)
	bc_handler.add_boundary_condition(bc_east_ovb)
	bc_handler.add_boundary_condition(bc_west_salt)
	bc_handler.add_boundary_condition(bc_west_ovb)
	bc_handler.add_boundary_condition(bc_south_salt)
	bc_handler.add_boundary_condition(bc_south_ovb)
	bc_handler.add_boundary_condition(bc_north_salt)
	bc_handler.add_boundary_condition(bc_north_ovb)

	T_gas = T_top
	h_conv = 5.0
	bc_cavern = heatBC.RobinBC("Cavern", nt*[T_gas], h_conv, time_values)
	bc_handler.add_boundary_condition(bc_cavern)

	heat_eq.set_boundary_conditions(bc_handler)





	# Set operation stage settings for momentum equation

	# Boundary conditions
	bc_west_salt = momBC.DirichletBC(boundary_name="West_salt", component=0, values=[0.0, 0.0], time_values=[0.0, tc_op.t_final])
	bc_west_ovb = momBC.DirichletBC(boundary_name = "West_ovb", component=0, values=[0.0, 0.0], time_values=[0.0, tc_op.t_final])

	bc_east_salt = momBC.DirichletBC(boundary_name="East_salt", component=0, values=[0.0, 0.0], time_values=[0.0, tc_op.t_final])
	bc_east_ovb = momBC.DirichletBC(boundary_name = "East_ovb", component=0, values=[0.0, 0.0], time_values=[0.0, tc_op.t_final])

	bc_bottom = momBC.DirichletBC(boundary_name="Bottom", component=2, values=[0.0, 0.0], time_values=[0.0, tc_op.t_final])

	bc_south_salt = momBC.DirichletBC(boundary_name="South_salt", component=1, values=[0.0, 0.0], time_values=[0.0, tc_op.t_final])
	bc_south_ovb = momBC.DirichletBC(boundary_name="South_ovb", component=1, values=[0.0, 0.0], time_values=[0.0, tc_op.t_final])

	bc_north_salt = momBC.DirichletBC(boundary_name="North_salt", component=1, values=[0.0, 0.0], time_values=[0.0, tc_op.t_final])
	bc_north_ovb = momBC.DirichletBC(boundary_name="North_ovb", component=1, values=[0.0, 0.0], time_values=[0.0, tc_op.t_final])

	bc_top = momBC.NeumannBC(boundary_name="Top", direction=2, density=0.0, ref_pos=z_surface, values=[0, 0], time_values=[0, tc_op.t_final], g=g_vec[2])

	p_values = 3*[0.8*p_roof, 0.8*p_roof, 0.2*p_roof, 0.2*p_roof] + [0.8*p_roof]
	t_values = [20*day*i for i in range(13)]

	bc_cavern = momBC.NeumannBC(boundary_name = "Cavern",
						direction = 2,
						density = gas_density,
						ref_pos = cavern_roof,
						values = p_values,
						time_values = t_values,
						g = g_vec[2])


	bc_operation = momBC.BcHandler(mom_eq)
	bc_operation.add_boundary_condition(bc_west_salt)
	bc_operation.add_boundary_condition(bc_west_ovb)
	bc_operation.add_boundary_condition(bc_east_salt)
	bc_operation.add_boundary_condition(bc_east_ovb)
	bc_operation.add_boundary_condition(bc_bottom)
	bc_operation.add_boundary_condition(bc_south_salt)
	bc_operation.add_boundary_condition(bc_south_ovb)
	bc_operation.add_boundary_condition(bc_north_salt)
	bc_operation.add_boundary_condition(bc_north_ovb)
	bc_operation.add_boundary_condition(bc_top)
	bc_operation.add_boundary_condition(bc_cavern)

	# Set boundary conditions
	mom_eq.set_boundary_conditions(bc_operation)

	# Define output folder
	output_folder_operation = os.path.join(output_folder, "operation")

	# Print output folder
	if MPI.COMM_WORLD.rank == 0:
		print(output_folder_operation)
		sys.stdout.flush()

	# Create output handlers
	output_mom = sf.SaveFields(mom_eq)
	output_mom.set_output_folder(output_folder_operation)
	output_mom.add_output_field("u", "Displacement (m)")
	output_mom.add_output_field("p_elems", "Mean stress (Pa)")
	output_mom.add_output_field("q_elems", "Von Mises stress (Pa)")

	output_heat = sf.SaveFields(heat_eq)
	output_heat.set_output_folder(output_folder_operation)
	output_heat.add_output_field("T", "Temperature (K)")

	outputs = [output_mom, output_heat]

	# Define simulator
	sim = sf.Simulator_TM(mom_eq, heat_eq, tc_op, outputs, False)
	sim.run()


if __name__ == '__main__':
	main()