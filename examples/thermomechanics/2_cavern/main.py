import safeincave as sf
import safeincave.Utils as ut
import safeincave.HeatBC as heatBC
import safeincave.MomentumBC as momBC
import os
import sys
import torch as to
from petsc4py import PETSc


def main():
	# Read grid
	grid_path = os.path.join("..", "..", "..", "grids", "cavern_irregular")
	grid = sf.GridHandlerGMSH("geom", grid_path)

	# Define output folder
	output_folder = os.path.join("output", "case_0")

	# Define momentum equation
	mom_eq = sf.LinearMomentum(grid, theta=0.5)

	# Define solver
	mom_solver = PETSc.KSP().create(grid.mesh.comm)
	mom_solver.setType("bicg")
	mom_solver.getPC().setType("asm")
	mom_solver.setTolerances(rtol=1e-12, max_it=100)
	mom_eq.set_solver(mom_solver)

	# Define material properties
	mat = sf.Material(mom_eq.n_elems)

	# Set material density
	salt_density = 2000
	rho = salt_density*to.ones(mom_eq.n_elems, dtype=to.float64)
	mat.set_density(rho)

	# Constitutive model
	E0 = 102*ut.GPa*to.ones(mom_eq.n_elems)
	nu0 = 0.3*to.ones(mom_eq.n_elems)
	spring_0 = sf.Spring(E0, nu0, "spring")

	# Create Kelvin-Voigt viscoelastic element
	eta = 105e11*to.ones(mom_eq.n_elems)
	E1 = 10*ut.GPa*to.ones(mom_eq.n_elems)
	nu1 = 0.32*to.ones(mom_eq.n_elems)
	kelvin = sf.Viscoelastic(eta, E1, nu1, "kelvin")

	# Create creep
	A = 1.9e-20*to.ones(mom_eq.n_elems)
	Q = 51600*to.ones(mom_eq.n_elems)
	n = 3.0*to.ones(mom_eq.n_elems)
	creep_0 = sf.DislocationCreep(A, Q, n, "creep")

	# Thermo-elastic element
	alpha = 44e-6*to.ones(mom_eq.n_elems)
	thermo = sf.Thermoelastic(alpha, "thermo")

	# Create constitutive model
	mat.add_to_elastic(spring_0)
	mat.add_to_thermoelastic(thermo)
	mat.add_to_non_elastic(kelvin)
	mat.add_to_non_elastic(creep_0)

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
	T0_field_elems = ut.create_field_elems(grid, T_field_fun)
	mom_eq.set_T0(T0_field_elems)
	mom_eq.set_T(T0_field_elems)

	# Time settings for equilibrium stage
	tc_equilibrium = sf.TimeController(dt=0.5, initial_time=0.0, final_time=10, time_unit="hour")

	# Boundary conditions
	time_values = [0*ut.hour,  1*ut.hour]
	nt = len(time_values)

	bc_west = momBC.DirichletBC(boundary_name = "West", 
					 		component = 0,
							values = [0.0, 0.0],
							time_values = [0.0, tc_equilibrium.t_final])

	bc_bottom = momBC.DirichletBC(boundary_name = "Bottom", 
					 	  component = 2,
					 	  values = [0.0, 0.0],
					 	  time_values = [0.0, tc_equilibrium.t_final])

	bc_south = momBC.DirichletBC(boundary_name = "South", 
					 	  component = 1,
					 	  values = [0.0, 0.0],
					 	  time_values = [0.0, tc_equilibrium.t_final])

	side_burden = 10.0*ut.MPa
	bc_east = momBC.NeumannBC(boundary_name = "East",
						direction = 2,
						density = salt_density,
						ref_pos = 660.0,
						values =      [side_burden, side_burden],
						time_values = [0.0, tc_equilibrium.t_final],
						g = g_vec[2])

	bc_north = momBC.NeumannBC(boundary_name = "North",
						direction = 2,
						density = salt_density,
						ref_pos = 660.0,
						values =      [side_burden, side_burden],
						time_values = [0.0, tc_equilibrium.t_final],
						g = g_vec[2])

	over_burden = 10.0*ut.MPa
	bc_top = momBC.NeumannBC(boundary_name = "Top",
						direction = 2,
						density = 0.0,
						ref_pos = 0.0,
						values =      [over_burden, over_burden],
						time_values = [0.0, tc_equilibrium.t_final],
						g = g_vec[2])

	gas_density = 0.082
	p_gas = 10.0*ut.MPa
	bc_cavern = momBC.NeumannBC(boundary_name = "Cavern",
						direction = 2,
						density = gas_density,
						# density = 0.082,
						ref_pos = 430.0,
						values =      [p_gas, p_gas],
						time_values = [0.0,            tc_equilibrium.t_final],
						g = g_vec[2])

	bc_equilibrium = momBC.BcHandler(mom_eq)
	bc_equilibrium.add_boundary_condition(bc_west)
	bc_equilibrium.add_boundary_condition(bc_bottom)
	bc_equilibrium.add_boundary_condition(bc_south)
	bc_equilibrium.add_boundary_condition(bc_east)
	bc_equilibrium.add_boundary_condition(bc_north)
	bc_equilibrium.add_boundary_condition(bc_top)
	bc_equilibrium.add_boundary_condition(bc_cavern)

	# Set boundary conditions
	mom_eq.set_boundary_conditions(bc_equilibrium)

	# Equilibrium output folder
	ouput_folder_equilibrium = os.path.join(output_folder, "equilibrium")

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
	sim = sf.Simulator_M(mom_eq, tc_equilibrium, outputs, True)
	sim.run()






	# Time settings for operation stage
	tc_operation = sf.TimeController(dt=1, initial_time=0.0, final_time=240, time_unit="day")

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
	T0_field_nodes = ut.create_field_nodes(grid, T_field_fun)
	heat_eq.set_initial_T(T0_field_nodes)

	# Define boundary conditions for heat diffusion
	time_values = [tc_operation.t_initial, tc_operation.t_final]
	nt = len(time_values)

	bc_handler = heatBC.BcHandler(heat_eq)

	bc_top = heatBC.DirichletBC("Top", nt*[T_top], time_values)
	bc_handler.add_boundary_condition(bc_top)

	bc_bottom = heatBC.NeumannBC("Bottom", nt*[dTdZ], time_values)
	bc_handler.add_boundary_condition(bc_bottom)

	bc_east = heatBC.NeumannBC("East", nt*[0.0], time_values)
	bc_handler.add_boundary_condition(bc_east)

	bc_west = heatBC.NeumannBC("West", nt*[0.0], time_values)
	bc_handler.add_boundary_condition(bc_west)

	bc_south = heatBC.NeumannBC("South", nt*[0.0], time_values)
	bc_handler.add_boundary_condition(bc_south)

	bc_north = heatBC.NeumannBC("North", nt*[0.0], time_values)
	bc_handler.add_boundary_condition(bc_north)

	T_gas = T_top
	h_conv = 5.0
	bc_cavern = heatBC.RobinBC("Cavern", nt*[T_gas], h_conv, time_values)
	bc_handler.add_boundary_condition(bc_cavern)

	heat_eq.set_boundary_conditions(bc_handler)





	# Set operation stage settings for momentum equation

	# Boundary conditions
	bc_west = momBC.DirichletBC("West", 0, [0.0, 0.0], [0.0, tc_operation.t_final])
	bc_bottom = momBC.DirichletBC("Bottom", 2, [0.0, 0.0], [0.0, tc_operation.t_final])
	bc_south = momBC.DirichletBC("South", 1, [0.0, 0.0], [0.0, tc_operation.t_final])
	bc_east = momBC.NeumannBC("East", 2, salt_density, 660.0, [side_burden, side_burden], [0.0, tc_operation.t_final], g_vec[2])
	bc_north = momBC.NeumannBC("North", 2, salt_density, 660.0, [side_burden, side_burden], [0.0, tc_operation.t_final], g_vec[2])
	bc_top = momBC.NeumannBC("Top", 2, 0.0, 0.0, [over_burden, over_burden], [0.0, tc_operation.t_final], g_vec[2])

	time_list = [0.0, 2.0*ut.hour, 14*ut.hour, 16*ut.hour, 24*ut.hour]
	p_list = [10.0*ut.MPa, 7.0*ut.MPa, 7.0*ut.MPa, 10.0*ut.MPa, 10.0*ut.MPa]
	bc_cavern = momBC.NeumannBC("Cavern", 2, gas_density, 430.0, p_list, time_list, g_vec[2])

	bc_operation = momBC.BcHandler(mom_eq)
	bc_operation.add_boundary_condition(bc_west)
	bc_operation.add_boundary_condition(bc_bottom)
	bc_operation.add_boundary_condition(bc_south)
	bc_operation.add_boundary_condition(bc_east)
	bc_operation.add_boundary_condition(bc_north)
	bc_operation.add_boundary_condition(bc_top)
	bc_operation.add_boundary_condition(bc_cavern)

	# Set boundary conditions
	mom_eq.set_boundary_conditions(bc_operation)

	# Define output folder
	output_folder_operation = os.path.join(output_folder, "operation")

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
	sim = sf.Simulator_TM(mom_eq, heat_eq, tc_operation, outputs, False)
	sim.run()


if __name__ == '__main__':
	main()