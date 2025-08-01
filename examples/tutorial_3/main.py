import os
import sys
sys.path.append(os.path.join("..", "..", "safeincave"))
from Grid2 import GridHandlerGMSH
from mpi4py import MPI
import ufl
import dolfinx as do
import torch as to
import numpy as np
from petsc4py import PETSc
import Utils as utils
from MaterialProps import *
from MomentumEquation import LinearMomentum
import MomentumBC as momBC
from OutputHandler import SaveFields
from Simulators import Simulator_M
from TimeHandler import TimeController, TimeControllerParabolic
import time

GPa = utils.GPa
day = utils.day

def get_geometry_parameters(path_to_grid):
	f = open(os.path.join(path_to_grid, "geom.geo"), "r")
	data = f.readlines()
	ovb_thickness = float(data[10][len("ovb_thickness = "):-2])
	salt_thickness = float(data[11][len("salt_thickness = "):-2])
	hanging_wall = float(data[12][len("hanging_wall = "):-2])
	return ovb_thickness, salt_thickness, hanging_wall





def main():
	comm = MPI.COMM_WORLD
	comm.Barrier()
	if MPI.COMM_WORLD.rank == 0:
	    start_time = MPI.Wtime()

	# Read grid
	grid_path = os.path.join("..", "..", "grids", "cavern_overburden_coarse")
	# grid_path = os.path.join("..", "..", "grids", "cavern_overburden")
	grid = GridHandlerGMSH("geom", grid_path)

	# Define output folder
	output_folder = os.path.join("output", "case_0")

	# Define momentum equation
	mom_eq = LinearMomentum(grid, theta=0.0)

	# Define solver
	mom_solver = PETSc.KSP().create(grid.mesh.comm)
	mom_solver.setType("bicg")
	mom_solver.getPC().setType("asm")
	mom_solver.setTolerances(rtol=1e-12, max_it=100)
	mom_eq.set_solver(mom_solver)

	# Define material properties
	mat = Material(mom_eq.n_elems)

	# Extract region indices
	ind_salt = grid.region_indices["Salt"]
	ind_ovb = grid.region_indices["Overburden"]

	# Set material density
	salt_density = 2200
	ovb_density = 2800
	gas_density = 10
	rho = to.zeros(mom_eq.n_elems, dtype=to.float64)
	rho[ind_salt] = salt_density
	rho[ind_ovb] = ovb_density
	mat.set_density(rho)

	# Constitutive model
	E0 = to.zeros(mom_eq.n_elems)
	E0[ind_salt] = 102*GPa
	E0[ind_ovb] = 180*GPa
	nu0 = 0.3*to.ones(mom_eq.n_elems)
	spring_0 = Spring(E0, nu0, "spring")

	# Create Kelvin-Voigt viscoelastic element
	eta = 105e11*to.ones(mom_eq.n_elems)
	E1 = 10*utils.GPa*to.ones(mom_eq.n_elems)
	nu1 = 0.32*to.ones(mom_eq.n_elems)
	kelvin = Viscoelastic(eta, E1, nu1, "kelvin")

	# Create creep
	A = to.zeros(mom_eq.n_elems)
	A[ind_salt] = 1.9e-20
	A[ind_ovb] = 0.0
	Q = 51600*to.ones(mom_eq.n_elems)
	n = 3.0*to.ones(mom_eq.n_elems)
	creep_0 = DislocationCreep(A, Q, n, "creep")

	# Create constitutive model
	mat.add_to_elastic(spring_0)
	# mat.add_to_non_elastic(kelvin)
	mat.add_to_non_elastic(creep_0)

	# Set constitutive model
	mom_eq.set_material(mat)

	# Set body forces
	g = -9.81
	g_vec = [0.0, 0.0, g]
	mom_eq.build_body_force(g_vec)

	# Set initial temperature field
	def T_field_fun(x,y,z):
		km = 1000
		dTdZ = 27/km
		T_surface = 20 + 273
		return T_surface - dTdZ*z
	T0_field = utils.create_field_elems(grid, T_field_fun)
	mom_eq.set_T0(T0_field)
	mom_eq.set_T(T0_field)

	# Time settings for equilibrium stage
	tc_eq = TimeControllerParabolic(final_time=5, initial_time=0.0, n_time_steps=50, time_unit="day")
	# tc_eq = TimeController(time_step=0.1, final_time=5, initial_time=0.0, time_unit="day")

	# Boundary conditions
	time_values = [0*utils.hour,  1*utils.hour]
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
						time_values = [0*day,  10*day],
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
	output_mom = SaveFields(mom_eq)
	output_mom.set_output_folder(ouput_folder_equilibrium)
	output_mom.add_output_field("u", "Displacement (m)")
	# output_mom.add_output_field("Temp", "Temperature (K)")
	output_mom.add_output_field("eps_tot", "Total strain (-)")
	output_mom.add_output_field("p_elems", "Mean stress (Pa)")
	output_mom.add_output_field("q_elems", "Von Mises stress (Pa)")
	output_mom.add_output_field("p_nodes", "Mean stress (Pa)")
	output_mom.add_output_field("q_nodes", "Von Mises stress (Pa)")
	outputs = [output_mom]

	# Define simulator
	sim = Simulator_M(mom_eq, tc_eq, outputs, True)
	sim.run()

	# Print time
	if MPI.COMM_WORLD.rank == 0:
		end_time = MPI.Wtime()
		elaspsed_time = end_time - start_time
		formatted_time = time.strftime("%H:%M:%S", time.gmtime(elaspsed_time))
		print(f"Time: {formatted_time} ({elaspsed_time} seconds)\n")








	# Time settings for equilibrium stage
	tc_op = TimeController(time_step=2, final_time=240, initial_time=0.0, time_unit="hour")

	# # Boundary conditions
	bc_west_salt = momBC.DirichletBC(boundary_name="West_salt", component=0, values=[0.0, 0.0], time_values=[0.0, tc_op.t_final])
	bc_west_ovb = momBC.DirichletBC(boundary_name = "West_ovb", component=0, values=[0.0, 0.0], time_values=[0.0, tc_op.t_final])

	bc_east_salt = momBC.DirichletBC(boundary_name="East_salt", component=0, values=[0.0, 0.0], time_values=[0.0, tc_op.t_final])
	bc_east_ovb = momBC.DirichletBC(boundary_name = "East_ovb", component=0, values=[0.0, 0.0], time_values=[0.0, tc_op.t_final])

	bc_bottom = momBC.DirichletBC(boundary_name="Bottom", component=2, values=[0.0, 0.0], time_values=[0.0, tc_op.t_final])

	bc_south_salt = momBC.DirichletBC(boundary_name="South_salt", component=1, values=[0.0, 0.0], time_values=[0.0, tc_op.t_final])
	bc_south_ovb = momBC.DirichletBC(boundary_name="South_ovb", component=1, values=[0.0, 0.0], time_values=[0.0, tc_op.t_final])

	bc_north_salt = momBC.DirichletBC(boundary_name="North_salt", component=1, values=[0.0, 0.0], time_values=[0.0, tc_op.t_final])
	bc_north_ovb = momBC.DirichletBC(boundary_name="North_ovb", component=1, values=[0.0, 0.0], time_values=[0.0, tc_op.t_final])

	bc_cavern = momBC.NeumannBC(boundary_name = "Cavern",
						direction = 2,
						density = gas_density,
						ref_pos = cavern_roof,
						values = [0.8*p_roof, 0.2*p_roof, 0.2*p_roof, 0.8*p_roof, 0.8*p_roof],
						time_values = [0*day,  2*day,  6*day, 8*day, 10*day],
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

	# Create output handlers
	output_mom = SaveFields(mom_eq)
	output_mom.set_output_folder(output_folder_operation)
	output_mom.add_output_field("u", "Displacement (m)")
	# output_mom.add_output_field("Temp", "Temperature (K)")
	output_mom.add_output_field("eps_tot", "Total strain (-)")
	output_mom.add_output_field("p_elems", "Mean stress (Pa)")
	output_mom.add_output_field("q_elems", "Von Mises stress (Pa)")
	output_mom.add_output_field("p_nodes", "Mean stress (Pa)")
	output_mom.add_output_field("q_nodes", "Von Mises stress (Pa)")
	outputs = [output_mom]

	# Define simulator
	sim = Simulator_M(mom_eq, tc_op, outputs, False)
	sim.run()

	# Print time
	if MPI.COMM_WORLD.rank == 0:
		end_time = MPI.Wtime()
		elaspsed_time = end_time - start_time
		formatted_time = time.strftime("%H:%M:%S", time.gmtime(elaspsed_time))
		print(f"Time: {formatted_time} ({elaspsed_time} seconds)\n")


if __name__ == '__main__':
	main()