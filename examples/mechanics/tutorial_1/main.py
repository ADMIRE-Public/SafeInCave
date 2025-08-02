import os
import sys
sys.path.append(os.path.join("..", "..", "..", "safeincave"))
from Grid import GridHandlerGMSH, GridHandlerFEniCS
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
from TimeHandler import TimeController
import time



def main():
	comm = MPI.COMM_WORLD
	comm.Barrier()
	if MPI.COMM_WORLD.rank == 0:
	    start_time = MPI.Wtime()

	# Read grid
	grid_path = os.path.join("..", "..", "grids", "cube")
	grid = GridHandlerGMSH("geom", grid_path)

	# Time settings for equilibrium stage
	unit = "hour"
	t_0 = 0.0
	dt = 0.01
	t_final = 1
	t_control = TimeController(time_step=dt, final_time=t_final, initial_time=t_0, time_unit=unit)

	# Define momentum equation
	mom_eq = LinearMomentum(grid, theta=0.5)

	# Define solver
	mom_solver = PETSc.KSP().create(grid.mesh.comm)
	mom_solver.setType("bicg")
	mom_solver.getPC().setType("asm")
	mom_solver.setTolerances(rtol=1e-12, max_it=100)
	mom_eq.set_solver(mom_solver)

	# Define material properties
	mat = Material(mom_eq.n_elems)

	# Set material density
	rho = 0.0*to.ones(mom_eq.n_elems, dtype=to.float64)
	mat.set_density(rho)

	# Extract region indices
	omega_A = grid.region_indices["OMEGA_A"]
	omega_B = grid.region_indices["OMEGA_B"]

	# Constitutive model
	E0 = to.zeros(mom_eq.n_elems)
	nu0 = to.zeros(mom_eq.n_elems)
	E0[omega_A] = 8*utils.GPa
	E0[omega_B] = 10*utils.GPa
	nu0[omega_A] = 0.2
	nu0[omega_B] = 0.3
	spring_0 = Spring(E0, nu0, "spring")

	# Create Kelvin-Voigt viscoelastic element
	eta = to.zeros(mom_eq.n_elems)
	E1 = to.zeros(mom_eq.n_elems)
	nu1 = to.zeros(mom_eq.n_elems)
	eta[omega_A] = 105e11
	eta[omega_B] = 38e11
	E1[omega_A] = 8*utils.GPa
	E1[omega_B] = 5*utils.GPa
	nu1[omega_A] = 0.35
	nu1[omega_B] = 0.28
	kelvin = Viscoelastic(eta, E1, nu1, "kelvin")

	# Create constitutive model
	mat.add_to_elastic(spring_0)
	mat.add_to_non_elastic(kelvin)

	# Set constitutive model
	mom_eq.set_material(mat)

	# Set body forces
	g = -9.81
	g_vec = [0.0, 0.0, g]
	mom_eq.build_body_force(g_vec)

	# Set initial temperature field
	T0_field = 298*to.ones(mom_eq.n_elems)
	mom_eq.set_T0(T0_field)
	mom_eq.set_T(T0_field)

	# Boundary conditions
	time_values = [0*utils.hour,  1*utils.hour]
	nt = len(time_values)

	bc_west = momBC.DirichletBC(boundary_name = "WEST", 
					 		component = 0,
							values = [0.0, 0.0],
							time_values = [0.0, t_control.t_final])

	bc_bottom = momBC.DirichletBC(boundary_name = "BOTTOM", 
					 	  component = 2,
					 	  values = [0.0, 0.0],
					 	  time_values = [0.0, t_control.t_final])

	bc_south = momBC.DirichletBC(boundary_name = "SOUTH", 
					 	  component = 1,
					 	  values = [0.0, 0.0],
					 	  time_values = [0.0, t_control.t_final])

	bc_east = momBC.NeumannBC(boundary_name = "EAST",
						direction = 2,
						density = 0.0,
						ref_pos = 0.0,
						values =      [5.0*utils.MPa, 5.0*utils.MPa],
						time_values = [0.0,           t_control.t_final],
						g = g_vec[2])

	bc_north = momBC.NeumannBC(boundary_name = "NORTH",
						direction = 2,
						density = 0.0,
						ref_pos = 0.0,
						values =      [5.0*utils.MPa, 5.0*utils.MPa],
						time_values = [0.0,           t_control.t_final],
						g = g_vec[2])

	bc_top = momBC.NeumannBC(boundary_name = "TOP",
						direction = 2,
						density = 0.0,
						ref_pos = 0.0,
						values =      [8.0*utils.MPa, 8.0*utils.MPa],
						time_values = [0.0,           t_control.t_final],
						g = g_vec[2])

	bc_handler = momBC.BcHandler(mom_eq)
	bc_handler.add_boundary_condition(bc_west)
	bc_handler.add_boundary_condition(bc_bottom)
	bc_handler.add_boundary_condition(bc_south)
	bc_handler.add_boundary_condition(bc_east)
	bc_handler.add_boundary_condition(bc_north)
	bc_handler.add_boundary_condition(bc_top)

	# Set boundary conditions
	mom_eq.set_boundary_conditions(bc_handler)

	# Define output folder
	output_folder = os.path.join("output", "case_0")

	# Create output handlers
	output_mom = SaveFields(mom_eq)
	output_mom.set_output_folder(output_folder)
	output_mom.add_output_field("u", "Displacement (m)")
	output_mom.add_output_field("eps_tot", "Total strain (-)")
	output_mom.add_output_field("sig", "Stress (Pa)")
	output_mom.add_output_field("p_elems", "Mean stress (Pa)")
	output_mom.add_output_field("q_elems", "Von Mises stress (Pa)")
	outputs = [output_mom]

	# Print output folder
	if MPI.COMM_WORLD.rank == 0:
		print(output_folder)

	# Define simulator
	sim = Simulator_M(mom_eq, t_control, outputs, True)
	sim.run()

	# Print time
	if MPI.COMM_WORLD.rank == 0:
		end_time = MPI.Wtime()
		elaspsed_time = end_time - start_time
		formatted_time = time.strftime("%H:%M:%S", time.gmtime(elaspsed_time))
		print(f"Time: {formatted_time} ({elaspsed_time} seconds)\n")



if __name__ == '__main__':
	main()