from safeincave import *
import safeincave.Utils as ut
import safeincave.MomentumBC as momBC
from mpi4py import MPI
from petsc4py import PETSc
import torch as to
import os
import time


def main():
	# Read grid
	grid_path = os.path.join("..", "..", "..", "grids", "cube")
	grid = GridHandlerGMSH("geom", grid_path)

	# Time settings for equilibrium stage
	unit = "hour"
	t_0 = 0.0
	dt = 0.01
	t_final = 1
	t_control = TimeController(dt=dt, initial_time=t_0, final_time=t_final, time_unit=unit)

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
	E0[omega_A] = 8*ut.GPa
	E0[omega_B] = 10*ut.GPa
	nu0[omega_A] = 0.2
	nu0[omega_B] = 0.3
	spring_0 = Spring(E0, nu0, "spring")

	# Create Kelvin-Voigt viscoelastic element
	eta = to.zeros(mom_eq.n_elems)
	E1 = to.zeros(mom_eq.n_elems)
	nu1 = to.zeros(mom_eq.n_elems)
	eta[omega_A] = 105e11
	eta[omega_B] = 38e11
	E1[omega_A] = 8*ut.GPa
	E1[omega_B] = 5*ut.GPa
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
	time_values = [0*ut.hour,  1*ut.hour]
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
						values =      [5.0*ut.MPa, 5.0*ut.MPa],
						time_values = [0.0,           t_control.t_final],
						g = g_vec[2])

	bc_north = momBC.NeumannBC(boundary_name = "NORTH",
						direction = 2,
						density = 0.0,
						ref_pos = 0.0,
						values =      [5.0*ut.MPa, 5.0*ut.MPa],
						time_values = [0.0,           t_control.t_final],
						g = g_vec[2])

	bc_top = momBC.NeumannBC(boundary_name = "TOP",
						direction = 2,
						density = 0.0,
						ref_pos = 0.0,
						values =      [8.0*ut.MPa, 8.0*ut.MPa],
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



if __name__ == '__main__':
	main()