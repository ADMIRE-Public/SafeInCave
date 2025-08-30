import safeincave as sf
import safeincave.Utils as ut
import safeincave.HeatBC as heatBC
import safeincave.MomentumBC as momBC
from mpi4py import MPI
import dolfinx as do
import os
import sys
import ufl
import torch as to
import numpy as np
from petsc4py import PETSc
import time



def main():
	# Read grid
	grid_path = os.path.join("..", "..", "..", "grids", "cube")
	grid = sf.GridHandlerGMSH("geom", grid_path)

	# Define output folder
	output_folder = os.path.join("output", "case_0")

	# Time settings for equilibrium stage
	t_control = sf.TimeControllerParabolic(n_time_steps=100, initial_time=0.0, final_time=10, time_unit="day")

	# Define equation
	heat_eq = sf.HeatDiffusion(grid)

	# Define solver
	solver_heat = PETSc.KSP().create(grid.mesh.comm)
	solver_heat.setType("cg")
	solver_heat.getPC().setType("asm")
	solver_heat.setTolerances(rtol=1e-12, max_it=100)
	heat_eq.set_solver(solver_heat)

	# Build material properties
	mat = sf.Material(heat_eq.n_elems)

	# Set material density
	rho = 2000.0*to.ones(heat_eq.n_elems, dtype=to.float64)
	mat.set_density(rho)

	# Set specific heat capacity
	cp = 850*to.ones(heat_eq.n_elems, dtype=to.float64)
	mat.set_specific_heat_capacity(cp)

	# Set thermal conductivity
	k = 7*to.ones(heat_eq.n_elems, dtype=to.float64)
	mat.set_thermal_conductivity(k)

	# Set material properties to heat_equation
	heat_eq.set_material(mat)

	# Define boundary conditions for heat diffusion
	time_values = [t_control.t_initial, t_control.t_final]
	nt = len(time_values)

	bc_east = heatBC.DirichletBC(boundary_name = "EAST", 
							values = nt*[273 + 1.0],
							time_values = time_values)

	bc_handler = heatBC.BcHandler(heat_eq)
	bc_handler.add_boundary_condition(bc_east)
	heat_eq.set_boundary_conditions(bc_handler)

	# Set initial temperature field
	fun = lambda x, y, z: 273 + 20
	T0_field = ut.create_field_nodes(heat_eq.grid, fun)
	heat_eq.set_initial_T(T0_field)



	# Define momentum equation
	mom_eq = sf.LinearMomentum(grid, 0.5)

	# Define solver
	mom_solver = PETSc.KSP().create(grid.mesh.comm)
	mom_solver.setType("bicg")
	mom_solver.getPC().setType("asm")
	mom_solver.setTolerances(rtol=1e-12, max_it=100)
	mom_eq.set_solver(mom_solver)

	# Constitutive model
	E = 102e9*to.ones(mom_eq.n_elems)
	nu = 0.3*to.ones(mom_eq.n_elems)
	spring_0 = sf.Spring(E, nu, "spring")

	# Extract region indices
	omega_A = grid.region_indices["OMEGA_A"]
	omega_B = grid.region_indices["OMEGA_B"]

	# Thermo-elastic element
	alpha = to.zeros(mom_eq.n_elems)
	alpha[omega_A] = 44e-6
	alpha[omega_B] = 74e-6
	# alpha = 44e-6*to.ones(mom_eq.n_elems)
	thermo = sf.Thermoelastic(alpha, "thermo")

	# Create constitutive model
	mat.add_to_elastic(spring_0)
	mat.add_to_thermoelastic(thermo)

	# Set constitutive model
	mom_eq.set_material(mat)

	# Set body forces
	g = -9.81
	g_vec = [0.0, 0.0, g]
	mom_eq.build_body_force(g_vec)

	# Boundary conditions
	bc_west_2 = momBC.DirichletBC(boundary_name = "WEST", 
					 		component = 2,
							values = nt*[0.0],
							time_values = time_values)

	bc_west_1 = momBC.DirichletBC(boundary_name = "WEST", 
					 		component = 1,
							values = nt*[0.0],
							time_values = time_values)

	bc_west_0 = momBC.DirichletBC(boundary_name = "WEST", 
					 	  component = 0,
					 	  values = nt*[0.0],
					 	  time_values = time_values)

	bc_bottom = momBC.DirichletBC(boundary_name = "BOTTOM", 
					 	  component = 2,
					 	  values = nt*[0.0],
					 	  time_values = time_values)

	bc_handler = momBC.BcHandler(mom_eq)
	bc_handler.add_boundary_condition(bc_west_0)
	bc_handler.add_boundary_condition(bc_west_1)
	bc_handler.add_boundary_condition(bc_west_2)
	bc_handler.add_boundary_condition(bc_bottom)

	# Set boundary conditions
	mom_eq.set_boundary_conditions(bc_handler)

	# Create output handlers
	output_mom = sf.SaveFields(mom_eq)
	output_mom.set_output_folder(output_folder)
	output_mom.add_output_field("u", "Displacement (m)")
	output_mom.add_output_field("sig", "Stress (Pa)")
	output_mom.add_output_field("p_nodes", "Mean stress (Pa)")
	output_mom.add_output_field("p_elems", "Mean stress (Pa)")
	output_mom.add_output_field("q_nodes", "Von Mises stress (Pa)")
	output_mom.add_output_field("q_elems", "Von Mises stress (Pa)")

	output_heat = sf.SaveFields(heat_eq)
	output_heat.set_output_folder(output_folder)
	output_heat.add_output_field("T", "Temperature (K)")

	outputs = [output_mom, output_heat]

	# Print output folder
	if MPI.COMM_WORLD.rank == 0:
		print(output_folder)
		sys.stdout.flush()

	# Define simulator
	sim = sf.Simulator_TM(mom_eq, heat_eq, t_control, outputs, True)
	sim.run()



if __name__ == '__main__':
	main()