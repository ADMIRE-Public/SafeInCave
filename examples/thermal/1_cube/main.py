import safeincave as sf
import safeincave.Utils as ut
import safeincave.HeatBC as heatBC
from petsc4py import PETSc
from mpi4py import MPI
import torch as to
import os
import sys


def main():
	# Read grid
	grid_path = os.path.join("..", "..", "..", "grids", "cube")
	grid = sf.GridHandlerGMSH("geom", grid_path)

	# Define output folder
	output_folder = os.path.join("output", "case_0")

	# Time settings for equilibrium stage
	t_control = sf.TimeControllerParabolic(n_time_steps=50, initial_time=0.0, final_time=5, time_unit="day")

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
							values = [273, 273],
							time_values = [t_control.t_initial, t_control.t_final])

	bc_west = heatBC.RobinBC(boundary_name = "WEST", 
							values = [273, 273],
							h = 5.0,
							time_values = [t_control.t_initial, t_control.t_final])

	bc_handler = heatBC.BcHandler(heat_eq)
	bc_handler.add_boundary_condition(bc_east)
	bc_handler.add_boundary_condition(bc_west)

	# Add boundary condition to heat equation
	heat_eq.set_boundary_conditions(bc_handler)

	# Set initial temperature field
	fun = lambda x, y, z: 293
	T0_field = ut.create_field_nodes(heat_eq.grid, fun)
	heat_eq.set_initial_T(T0_field)

	# Create output handlers
	output_heat = sf.SaveFields(heat_eq)
	output_heat.set_output_folder(output_folder)
	output_heat.add_output_field("T", "Temperature (K)")
	outputs = [output_heat]

	# Print output folder
	if MPI.COMM_WORLD.rank == 0:
		print(output_folder)
		sys.stdout.flush()

	# Define simulator
	sim = sf.Simulator_T(heat_eq, t_control, outputs, True)
	sim.run()



if __name__ == '__main__':
	main()