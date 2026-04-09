import safeincave as sf
import safeincave.Utils as ut
import safeincave.HeatBC as heatBC
from petsc4py import PETSc
import torch as to
import os

def main():
	# Read grid
	grid_path = os.path.join("..", "..", "..", "grids", "cavern_regular")
	grid = sf.GridHandlerGMSH("geom", grid_path)

	# Define output folder
	output_folder = os.path.join("output", "case_0")

	# Time settings for equilibrium stage
	t_control = sf.TimeControllerParabolic(n_time_steps=100, initial_time=0, final_time=5, time_unit="year")

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

	km = 1000
	dTdZ = 27/km
	T_top = 273 + 20
	T_gas = 273 + 10
	h_conv = 5.0

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

	bc_cavern = heatBC.RobinBC("Cavern", nt*[T_gas], h_conv, time_values)
	bc_handler.add_boundary_condition(bc_cavern)

	heat_eq.set_boundary_conditions(bc_handler)

	# Set initial temperature field
	fun = lambda x, y, z: T_top - dTdZ*(z - 660)
	T0_field = ut.create_field_nodes(heat_eq.grid, fun)
	heat_eq.set_initial_T(T0_field)

	# Create output handlers
	output_heat = sf.SaveFields(heat_eq)
	output_heat.set_output_folder(output_folder)
	output_heat.add_output_field("T", "Temperature (K)")
	outputs = [output_heat]

	# Define simulator
	sim = sf.Simulator_T(heat_eq, t_control, outputs, True)
	sim.run()


if __name__ == '__main__':
	main()