import safeincave as sf
import safeincave.MassBC as massBC
from petsc4py import PETSc
from mpi4py import MPI
import torch as to
import os
import sys


class Permeability(sf.PermeabilityBase):
	def __init__(self, n_elems: int):
		super().__init__(n_elems)

	def compute(self) -> to.Tensor:
		k = 1e-12*to.ones(self.n_elems, dtype=to.float64)
		return k


def main():
	# Read grid
	grid_path = os.path.join("..", "..", "..", "grids", "cube")
	grid = sf.GridHandlerGMSH("geom", grid_path)

	# Define output folder
	output_folder = os.path.join("output", "case_0", "linear")

	# Define equation
	mass_eq = sf.MassPorousMedia(grid, is_linear=True)

	# Define solver
	solver_heat = PETSc.KSP().create(grid.mesh.comm)
	solver_heat.setType("cg")
	solver_heat.getPC().setType("asm")
	solver_heat.setTolerances(rtol=1e-12, max_it=100)
	mass_eq.set_solver(solver_heat)

	# Build material properties
	mat = sf.Material(mass_eq.n_elems)

	# Set solid grain compressibility
	cs = 3.59e-11*to.ones(mass_eq.n_elems, dtype=to.float64)
	mat.set_cs(cs)

	# Set fluid compressibility
	cf = 4.59e-4*to.ones(mass_eq.n_elems, dtype=to.float64)
	mat.set_cf(cf)

	# Set porosity
	phi = 0.03*to.ones(mass_eq.n_elems, dtype=to.float64)
	mat.set_porosity(phi)

	# Set permeability
	kappa_model = Permeability(mass_eq.n_elems)
	mat.set_permeability_model(kappa_model)

	# Set material properties to mass equation
	mass_eq.set_material(mat)

	# Time settings for equilibrium stage
	t_control = sf.TimeControllerParabolic(n_time_steps=50, initial_time=0.0, final_time=5, time_unit="day")

	# Define boundary conditions for mass equation
	bc_east = massBC.DirichletBC(boundary_name = "EAST", 
							values = [1.0, 1.0],
							time_values = [t_control.t_initial, t_control.t_final])
	
	bc_west = massBC.DirichletBC(boundary_name = "WEST", 
							values = [0.0, 0.0],
							time_values = [t_control.t_initial, t_control.t_final])
    
	bc_handler = massBC.BcHandler(mass_eq)
	bc_handler.add_boundary_condition(bc_east)
	bc_handler.add_boundary_condition(bc_west)
	mass_eq.set_boundary_conditions(bc_handler)


	# Set initial temperature field
	P0_field = 0.0*to.ones(mass_eq.n_nodes, dtype=to.float64)
	mass_eq.set_initial_P(P0_field)

	# Create output handlers
	output_mass = sf.SaveFields(mass_eq)
	output_mass.set_output_folder(output_folder)
	output_mass.add_output_field("P", "Pressure (Pa)")
	outputs = [output_mass]

	# Print output folder
	if MPI.COMM_WORLD.rank == 0:
		print(output_folder)
		sys.stdout.flush()

	# Define simulator
	sim = sf.Simulator_H(mass_eq, t_control, outputs)
	sim.run()



if __name__ == '__main__':
	main()