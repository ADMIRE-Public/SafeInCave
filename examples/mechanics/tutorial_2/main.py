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
from TimeHandler import TimeController, TimeControllerParabolic
import time


class LinearMomentumMod(LinearMomentum):
	def __init__(self, grid, theta):
		super().__init__(grid, theta)

	def initialize(self) -> None:
		self.C.x.array[:] = to.flatten(self.mat.C)
		self.Fvp = do.fem.Function(self.DG0_1)
		self.alpha = do.fem.Function(self.DG0_1)
		self.eps_vp = do.fem.Function(self.DG0_3x3)

	def run_after_solve(self):
		try:
			self.eps_vp.x.array[:] = to.flatten(self.mat.elems_ne[-1].eps_ne_k)
			self.Fvp.x.array[:] = self.mat.elems_ne[-1].Fvp
			self.alpha.x.array[:] = self.mat.elems_ne[-1].alpha
		except:
			pass



def main():
	comm = MPI.COMM_WORLD
	comm.Barrier()
	if MPI.COMM_WORLD.rank == 0:
	    start_time = MPI.Wtime()

	# Read grid
	grid_path = os.path.join("..", "..", "grids", "cavern_irregular")
	grid = GridHandlerGMSH("geom", grid_path)

	# Define output folder
	output_folder = os.path.join("output", "case_2")

	# Define momentum equation
	mom_eq = LinearMomentumMod(grid, theta=0.5)

	# Define solver
	mom_solver = PETSc.KSP().create(grid.mesh.comm)
	mom_solver.setType("bicg")
	mom_solver.getPC().setType("asm")
	mom_solver.setTolerances(rtol=1e-12, max_it=100)
	mom_eq.set_solver(mom_solver)

	# Define material properties
	mat = Material(mom_eq.n_elems)

	# Set material density
	salt_density = 2000
	rho = salt_density*to.ones(mom_eq.n_elems, dtype=to.float64)
	mat.set_density(rho)

	# Constitutive model
	E0 = 102*utils.GPa*to.ones(mom_eq.n_elems)
	nu0 = 0.3*to.ones(mom_eq.n_elems)
	spring_0 = Spring(E0, nu0, "spring")

	# Create Kelvin-Voigt viscoelastic element
	eta = 105e11*to.ones(mom_eq.n_elems)
	E1 = 10*utils.GPa*to.ones(mom_eq.n_elems)
	nu1 = 0.32*to.ones(mom_eq.n_elems)
	kelvin = Viscoelastic(eta, E1, nu1, "kelvin")

	# Create creep
	A = 1.9e-20*to.ones(mom_eq.n_elems)
	Q = 51600*to.ones(mom_eq.n_elems)
	n = 3.0*to.ones(mom_eq.n_elems)
	creep_0 = DislocationCreep(A, Q, n, "creep")

	# Create constitutive model
	mat.add_to_elastic(spring_0)
	mat.add_to_non_elastic(kelvin)
	mat.add_to_non_elastic(creep_0)

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

	# Time settings for equilibrium stage
	# tc_equilibrium = TimeControllerParabolic(final_time=10, initial_time=0.0, n_time_steps=50, time_unit="day")
	tc_equilibrium = TimeController(time_step=0.5, final_time=10, initial_time=0.0, time_unit="hour")

	# Boundary conditions
	time_values = [0*utils.hour,  1*utils.hour]
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

	
	# bc_east = momBC.DirichletBC(boundary_name = "East", 
	# 				 		component = 0,
	# 						values = [0.0, 0.0],
	# 						time_values = [0.0, tc_equilibrium.t_final])

	# bc_north = momBC.DirichletBC(boundary_name = "North", 
	# 				 	  component = 1,
	# 				 	  values = [0.0, 0.0],
	# 				 	  time_values = [0.0, tc_equilibrium.t_final])

	side_burden = 10.0*utils.MPa
	bc_east = momBC.NeumannBC(boundary_name = "East",
						direction = 2,
						density = salt_density,
						ref_pos = 660.0,
						values =      [side_burden, side_burden],
						time_values = [0.0,            tc_equilibrium.t_final],
						g = g_vec[2])

	bc_north = momBC.NeumannBC(boundary_name = "North",
						direction = 2,
						density = salt_density,
						ref_pos = 660.0,
						values =      [side_burden, side_burden],
						time_values = [0.0,            tc_equilibrium.t_final],
						g = g_vec[2])

	over_burden = 10.0*utils.MPa
	bc_top = momBC.NeumannBC(boundary_name = "Top",
						direction = 2,
						density = 0.0,
						ref_pos = 0.0,
						values =      [over_burden, over_burden],
						time_values = [0.0,            tc_equilibrium.t_final],
						g = g_vec[2])

	gas_density = 0.082
	p_gas = 10.0*utils.MPa
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

	# Print output folder
	if MPI.COMM_WORLD.rank == 0:
		print(ouput_folder_equilibrium)

	# Create output handlers
	output_mom = SaveFields(mom_eq)
	output_mom.set_output_folder(ouput_folder_equilibrium)
	output_mom.add_output_field("u", "Displacement (m)")
	output_mom.add_output_field("eps_tot", "Total strain (-)")
	output_mom.add_output_field("sig", "Stress (Pa)")
	output_mom.add_output_field("p_elems", "Mean stress (Pa)")
	output_mom.add_output_field("q_elems", "Von Mises stress (Pa)")
	outputs = [output_mom]

	# Define simulator
	sim = Simulator_M(mom_eq, tc_equilibrium, outputs, True)
	sim.run()

	# Print time
	if MPI.COMM_WORLD.rank == 0:
		end_time = MPI.Wtime()
		elaspsed_time = end_time - start_time
		formatted_time = time.strftime("%H:%M:%S", time.gmtime(elaspsed_time))
		print(f"Time: {formatted_time} ({elaspsed_time} seconds)\n")








	# Create Desai's viscoplastic model
	mu_1 = 5.3665857009859815e-11*to.ones(mom_eq.n_elems)
	N_1 = 3.1*to.ones(mom_eq.n_elems)
	n = 3.0*to.ones(mom_eq.n_elems)
	a_1 = 1.965018496922832e-05*to.ones(mom_eq.n_elems)
	eta = 0.8275682807874163*to.ones(mom_eq.n_elems)
	beta_1 = 0.0048*to.ones(mom_eq.n_elems)
	beta = 0.995*to.ones(mom_eq.n_elems)
	m = -0.5*to.ones(mom_eq.n_elems)
	gamma = 0.095*to.ones(mom_eq.n_elems)
	alpha_0 = 0.0022*to.ones(mom_eq.n_elems)
	sigma_t = 5.0*to.ones(mom_eq.n_elems)
	desai = ViscoplasticDesai(mu_1, N_1, a_1, eta, n, beta_1, beta, m, gamma, sigma_t, alpha_0, "desai")

	# Compute initial hardening parameter
	stress_to = utils.numpy2torch(mom_eq.sig.x.array.reshape((mom_eq.n_elems, 3, 3)))
	desai.compute_initial_hardening(stress_to, Fvp_0=0.0)

	# Add viscoplastic element to constitutive model
	# mom_eq.mat.add_to_non_elastic(desai)

	# Time settings for equilibrium stage
	tc_operation = TimeController(time_step=0.1, final_time=24, initial_time=0.0, time_unit="hour")

	# Boundary conditions
	bc_west = momBC.DirichletBC(boundary_name = "West", 
					 		component = 0,
							values = [0.0, 0.0],
							time_values = [0.0, tc_operation.t_final])

	bc_bottom = momBC.DirichletBC(boundary_name = "Bottom", 
					 	  component = 2,
					 	  values = [0.0, 0.0],
					 	  time_values = [0.0, tc_operation.t_final])

	bc_south = momBC.DirichletBC(boundary_name = "South", 
					 	  component = 1,
					 	  values = [0.0, 0.0],
					 	  time_values = [0.0, tc_operation.t_final])

	# bc_east = momBC.DirichletBC(boundary_name = "East", 
	# 				 		component = 0,
	# 						values = [0.0, 0.0],
	# 						time_values = [0.0, tc_operation.t_final])

	# bc_north = momBC.DirichletBC(boundary_name = "North", 
	# 				 	  component = 1,
	# 				 	  values = [0.0, 0.0],
	# 				 	  time_values = [0.0, tc_operation.t_final])

	bc_east = momBC.NeumannBC(boundary_name = "East",
						direction = 2,
						density = salt_density,
						ref_pos = 660.0,
						values =      [side_burden, side_burden],
						time_values = [0.0,            tc_operation.t_final],
						g = g_vec[2])

	bc_north = momBC.NeumannBC(boundary_name = "North",
						direction = 2,
						density = salt_density,
						ref_pos = 660.0,
						values =      [side_burden, side_burden],
						time_values = [0.0,            tc_operation.t_final],
						g = g_vec[2])

	bc_top = momBC.NeumannBC(boundary_name = "Top",
						direction = 2,
						density = 0.0,
						ref_pos = 0.0,
						values =      [over_burden, over_burden],
						time_values = [0.0,            tc_operation.t_final],
						g = g_vec[2])

	bc_cavern = momBC.NeumannBC(boundary_name = "Cavern",
						direction = 2,
						density = gas_density,
						ref_pos = 430.0,
						values =      [10.0*utils.MPa, 7.0*utils.MPa, 7.0*utils.MPa, 10.0*utils.MPa, 10.0*utils.MPa],
						time_values = [0.0, 2.0*utils.hour, 14*utils.hour, 16*utils.hour, 24*utils.hour],
						g = g_vec[2])

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

	# Print output folder
	if MPI.COMM_WORLD.rank == 0:
		print(output_folder_operation)

	# Create output handlers
	output_mom = SaveFields(mom_eq)
	output_mom.set_output_folder(output_folder_operation)
	output_mom.add_output_field("u", "Displacement (m)")
	output_mom.add_output_field("eps_tot", "Total strain (-)")
	output_mom.add_output_field("eps_vp", "Viscoplastic strain (-)")
	output_mom.add_output_field("alpha", "Hardening parameter (-)")
	output_mom.add_output_field("Fvp", "Yield function (-)")
	output_mom.add_output_field("p_elems", "Mean stress (Pa)")
	output_mom.add_output_field("q_elems", "Von Mises stress (Pa)")
	outputs = [output_mom]

	# Define simulator
	sim = Simulator_M(mom_eq, tc_operation, outputs, False)
	sim.run()

	# Print time
	if MPI.COMM_WORLD.rank == 0:
		end_time = MPI.Wtime()
		elaspsed_time = end_time - start_time
		formatted_time = time.strftime("%H:%M:%S", time.gmtime(elaspsed_time))
		print(f"Time: {formatted_time} ({elaspsed_time} seconds)\n")


if __name__ == '__main__':
	main()