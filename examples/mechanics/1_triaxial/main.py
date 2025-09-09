import dolfinx as do
import os
import torch as to
from petsc4py import PETSc
import safeincave as sf
import safeincave.Utils as ut
import safeincave.MomentumBC as momBC

class LinearMomentumMod(sf.LinearMomentum):
	def __init__(self, grid, theta):
		super().__init__(grid, theta)

	def initialize(self) -> None:
		self.C.x.array[:] = to.flatten(self.mat.C)
		self.Fvp = do.fem.Function(self.DG0_1)
		self.eps_ve = do.fem.Function(self.DG0_3x3)
		self.eps_cr = do.fem.Function(self.DG0_3x3)
		self.eps_vp = do.fem.Function(self.DG0_3x3)

	def run_after_solve(self):
		self.eps_ve.x.array[:] = to.flatten(self.mat.elems_ne[0].eps_ne_k)
		self.eps_cr.x.array[:] = to.flatten(self.mat.elems_ne[1].eps_ne_k)
		self.eps_vp.x.array[:] = to.flatten(self.mat.elems_ne[2].eps_ne_k)
		self.Fvp.x.array[:] = self.mat.elems_ne[2].Fvp



def main():
	# Read grid
	grid_path = os.path.join("..", "..", "..", "grids", "cube")
	grid = sf.GridHandlerGMSH("geom", grid_path)

	# Define output folder
	output_folder = os.path.join("output", "case_0")

	# Time settings for equilibrium stage
	unit = "hour"
	t_0 = 0.0
	dt = 0.5
	t_final = 24
	t_control = sf.TimeController(dt=dt, initial_time=t_0, final_time=t_final, time_unit=unit)

	# Define momentum equation
	mom_eq = LinearMomentumMod(grid, theta=0.5)

	# Define solver
	mom_solver = PETSc.KSP().create(grid.mesh.comm)
	mom_solver.setType("cg")
	mom_solver.getPC().setType("asm")
	mom_solver.setTolerances(rtol=1e-12, max_it=100)
	mom_eq.set_solver(mom_solver)

	# Define material properties
	mat = sf.Material(mom_eq.n_elems)

	# Set material density
	rho = 2000.0*to.ones(mom_eq.n_elems, dtype=to.float64)
	mat.set_density(rho)

	# Constitutive model
	E = 102e9*to.ones(mom_eq.n_elems)
	nu = 0.3*to.ones(mom_eq.n_elems)
	spring_0 = sf.Spring(E, nu, "spring")

	# Create Kelvin-Voigt viscoelastic element
	eta = 105e11*to.ones(mom_eq.n_elems)
	E = 10e9*to.ones(mom_eq.n_elems)
	nu = 0.32*to.ones(mom_eq.n_elems)
	kelvin = sf.Viscoelastic(eta, E, nu, "kelvin")

	# Create creep
	A = 1.9e-20*to.ones(mom_eq.n_elems)
	Q = 51600*to.ones(mom_eq.n_elems)
	n = 3.0*to.ones(mom_eq.n_elems)
	creep_0 = sf.DislocationCreep(A, Q, n, "creep")

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
	desai = sf.ViscoplasticDesai(mu_1, N_1, a_1, eta, n, beta_1, beta, m, gamma, sigma_t, alpha_0, "desai")

	# Create constitutive model
	mat.add_to_elastic(spring_0)
	mat.add_to_non_elastic(kelvin)
	mat.add_to_non_elastic(creep_0)
	mat.add_to_non_elastic(desai)

	# Set constitutive model
	mom_eq.set_material(mat)

	# Set body forces
	g_vec = [0.0, 0.0, 0.0]
	mom_eq.build_body_force(g_vec)

	# Set initial temperature field
	T0_field = 293*to.ones(mom_eq.n_elems)
	mom_eq.set_T0(T0_field)
	mom_eq.set_T(T0_field)

	# Boundary conditions
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
						values =      [4.0*ut.MPa, 4.0*ut.MPa],
						time_values = [0.0, t_control.t_final],
						g = g_vec[2])

	bc_north = momBC.NeumannBC(boundary_name = "NORTH",
						direction = 2,
						density = 0.0,
						ref_pos = 0.0,
						values =      [4.0*ut.MPa, 4.0*ut.MPa],
						time_values = [0.0, t_control.t_final],
						g = g_vec[2])

	bc_top = momBC.NeumannBC(boundary_name = "TOP",
						direction = 2,
						density = 0.0,
						ref_pos = 0.0,
						values =      [4.1*ut.MPa, 16*ut.MPa, 16*ut.MPa,  6*ut.MPa,   6*ut.MPa],
						time_values = [0*ut.hour,  2*ut.hour, 14*ut.hour, 16*ut.hour, 24*ut.hour],
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

	# Create output handlers
	output_mom = sf.SaveFields(mom_eq)
	output_mom.set_output_folder(output_folder)
	output_mom.add_output_field("u", "Displacement (m)")
	output_mom.add_output_field("eps_tot", "Total strain (-)")
	output_mom.add_output_field("eps_ve", "Viscoelastic strain (-)")
	output_mom.add_output_field("eps_cr", "Creep strain (-)")
	output_mom.add_output_field("eps_vp", "Viscoplastic strain (-)")
	output_mom.add_output_field("Fvp", "Yield function (-)")
	outputs = [output_mom]

	# Define simulator
	sim = sf.Simulator_M(mom_eq, t_control, outputs, compute_elastic_response=True)
	sim.run()



if __name__ == '__main__':
	main()