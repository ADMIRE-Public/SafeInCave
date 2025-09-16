import safeincave as sf
import safeincave.Utils as ut
# from petsc4py import PETSc
import safeincave.MomentumBC as momBC
# import numpy as np
# import torch as to
# from mpi4py import MPI
# import os



# class Simulator_GUI():
# 	def __init__(self, input_file):
# 		self.input_file = input_file

# 		self.output_folder = self.input_file["output"]["path"]

# 		self.build_grid()
# 		self.initialize_equation()
# 		self.build_solver()
# 		self.initialize_material()
# 		self.set_gravity()


# 	def build_grid(self):
# 		grid_path = self.input_file["grid"]["path"]
# 		grid_name = self.input_file["grid"]["name"]
# 		self.grid = sf.GridHandlerGMSH(grid_name, grid_path)

# 	def initialize_equation(self):
# 		theta = self.input_file["time_settings"]["theta"]
# 		self.mom_eq = sf.LinearMomentum(self.grid, theta=theta)

# 	def set_gravity(self):
# 		g_vec = [0.0, 0.0, 0.0]
# 		i = self.input_file["body_force"]["direction"]
# 		self.g = self.input_file["body_force"]["gravity"]
# 		g_vec[i] = self.g
# 		self.mom_eq.build_body_force(g_vec)

# 	def initialize_material(self):
# 		self.mat = sf.Material(self.grid.n_elems)
# 		density = self.grid.get_parameter(self.input_file["body_force"]["density"])
# 		self.mat.set_density(density)

# 		for elem_name in self.input_file["constitutive_model"]["elastic"].keys():
# 			E = self.grid.get_parameter(self.input_file["constitutive_model"]["elastic"][elem_name]["parameters"]["E"])
# 			nu = self.grid.get_parameter(self.input_file["constitutive_model"]["elastic"][elem_name]["parameters"]["nu"])
# 			spring_0 = sf.Spring(E, nu, elem_name)
# 			self.mat.add_to_elastic(spring_0)

# 		self.mom_eq.set_material(self.mat)

# 	def build_solver(self):
# 		solver = PETSc.KSP().create(self.grid.mesh.comm)
# 		solver.setType("cg")
# 		solver.getPC().setType("ilu")
# 		solver.setTolerances(rtol=1e-12, max_it=100)
# 		self.mom_eq.set_solver(solver)


# 	def run_equilibrium(self):
# 		# Build material: non-elastic element
# 		for elem_name in self.input_file["constitutive_model"]["nonelastic"].keys():
# 			if self.input_file["constitutive_model"]["nonelastic"][elem_name]["active"]:
# 				if self.input_file["constitutive_model"]["nonelastic"][elem_name]["equilibrium"]:
# 					if self.input_file["constitutive_model"]["nonelastic"][elem_name]["type"] == "KelvinVoigt":
# 						E = self.grid.get_parameter(self.input_file["constitutive_model"]["nonelastic"][elem_name]["parameters"]["E"])
# 						nu = self.grid.get_parameter(self.input_file["constitutive_model"]["nonelastic"][elem_name]["parameters"]["nu"])
# 						eta = self.grid.get_parameter(self.input_file["constitutive_model"]["nonelastic"][elem_name]["parameters"]["eta"])
# 						kelvin = sf.Viscoelastic(eta, E, nu, elem_name)
# 						self.mom_eq.mat.add_to_non_elastic(kelvin)
# 					elif self.input_file["constitutive_model"]["nonelastic"][elem_name]["type"] == "DislocationCreep":
# 						A = self.grid.get_parameter(self.input_file["constitutive_model"]["nonelastic"][elem_name]["parameters"]["A"])
# 						n = self.grid.get_parameter(self.input_file["constitutive_model"]["nonelastic"][elem_name]["parameters"]["n"])
# 						Q = self.grid.get_parameter(self.input_file["constitutive_model"]["nonelastic"][elem_name]["parameters"]["Q"])
# 						creep_0 = sf.DislocationCreep(A, Q, n, elem_name)
# 						self.mom_eq.mat.add_to_non_elastic(creep_0)
# 						T = self.grid.get_parameter(self.input_file["constitutive_model"]["nonelastic"][elem_name]["parameters"]["T"])
# 						self.mom_eq.set_T0(T)
# 						self.mom_eq.set_T(T)
# 					elif self.input_file["constitutive_model"]["nonelastic"][elem_name]["type"] == "ViscoplasticDesai":
# 						mu_1 = self.grid.get_parameter(self.input_file["constitutive_model"]["nonelastic"][elem_name]["parameters"]["mu_1"])
# 						N_1 = self.grid.get_parameter(self.input_file["constitutive_model"]["nonelastic"][elem_name]["parameters"]["N_1"])
# 						n = self.grid.get_parameter(self.input_file["constitutive_model"]["nonelastic"][elem_name]["parameters"]["n"])
# 						a_1 = self.grid.get_parameter(self.input_file["constitutive_model"]["nonelastic"][elem_name]["parameters"]["a_1"])
# 						eta = self.grid.get_parameter(self.input_file["constitutive_model"]["nonelastic"][elem_name]["parameters"]["eta"])
# 						beta_1 = self.grid.get_parameter(self.input_file["constitutive_model"]["nonelastic"][elem_name]["parameters"]["beta_1"])
# 						beta = self.grid.get_parameter(self.input_file["constitutive_model"]["nonelastic"][elem_name]["parameters"]["beta"])
# 						m = self.grid.get_parameter(self.input_file["constitutive_model"]["nonelastic"][elem_name]["parameters"]["m"])
# 						gamma = self.grid.get_parameter(self.input_file["constitutive_model"]["nonelastic"][elem_name]["parameters"]["gamma"])
# 						alpha_0 = self.grid.get_parameter(self.input_file["constitutive_model"]["nonelastic"][elem_name]["parameters"]["alpha_0"])
# 						sigma_t = self.grid.get_parameter(self.input_file["constitutive_model"]["nonelastic"][elem_name]["parameters"]["sigma_t"])
# 						desai = sf.ViscoplasticDesai(mu_1, N_1, a_1, eta, n, beta_1, beta, m, gamma, sigma_t, alpha_0, elem_name)
# 						self.mom_eq.mat.add_to_non_elastic(desai)
# 					else:
# 						elem_type = self.input_file["constitutive_model"]["nonelastic"][elem_name]["type"]
# 						raise Exception(f"Element type {elem_type} not supported.")


# 		# Time settings for equilibrium stage
# 		dt = self.input_file["simulation_settings"]["equilibrium"]["dt_max"]
# 		tf = self.input_file["simulation_settings"]["equilibrium"]["ite_max"]*dt
# 		tc_equilibrium = sf.TimeController(dt=dt, initial_time=0.0, final_time=tf, time_unit="second")

# 		# Loop over boundaries
# 		bc_equilibrium = momBC.BcHandler(self.mom_eq)
# 		t_values = [0.0, tc_equilibrium.t_final]
# 		for b_name in self.input_file["boundary_conditions"].keys():
# 			bc_value = self.input_file["boundary_conditions"][b_name]["values"][0]
# 			bc_values = bc_value*np.ones(len(t_values))
# 			if self.input_file["boundary_conditions"][b_name]["type"] == "neumann":
# 				bc = momBC.NeumannBC(boundary_name = b_name,
# 									 direction = self.input_file["boundary_conditions"][b_name]["direction"],
# 									 density = self.input_file["boundary_conditions"][b_name]["density"],
# 									 ref_pos = self.input_file["boundary_conditions"][b_name]["reference_position"],
# 									 values = bc_values,
# 									 time_values = t_values,
# 									 g = self.g)
# 			elif self.input_file["boundary_conditions"][b_name]["type"] == "dirichlet":
# 				bc = momBC.DirichletBC(boundary_name = b_name, 
# 								 	   component = self.input_file["boundary_conditions"][b_name]["component"],
# 									   values = bc_values,
# 									   time_values = t_values)
# 			else:
# 				b_type = self.input_file["boundary_conditions"][b_name]["type"]
# 				raise Exception(f"Boundary condition type {b_type} not supported.")
# 			bc_equilibrium.add_boundary_condition(bc)

# 		# Set boundary conditions
# 		self.mom_eq.set_boundary_conditions(bc_equilibrium)

# 		# Create output handlers
# 		output_mom = sf.SaveFields(self.mom_eq)
# 		output_mom.set_output_folder(os.path.join(self.output_folder, "equilibrium"))
# 		output_mom.add_output_field("u", "Displacement (m)")
# 		output_mom.add_output_field("p_elems", "Mean Stress (MPa)")
# 		outputs = [output_mom]

# 		# Define simulator
# 		sim = sf.Simulator_M(self.mom_eq, tc_equilibrium, outputs, compute_elastic_response=True)
# 		sim.run()

# 	def run_operation(self):
# 		# Build material: non-elastic element
# 		for elem_name in self.input_file["constitutive_model"]["nonelastic"].keys():
# 			if self.input_file["constitutive_model"]["nonelastic"][elem_name]["active"]:
# 				if self.input_file["constitutive_model"]["nonelastic"][elem_name]["equilibrium"]:
# 					if self.input_file["constitutive_model"]["nonelastic"][elem_name]["type"] == "KelvinVoigt":
# 						E = self.grid.get_parameter(self.input_file["constitutive_model"]["nonelastic"][elem_name]["parameters"]["E"])
# 						nu = self.grid.get_parameter(self.input_file["constitutive_model"]["nonelastic"][elem_name]["parameters"]["nu"])
# 						eta = self.grid.get_parameter(self.input_file["constitutive_model"]["nonelastic"][elem_name]["parameters"]["eta"])
# 						kelvin = sf.Viscoelastic(eta, E, nu, elem_name)
# 						self.mom_eq.mat.add_to_non_elastic(kelvin)
# 					elif self.input_file["constitutive_model"]["nonelastic"][elem_name]["type"] == "DislocationCreep":
# 						A = self.grid.get_parameter(self.input_file["constitutive_model"]["nonelastic"][elem_name]["parameters"]["A"])
# 						n = self.grid.get_parameter(self.input_file["constitutive_model"]["nonelastic"][elem_name]["parameters"]["n"])
# 						Q = self.grid.get_parameter(self.input_file["constitutive_model"]["nonelastic"][elem_name]["parameters"]["Q"])
# 						creep_0 = sf.DislocationCreep(A, Q, n, elem_name)
# 						self.mom_eq.mat.add_to_non_elastic(creep_0)
# 						T = self.grid.get_parameter(self.input_file["constitutive_model"]["nonelastic"][elem_name]["parameters"]["T"])
# 						self.mom_eq.set_T0(T)
# 						self.mom_eq.set_T(T)
# 					elif self.input_file["constitutive_model"]["nonelastic"][elem_name]["type"] == "ViscoplasticDesai":
# 						mu_1 = self.grid.get_parameter(self.input_file["constitutive_model"]["nonelastic"][elem_name]["parameters"]["mu_1"])
# 						N_1 = self.grid.get_parameter(self.input_file["constitutive_model"]["nonelastic"][elem_name]["parameters"]["N_1"])
# 						n = self.grid.get_parameter(self.input_file["constitutive_model"]["nonelastic"][elem_name]["parameters"]["n"])
# 						a_1 = self.grid.get_parameter(self.input_file["constitutive_model"]["nonelastic"][elem_name]["parameters"]["a_1"])
# 						eta = self.grid.get_parameter(self.input_file["constitutive_model"]["nonelastic"][elem_name]["parameters"]["eta"])
# 						beta_1 = self.grid.get_parameter(self.input_file["constitutive_model"]["nonelastic"][elem_name]["parameters"]["beta_1"])
# 						beta = self.grid.get_parameter(self.input_file["constitutive_model"]["nonelastic"][elem_name]["parameters"]["beta"])
# 						m = self.grid.get_parameter(self.input_file["constitutive_model"]["nonelastic"][elem_name]["parameters"]["m"])
# 						gamma = self.grid.get_parameter(self.input_file["constitutive_model"]["nonelastic"][elem_name]["parameters"]["gamma"])
# 						alpha_0 = self.grid.get_parameter(self.input_file["constitutive_model"]["nonelastic"][elem_name]["parameters"]["alpha_0"])
# 						sigma_t = self.grid.get_parameter(self.input_file["constitutive_model"]["nonelastic"][elem_name]["parameters"]["sigma_t"])
# 						desai = sf.ViscoplasticDesai(mu_1, N_1, a_1, eta, n, beta_1, beta, m, gamma, sigma_t, alpha_0, elem_name)
# 						self.mom_eq.mat.add_to_non_elastic(desai)
# 					else:
# 						elem_type = self.input_file["constitutive_model"]["nonelastic"][elem_name]["type"]
# 						raise Exception(f"Element type {elem_type} not supported.")


# 		# Time settings for operation stage
# 		t_values = self.input_file["time_settings"]["time_list"]
# 		dt = self.input_file["simulation_settings"]["operation"]["dt_max"]
# 		tf = t_values[-1]
# 		tc_op = sf.TimeController(dt=dt, initial_time=0.0, final_time=tf, time_unit="second")

# 		# Loop over boundaries
# 		bc_op = momBC.BcHandler(self.mom_eq)
# 		for b_name in self.input_file["boundary_conditions"].keys():
# 			bc_values = self.input_file["boundary_conditions"][b_name]["values"]
# 			if self.input_file["boundary_conditions"][b_name]["type"] == "neumann":
# 				bc = momBC.NeumannBC(boundary_name = b_name,
# 									 direction = self.input_file["boundary_conditions"][b_name]["direction"],
# 									 density = self.input_file["boundary_conditions"][b_name]["density"],
# 									 ref_pos = self.input_file["boundary_conditions"][b_name]["reference_position"],
# 									 values = bc_values,
# 									 time_values = t_values,
# 									 g = self.g)
# 			elif self.input_file["boundary_conditions"][b_name]["type"] == "dirichlet":
# 				bc = momBC.DirichletBC(boundary_name = b_name, 
# 								 	   component = self.input_file["boundary_conditions"][b_name]["component"],
# 									   values = bc_values,
# 									   time_values = t_values)
# 			else:
# 				b_type = self.input_file["boundary_conditions"][b_name]["type"]
# 				raise Exception(f"Boundary condition type {b_type} not supported.")
# 			bc_op.add_boundary_condition(bc)

# 		# Set boundary conditions
# 		self.mom_eq.set_boundary_conditions(bc_op)

# 		# Create output handlers
# 		output_mom = sf.SaveFields(self.mom_eq)
# 		output_mom.set_output_folder(os.path.join(self.output_folder, "operation"))
# 		output_mom.add_output_field("u", "Displacement (m)")
# 		output_mom.add_output_field("p_elems", "Mean Stress (MPa)")
# 		outputs = [output_mom]

# 		# Define simulator
# 		compute_elastic_response = True
# 		if self.input_file["simulation_settings"]["equilibrium"]["active"] == True:
# 			compute_elastic_response = False
# 		sim = sf.Simulator_M(self.mom_eq, tc_op, outputs, compute_elastic_response=compute_elastic_response)
# 		sim.run()


# 	def run(self):
# 		if self.input_file["simulation_settings"]["equilibrium"]["active"] == True:
# 			self.run_equilibrium()
# 		self.run_operation()

def main():
	input_file = ut.read_json("input_file.json")

	sim = sf.Simulator_GUI(input_file)
	sim.run()

if __name__ == '__main__':
	main()