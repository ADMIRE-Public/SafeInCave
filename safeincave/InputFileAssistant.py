import json
import numpy as np

class BuildInputFile():
	def __init__(self):
		self.input_file = {}

	def save_input_file(self, input_file_name):
		with open(input_file_name, "w") as f:
		    json.dump(self.input_file, f, indent=4)

	def section_input_grid(self, path_to_grid, grid_name="geom"):
		self.input_file["grid"] = {}
		self.input_file["grid"]["path"] = path_to_grid
		self.input_file["grid"]["name"] = grid_name

		from Grid import GridHandlerGMSH
		grid = GridHandlerGMSH(grid_name, path_to_grid)
		self.list_of_boundary_names = list(grid.get_boundary_names())
		self.list_of_subdomain_names = list(grid.get_subdomain_names())
		self.n_elems = grid.mesh.num_cells()

	def section_output(self, output_folder):
		self.input_file["output"] = {}
		self.input_file["output"]["path"] = output_folder

	def section_solver(self, solver_settings):
		self.input_file["solver_settings"] = solver_settings

	def section_simulation(self, simulation_settings):
		self.input_file["simulation_settings"] = simulation_settings

	def section_time(self, time_list, theta=0.5):
		self.input_file["time_settings"] = {}
		self.input_file["time_settings"]["theta"] = theta
		self.input_file["time_settings"]["time_list"] = time_list

	def section_boundary_conditions(self):
		self.input_file["boundary_conditions"] = {}
		for boundary_name in self.list_of_boundary_names:
			self.input_file["boundary_conditions"][boundary_name] = {}
			self.input_file["boundary_conditions"][boundary_name]["type"] = "neumann"
			self.input_file["boundary_conditions"][boundary_name]["component"] = 0
			self.input_file["boundary_conditions"][boundary_name]["values"] = list(np.zeros(len(self.input_file["time_settings"]["time_list"])))

	def add_boundary_condition(self, boundary_name, bc_data):
		bc_type = bc_data["type"]
		assert boundary_name in self.list_of_boundary_names, f"{boundary_name} is not in {self.list_of_boundary_names}."
		assert bc_type in ["dirichlet", "neumann"], f"{bc_type} must be either 'dirichlet' or 'neumann'."
		self.input_file["boundary_conditions"][boundary_name] = bc_data

	def section_body_forces(self, value, direction):
		self.input_file["body_force"] = {
											"gravity": -9.81,
											"density": value,
											"direction": direction
		}

	def section_constitutive_model(self):
		self.input_file["constitutive_model"] = {
		    "Elastic": {},
		    "Viscoelastic": {},
		    "Inelastic": {}
		}

	def add_element(self, element_name, element_parameters, element_type="Elastic"):
		self.input_file["constitutive_model"][element_type][element_name] = element_parameters

	def add_elastic_element(self, element_name, element_parameters):
		self.add_element(element_name, element_parameters, element_type="Elastic")

	def add_viscoelastic_element(self, element_name, element_parameters):
		self.add_element(element_name, element_parameters, element_type="Viscoelastic")

	def add_inelastic_element(self, element_name, element_parameters):
		self.add_element(element_name, element_parameters, element_type="Inelastic")