"""
Useful class to assist building the input_file.json.
"""
# Copyright 2024 The safeincave community.
#
# This file is part of safeincave.
#
# Licensed under the GNU GENERAL PUBLIC LICENSE, Version 3 (the "License"); you may not
# use this file except in compliance with the License.  You may obtain a copy
# of the License at
#
#     https://spdx.org/licenses/GPL-3.0-or-later.html
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS, WITHOUT
# WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.  See the
# License for the specific language governing permissions and limitations under
# the License.

import json
import numpy as np
import torch as to
import dolfin as do

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
		self.grid = GridHandlerGMSH(grid_name, path_to_grid)
		self.list_of_boundary_names = list(self.grid.get_boundary_names())
		self.list_of_subdomain_names = list(self.grid.get_subdomain_names())
		self.n_elems = self.grid.mesh.num_cells()
		self.__load_region_indices()

		self.input_file["grid"]["regions"] = {value: key for key, value in self.tags_dict.items()}
		self.input_file["grid"]["boundaries"] = self.list_of_boundary_names



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
			self.input_file["boundary_conditions"][boundary_name]["direction"] = 0
			self.input_file["boundary_conditions"][boundary_name]["density"] = 0
			self.input_file["boundary_conditions"][boundary_name]["reference_position"] = 0.0
			self.input_file["boundary_conditions"][boundary_name]["values"] = list(np.zeros(len(self.input_file["time_settings"]["time_list"])))

	def add_neumann(self, name, values, direction=0, density=0.0, reference_position=0.0):
		assert name in self.list_of_boundary_names, f"{name} is not in {self.list_of_boundary_names}."
		self.input_file["boundary_conditions"][name] = {
			"type": "neumann",
			"direction": direction,
			"density": density,
			"reference_position": reference_position,
			"values": values
		}

	def add_dirichlet(self, name, values, component):
		assert name in self.list_of_boundary_names, f"{name} is not in {self.list_of_boundary_names}."
		self.input_file["boundary_conditions"][name] = {
			"type": "dirichlet",
			"component": component,
			"values": values
		}


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

	def __correct_data_type(self, data):
		if type(data)==np.ndarray or type(data)==to.Tensor:
			data = data.tolist()
		return data

	def add_element(self, element_name, element_parameters, element_type="Elastic"):
		self.input_file["constitutive_model"][element_type][element_name] = element_parameters

	def add_elastic_element(self, name, E, nu, active=True):
		self.input_file["constitutive_model"]["Elastic"][name] = {
			"type": "Spring",
			"active": active,
			"parameters": {
				"E": self.__correct_data_type(E),
				"nu": self.__correct_data_type(nu)
			}
		}

	def add_viscoelastic_element(self, name, E, nu, eta, active=True):
		self.input_file["constitutive_model"]["Viscoelastic"][name] = {
			"type": "KelvinVoigt",
			"active": active,
			"parameters": {
				"E": self.__correct_data_type(E),
				"nu": self.__correct_data_type(nu),
				"eta": self.__correct_data_type(eta),
			}
		}

	def add_dislocation_creep_element(self, name, A, n, Q, T, active=True):
		self.input_file["constitutive_model"]["Inelastic"][name] = {
			"type": "DislocationCreep",
			"active": active,
			"parameters": {
				"A": self.__correct_data_type(A),
				"n": self.__correct_data_type(n),
				"Q": self.__correct_data_type(Q),
				"T": self.__correct_data_type(T),
			}
		}

	def add_desai_element(self, name, mu_1, N_1, n, a_1, eta, beta_1, beta, m, gamma, alpha_0, sigma_t, active=True):
		self.input_file["constitutive_model"]["Inelastic"][name] = {
			"type": "ViscoplasticDesai",
			"active": active,
			"parameters": {
				"mu_1": 	self.__correct_data_type(mu_1),
				"N_1": 		self.__correct_data_type(N_1),
				"n": 		self.__correct_data_type(n),
				"a_1":		self.__correct_data_type(a_1),
				"eta": 		self.__correct_data_type(eta),
				"beta_1": 	self.__correct_data_type(beta_1),
				"beta": 	self.__correct_data_type(beta),
				"m": 		self.__correct_data_type(m),
				"gamma": 	self.__correct_data_type(gamma),
				"alpha_0": 	self.__correct_data_type(alpha_0),
				"sigma_t": 	self.__correct_data_type(sigma_t)
			}
		}

	# def add_viscoelastic_element(self, element_name, element_parameters):
	# 	self.add_element(element_name, element_parameters, element_type="Viscoelastic")

	def add_inelastic_element(self, element_name, element_parameters):
		self.add_element(element_name, element_parameters, element_type="Inelastic")

	def __load_region_indices(self):
		region_names = list(self.grid.get_subdomain_names())
		self.region_indices = {}
		self.tags_dict = {}
		for i in range(len(region_names)):
			self.region_indices[region_names[i]] = []
			tag = self.grid.get_subdomain_tags(region_names[i])
			self.tags_dict[tag] = region_names[i]

		for cell in do.cells(self.grid.mesh):
			region_marker = self.grid.subdomains[cell]
			self.region_indices[self.tags_dict[region_marker]].append(cell.index())

