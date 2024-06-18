import os
import json
import sys
import numpy as np
sys.path.append(os.path.join("..", "..", "libs"))
from Grid import GridHandlerGMSH

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

		grid = GridHandlerGMSH(grid_name, path_to_grid)
		self.list_of_boundary_names = list(grid.get_boundary_names())
		self.list_of_subdomain_names = list(grid.get_subdomain_names())
		self.n_elems = grid.mesh.num_cells()

	def section_time_settings(self, time_list, theta=0.5):
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

	def add_boundary_condition(self, boundary_name, bc_values, bc_type="dirichlet", component=0.):
		assert boundary_name in self.list_of_boundary_names, f"{boundary_name} is not in {self.list_of_boundary_names}."
		assert bc_type in ["dirichlet", "neumann"], f"{bc_type} must be either 'dirichlet' or 'neumann'."
		self.input_file["boundary_conditions"][boundary_name] = {}
		self.input_file["boundary_conditions"][boundary_name]["values"] = bc_values
		self.input_file["boundary_conditions"][boundary_name]["component"] = component
		self.input_file["boundary_conditions"][boundary_name]["type"] = bc_type

	def section_material_properties(self):
		self.input_file["material_properties"] = {
		    "Elastic": {},
		    "Viscoelastic": {},
		    "Inelastic": {}
		}

	def add_elastic_element(self, element_name, element_parameters):
		self.input_file["material_properties"]["Elastic"][element_name] = element_parameters

	def add_viscoelastic_element(self, element_name, element_parameters):
		self.input_file["material_properties"]["Viscoelastic"][element_name] = element_parameters

	def add_inelastic_element(self, element_name, element_parameters):
		self.input_file["material_properties"]["Inelastic"][element_name] = element_parameters





if __name__ == '__main__':
	# Useful units
	hour = 60*60
	day = 24*hour
	MPa = 1e6

	# Initialize input file object
	bif = BuildInputFile()

	# Create input_grid section
	path_to_grid = os.path.join("..", "..", "grids", "cube_0")
	bif.section_input_grid(path_to_grid, "geom")

	# Create time_settings section
	time_list = [0*hour,  2*hour,  10*hour, 12*hour, 14*hour, 16*hour, 20*hour, 22*hour, 24*hour]
	bif.section_time_settings(time_list, theta=0.5)

	# # Create boundary_conditions section
	dirichlet_zeros = list(np.zeros(len(time_list)))
	bif.section_boundary_conditions()
	bif.add_boundary_condition("EAST", dirichlet_zeros, "dirichlet", 0)
	bif.add_boundary_condition("SOUTH", dirichlet_zeros, "dirichlet", 1)
	bif.add_boundary_condition("BOTTOM", dirichlet_zeros, "dirichlet", 2)

	stress_horizontal = [5*MPa, 5*MPa, 5*MPa, 5*MPa, 5*MPa, 5*MPa, 5*MPa, 5*MPa, 5*MPa]
	stress_vertical = [6*MPa, 10*MPa, 10*MPa, 6*MPa, 6*MPa, 12*MPa, 12*MPa, 6*MPa, 6*MPa]
	bif.add_boundary_condition("WEST", stress_horizontal, "neumann")
	bif.add_boundary_condition("NORTH", stress_horizontal, "neumann")
	bif.add_boundary_condition("TOP", stress_vertical, "neumann")

	# Assign material properties
	bif.section_material_properties()

	# Add elastic properties
	elastic_parameters = {
        "type": "Spring",
        "active": True,
        "parameters": {
            "E": 	list(102e9*np.ones(bif.n_elems)),
            "nu": 	list(0.3*np.ones(bif.n_elems))
        }
    }
	bif.add_elastic_element("Spring_0", elastic_parameters)

	# Add viscoelastic properties
	viscoelastic_parameters = {
        "type": "KelvinVoigt",
        "active": True,
        "parameters": {
            "E": 	list(102e9*np.ones(bif.n_elems)),
            "nu": 	list(0.3*np.ones(bif.n_elems)),
            "eta": 	list(10.5e12*np.ones(bif.n_elems))
		}
    }
	bif.add_viscoelastic_element("KelvinVoigt_0", viscoelastic_parameters)

	# Add viscoplastic parameters
	viscoplastic_parameters = {
        "type": "ViscoplasticDesai",
        "active": False,
        "parameters": {
            "F_0": 		list(1.0*np.ones(bif.n_elems)),
            "mu_1": 	list(5.3665857009859815e-11*np.ones(bif.n_elems)),
            "N_1": 		list(3.1*np.ones(bif.n_elems)),
            "n": 		list(3.0*np.ones(bif.n_elems)),
            "a_1":		list(1.965018496922832e-05*np.ones(bif.n_elems)),
            "eta": 		list(0.8275682807874163*np.ones(bif.n_elems)),
            "beta_1": 	list(0.0048*np.ones(bif.n_elems)),
            "beta": 	list(0.995*np.ones(bif.n_elems)),
            "m": 		list(-0.5*np.ones(bif.n_elems)),
            "gamma": 	list(0.095*np.ones(bif.n_elems)),
            "alpha_1": 	list(0.0022*np.ones(bif.n_elems)),
            "alpha_0": 	list(0.0040715714049800586*np.ones(bif.n_elems)),
            "k_v": 		list(0.0*np.ones(bif.n_elems)),
            "sigma_t": 	list(5.0*np.ones(bif.n_elems))
		}
    }
	bif.add_inelastic_element("desai", viscoplastic_parameters)

	# Add dislocation creep parameters
	creep_parameters = {
        "type": "DislocationCreep",
        "active": True,
        "parameters": {
            "A": 	list(1.9e-20*np.ones(bif.n_elems)),
            "n": 	list(3.0*np.ones(bif.n_elems)),
            "T": 	list(298*np.ones(bif.n_elems)),
            "Q": 	list(51600*np.ones(bif.n_elems)),
            "R": 	list(8.32*np.ones(bif.n_elems))
		}
    }
	bif.add_inelastic_element("creep", creep_parameters)




	# Save input_file.json
	bif.save_input_file("input_file_0.json")