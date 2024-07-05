import os
import sys
sys.path.append(os.path.join("..", "..", "libs"))
from ConstitutiveModel import *
from Grid import GridHandlerGMSH
from Utils import *
from dolfin import *
from equilibrium import compute_equilibrium
from operation import run_simulation
import copy
import time

def build_constitutive_model(input_model, n_elems):
	theta = 0.5
	m = ConstitutiveModelHandler(theta, n_elems)

	# Add elastic element(s)
	elems_e = get_list_of_elements(input_model, n_elems, element_class="Elastic")
	for elem_e in elems_e:
		m.add_elastic_element(elem_e)

	# Add viscoelastic element
	elems_ve = get_list_of_elements(input_model, n_elems, element_class="Viscoelastic")
	for elem_ve in elems_ve:
		m.add_viscoelastic_element(elem_ve)

	# Add viscoelastic element
	elems_ie = get_list_of_elements(input_model, n_elems, element_class="Inelastic")
	for elem_ie in elems_ie:
		m.add_inelastic_element(elem_ie)
	
	# Initialize constitutive model
	m.initialize()
	return m

theta_method_dict = {
	"IMP": 0.0,
	"CN": 0.5,
	"EXP": 1.0
}

def build_model_name(input_model):
	model_name = "model"
	for key in input_model["Elastic"].keys():
		if input_model["Elastic"][key]["active"] == True:
			model_name += "_e"
	for key in input_model["Viscoelastic"].keys():
		if input_model["Viscoelastic"][key]["active"] == True and input_model["Viscoelastic"][key]["parameters"]["eta"] > 1.0:
			model_name += "_ve"
	for key in input_model["Inelastic"].keys():
		if input_model["Inelastic"][key]["active"] == True:
			if input_model["Inelastic"][key]["type"] == "ViscoplasticDesai":
				model_name += "_vp"
			elif input_model["Inelastic"][key]["type"] == "DislocationCreep":
				model_name += "_cr"
	return model_name



def main():
	start_0 = time.time()

	# Input mesh
	mesh_name = "cavern_regular"

	# Choose integration method
	theta_name = "EXP"

	# Create mesh
	mesh_folder = os.path.join("..", "..", "grids", mesh_name)
	g = GridHandlerGMSH("geom", mesh_folder)
	n_elems = g.mesh.num_cells()
	
	# Read input model
	input_model = read_json("input_model.json")
	model_name = build_model_name(input_model)

	# This is input_model to be saved
	input_model_to_be_saved = copy.deepcopy(input_model)
	
	# Read input_bc
	input_bc = read_json(f"input_bc.json")

	# Define constitutive model
	m = build_constitutive_model(input_model, n_elems)

	# Choose case name
	case_name = os.path.join(mesh_name, model_name)

	# Compose output name
	output_folder = os.path.join("output", mesh_name, model_name, theta_name)
	print(output_folder)

	# Run initial simulation with elastic and viscoelastic elements
	u_0 = compute_equilibrium(m, g, input_bc, output_folder)

	# Run simulation with full constitutive model
	m.theta = theta_method_dict[theta_name]
	run_simulation(m, g, u_0, input_bc, output_folder)

	# Save inputs
	save_json(input_bc, os.path.join(output_folder, "input_bc.json"))
	save_json(input_model_to_be_saved, os.path.join(output_folder, "input_model.json"))

	# Print simulation CPU time
	final = time.time()
	print(f"Time: {final - start_0} seconds.")


if __name__ == '__main__':
	main()