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

def define_solver(input_file):
	solver_type = input_file["solver_settings"]["type"]
	if solver_type == "KrylovSolver":
		solver = KrylovSolver(
		                      	method = input_file["solver_settings"]["method"],
		                      	preconditioner = input_file["solver_settings"]["preconditioner"]
		                      )
		solver.parameters["relative_tolerance"] = input_file["solver_settings"]["relative_tolerance"]
	elif solver_type == "LU":
		solver = LUSolver(input_file["solver_settings"]["method"])
		solver.parameters["symmetric"] = input_file["solver_settings"]["symmetric"]
	else:
		raise Exception(f"Solver type {solver_type} not supported. Choose between KrylovSolver and LU.")
	return solver


def main_0():
	start_0 = time.time()

	# Read input file
	input_file = read_json("input_file_0.json")

	# This is input_file to be saved
	input_file_to_be_saved = copy.deepcopy(input_file)

	# Output folder
	output_folder = input_file["output"]["path"]
	
	# Create mesh
	g = GridHandlerGMSH(input_file["grid"]["name"], input_file["grid"]["path"])

	# Define constitutive model
	m = ConstitutiveModelHandler(g, input_file)

	# Define solver
	solver = define_solver(input_file)

	# Save input file
	filename = os.path.join(os.path.join(input_file["output"]["path"], "input_file.json"))
	os.makedirs(os.path.dirname(filename), exist_ok=True)
	save_json(input_file_to_be_saved, filename)

	# Compute equilibrium
	u_0 = compute_equilibrium(g, m, solver, input_file)

	# Run (transient) simulation
	run_simulation(g, m, solver, u_0, input_file)

	# Print simulation CPU time
	final = time.time()
	print(f"Time: {final - start_0} seconds.")






if __name__ == '__main__':
	main_0()