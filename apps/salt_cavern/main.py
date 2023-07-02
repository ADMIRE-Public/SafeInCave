import os
import sys
import numpy as np
import time
import shutil
sys.path.append(os.path.join("..", "..", "libs"))
from Simulators import H2CaveSimulator
from Utils import *

look_up_table = np.array([
	[0.0*day, 20*MPa],
	[2.0*day, 10*MPa],
	[3.0*day, 15*MPa],
	[4.0*day, 10*MPa],
	[5.0*day, 15*MPa],
	[6.0*day, 10*MPa],
	[7.0*day, 15*MPa],
	[8.0*day, 10*MPa],
	[9.0*day, 15*MPa],
	[10.0*day, 10*MPa],
	[11.0*day, 15*MPa],
	[12.0*day, 10*MPa],
	[13.0*day, 15*MPa],
	[14.0*day, 10*MPa],
	[15.0*day, 15*MPa],
	[16.0*day, 10*MPa],
	[17.0*day, 15*MPa],
	[18.0*day, 10*MPa],
	[19.0*day, 15*MPa],
	[20.0*day, 10*MPa],
])

def pressure(t_list):
	p_list = []
	for t in t_list:
		p = np.interp(t, look_up_table[:,0], look_up_table[:,1])
		p_list.append(p)
	return p_list

def write_settings(input_bc):
	# Define time levels
	n_steps = 75
	t_f = 5*day
	time_list = list(np.linspace(0, t_f, n_steps))
	input_bc["Time"]["timeList"] = time_list

	# Define boundary conditions
	input_bc["BoundaryConditions"]["u_x"]["SIDE_X"]["type"] = "DIRICHLET"
	input_bc["BoundaryConditions"]["u_x"]["SIDE_X"]["value"] = list(np.repeat(0.0, len(time_list)))
	input_bc["BoundaryConditions"]["u_x"]["SIDE_Y"]["type"] = "NEUMANN"
	input_bc["BoundaryConditions"]["u_x"]["SIDE_Y"]["value"] = list(np.repeat(0.0, len(time_list)))
	input_bc["BoundaryConditions"]["u_x"]["OUTER"]["type"] = "DIRICHLET"
	input_bc["BoundaryConditions"]["u_x"]["OUTER"]["value"] = list(np.repeat(0.0, len(time_list)))
	input_bc["BoundaryConditions"]["u_x"]["BOTTOM"]["type"] = "NEUMANN"
	input_bc["BoundaryConditions"]["u_x"]["BOTTOM"]["value"] = list(np.repeat(0.0, len(time_list)))
	input_bc["BoundaryConditions"]["u_x"]["TOP"]["type"] = "NEUMANN"
	input_bc["BoundaryConditions"]["u_x"]["TOP"]["value"] = list(np.repeat(0.0, len(time_list)))
	input_bc["BoundaryConditions"]["u_x"]["WALL"]["type"] = "NEUMANN"
	input_bc["BoundaryConditions"]["u_x"]["WALL"]["value"] = list(np.repeat(0.0, len(time_list)))

	input_bc["BoundaryConditions"]["u_y"]["SIDE_X"]["type"] = "NEUMANN"
	input_bc["BoundaryConditions"]["u_y"]["SIDE_X"]["value"] = list(np.repeat(0.0, len(time_list)))
	input_bc["BoundaryConditions"]["u_y"]["SIDE_Y"]["type"] = "DIRICHLET"
	input_bc["BoundaryConditions"]["u_y"]["SIDE_Y"]["value"] = list(np.repeat(0.0, len(time_list)))
	input_bc["BoundaryConditions"]["u_y"]["OUTER"]["type"] = "DIRICHLET"
	input_bc["BoundaryConditions"]["u_y"]["OUTER"]["value"] = list(np.repeat(0.0, len(time_list)))
	input_bc["BoundaryConditions"]["u_y"]["BOTTOM"]["type"] = "NEUMANN"
	input_bc["BoundaryConditions"]["u_y"]["BOTTOM"]["value"] = list(np.repeat(0.0, len(time_list)))
	input_bc["BoundaryConditions"]["u_y"]["TOP"]["type"] = "NEUMANN"
	input_bc["BoundaryConditions"]["u_y"]["TOP"]["value"] = list(np.repeat(0.0, len(time_list)))
	input_bc["BoundaryConditions"]["u_y"]["WALL"]["type"] = "NEUMANN"
	input_bc["BoundaryConditions"]["u_y"]["WALL"]["value"] = list(np.repeat(0.0, len(time_list)))

	input_bc["BoundaryConditions"]["u_z"]["SIDE_X"]["type"] = "NEUMANN"
	input_bc["BoundaryConditions"]["u_z"]["SIDE_X"]["value"] = list(np.repeat(0.0, len(time_list)))
	input_bc["BoundaryConditions"]["u_z"]["SIDE_Y"]["type"] = "NEUMANN"
	input_bc["BoundaryConditions"]["u_z"]["SIDE_Y"]["value"] = list(np.repeat(0.0, len(time_list)))
	input_bc["BoundaryConditions"]["u_z"]["OUTER"]["type"] = "NEUMANN"
	input_bc["BoundaryConditions"]["u_z"]["OUTER"]["value"] = list(np.repeat(0.0, len(time_list)))
	input_bc["BoundaryConditions"]["u_z"]["BOTTOM"]["type"] = "DIRICHLET"
	input_bc["BoundaryConditions"]["u_z"]["BOTTOM"]["value"] = list(np.repeat(0.0, len(time_list)))
	input_bc["BoundaryConditions"]["u_z"]["TOP"]["type"] = "NEUMANN"
	input_bc["BoundaryConditions"]["u_z"]["TOP"]["value"] = list(np.repeat(0.0, len(time_list)))
	input_bc["BoundaryConditions"]["u_z"]["WALL"]["type"] = "NEUMANN"
	input_bc["BoundaryConditions"]["u_z"]["WALL"]["value"] = pressure(time_list)

	return input_bc


def main():
	# Read input_model
	input_model = read_json("input_model.json")

	# Read input_bc
	input_bc = read_json("input_bc_fem.json")
	# input_bc = write_settings(input_bc)
	# save_json(input_bc, "input_bc_fem.json")

	# Build simulation and run
	H2CaveSimulator(input_model, input_bc)


if __name__ == "__main__":
	main()