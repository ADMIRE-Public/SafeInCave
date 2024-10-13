import unittest
import os
import sys
sys.path.append(os.path.join("..", "..", "..", "safeincave"))
import torch as to
import numpy as np
from Equations import LinearMomentum
from ConstitutiveModel import ConstitutiveModel
from Grid import GridHandlerGMSH
from Utils import read_json, save_json
import dolfin as do
import ufl as ufl
from Utils import numpy2torch

def main():
	grid = GridHandlerGMSH("geom", os.path.join("..", "cube_coarse"))
	n_elems = grid.mesh.num_cells()
	input_file = read_json(os.path.join("..", "cube_coarse", "input_file.json"))
	theta = input_file["time_settings"]["theta"]
	eq = LinearMomentum(grid, theta, input_file)

	eq.initialize()

	data = {}
	data["u_0"] = list(eq.u.vector()[:])
	data["eps_tot_0"] = list(eq.eps_tot.vector()[:])
	data["alpha_0"] = eq.m.elems_ie[0].alpha.tolist()
		
	t = 0
	dt = 3600.
	eq.solve(t, dt)

	data["u_1"] = list(eq.u.vector()[:])
	data["eps_tot_1"] = list(eq.eps_tot.vector()[:])
	data["alpha_1"] = eq.m.elems_ie[0].alpha.tolist()

	data["alpha_2"] = [	 0.00537296, 0.00537268, 0.00537241, 0.00537275, 0.00537117, 0.00537477,
						 0.00537461, 0.00537025, 0.00537573, 0.00536915, 0.00537585, 0.00536903,
						 0.00537496, 0.00537498, 0.00537523, 0.00536961, 0.00537512, 0.00536995,
						 0.00537045, 0.00537026, 0.00537073, 0.00537447, 0.00537423, 0.00537137,
						 0.00537221, 0.00537219, 0.00537481, 0.0053726,  0.00537006, 0.00537261,
						 0.00537481, 0.00536994, 0.00537569, 0.00537575, 0.00536857, 0.00536851,
						 0.00537505, 0.00536931, 0.00537491, 0.00536944, 0.005375,   0.00536978,
						 0.0053697,  0.00537494, 0.00537097, 0.00537018, 0.00537411, 0.00537435]


	save_json(data, "expected_values.json")


if __name__ == '__main__':
	main()