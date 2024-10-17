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

	data["u_equilibrium"] = [-1.55961103e-04, -1.56019017e-04, -1.71595569e-04, -1.55758388e-04,
							 -1.55723575e-04,  6.19815414e-17, -7.79219424e-05, -1.55837295e-04,
							 -8.60224137e-05, -1.55852532e-04, -1.16886505e-04, -8.60167910e-05,
							 -1.55971835e-04, -7.79898182e-05, -1.71694761e-04, -7.79756695e-05,
							 -1.16992317e-04, -1.71747930e-04,  3.58340330e-17, -1.56012560e-04,
							 -1.71665906e-04,  3.58340330e-17, -1.55702987e-04,  6.19815414e-17,
							  3.58340330e-17, -1.16877730e-04, -8.60633951e-05,  6.33216234e-18,
							 -7.79849661e-05, -1.71771063e-04, -7.79220932e-05, -7.79166225e-05,
							 -8.60916280e-05, -1.55747454e-04, -7.78631946e-05,  6.19815414e-17,
							 -7.78788398e-05, -1.16782696e-04,  6.19815414e-17, -7.79785887e-05,
							 -3.89897381e-05, -1.71794617e-04, -1.55974529e-04,  3.18383255e-17,
							 -1.71709918e-04, -1.55845407e-04, -3.89608475e-05, -8.60395187e-05,
							  3.58340330e-17, -7.78593401e-05,  6.19815414e-17,  6.33216234e-18,
							 -3.89592211e-05, -8.60868523e-05,  6.33216234e-18,  3.18383255e-17,
							 -1.71787692e-04, -7.79195523e-05,  3.18383255e-17, -8.60965664e-05,
							 -1.55736980e-04,  3.18383255e-17,  2.24066759e-17, -7.78719207e-05,
							 -3.89324158e-05,  6.19815414e-17,  6.33216234e-18,  3.18383255e-17,
							  2.24066759e-17]


	save_json(data, "expected_values.json")


if __name__ == '__main__':
	main()