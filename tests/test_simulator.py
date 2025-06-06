import unittest
import os
import sys
sys.path.append(os.path.join("..", "safeincave"))
import torch as to
import numpy as np
from Simulator import Simulator
import Simulator as s1
from Utils import read_json

from Gridx import GridHandlerGMSH

class Test1(unittest.TestCase):

	def setUp(self):
		input_file1 = read_json(os.path.join("files", "cube_coarse", "input_file.json"))
		input_file1["grid"]["path"] = os.path.join("files", "cube_coarse")
		self.sim = s1.Simulator(input_file1)
		self.load_expected_values()

	def load_expected_values(self):
		self.true_data = read_json(os.path.join("files", "expected_values_equations", "expected_values.json"))
		self.true_data["eps_tot_0"] = np.array(self.true_data["eps_tot_0"])
		self.true_data["u_0"] = np.array(self.true_data["u_0"])
		self.true_data["alpha_0"] = np.array(self.true_data["alpha_0"])
		self.true_data["eps_tot_1"] = np.array(self.true_data["eps_tot_1"])
		self.true_data["u_1"] = np.array(self.true_data["u_1"])
		self.true_data["alpha_1"] = np.array(self.true_data["alpha_1"])
		self.true_data["u_equilibrium"] = np.array(self.true_data["u_equilibrium"])

	def test_0(self):
		self.sim.run_simulation(solve_equilibrium=False, verbose=True)

		# print()
		# print(self.sim.eq_mom.m.elems_ie[0].alpha.numpy())
		# print()
		# print(self.sim.eq_mom.u.x.array)
		# print()
		np.testing.assert_allclose(self.sim.eq_mom.m.elems_ie[0].alpha.numpy(), np.array(self.true_data["alpha_2"]), rtol=1e-8, atol=1e-8)
		np.testing.assert_allclose(self.sim.eq_mom.u.x.array, np.array(self.true_data["u_operation"]), rtol=1e-8, atol=1e-8)
