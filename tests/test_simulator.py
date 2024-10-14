import unittest
import os
import sys
sys.path.append(os.path.join("..", "safeincave"))
import torch as to
import numpy as np
from Simulator import Simulator
from Utils import read_json

class Test1(unittest.TestCase):
	def setUp(self):
		self.input_file = read_json(os.path.join("files", "cube_coarse", "input_file.json"))
		self.sim = Simulator(self.input_file)
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
		self.sim.run_simulation(verbose=False)
		np.testing.assert_allclose(self.sim.eq_mom.m.elems_ie[0].alpha.numpy(), np.array(self.true_data["alpha_2"]), rtol=1e-8, atol=1e-8)

	def test_1(self):
		self.sim.run_equilibrium(verbose=False)
		np.testing.assert_allclose(self.sim.eq_mom.u.vector()[:], np.array(self.true_data["u_equilibrium"]), rtol=1e-8, atol=1e-8)

