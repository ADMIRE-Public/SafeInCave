import unittest
import os
import sys
sys.path.append(os.path.join("..", "safeincave"))
import numpy as np
from Equations import LinearMomentum
from Gridx import GridHandlerGMSH
from Utils import read_json

class Test1(unittest.TestCase):

	def setUp(self):
		self.grid = GridHandlerGMSH("geom", os.path.join("files", "cube_coarse"))
		self.input_file = read_json(os.path.join("files", "cube_coarse", "input_file.json"))
		self.theta = self.input_file["time_settings"]["theta"]
		self.eq = LinearMomentum(self.grid, self.theta, self.input_file)
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

	def test_equilibrium_0(self):
		self.eq.solve_equilibrium(verbose=False, save_results=False)
		np.testing.assert_allclose(self.eq.u.vector()[:], np.array(self.true_data["u_equilibrium"]), rtol=1e-8, atol=1e-8)

	def test_equilibrium_1(self):
		self.eq.initialize(solve_equilibrium=True, verbose=False)
		np.testing.assert_allclose(self.eq.u.vector()[:], np.array(self.true_data["u_equilibrium"]), rtol=1e-8, atol=1e-8)

	def test_full(self):
		self.eq.initialize(solve_equilibrium=False, verbose=False)

		self.assertEqual(len(self.eq.integral_neumann), 3)
		self.assertEqual(len(self.eq.bcs), 3)

		for bc_neumann in self.eq.integral_neumann:
			self.assertIsInstance(bc_neumann, ufl.form.Form)
		for bc_dirichlet in self.eq.bcs:
			self.assertIsInstance(bc_dirichlet, do.fem.dirichletbc.DirichletBC)

		np.testing.assert_allclose(self.eq.eps_tot.vector()[:], self.true_data["eps_tot_0"], rtol=1e-8, atol=1e-8)
		np.testing.assert_allclose(self.eq.u.vector()[:], self.true_data["u_0"], rtol=1e-8, atol=1e-8)
		np.testing.assert_allclose(self.eq.m.elems_ie[0].alpha.numpy(), self.true_data["alpha_0"], rtol=1e-8, atol=1e-8)
		
		t = 0
		dt = 3600.
		self.eq.solve(t, dt)

		np.testing.assert_allclose(self.eq.eps_tot.vector()[:], self.true_data["eps_tot_1"], rtol=1e-8, atol=1e-8)
		np.testing.assert_allclose(self.eq.u.vector()[:], self.true_data["u_1"], rtol=1e-8, atol=1e-8)
		np.testing.assert_allclose(self.eq.m.elems_ie[0].alpha.numpy(), self.true_data["alpha_1"], rtol=1e-8, atol=1e-8)

