import unittest
import os
import sys
sys.path.append(os.path.join("..", "safeincave"))
import torch as to
import numpy as np
from ConstitutiveModel import ConstitutiveModel
from Elements import Spring, Viscoelastic, DislocationCreep, ViscoplasticDesai
import dolfin as do

class Test1(unittest.TestCase):
	def setUp(self):
		self.n_elems = 1
		input_cm = {
			"Elastic": {
				"Spring0": {
					"type": "Spring",
					"active": True,
					"parameters": {
						"E":  list(102e9*np.ones(self.n_elems)),
						"nu": list(0.3*np.ones(self.n_elems))
					}
				},
				"Spring1": {
					"type": "Spring",
					"active": True,
					"parameters": {
						"E":  list(72e9*np.ones(self.n_elems)),
						"nu": list(0.4*np.ones(self.n_elems))
					}
				}
			},
			"Viscoelastic": {
				"KelvinVoigt1": {
					"type": "KelvinVoigt",
					"active": True,
					"parameters": {
						"E":   list(10e9*np.ones(self.n_elems)),
						"nu":  list(0.32*np.ones(self.n_elems)),
						"eta": list(105e11*np.ones(self.n_elems))
					}
				},
				"KelvinVoigt2": {
					"type": "KelvinVoigt",
					"active": True,
					"parameters": {
						"E":   list(6e9*np.ones(self.n_elems)),
						"nu":  list(0.22*np.ones(self.n_elems)),
						"eta": list(15e11*np.ones(self.n_elems))
					}
				}
			},
			"Inelastic": {
				"ViscPlastDesai": {
					"type": "ViscoplasticDesai",
					"active": True,
					"parameters": {
						"mu_1": 	list(5.3665857009859815e-11*np.ones(self.n_elems)),
						"N_1": 		list(3.1*np.ones(self.n_elems)),
						"n": 		list(3.0*np.ones(self.n_elems)),
						"a_1":		list(1.965018496922832e-05*np.ones(self.n_elems)),
						"eta": 		list(0.8275682807874163*np.ones(self.n_elems)),
						"beta_1": 	list(0.0048*np.ones(self.n_elems)),
						"beta": 	list(0.995*np.ones(self.n_elems)),
						"m": 		list(-0.5*np.ones(self.n_elems)),
						"gamma": 	list(0.095*np.ones(self.n_elems)),
						"alpha_0": 	list(0.0022*np.ones(self.n_elems)),
						"k_v": 		list(0.0*np.ones(self.n_elems)),
						"sigma_t": 	list(5.0*np.ones(self.n_elems))
					}
				},
				"DisCreep": {
					"type": "DislocationCreep",
					"active": True,
					"parameters": {
						"A": list(1.9e-20*np.ones(self.n_elems)),
						"n": list(3.0*np.ones(self.n_elems)),
						"T": list(298*np.ones(self.n_elems)),
						"Q": list(51600*np.ones(self.n_elems)),
						"R": list(8.32*np.ones(self.n_elems))
					}
				}
			}
		}

		self.cm = ConstitutiveModel(self.n_elems, input_cm)

		self.true_C0_inv = to.tensor([[	[ 2.3693e-11, -8.4967e-12, -8.4967e-12,  0.0000e+00,  0.0000e+00, 0.0000e+00],
										[-8.4967e-12,  2.3693e-11, -8.4967e-12,  0.0000e+00,  0.0000e+00, 0.0000e+00],
										[-8.4967e-12, -8.4967e-12,  2.3693e-11,  0.0000e+00,  0.0000e+00, 0.0000e+00],
										[ 0.0000e+00,  0.0000e+00,  0.0000e+00,  3.2190e-11,  0.0000e+00, 0.0000e+00],
										[ 0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  3.2190e-11, 0.0000e+00],
										[ 0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00, 3.2190e-11]]], dtype=to.float64)

	def test_C0_inv(self):
		to.testing.assert_close(self.cm.C0_inv, self.true_C0_inv, rtol=1e-4, atol=1e-9)

	def test_elems(self):
		self.assertEqual(len(self.cm.elems_e), 2)
		self.assertIsInstance(self.cm.elems_e[0], Spring)
		self.assertIsInstance(self.cm.elems_e[1], Spring)

		self.assertEqual(len(self.cm.elems_ve), 2)
		self.assertIsInstance(self.cm.elems_ve[0], Viscoelastic)
		self.assertIsInstance(self.cm.elems_ve[1], Viscoelastic)

		self.assertEqual(len(self.cm.elems_ie), 2)
		self.assertIsInstance(self.cm.elems_ie[0], ViscoplasticDesai)
		self.assertIsInstance(self.cm.elems_ie[1], DislocationCreep)

	def test_n_elems(self):
		self.assertEqual(self.n_elems, self.cm.n_elems)

