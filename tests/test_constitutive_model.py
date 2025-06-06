import unittest
import os
import sys
sys.path.append(os.path.join("..", "safeincave"))
import torch as to
import numpy as np
from Grid import GridHandlerFEniCS
from ConstitutiveModel import ConstitutiveModel
from Elements import Spring, Viscoelastic, DislocationCreep, ViscoplasticDesai
import dolfinx as do
from mpi4py import MPI

class Test1(unittest.TestCase):
	def setUp(self):
		mesh = do.mesh.create_box(	MPI.COMM_WORLD,
									[np.array([0., 0., 0.]), np.array([1., 1., 1.])],
									[1, 1, 1],
									cell_type = do.mesh.CellType.tetrahedron)
		self.grid = GridHandlerFEniCS(mesh)
		self.n_elems = self.grid.n_elems
		input_cm = {
			"elastic": {
				"Spring0": {
					"type": "Spring",
					"active": True,
					"parameters": {
						"E":  102e9,
						"nu": 0.3
					}
				},
				"Spring1": {
					"type": "Spring",
					"active": True,
					"parameters": {
						"E":  72e9,
						"nu": 0.4
					}
				}
			},
			"viscoelastic": {
				"KelvinVoigt1": {
					"type": "KelvinVoigt",
					"active": True,
					"parameters": {
						"E":   10e9,
						"nu":  0.32,
						"eta": 105e11
					}
				},
				"KelvinVoigt2": {
					"type": "KelvinVoigt",
					"active": True,
					"parameters": {
						"E":   6e9,
						"nu":  0.22,
						"eta": 15e11
					}
				}
			},
			"inelastic": {
				"ViscPlastDesai": {
					"type": "ViscoplasticDesai",
					"active": True,
					"parameters": {
						"mu_1": 	5.3665857009859815e-11,
						"N_1": 		3.1,
						"n": 		3.0,
						"a_1":		1.965018496922832e-05,
						"eta": 		0.8275682807874163,
						"beta_1": 	0.0048,
						"beta": 	0.995,
						"m": 		-0.5,
						"gamma": 	0.095,
						"alpha_0": 	0.0022,
						"k_v": 		0.0,
						"sigma_t": 	5.0
					}
				},
				"DisCreep": {
					"type": "DislocationCreep",
					"active": True,
					"parameters": {
						"A": 1.9e-20,
						"n": 3.0,
						"T": 298,
						"Q": 51600,
						"R": 8.32
					}
				}
			}
		}

		self.cm = ConstitutiveModel(self.grid, input_cm)

		self.true_C0_inv = to.tensor([ 	[ 2.3693e-11, -8.4967e-12, -8.4967e-12,  0.0000e+00,  0.0000e+00, 0.0000e+00],
										[-8.4967e-12,  2.3693e-11, -8.4967e-12,  0.0000e+00,  0.0000e+00, 0.0000e+00],
										[-8.4967e-12, -8.4967e-12,  2.3693e-11,  0.0000e+00,  0.0000e+00, 0.0000e+00],
										[ 0.0000e+00,  0.0000e+00,  0.0000e+00,  3.2190e-11,  0.0000e+00, 0.0000e+00],
										[ 0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  3.2190e-11, 0.0000e+00],
										[ 0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00, 3.2190e-11]], dtype=to.float64)

	def test_C0_inv(self):
		# print(self.cm.C0_inv)
		to.testing.assert_close(self.cm.C0_inv[0], self.true_C0_inv, rtol=1e-4, atol=1e-9)

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
		self.assertEqual(self.n_elems, self.cm.grid.n_elems)

