import unittest
import os
import sys
sys.path.append(os.path.join("..", "safeincave"))
import torch as to
import numpy as np
from Utils import dotdot2, compute_C, dotdot, tensor2voigt
import dolfinx as do
from mpi4py import MPI
from ufl.tensors import ListTensor

class Test1(unittest.TestCase):
	def setUp(self):
		self.n_elems = 2
		self.nu = to.tensor([0.2, 0.3], dtype=to.float64)
		self.E = to.tensor([1e9, 2e9], dtype=to.float64)
		self.expected_C = to.tensor([[ 	[1.1111e+09, 2.7778e+08, 2.7778e+08, 0.0000e+00, 0.0000e+00,0.0000e+00],
										[2.7778e+08, 1.1111e+09, 2.7778e+08, 0.0000e+00, 0.0000e+00,0.0000e+00],
										[2.7778e+08, 2.7778e+08, 1.1111e+09, 0.0000e+00, 0.0000e+00,0.0000e+00],
										[0.0000e+00, 0.0000e+00, 0.0000e+00, 8.3333e+08, 0.0000e+00,0.0000e+00],
										[0.0000e+00, 0.0000e+00, 0.0000e+00, 0.0000e+00, 8.3333e+08,0.0000e+00],
										[0.0000e+00, 0.0000e+00, 0.0000e+00, 0.0000e+00, 0.0000e+00,8.3333e+08] ],
						          	 [  [2.6923e+09, 1.1538e+09, 1.1538e+09, 0.0000e+00, 0.0000e+00,0.0000e+00],
										[1.1538e+09, 2.6923e+09, 1.1538e+09, 0.0000e+00, 0.0000e+00,0.0000e+00],
										[1.1538e+09, 1.1538e+09, 2.6923e+09, 0.0000e+00, 0.0000e+00,0.0000e+00],
										[0.0000e+00, 0.0000e+00, 0.0000e+00, 1.5385e+09, 0.0000e+00,0.0000e+00],
										[0.0000e+00, 0.0000e+00, 0.0000e+00, 0.0000e+00, 1.5385e+09,0.0000e+00],
										[0.0000e+00, 0.0000e+00, 0.0000e+00, 0.0000e+00, 0.0000e+00,1.5385e+09] ]], dtype=to.float64)
		self.epsilon = to.tensor([[ [1., 4., 5.],
                                 	[4., 2., 6.],
                                 	[5., 6., 3.] ],
                                  [ [6., 1., 2.],
                                 	[1., 5., 3.],
                                 	[2., 3., 4.] ]], dtype=to.float64)
		self.expected_sigma = to.tensor([[	 [2.5000e+09, 3.3333e+09, 4.1666e+09],
									         [3.3333e+09, 3.3333e+09, 5.0000e+09],
									         [4.1666e+09, 5.0000e+09, 4.1666e+09]],
									        [[2.6538e+10, 1.5385e+09, 3.0770e+09],
									         [1.5385e+09, 2.5000e+10, 4.6155e+09],
									         [3.0770e+09, 4.6155e+09, 2.3461e+10]]], dtype=to.float64)

	def test_compute_C(self):
		C = compute_C(self.n_elems, self.nu, self.E)
		to.testing.assert_close(C, self.expected_C, rtol=1e-4, atol=1e-9)

	def test_dotdot2(self):
		sigma = dotdot2(self.expected_C, self.epsilon)
		to.testing.assert_close(sigma, self.expected_sigma, rtol=1e-4, atol=1e-9)


class Test2(unittest.TestCase):
	def setUp(self):
		# self.mesh = do.UnitCubeMesh(1, 1, 1)
		self.mesh = do.mesh.create_box(	MPI.COMM_WORLD,
										[np.array([0., 0., 0.]), np.array([1., 1., 1.])],
										[1, 1, 1],
										cell_type = do.mesh.CellType.tetrahedron)
		self.n_elems = self.mesh.topology.index_map(3).size_local

		self.DG_6x6 = do.fem.functionspace(self.mesh, ("DG", 0, (6, 6)))
		self.DG_3x3 = do.fem.functionspace(self.mesh, ("DG", 0, (3, 3)))

		self.C = do.fem.Function(self.DG_6x6)
		self.eps = do.fem.Function(self.DG_3x3)

		eps_template = np.array([  	[1., 4., 5.],
									[4., 2., 6.],
									[5., 6., 3.] ])
		epsilon_vec = np.tile(eps_template, (self.n_elems, 1, 1))
		self.eps.x.array[:] = epsilon_vec.flatten()

		C_template = np.array([	[1.1111e+09, 2.7778e+08, 2.7778e+08, 0.0000e+00, 0.0000e+00,0.0000e+00],
								[2.7778e+08, 1.1111e+09, 2.7778e+08, 0.0000e+00, 0.0000e+00,0.0000e+00],
								[2.7778e+08, 2.7778e+08, 1.1111e+09, 0.0000e+00, 0.0000e+00,0.0000e+00],
								[0.0000e+00, 0.0000e+00, 0.0000e+00, 8.3333e+08, 0.0000e+00,0.0000e+00],
								[0.0000e+00, 0.0000e+00, 0.0000e+00, 0.0000e+00, 8.3333e+08,0.0000e+00],
								[0.0000e+00, 0.0000e+00, 0.0000e+00, 0.0000e+00, 0.0000e+00,8.3333e+08] ])
		C_vec = np.tile(C_template, (self.n_elems, 1, 1))
		self.C.x.array[:] = C_vec.flatten()
		# self.epsilon = to.from_numpy(np.tile(A, (self.mesh.num_cells(), 1, 1)))

	def test_dotdot(self):
		sigma = dotdot(self.C, self.eps)
		self.assertIsInstance(sigma, ListTensor)
		self.assertEqual(sigma.ufl_shape, (3,3))
