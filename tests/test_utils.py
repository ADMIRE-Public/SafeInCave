# import unittest
# import os
# import sys
# sys.path.append(os.path.join("..", "safeincave"))
# import torch as to
# import numpy as np
# import numpy.typing as NDArray
# import Utils as ut
# import dolfinx as do
# import ufl
# from mpi4py import MPI
# from Grid import GridHandlerGMSH

import os
from safeincave import GridHandlerGMSH
import safeincave.Utils as ut
import torch as to
import numpy as np
import unittest

from mpi4py import MPI
import dolfinx as do
import ufl


class Test1(unittest.TestCase):
	def setUp(self):
		# self.mesh = do.UnitCubeMesh(1, 1, 1)
		self.mesh = do.mesh.create_box(	MPI.COMM_WORLD,
										[np.array([0., 0., 0.]), np.array([1., 1., 1.])],
										[1, 1, 1],
										cell_type = do.mesh.CellType.tetrahedron)
		self.n_elems = self.mesh.topology.index_map(3).size_local

		self.DG_6x6 = do.fem.functionspace(self.mesh, ("DG", 0, (6, 6)))
		self.DG_3x3 = do.fem.functionspace(self.mesh, ("DG", 0, (3, 3)))
		self.CG1_3x1 = do.fem.functionspace(self.mesh, ("Lagrange", 1, (self.mesh.topology.dim, )))

		self.C = do.fem.Function(self.DG_6x6)
		self.eps = do.fem.Function(self.DG_3x3)
		self.u = do.fem.Function(self.CG1_3x1)

		# Populate u
		def initial_u(x):
		    val = x[0]**2 + x[1]**2 + x[2]**2  # x² + y² + z²
		    return np.array([val, val, val], dtype=np.float64)
		self.u.interpolate(initial_u)

		self.eps_tensor = to.tensor([[ [1., 4., 5.],
                            	   	   [4., 2., 6.],
                            	   	   [5., 6., 3.] ],
                            	     [ [6., 1., 2.],
                            	   	   [1., 5., 3.],
                            	   	   [2., 3., 4.] ],
                            	     [ [1., 4., 5.],
                            	   	   [4., 2., 6.],
                            	   	   [5., 6., 3.] ],
                            	     [ [6., 1., 2.],
                            	   	   [1., 5., 3.],
                            	   	   [2., 3., 4.] ],
                            	     [ [1., 4., 5.],
                            	   	   [4., 2., 6.],
                            	   	   [5., 6., 3.] ],
                            	     [ [6., 1., 2.],
                            	   	   [1., 5., 3.],
                            	   	   [2., 3., 4.] ]
                            	    ], dtype=to.float64)
		self.eps.x.array[:] = to.flatten(self.eps_tensor)

		self.C_tensor = to.tensor([[ 	[1.1111e+09, 2.7778e+08, 2.7778e+08, 0.0000e+00, 0.0000e+00,0.0000e+00],
										[2.7778e+08, 1.1111e+09, 2.7778e+08, 0.0000e+00, 0.0000e+00,0.0000e+00],
										[2.7778e+08, 2.7778e+08, 1.1111e+09, 0.0000e+00, 0.0000e+00,0.0000e+00],
										[0.0000e+00, 0.0000e+00, 0.0000e+00, 8.3333e+08, 0.0000e+00,0.0000e+00],
										[0.0000e+00, 0.0000e+00, 0.0000e+00, 0.0000e+00, 8.3333e+08,0.0000e+00],
										[0.0000e+00, 0.0000e+00, 0.0000e+00, 0.0000e+00, 0.0000e+00,8.3333e+08] ],
					         	   [   [2.6923e+09, 1.1538e+09, 1.1538e+09, 0.0000e+00, 0.0000e+00,0.0000e+00],
										[1.1538e+09, 2.6923e+09, 1.1538e+09, 0.0000e+00, 0.0000e+00,0.0000e+00],
										[1.1538e+09, 1.1538e+09, 2.6923e+09, 0.0000e+00, 0.0000e+00,0.0000e+00],
										[0.0000e+00, 0.0000e+00, 0.0000e+00, 1.5385e+09, 0.0000e+00,0.0000e+00],
										[0.0000e+00, 0.0000e+00, 0.0000e+00, 0.0000e+00, 1.5385e+09,0.0000e+00],
										[0.0000e+00, 0.0000e+00, 0.0000e+00, 0.0000e+00, 0.0000e+00,1.5385e+09] ],
									[ 	[1.1111e+09, 2.7778e+08, 2.7778e+08, 0.0000e+00, 0.0000e+00,0.0000e+00],
										[2.7778e+08, 1.1111e+09, 2.7778e+08, 0.0000e+00, 0.0000e+00,0.0000e+00],
										[2.7778e+08, 2.7778e+08, 1.1111e+09, 0.0000e+00, 0.0000e+00,0.0000e+00],
										[0.0000e+00, 0.0000e+00, 0.0000e+00, 8.3333e+08, 0.0000e+00,0.0000e+00],
										[0.0000e+00, 0.0000e+00, 0.0000e+00, 0.0000e+00, 8.3333e+08,0.0000e+00],
										[0.0000e+00, 0.0000e+00, 0.0000e+00, 0.0000e+00, 0.0000e+00,8.3333e+08] ],
					         	 	[   [2.6923e+09, 1.1538e+09, 1.1538e+09, 0.0000e+00, 0.0000e+00,0.0000e+00],
										[1.1538e+09, 2.6923e+09, 1.1538e+09, 0.0000e+00, 0.0000e+00,0.0000e+00],
										[1.1538e+09, 1.1538e+09, 2.6923e+09, 0.0000e+00, 0.0000e+00,0.0000e+00],
										[0.0000e+00, 0.0000e+00, 0.0000e+00, 1.5385e+09, 0.0000e+00,0.0000e+00],
										[0.0000e+00, 0.0000e+00, 0.0000e+00, 0.0000e+00, 1.5385e+09,0.0000e+00],
										[0.0000e+00, 0.0000e+00, 0.0000e+00, 0.0000e+00, 0.0000e+00,1.5385e+09] ],
									[ 	[1.1111e+09, 2.7778e+08, 2.7778e+08, 0.0000e+00, 0.0000e+00,0.0000e+00],
										[2.7778e+08, 1.1111e+09, 2.7778e+08, 0.0000e+00, 0.0000e+00,0.0000e+00],
										[2.7778e+08, 2.7778e+08, 1.1111e+09, 0.0000e+00, 0.0000e+00,0.0000e+00],
										[0.0000e+00, 0.0000e+00, 0.0000e+00, 8.3333e+08, 0.0000e+00,0.0000e+00],
										[0.0000e+00, 0.0000e+00, 0.0000e+00, 0.0000e+00, 8.3333e+08,0.0000e+00],
										[0.0000e+00, 0.0000e+00, 0.0000e+00, 0.0000e+00, 0.0000e+00,8.3333e+08] ],
					         	 	[   [2.6923e+09, 1.1538e+09, 1.1538e+09, 0.0000e+00, 0.0000e+00,0.0000e+00],
										[1.1538e+09, 2.6923e+09, 1.1538e+09, 0.0000e+00, 0.0000e+00,0.0000e+00],
										[1.1538e+09, 1.1538e+09, 2.6923e+09, 0.0000e+00, 0.0000e+00,0.0000e+00],
										[0.0000e+00, 0.0000e+00, 0.0000e+00, 1.5385e+09, 0.0000e+00,0.0000e+00],
										[0.0000e+00, 0.0000e+00, 0.0000e+00, 0.0000e+00, 1.5385e+09,0.0000e+00],
										[0.0000e+00, 0.0000e+00, 0.0000e+00, 0.0000e+00, 0.0000e+00,1.5385e+09] ]
									], dtype=to.float64)
		self.C.x.array[:] = to.flatten(self.C_tensor)

		self.expected_sigma = to.tensor([   [[2.5000e+09, 3.3333e+09, 4.1666e+09],
								             [3.3333e+09, 3.3333e+09, 5.0000e+09],
								             [4.1666e+09, 5.0000e+09, 4.1666e+09]],
								            [[2.6538e+10, 1.5385e+09, 3.0770e+09],
								             [1.5385e+09, 2.5000e+10, 4.6155e+09],
								             [3.0770e+09, 4.6155e+09, 2.3461e+10]],
								            [[2.5000e+09, 3.3333e+09, 4.1666e+09],
								             [3.3333e+09, 3.3333e+09, 5.0000e+09],
								             [4.1666e+09, 5.0000e+09, 4.1666e+09]],
								            [[2.6538e+10, 1.5385e+09, 3.0770e+09],
								             [1.5385e+09, 2.5000e+10, 4.6155e+09],
								             [3.0770e+09, 4.6155e+09, 2.3461e+10]],
								            [[2.5000e+09, 3.3333e+09, 4.1666e+09],
								             [3.3333e+09, 3.3333e+09, 5.0000e+09],
								             [4.1666e+09, 5.0000e+09, 4.1666e+09]],
								            [[2.6538e+10, 1.5385e+09, 3.0770e+09],
								             [1.5385e+09, 2.5000e+10, 4.6155e+09],
								             [3.0770e+09, 4.6155e+09, 2.3461e+10]]
							             ], dtype=to.float64)


	def test_dotdot_ufl(self):
		sigma = ut.dotdot_ufl(self.C, self.eps)
		self.assertIsInstance(sigma, ufl.core.expr.Expr)
		self.assertEqual(sigma.ufl_shape, (3,3))

	def test_dotdot_torch(self):
		sigma = ut.dotdot_torch(self.C_tensor, self.eps_tensor)
		self.assertIsInstance(sigma, to.Tensor)
		self.assertEqual(sigma.shape, (self.n_elems,3,3))
		to.testing.assert_close(sigma, self.expected_sigma, rtol=1e-4, atol=1e-9)

	def test_epsilon(self):
		eps = ut.epsilon(self.u)
		self.assertIsInstance(eps, ufl.core.expr.Expr)
		self.assertEqual(eps.ufl_shape, (3,3))
		self.assertEqual(ufl.rank(eps), 2)

	def test_tensor2voigt(self):
		eps_voigt = ut.tensor2voigt(self.eps)
		self.assertIsInstance(eps_voigt, ufl.core.expr.Expr)
		self.assertEqual(eps_voigt.ufl_shape, (6,))

	def test_voigt2tensor(self):
		eps_voigt = ut.tensor2voigt(self.eps)
		eps = ut.voigt2tensor(eps_voigt)
		self.assertIsInstance(eps, ufl.core.expr.Expr)
		self.assertEqual(eps.ufl_shape, (3,3))

	def test_numpy2torch(self):
	    a = np.array([[1, 2, 3], [4, 5, 6]], dtype=np.float32)  # not float64 on purpose
	    t = ut.numpy2torch(a)
	    assert isinstance(t, to.Tensor)
	    assert t.dtype == to.float64
	    assert t.shape == (2, 3)
	    np.testing.assert_allclose(t.numpy(), a.astype(np.float64))


class Test2(unittest.TestCase):
	def setUp(self):
		self.grid = GridHandlerGMSH("geom", os.path.join("..", "grids", "cube"))

	def test_fields(self):
		fun = lambda x,y,z: x**2 + y**2 + z**2
		field_nodes = ut.create_field_nodes(self.grid, fun)
		field_elems = ut.create_field_elems(self.grid, fun)
		data = ut.read_json(os.path.join("files", "expected_values_equations", "field_nodes_elems.json"))
		expected_field_nodes = to.tensor(data["field_nodes"], dtype=to.float64)
		expected_field_elems = to.tensor(data["field_elems"], dtype=to.float64)

		to.testing.assert_close(field_nodes, expected_field_nodes, rtol=1e-4, atol=1e-9)
		to.testing.assert_close(field_elems, expected_field_elems, rtol=1e-4, atol=1e-9)
