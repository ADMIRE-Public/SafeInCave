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
        self.mesh = do.mesh.create_box(
            MPI.COMM_WORLD,
            [np.array([0.0, 0.0, 0.0]), np.array([1.0, 1.0, 1.0])],
            [1, 1, 1],
            cell_type=do.mesh.CellType.tetrahedron,
        )
        self.n_elems = self.mesh.topology.index_map(3).size_local

        self.DG_6x6 = do.fem.functionspace(self.mesh, ("DG", 0, (6, 6)))
        self.DG_3x3 = do.fem.functionspace(self.mesh, ("DG", 0, (3, 3)))
        self.CG1_3x1 = do.fem.functionspace(
            self.mesh, ("Lagrange", 1, (self.mesh.topology.dim,))
        )

        self.C = do.fem.Function(self.DG_6x6)
        self.eps = do.fem.Function(self.DG_3x3)
        self.u = do.fem.Function(self.CG1_3x1)

        # Populate u
        def initial_u(x):
            val = x[0] ** 2 + x[1] ** 2 + x[2] ** 2  # x² + y² + z²
            return np.array([val, val, val], dtype=np.float64)

        self.u.interpolate(initial_u)

        self.eps_tensor = to.tensor(
            [
                [[1.0, 4.0, 5.0], [4.0, 2.0, 6.0], [5.0, 6.0, 3.0]],
                [[6.0, 1.0, 2.0], [1.0, 5.0, 3.0], [2.0, 3.0, 4.0]],
                [[1.0, 4.0, 5.0], [4.0, 2.0, 6.0], [5.0, 6.0, 3.0]],
                [[6.0, 1.0, 2.0], [1.0, 5.0, 3.0], [2.0, 3.0, 4.0]],
                [[1.0, 4.0, 5.0], [4.0, 2.0, 6.0], [5.0, 6.0, 3.0]],
                [[6.0, 1.0, 2.0], [1.0, 5.0, 3.0], [2.0, 3.0, 4.0]],
            ],
            dtype=to.float64,
        )
        self.eps.x.array[:] = to.flatten(self.eps_tensor)

        self.C_tensor = to.tensor(
            [
                [
                    [1.1111e09, 2.7778e08, 2.7778e08, 0.0000e00, 0.0000e00, 0.0000e00],
                    [2.7778e08, 1.1111e09, 2.7778e08, 0.0000e00, 0.0000e00, 0.0000e00],
                    [2.7778e08, 2.7778e08, 1.1111e09, 0.0000e00, 0.0000e00, 0.0000e00],
                    [0.0000e00, 0.0000e00, 0.0000e00, 8.3333e08, 0.0000e00, 0.0000e00],
                    [0.0000e00, 0.0000e00, 0.0000e00, 0.0000e00, 8.3333e08, 0.0000e00],
                    [0.0000e00, 0.0000e00, 0.0000e00, 0.0000e00, 0.0000e00, 8.3333e08],
                ],
                [
                    [2.6923e09, 1.1538e09, 1.1538e09, 0.0000e00, 0.0000e00, 0.0000e00],
                    [1.1538e09, 2.6923e09, 1.1538e09, 0.0000e00, 0.0000e00, 0.0000e00],
                    [1.1538e09, 1.1538e09, 2.6923e09, 0.0000e00, 0.0000e00, 0.0000e00],
                    [0.0000e00, 0.0000e00, 0.0000e00, 1.5385e09, 0.0000e00, 0.0000e00],
                    [0.0000e00, 0.0000e00, 0.0000e00, 0.0000e00, 1.5385e09, 0.0000e00],
                    [0.0000e00, 0.0000e00, 0.0000e00, 0.0000e00, 0.0000e00, 1.5385e09],
                ],
                [
                    [1.1111e09, 2.7778e08, 2.7778e08, 0.0000e00, 0.0000e00, 0.0000e00],
                    [2.7778e08, 1.1111e09, 2.7778e08, 0.0000e00, 0.0000e00, 0.0000e00],
                    [2.7778e08, 2.7778e08, 1.1111e09, 0.0000e00, 0.0000e00, 0.0000e00],
                    [0.0000e00, 0.0000e00, 0.0000e00, 8.3333e08, 0.0000e00, 0.0000e00],
                    [0.0000e00, 0.0000e00, 0.0000e00, 0.0000e00, 8.3333e08, 0.0000e00],
                    [0.0000e00, 0.0000e00, 0.0000e00, 0.0000e00, 0.0000e00, 8.3333e08],
                ],
                [
                    [2.6923e09, 1.1538e09, 1.1538e09, 0.0000e00, 0.0000e00, 0.0000e00],
                    [1.1538e09, 2.6923e09, 1.1538e09, 0.0000e00, 0.0000e00, 0.0000e00],
                    [1.1538e09, 1.1538e09, 2.6923e09, 0.0000e00, 0.0000e00, 0.0000e00],
                    [0.0000e00, 0.0000e00, 0.0000e00, 1.5385e09, 0.0000e00, 0.0000e00],
                    [0.0000e00, 0.0000e00, 0.0000e00, 0.0000e00, 1.5385e09, 0.0000e00],
                    [0.0000e00, 0.0000e00, 0.0000e00, 0.0000e00, 0.0000e00, 1.5385e09],
                ],
                [
                    [1.1111e09, 2.7778e08, 2.7778e08, 0.0000e00, 0.0000e00, 0.0000e00],
                    [2.7778e08, 1.1111e09, 2.7778e08, 0.0000e00, 0.0000e00, 0.0000e00],
                    [2.7778e08, 2.7778e08, 1.1111e09, 0.0000e00, 0.0000e00, 0.0000e00],
                    [0.0000e00, 0.0000e00, 0.0000e00, 8.3333e08, 0.0000e00, 0.0000e00],
                    [0.0000e00, 0.0000e00, 0.0000e00, 0.0000e00, 8.3333e08, 0.0000e00],
                    [0.0000e00, 0.0000e00, 0.0000e00, 0.0000e00, 0.0000e00, 8.3333e08],
                ],
                [
                    [2.6923e09, 1.1538e09, 1.1538e09, 0.0000e00, 0.0000e00, 0.0000e00],
                    [1.1538e09, 2.6923e09, 1.1538e09, 0.0000e00, 0.0000e00, 0.0000e00],
                    [1.1538e09, 1.1538e09, 2.6923e09, 0.0000e00, 0.0000e00, 0.0000e00],
                    [0.0000e00, 0.0000e00, 0.0000e00, 1.5385e09, 0.0000e00, 0.0000e00],
                    [0.0000e00, 0.0000e00, 0.0000e00, 0.0000e00, 1.5385e09, 0.0000e00],
                    [0.0000e00, 0.0000e00, 0.0000e00, 0.0000e00, 0.0000e00, 1.5385e09],
                ],
            ],
            dtype=to.float64,
        )
        self.C.x.array[:] = to.flatten(self.C_tensor)

        self.expected_sigma = to.tensor(
            [
                [
                    [2.5000e09, 3.3333e09, 4.1666e09],
                    [3.3333e09, 3.3333e09, 5.0000e09],
                    [4.1666e09, 5.0000e09, 4.1666e09],
                ],
                [
                    [2.6538e10, 1.5385e09, 3.0770e09],
                    [1.5385e09, 2.5000e10, 4.6155e09],
                    [3.0770e09, 4.6155e09, 2.3461e10],
                ],
                [
                    [2.5000e09, 3.3333e09, 4.1666e09],
                    [3.3333e09, 3.3333e09, 5.0000e09],
                    [4.1666e09, 5.0000e09, 4.1666e09],
                ],
                [
                    [2.6538e10, 1.5385e09, 3.0770e09],
                    [1.5385e09, 2.5000e10, 4.6155e09],
                    [3.0770e09, 4.6155e09, 2.3461e10],
                ],
                [
                    [2.5000e09, 3.3333e09, 4.1666e09],
                    [3.3333e09, 3.3333e09, 5.0000e09],
                    [4.1666e09, 5.0000e09, 4.1666e09],
                ],
                [
                    [2.6538e10, 1.5385e09, 3.0770e09],
                    [1.5385e09, 2.5000e10, 4.6155e09],
                    [3.0770e09, 4.6155e09, 2.3461e10],
                ],
            ],
            dtype=to.float64,
        )

    def test_dotdot_ufl(self):
        sigma = ut.dotdot_ufl(self.C, self.eps)
        self.assertIsInstance(sigma, ufl.core.expr.Expr)
        self.assertEqual(sigma.ufl_shape, (3, 3))

    def test_dotdot_torch(self):
        sigma = ut.dotdot_torch(self.C_tensor, self.eps_tensor)
        self.assertIsInstance(sigma, to.Tensor)
        self.assertEqual(sigma.shape, (self.n_elems, 3, 3))
        to.testing.assert_close(sigma, self.expected_sigma, rtol=1e-4, atol=1e-9)

    def test_epsilon(self):
        eps = ut.epsilon(self.u)
        self.assertIsInstance(eps, ufl.core.expr.Expr)
        self.assertEqual(eps.ufl_shape, (3, 3))
        self.assertEqual(ufl.rank(eps), 2)

    def test_tensor2voigt(self):
        eps_voigt = ut.tensor2voigt(self.eps)
        self.assertIsInstance(eps_voigt, ufl.core.expr.Expr)
        self.assertEqual(eps_voigt.ufl_shape, (6,))

    def test_voigt2tensor(self):
        eps_voigt = ut.tensor2voigt(self.eps)
        eps = ut.voigt2tensor(eps_voigt)
        self.assertIsInstance(eps, ufl.core.expr.Expr)
        self.assertEqual(eps.ufl_shape, (3, 3))

    def test_numpy2torch(self):
        a = np.array([[1, 2, 3], [4, 5, 6]], dtype=np.float32)  # not float64 on purpose
        t = ut.numpy2torch(a)
        assert isinstance(t, to.Tensor)
        assert t.dtype == to.float64
        assert t.shape == (2, 3)
        np.testing.assert_allclose(t.numpy(), a.astype(np.float64))


class Test2(unittest.TestCase):
    def setUp(self):
        self.grid = GridHandlerGMSH("geom", os.path.join("..", "grids", "cube_regions"))

    def test_fields(self):
        def fun_1x1(x, y, z):
            return x**2 + y**2 + z**2

        def fun_2x2(x, y, z):
            return to.tensor([[x, y], [y, z]], dtype=to.float64)

        def fun_3x3(x, y, z):
            return to.tensor([[x, y, z], [y, z, x], [z, x, y]], dtype=to.float64)

        field_nodes = ut.create_field_nodes(self.grid, fun_1x1)
        field_elems = ut.create_field_elems(self.grid, fun_1x1)
        data = ut.read_json(
            os.path.join("files", "expected_values_equations", "field_nodes_elems.json")
        )
        expected_field_nodes = to.tensor(data["field_nodes"], dtype=to.float64)
        expected_field_elems = to.tensor(data["field_elems"], dtype=to.float64)
        to.testing.assert_close(field_nodes, expected_field_nodes, rtol=1e-4, atol=1e-9)
        to.testing.assert_close(field_elems, expected_field_elems, rtol=1e-4, atol=1e-9)

        field_elems = ut.create_field_elems(self.grid, fun_2x2)
        assert field_elems.shape == (self.grid.n_elems, 2, 2)

        field_elems = ut.create_field_elems(self.grid, fun_3x3)
        assert field_elems.shape == (self.grid.n_elems, 3, 3)

        field_nodes = ut.create_field_nodes(self.grid, fun_2x2)
        assert field_nodes.shape == (self.grid.n_nodes, 2, 2)

        field_nodes = ut.create_field_nodes(self.grid, fun_3x3)
        assert field_nodes.shape == (self.grid.n_nodes, 3, 3)
