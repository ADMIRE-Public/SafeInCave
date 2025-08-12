import unittest
import os
import sys
sys.path.append(os.path.join("..", "safeincave"))
import torch as to
import numpy as np
from Grid import GridHandlerGMSH, GridHandlerFEniCS
import dolfinx as do
from mpi4py import MPI

class Test1(unittest.TestCase):
	def setUp(self):
		self.grid = GridHandlerGMSH("geom", os.path.join("..", "grids", "cube"))
		self.expected_bNames = ["NORTH", "SOUTH", "WEST", "EAST", "BOTTOM", "TOP"]
		self.expected_btags_20 = np.array([26, 26, 26, 26, 26, 24, 24, 21, 21, 24, 26, 24, 26, 24, 24, 24, 21, 21, 24, 26])

		self.expected_dNames = ["OMEGA_A", "OMEGA_B"]
		self.expected_dtags_10b = np.array(10*[27], dtype=np.int32)
		self.expected_dtags_10f = np.array(10*[28], dtype=np.int32)

	def test_mesh_type(self):
		self.assertIsInstance(self.grid.mesh, do.mesh.Mesh)

	def test_boundaries(self):
		self.assertIsInstance(self.grid.get_boundaries(), do.mesh.MeshTags)

		boundary_names = self.grid.get_boundary_names()
		print(boundary_names)

		self.assertEqual(boundary_names[0], self.expected_bNames[0])
		self.assertEqual(boundary_names[1], self.expected_bNames[1])
		self.assertEqual(boundary_names[2], self.expected_bNames[2])
		self.assertEqual(boundary_names[3], self.expected_bNames[3])
		self.assertEqual(boundary_names[4], self.expected_bNames[4])
		self.assertEqual(boundary_names[5], self.expected_bNames[5])

		self.assertEqual(self.grid.get_boundary_tag(boundary_names[0]), 21)
		self.assertEqual(self.grid.get_boundary_tag(boundary_names[1]), 22)
		self.assertEqual(self.grid.get_boundary_tag(boundary_names[2]), 23)
		self.assertEqual(self.grid.get_boundary_tag(boundary_names[3]), 24)
		self.assertEqual(self.grid.get_boundary_tag(boundary_names[4]), 25)
		self.assertEqual(self.grid.get_boundary_tag(boundary_names[5]), 26)

		boundaries = self.grid.get_boundaries()
		np.testing.assert_allclose(boundaries.values[:20], self.expected_btags_20)
		

	def test_subdomains(self):
		subdomains = self.grid.get_subdomains()
		self.assertIsInstance(subdomains, do.mesh.MeshTags)

		subdomain_names = self.grid.get_subdomain_names()
		self.assertEqual(subdomain_names[0], self.expected_dNames[0])
		self.assertEqual(subdomain_names[1], self.expected_dNames[1])
		self.assertEqual(self.grid.get_subdomain_tag(subdomain_names[0]), 27)
		self.assertEqual(self.grid.get_subdomain_tag(subdomain_names[1]), 28)
		
		np.testing.assert_allclose(subdomains.values[:10], self.expected_dtags_10f)
		np.testing.assert_allclose(subdomains.values[-10:], self.expected_dtags_10b)

# class Test2(unittest.TestCase):
# 	def setUp(self):
# 		mesh = do.mesh.create_box(	MPI.COMM_WORLD,
# 									[np.array([0., 0., 0.]), np.array([1., 1., 1.])],
# 									[5, 5, 5],
# 									cell_type = do.mesh.CellType.tetrahedron)
# 		self.grid = GridHandlerFEniCS(mesh)
# 		self.expected_bNames = ["WEST", "EAST", "SOUTH", "NORTH", "BOTTOM", "TOP"]
# 		self.expected_btags_20 = np.array([5, 3, 5, 3, 5, 3, 2, 2, 2, 5, 2, 2, 2, 2, 3, 2, 5, 5, 5, 3])
# 		self.expected_dNames = ["BODY"]
# 		self.expected_dtags_10f = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
# 		self.expected_dtags_10b = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0])

# 	def test_mesh_type(self):
# 		self.assertIsInstance(self.grid.mesh, do.mesh.Mesh)

# 	def test_boundaries(self):
# 		for i, b_name in enumerate(self.grid.get_boundary_names()):
# 			self.assertEqual(b_name, self.expected_bNames[i])
# 			self.assertEqual(self.grid.get_boundary_tags(b_name), i+1)
# 		bounds = self.grid.get_boundaries()
# 		np.testing.assert_allclose(bounds.values[:20], self.expected_btags_20)

# 	def test_subdomains(self):
# 		subdomains = self.grid.get_subdomains()
# 		self.assertIsInstance(subdomains, do.mesh.MeshTags)

# 		for i, domain_name in enumerate(self.grid.get_subdomain_names()):
# 			self.assertEqual(domain_name, self.expected_dNames[i])
# 			self.assertEqual(self.grid.get_subdomain_tag(domain_name), i+1)
