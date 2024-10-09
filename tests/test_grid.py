import unittest
import os
import sys
sys.path.append(os.path.join("..", "safeincave"))
import torch as to
import numpy as np
from Grid import GridHandlerGMSH, GridHandlerFEniCS
import dolfin as do

class Test1(unittest.TestCase):
	def setUp(self):
		self.grid = GridHandlerGMSH("geom", os.path.join("..", "grids", "cube"))
		self.expected_bNames = ["NORTH", "SOUTH", "WEST", "EAST", "BOTTOM", "TOP"]
		self.expected_btags_20 = np.array([5, 3, 5, 2, 3, 2, 0, 0, 0, 5, 4, 5, 2, 4, 2, 0, 0, 0, 4, 5])
		self.expected_dNames = ["OMEGA_A", "OMEGA_B"]
		self.expected_dtags_10f = np.array([2, 2, 2, 2, 2, 2, 2, 2, 2, 2])
		self.expected_dtags_10b = np.array([1, 1, 1, 1, 1, 1, 1, 1, 1, 1])

	def test_mesh_type(self):
		self.assertIsInstance(self.grid.mesh, do.Mesh)

	def test_boundaries(self):
		self.assertIsInstance(self.grid.get_boundaries(), do.cpp.mesh.MeshFunctionSizet)
		for i, b_name in enumerate(self.grid.get_boundary_names()):
			self.assertEqual(b_name, self.expected_bNames[i])
			self.assertEqual(self.grid.get_boundary_tags(b_name), i+1)
		bounds = self.grid.get_boundaries()
		np.testing.assert_allclose(bounds.array()[:20], self.expected_btags_20)

	def test_subdomains(self):
		self.assertIsInstance(self.grid.get_subdomains(), do.cpp.mesh.MeshFunctionSizet)
		for i, domain_name in enumerate(self.grid.get_subdomain_names()):
			self.assertEqual(domain_name, self.expected_dNames[i])
			self.assertEqual(self.grid.get_subdomain_tags(domain_name), i+1)
		domains = self.grid.get_subdomains()
		np.testing.assert_allclose(domains.array()[:10], self.expected_dtags_10f)
		np.testing.assert_allclose(domains.array()[-10:], self.expected_dtags_10b)


class Test2(unittest.TestCase):
	def setUp(self):
		mesh = do.UnitCubeMesh(5, 5, 5)
		self.grid = GridHandlerFEniCS(mesh)
		self.expected_bNames = ["WEST", "EAST", "SOUTH", "NORTH", "BOTTOM", "TOP"]
		self.expected_btags_20 = np.array([5, 3, 0, 5, 1, 0, 0, 3, 1, 0, 0, 0, 5, 3, 0, 5, 0, 0, 0, 3])
		self.expected_dNames = ["BODY"]
		self.expected_dtags_10f = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
		self.expected_dtags_10b = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0])

	def test_mesh_type(self):
		self.assertIsInstance(self.grid.mesh, do.Mesh)

	def test_boundaries(self):
		for i, b_name in enumerate(self.grid.get_boundary_names()):
			self.assertEqual(b_name, self.expected_bNames[i])
			self.assertEqual(self.grid.get_boundary_tags(b_name), i+1)
		bounds = self.grid.get_boundaries()
		np.testing.assert_allclose(bounds.array()[:20], self.expected_btags_20)

	def test_subdomains(self):
		self.assertIsInstance(self.grid.get_subdomains(), do.cpp.mesh.MeshFunctionSizet)
		for i, domain_name in enumerate(self.grid.get_subdomain_names()):
			self.assertEqual(domain_name, self.expected_dNames[i])
			self.assertEqual(self.grid.get_subdomain_tags(domain_name), i+1)
		domains = self.grid.get_subdomains()


