import unittest
import os
import sys
sys.path.append(os.path.join("..", "safeincave"))
import torch as to
import numpy as np
from Equations import LinearMomentum
from ConstitutiveModel import ConstitutiveModel
from Grid import GridHandlerGMSH
from Utils import read_json
import dolfin as do

class Test1(unittest.TestCase):
	def setUp(self):
		self.grid = GridHandlerGMSH("geom", os.path.join("..", "grids", "cube"))
		self.n_elems = self.grid.mesh.num_cells()
		self.input_file = read_json(os.path.join("files", "input_file.json"))
		self.theta = self.input_file["time_settings"]["theta"]

		self.eq = LinearMomentum(self.grid, self.theta, self.input_file)

		# # Define constitutive model
		# self.m = ConstitutiveModel(self.n_elems, input_file["constitutive_model"])

	def test_0(self):
		self.eq.initialize()
		
		# t = 0
		# dt = 3600.
		# self.eq.solve(t, dt)