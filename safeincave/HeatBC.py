# Copyright 2025 The safeincave community.
#
# This file is part of safeincave.
#
# Licensed under the GNU GENERAL PUBLIC LICENSE, Version 3 (the "License"); you may not
# use this file except in compliance with the License.  You may obtain a copy
# of the License at
#
#     https://spdx.org/licenses/GPL-3.0-or-later.html
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS, WITHOUT
# WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.  See the
# License for the specific language governing permissions and limitations under
# the License.
from abc import ABC
import numpy as np
import dolfinx as do
import ufl

class GeneralBC(ABC):
	def __init__(self, boundary_name : str, values : list, time_values : list):
		self.boundary_name = boundary_name
		self.values = values
		self.time_values = time_values
		self.type = None


class DirichletBC(GeneralBC):
	def __init__(self, boundary_name : str, values : list, time_values : list):
		super().__init__(boundary_name, values, time_values)
		self.type = "dirichlet"

class NeumannBC(GeneralBC):
	def __init__(self, boundary_name : str, values : list, time_values : list):
		super().__init__(boundary_name, values, time_values)
		self.type = "neumann"

class RobinBC(GeneralBC):
	def __init__(self, boundary_name : str, values : list, h : float, time_values : list):
		super().__init__(boundary_name, values, time_values)
		self.type = "robin"
		self.h = h



class BcHandler():
	def __init__(self, equation):
		self.eq = equation
		self.dirichlet_boundaries = []
		self.neumann_boundaries = []
		self.robin_boundaries = []

	def reset_boundary_conditions(self):
		self.dirichlet_boundaries = []
		self.neumann_boundaries = []
		self.robin_boundaries = []

	def add_boundary_condition(self, bc : GeneralBC):
		if bc.type == "dirichlet":
			self.dirichlet_boundaries.append(bc)
		elif bc.type == "neumann":
			self.neumann_boundaries.append(bc)
		elif bc.type == "robin":
			self.robin_boundaries.append(bc)
		else:
			raise Exception(f"Boundary type {bc.type} not supported.")

	def update_bcs(self, t):
		self.update_dirichlet(t)
		self.update_neumann(t)
		self.update_robin(t)

	def update_dirichlet(self, t):
		self.dirichlet_bcs = []
		for bc in self.dirichlet_boundaries:
			value = np.interp(t, bc.time_values, bc.values)
			dofs = do.fem.locate_dofs_topological(
				self.eq.V,
				self.eq.grid.boundary_dim,
				self.eq.grid.get_boundary_tags(bc.boundary_name)
			)
			self.dirichlet_bcs.append(
				do.fem.dirichletbc(
					do.default_scalar_type(value),
					dofs,
					self.eq.V
				)
			)

	def update_neumann(self, t):
		self.neumann_bcs = []
		for bc in self.neumann_boundaries:
			value = np.interp(t, bc.time_values, bc.values)
			self.neumann_bcs.append(value*self.eq.T_*self.eq.ds(self.eq.grid.get_boundary_tag(bc.boundary_name)))

	def update_robin(self, t):
		self.robin_bcs_a = []
		self.robin_bcs_b = []
		for bc in self.robin_boundaries:
			T_inf = np.interp(t, bc.time_values, bc.values)
			self.robin_bcs_a.append(bc.h*self.eq.dT*self.eq.T_*self.eq.ds(self.eq.grid.get_boundary_tag(bc.boundary_name)))
			self.robin_bcs_b.append(bc.h*T_inf*self.eq.T_*self.eq.ds(self.eq.grid.get_boundary_tag(bc.boundary_name)))

