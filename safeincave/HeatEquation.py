# Copyright 2024 The safeincave community.
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
from abc import ABC, abstractmethod
import dolfinx as do
import basix
import ufl
from petsc4py import PETSc
import torch as to
from Utils import dotdot2
from MaterialProps import Material
from HeatBC import BcHandler
import Utils as utils

class HeatDiffusion():
	def __init__(self, grid):
		self.grid = grid

		self.create_function_spaces()

		self.n_elems = self.DG0_1.dofmap.index_map.size_local + len(self.DG0_1.dofmap.index_map.ghosts)
		self.n_nodes = self.V.dofmap.index_map.size_local + len(self.V.dofmap.index_map.ghosts)

		self.create_trial_test_functions()
		self.create_ds_dx()
		self.create_fenicsx_fields()
		# self.create_pytorch_fields()

	def set_material(self, material : Material):
		self.mat = material
		self.initialize()

	def set_solver(self, solver : PETSc.KSP):
		self.solver = solver

	def set_boundary_conditions(self, bc : BcHandler):
		self.bc = bc

	def create_trial_test_functions(self):
		self.dT = ufl.TrialFunction(self.V)
		self.T_ = ufl.TestFunction(self.V)

	def create_function_spaces(self):
		self.DG0_1 = do.fem.functionspace(self.grid.mesh, ("DG", 0))
		self.V = do.fem.functionspace(self.grid.mesh, ("Lagrange", 1))

	def create_ds_dx(self):
		self.ds = ufl.Measure("ds", domain=self.grid.mesh, subdomain_data=self.grid.get_boundaries())
		self.dx = ufl.Measure("dx", domain=self.grid.mesh, subdomain_data=self.grid.get_subdomains())

	def create_fenicsx_fields(self):
		self.k = do.fem.Function(self.DG0_1)
		self.rho = do.fem.Function(self.DG0_1)
		self.cp = do.fem.Function(self.DG0_1)
		self.T_old = do.fem.Function(self.V)
		self.T = do.fem.Function(self.V)
		self.X = do.fem.Function(self.V)

	def initialize(self):
		self.k.x.array[:] = self.mat.k
		self.rho.x.array[:] = self.mat.density
		self.cp.x.array[:] = self.mat.cp

	def split_solution(self):
		self.T = self.X

	def update_T_old(self):
		self.T_old.x.array[:] = self.T.x.array

	def set_initial_T(self, T_field):
		self.T_old.x.array[:] = T_field
		self.T.x.array[:] = T_field

	def get_T_elems(self):
		T_elems = utils.project(self.T, self.DG0_1)
		return utils.numpy2torch(T_elems.x.array)


	def solve(self, t, dt):
		# Update boundary conditions
		self.bc.update_bcs(t)

		# Build bilinear form
		a = (self.rho*self.cp*self.dT*self.T_/dt + self.k*ufl.dot(ufl.grad(self.dT), ufl.grad(self.T_)))*self.dx
		a += sum(self.bc.robin_bcs_a)
		bilinear_form = do.fem.form(a)
		A = do.fem.petsc.assemble_matrix(bilinear_form, bcs=self.bc.dirichlet_bcs)
		A.assemble()

		# Build linear form
		L = (self.rho*self.cp*self.T_old*self.T_/dt)*self.dx + sum(self.bc.neumann_bcs) + sum(self.bc.robin_bcs_b)
		linear_form = do.fem.form(L)
		b = do.fem.petsc.assemble_vector(linear_form)
		do.fem.petsc.apply_lifting(b, [bilinear_form], [self.bc.dirichlet_bcs])
		b.ghostUpdate(addv=PETSc.InsertMode.ADD_VALUES, mode=PETSc.ScatterMode.REVERSE)
		do.fem.petsc.set_bc(b, self.bc.dirichlet_bcs)
		b.ghostUpdate(addv=PETSc.InsertMode.INSERT, mode=PETSc.ScatterMode.FORWARD)

		# Solve linear system
		self.solver.setOperators(A)
		self.solver.solve(b, self.X.x.petsc_vec)
		self.X.x.scatter_forward()
		self.split_solution()

		# Update old temperature field
		self.update_T_old()


