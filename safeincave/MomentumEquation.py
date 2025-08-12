"""
Discretization of the momentum balance equations
"""
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
from MomentumBC import BcHandler
import Utils as utils

class LinearMomentumBase(ABC):
	def __init__(self, grid, theta):
		self.grid = grid
		self.theta = theta

		self.create_function_spaces()
		self.create_ds_dx()

		self.n_elems = self.DG0_1.dofmap.index_map.size_local + len(self.DG0_1.dofmap.index_map.ghosts)
		self.n_nodes = self.CG1_1.dofmap.index_map.size_local + len(self.CG1_1.dofmap.index_map.ghosts)

		self.commom_fields()
		self.create_fenicsx_fields()
		self.create_pytorch_fields()

	def commom_fields(self):
		self.T0 = to.zeros(self.n_elems, dtype=to.float64)
		self.Temp = to.zeros(self.n_elems, dtype=to.float64)
		self.sig = do.fem.Function(self.DG0_3x3)
		self.eps_tot = do.fem.Function(self.DG0_3x3)
		self.u = do.fem.Function(self.CG1_3x1)
		self.q_elems = do.fem.Function(self.DG0_1)
		self.q_nodes = do.fem.Function(self.CG1_1)
		self.p_elems = do.fem.Function(self.DG0_1)
		self.p_nodes = do.fem.Function(self.CG1_1)

	def set_material(self, material : Material):
		self.mat = material
		self.initialize()

	def set_T(self, T):
		self.Temp = T

	def set_T0(self, T0):
		self.T0 = T0

	def set_solver(self, solver : PETSc.KSP):
		self.solver = solver

	def set_boundary_conditions(self, bc : BcHandler):
		self.bc = bc

	def create_function_spaces(self):
		self.CG1_3x1 = do.fem.functionspace(self.grid.mesh, ("Lagrange", 1, (self.grid.domain_dim, )))
		self.DG0_1 = do.fem.functionspace(self.grid.mesh, ("DG", 0))
		self.CG1_1 = do.fem.functionspace(self.grid.mesh, ("Lagrange", 1))
		self.DG0_3x3 = do.fem.functionspace(self.grid.mesh, ("DG", 0, (3, 3)))
		self.DG0_6x6 = do.fem.functionspace(self.grid.mesh, ("DG", 0, (6, 6)))

	def create_ds_dx(self):
		self.ds = ufl.Measure("ds", domain=self.grid.mesh, subdomain_data=self.grid.get_boundaries())
		self.dx = ufl.Measure("dx", domain=self.grid.mesh, subdomain_data=self.grid.get_subdomains())

	def create_normal(self):
		n = ufl.FacetNormal(self.grid.mesh)
		self.normal = ufl.dot(n, self.u_)

	def build_body_force(self, g : list):
		density = do.fem.Function(self.DG0_1)
		density.x.array[:] = self.mat.density
		body_force = density*do.fem.Constant(self.grid.mesh, do.default_scalar_type(tuple(g)))
		self.b_body = ufl.dot(body_force, self.u_)*self.dx

	# def compute_q_nodes(self) -> do.fem.Function:
	# 	dev = self.sig - (1/3)*ufl.tr(self.sig)*ufl.Identity(3)
	# 	q_form = ufl.sqrt((3/2)*ufl.inner(dev, dev))
	# 	self.q_nodes = utils.project(q_form, self.CG1_1)

	# def compute_q_elems(self) -> do.fem.Function:
	# 	dev = self.sig - (1/3)*ufl.tr(self.sig)*ufl.Identity(3)
	# 	q_form = ufl.sqrt((3/2)*ufl.inner(dev, dev))
	# 	self.q_elems = utils.project(q_form, self.DG0_1)

	def compute_q_nodes(self) -> do.fem.Function:
		stress = utils.numpy2torch(self.sig.x.array.reshape((self.n_elems, 3, 3)))
		I1 = stress[:,0,0] + stress[:,1,1] + stress[:,2,2]
		I2 = stress[:,0,0]*stress[:,1,1] + stress[:,1,1]*stress[:,2,2] + stress[:,0,0]*stress[:,2,2] - stress[:,0,1]**2 - stress[:,0,2]**2 - stress[:,1,2]**2
		J2 = (1/3)*I1**2 - I2
		q_to = to.sqrt(3*J2)
		self.q_nodes.x.array[:] = self.grid.A_csr.dot(q_to.numpy())

	def compute_q_elems(self) -> do.fem.Function:
		stress = utils.numpy2torch(self.sig.x.array.reshape((self.n_elems, 3, 3)))
		I1 = stress[:,0,0] + stress[:,1,1] + stress[:,2,2]
		I2 = stress[:,0,0]*stress[:,1,1] + stress[:,1,1]*stress[:,2,2] + stress[:,0,0]*stress[:,2,2] - stress[:,0,1]**2 - stress[:,0,2]**2 - stress[:,1,2]**2
		J2 = (1/3)*I1**2 - I2
		q_to = to.sqrt(3*J2)
		self.q_elems.x.array[:] = self.grid.smoother.dot(q_to.numpy())

	def compute_total_strain(self):
		self.eps_tot = utils.project(utils.epsilon(self.u), self.DG0_3x3)
		eps_to = utils.numpy2torch(self.eps_tot.x.array.reshape((self.n_elems, 3, 3)))
		return eps_to

	def compute_eps_th(self):
		eps_th = to.zeros((self.n_elems, 3, 3), dtype=to.float64)
		deltaT = self.Temp - self.T0
		for elem_th in self.mat.elems_th:
			elem_th.compute_eps_th(deltaT)
			eps_th += elem_th.eps_th
		return eps_th

	def compute_eps_ne_k(self, dt):
		eps_ne_k = to.zeros((self.n_elems, 3, 3), dtype=to.float64)
		for elem_ne in self.mat.elems_ne:
			elem_ne.compute_eps_ne_k(dt*self.theta, dt*(1 - self.theta))
			eps_ne_k += elem_ne.eps_ne_k
		return eps_ne_k

	def compute_eps_ne_rate(self, stress, dt):
		for elem_ne in self.mat.elems_ne:
			elem_ne.compute_eps_ne_rate(stress, dt*self.theta, self.Temp, return_eps_ne=False)

	def update_eps_ne_rate_old(self):
		for elem_ne in self.mat.elems_ne:
			elem_ne.update_eps_ne_rate_old()

	def update_eps_ne_old(self, stress, stress_k, dt):
		for elem_ne in self.mat.elems_ne:
			elem_ne.update_eps_ne_old(stress, stress_k, dt*(1-self.theta))

	def increment_internal_variables(self, stress, stress_k, dt):
		for elem_ne in self.mat.elems_ne:
			elem_ne.increment_internal_variables(stress, stress_k, dt)

	def update_internal_variables(self):
		for elem_ne in self.mat.elems_ne:
			elem_ne.update_internal_variables()

	def create_solution_vector(self):
		self.X = do.fem.Function(self.V)

	def run_after_solve(self):
		pass

	@abstractmethod
	def compute_CT(self, dt, stress_k):
		pass

	@abstractmethod
	def compute_eps_rhs(self, dt, stress_k, eps_k):
		pass

	@abstractmethod
	def compute_elastic_stress(self, eps_e):
		pass

	@abstractmethod
	def compute_stress(self, eps_tot, eps_rhs, p):
		pass

	@abstractmethod
	def create_fenicsx_fields(self):
		pass

	@abstractmethod
	def create_pytorch_fields(self):
		pass

	@abstractmethod
	def create_trial_test_functions(self):
		pass

	@abstractmethod
	def get_uV(self):
		"""Function space for displacement field"""
		pass

	@abstractmethod
	def initialize(self) -> None:
		pass

	@abstractmethod
	def split_solution(self):
		pass

	@abstractmethod
	def split_solution(self):
		pass

	@abstractmethod
	def compute_p_nodes(self):
		pass

	@abstractmethod
	def solve_elastic_response(self):
		pass

	@abstractmethod
	def solve(self):
		pass





class LinearMomentum(LinearMomentumBase):
	def __init__(self, grid, theta):
		super().__init__(grid, theta)
		self.V = self.CG1_3x1
		self.create_trial_test_functions()
		self.create_normal()
		self.create_solution_vector()

	def create_fenicsx_fields(self):
		self.C = do.fem.Function(self.DG0_6x6)
		self.CT = do.fem.Function(self.DG0_6x6)
		self.eps_rhs = do.fem.Function(self.DG0_3x3)

	def create_pytorch_fields(self):
		self.eps_rhs_to = to.zeros((self.n_elems, 3, 3))

	def create_trial_test_functions(self):
		self.du = ufl.TrialFunction(self.V)
		self.u_ = ufl.TestFunction(self.V)

	def get_uV(self):
		return self.V

	def initialize(self) -> None:
		self.C.x.array[:] = to.flatten(self.mat.C)

	def compute_CT(self, stress_k, dt):
		self.mat.compute_G_B(stress_k, dt, self.theta, self.Temp)
		self.mat.compute_CT(dt, self.theta)
		self.CT.x.array[:] = to.flatten(self.mat.CT)

	def compute_elastic_stress(self, eps_e : to.Tensor) -> to.Tensor:
		stress_to = dotdot2(self.mat.C, eps_e)
		self.sig.x.array[:] = to.flatten(stress_to)
		return stress_to

	def compute_stress(self, eps_tot_to : to.Tensor, *_) -> to.Tensor:
		stress_to = dotdot2(self.mat.CT, eps_tot_to - self.eps_rhs_to)
		self.sig.x.array[:] = to.flatten(stress_to)
		return stress_to

	def compute_eps_rhs(self, dt : float, stress_k : to.Tensor) -> None:
		eps_ne_k = self.compute_eps_ne_k(dt)
		eps_th = self.compute_eps_th()
		self.eps_rhs_to = eps_ne_k + eps_th - dt*(1 - self.theta)*(self.mat.B + dotdot2(self.mat.G, stress_k))
		self.eps_rhs.x.array[:] = to.flatten(self.eps_rhs_to)

	def solve_elastic_response(self):
		# Build bilinear form
		a = ufl.inner(utils.dotdot(self.C, utils.epsilon(self.du)), utils.epsilon(self.u_))*self.dx
		bilinear_form = do.fem.form(a)
		A = do.fem.petsc.assemble_matrix(bilinear_form, bcs=self.bc.dirichlet_bcs)
		A.assemble()

		# Build linear form
		linear_form = do.fem.form(self.b_body + sum(self.bc.neumann_bcs))
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

	def split_solution(self):
		self.u = self.X

	# def compute_p_nodes(self) -> do.fem.Function:
	# 	self.p_nodes = utils.project(ufl.tr(self.sig)/3, self.CG1_1)

	# def compute_p_elems(self) -> do.fem.Function:
	# 	# self.p_elems = utils.project(ufl.tr(self.sig)/3, self.DG0_1)
	# 	stress_to = utils.numpy2torch(self.sig.x.array.reshape((self.n_elems, 3, 3)))
	# 	p_to = to.einsum("kii->k", stress_to)
	# 	self.p_elems.x.array[:] = to.flatten(p_to)

	def compute_p_nodes(self) -> do.fem.Function:
		stress = utils.numpy2torch(self.sig.x.array.reshape((self.n_elems, 3, 3)))
		I1 = stress[:,0,0] + stress[:,1,1] + stress[:,2,2]
		p_to = I1/3
		self.p_nodes.x.array[:] = self.grid.A_csr.dot(p_to)

	def compute_p_elems(self) -> do.fem.Function:
		stress_to = utils.numpy2torch(self.sig.x.array.reshape((self.n_elems, 3, 3)))
		p_to = to.einsum("kii->k", stress_to)
		p_to = self.grid.smoother.dot(p_to.numpy())
		self.p_elems.x.array[:] = p_to

	def solve(self, stress_k_to, t, dt):

		# Compute consistent tangent matrix
		self.compute_CT(stress_k_to, dt)

		# Compute right-hand side epsilon
		self.compute_eps_rhs(dt, stress_k_to)

		# Build bilinear form
		a = ufl.inner(utils.dotdot(self.CT, utils.epsilon(self.du)), utils.epsilon(self.u_))*self.dx
		bilinear_form = do.fem.form(a)
		A = do.fem.petsc.assemble_matrix(bilinear_form, bcs=self.bc.dirichlet_bcs)
		A.assemble()

		# Build linear form
		b_rhs = ufl.inner(utils.dotdot(self.CT, self.eps_rhs), utils.epsilon(self.u_))*self.dx
		linear_form = do.fem.form(self.b_body + sum(self.bc.neumann_bcs) + b_rhs)
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

		self.run_after_solve()



