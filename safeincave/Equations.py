"""
Everything related to the solution of the linear momentum balance equation.
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

import os
import torch as to
import dolfin as do
import numpy as np
import Utils as utils
from ConstitutiveModel import ConstitutiveModel

class LinearMomentum():
	"""
	This class is intended to solve the linear momemum balance equation considering
	a general non-linear constitutive models.

	Parameters
	----------
	m : :class:`safeincave.ConstitutiveModel.ConstitutiveModel`
		Object containing the constitutive model data.
	theta : int
		Choice of the time integration method: explicit (theta=1.0), Crank-Nicolson (theta=0.5) or fully-implicit (theta=0.0).
	"""
	def __init__(self, grid, theta, input_file):
		self.grid = grid
		self.input_file = input_file

		# Output folder
		self.operation_output_folder = os.path.join(self.input_file["output"]["path"], "operation")

		# Get number of elements
		self.n_elems = self.grid.n_elems
		self.n_nodes = self.grid.n_nodes

		# Define constitutive model
		self.m = ConstitutiveModel(self.grid, self.input_file["constitutive_model"])

		# Time integration method
		self.theta = theta

		# Transient settings
		self.time_list = self.input_file["time_settings"]["time_list"]
		self.t0 = self.time_list[0]

		# Create function spaces
		self.CG_3x1 = do.VectorFunctionSpace(self.grid.mesh, "CG", 1)
		self.DG_1x1 = do.FunctionSpace(self.grid.mesh, "DG", 0)
		self.DG_3x3 = do.TensorFunctionSpace(self.grid.mesh, "DG", 0)
		self.DG_6x6 = do.TensorFunctionSpace(self.grid.mesh, "DG", 0, shape=(6, 6))

		# Create tensor fields
		self.C0 = do.Function(self.DG_6x6)
		self.C1 = do.Function(self.DG_6x6)
		self.CT = do.Function(self.DG_6x6)
		self.eps_tot = do.Function(self.DG_3x3)
		self.eps_rhs = do.Function(self.DG_3x3)

		self.sigma = do.Function(self.DG_3x3)
		self.sigma_0 = do.Function(self.DG_3x3)
		self.sigma_smooth = do.Function(self.DG_3x3)
		self.sigma_smooth.rename("Stress", "MPa")

		self.sigma_v = do.Function(self.DG_1x1)
		self.sigma_v.rename("Mean stress", "MPa")
		self.sigma_v_smooth = do.Function(self.DG_1x1)
		self.sigma_v_smooth.rename("Mean stress", "MPa")

		self.von_mises = do.Function(self.DG_1x1)
		self.von_mises.rename("Von Mises stress", "MPa")
		self.von_mises_smooth = do.Function(self.DG_1x1)
		self.von_mises_smooth.rename("Von Mises stress", "MPa")

		# Define variational problem
		self.du = do.TrialFunction(self.CG_3x1)
		self.v = do.TestFunction(self.CG_3x1)
		self.ds = do.Measure("ds", domain=self.grid.mesh, subdomain_data=self.grid.get_boundaries())
		self.dx = do.Measure("dx", domain=self.grid.mesh, subdomain_data=self.grid.get_subdomains())
		self.normal = do.dot(self.v, do.FacetNormal(self.grid.mesh))
		self.n = do.FacetNormal(self.grid.mesh)

		# Create displacement vector
		self.u = do.Function(self.CG_3x1)
		self.u.rename("Displacement", "m")

		# Define pytorch tensor quantities
		self.eps_e_torch = to.zeros((self.n_elems, 3, 3), dtype=to.float64)
		self.eps_ve_torch = to.zeros((self.n_elems, 3, 3), dtype=to.float64)
		self.eps_ie_torch = to.zeros((self.n_elems, 3, 3), dtype=to.float64)
		self.eps_bar_ve_torch = to.zeros((self.n_elems, 3, 3), dtype=to.float64)
		self.eps_bar_ie_torch = to.zeros((self.n_elems, 3, 3), dtype=to.float64)
		self.GT_ve_torch = to.zeros((self.n_elems, 6, 6), dtype=to.float64)
		self.GT_ie_torch = to.zeros((self.n_elems, 6, 6), dtype=to.float64)
		self.BT_ve_torch = to.zeros((self.n_elems, 3, 3), dtype=to.float64)
		self.BT_ie_torch = to.zeros((self.n_elems, 3, 3), dtype=to.float64)
		self.CT_torch = to.zeros((self.n_elems, 6, 6), dtype=to.float64)
		self.stress_torch = to.zeros((self.n_elems, 3, 3), dtype=to.float64)
		self.stress_k_torch = to.zeros((self.n_elems, 3, 3), dtype=to.float64)

		# Define salt specific weight
		self.direction = self.input_file["body_force"]["direction"]
		self.density = do.Function(self.DG_1x1)
		self.density.vector()[:] = self.grid.get_parameter(self.input_file["body_force"]["density"])
		self.gravity = self.input_file["body_force"]["gravity"]

		self.g = [0.0, 0.0, 0.0]
		self.g[self.direction] = self.gravity
		self.body_force = self.density*do.Constant(tuple(self.g))

		# Build body forces on RHS vector
		self.b_body = do.dot(self.body_force, self.v)*self.dx

		# Define linear solver
		self.solver = self.define_solver()

		# Create output file
		self.u_vtk = do.File(os.path.join(self.operation_output_folder, "vtk", "displacement", "displacement.pvd"))
		self.stress_vtk = do.File(os.path.join(self.operation_output_folder, "vtk", "stress", "stress.pvd"))
		self.stress_smooth_vtk = do.File(os.path.join(self.operation_output_folder, "vtk", "stress_smooth", "stress_smooth.pvd"))

		self.q_vtk = do.File(os.path.join(self.operation_output_folder, "vtk", "q", "q.pvd"))
		self.qs_vtk = do.File(os.path.join(self.operation_output_folder, "vtk", "q_smooth", "q_smooth.pvd"))
		self.p_vtk = do.File(os.path.join(self.operation_output_folder, "vtk", "p", "p.pvd"))
		self.ps_vtk = do.File(os.path.join(self.operation_output_folder, "vtk", "p_smooth", "p_smooth.pvd"))



	def save_solution(self, t):
		self.u_vtk << (self.u, t)
		self.stress_vtk << (self.sigma, t)

		stress_np = self.sigma.vector()[:].reshape((self.n_elems, 3, 3))
		# self.sigma_smooth.vector()[:] = self.apply_smoother(stress_np)
		self.sigma_smooth.vector()[:] = self.apply_smoother(self.stress_torch.numpy())
		self.stress_smooth_vtk << (self.sigma_smooth, t)

		MPa = 1e6
		sxx = self.stress_torch[:,0,0]/MPa
		syy = self.stress_torch[:,1,1]/MPa
		szz = self.stress_torch[:,2,2]/MPa
		sxy = self.stress_torch[:,0,1]/MPa
		sxz = self.stress_torch[:,0,2]/MPa
		syz = self.stress_torch[:,1,2]/MPa

		I1 = sxx + syy + szz
		I2 = sxx*syy + syy*szz + sxx*szz - sxy**2 - syz**2 - sxz**2
		I3 = sxx*syy*szz + 2*sxy*syz*sxz - szz*sxy**2 - sxx*syz**2 - syy*sxz**2
		J2 = (1/3)*I1**2 - I2
		q = np.sqrt(3*J2)
		p = I1/3

		self.von_mises_smooth.vector()[:] = self.grid.smoother.dot(q.numpy())
		self.qs_vtk << (self.von_mises_smooth, t)

		self.von_mises.vector()[:] = q
		self.q_vtk << (self.von_mises, t)

		self.sigma_v_smooth.vector()[:] = self.grid.smoother.dot(p.numpy())
		self.ps_vtk << (self.sigma_v_smooth, t)

		self.sigma_v.vector()[:] = p
		self.p_vtk << (self.sigma_v, t)



	def apply_smoother(self, field_np):
		field_copy = field_np.copy()
		field_copy[:,0,0] = self.grid.smoother.dot(field_copy[:,0,0])
		field_copy[:,1,1] = self.grid.smoother.dot(field_copy[:,1,1])
		field_copy[:,2,2] = self.grid.smoother.dot(field_copy[:,2,2])
		field_copy[:,0,1] = field_copy[:,1,0] = self.grid.smoother.dot(field_copy[:,0,1])
		field_copy[:,0,2] = field_copy[:,2,0] = self.grid.smoother.dot(field_copy[:,0,2])
		field_copy[:,1,2] = field_copy[:,2,1] = self.grid.smoother.dot(field_copy[:,1,2])
		return field_copy.flatten()



	def initialize(self, solve_equilibrium=False, verbose=True, save_results=False):
		# Apply Neumann boundary condition
		self.apply_neumann_bc(self.t0)

		# Define Dirichlet boundary conditions
		self.define_dirichlet_bc(self.t0)

		# Build RHS vector
		b = do.assemble(self.b_body + sum(self.integral_neumann))

		# Initialize elastic stiffness matrix, C0
		self.C0.vector()[:] = to.flatten(self.m.C0)

		# Build stiffness matrix
		a_form = do.inner(utils.dotdot(self.C0, utils.epsilon(self.du)), utils.epsilon(self.v))*self.dx
		A = do.assemble(a_form)

		# Solve linear system
		[bc.apply(A, b) for bc in self.bcs]
		self.solver.solve(A, self.u.vector(), b)

		# Compute total strain
		self.eps_tot.assign(utils.local_projection(utils.epsilon(self.u), self.DG_3x3))

		# Compute stress
		eps_tot_torch = utils.numpy2torch(self.eps_tot.vector()[:].reshape((self.n_elems, 3, 3)))
		self.compute_stress_C0(eps_tot_torch)

		if solve_equilibrium:
			self.solve_equilibrium(verbose, save_results)

		# Compute initial hardening
		self.compute_initial_hardening(verbose)

		# Compute old ielastic strain rates
		self.compute_eps_ie_rate()

		# Update inelastic strain rate (Warning! Do NOT update eps_ie here, because this is wrong!)
		self.update_eps_ie_rate_old()

		# Assign stress
		self.sigma.vector()[:] = self.stress_torch.flatten()


	def solve_equilibrium(self, verbose=True, save_results=False):
		# Get maximum time step size
		dt = self.input_file["simulation_settings"]["equilibrium"]["dt_max"]
		t = self.t0

		# Apply Neumann boundary condition
		self.apply_neumann_bc(self.t0)

		# Define Dirichlet boundary conditions
		self.define_dirichlet_bc(self.t0)

		# Build RHS vector
		b = do.assemble(self.b_body + sum(self.integral_neumann))

		# Initialize elastic stiffness matrix, C0
		self.C0.vector()[:] = to.flatten(self.m.C0)

		# Build stiffness matrix
		a_form = do.inner(utils.dotdot(self.C0, utils.epsilon(self.du)), utils.epsilon(self.v))*self.dx
		A = do.assemble(a_form)

		# Solve linear system
		[bc.apply(A, b) for bc in self.bcs]
		self.solver.solve(A, self.u.vector(), b)

		# Compute total strain
		self.eps_tot.assign(utils.local_projection(utils.epsilon(self.u), self.DG_3x3))

		# Compute stress
		eps_tot_torch = utils.numpy2torch(self.eps_tot.vector()[:].reshape((self.n_elems, 3, 3)))
		self.compute_stress_C0(eps_tot_torch)

		# Compute old viscoelastic strain rates
		self.compute_eps_ve_rate(0)

		# Update viscoelastic strain rate (Warning! Do NOT update eps_ie here, because this is wrong!)
		self.update_eps_ve_rate_old()

		if save_results:
			# Output folder
			equilibrium_output_folder = os.path.join(self.input_file["output"]["path"], "equilibrium")

			# Create output file
			u_vtk = do.File(os.path.join(equilibrium_output_folder, "vtk", "displacement", "displacement.pvd"))
			stress_vtk = do.File(os.path.join(equilibrium_output_folder, "vtk", "stress", "stress.pvd"))

			# Save output fields
			u_vtk << (self.u, t)
			stress_vtk << (self.sigma, t)

		# Initialize pseudo time control settings
		n_step = 1
		tol_time = self.input_file["simulation_settings"]["equilibrium"]["time_tol"]
		error_time = 2*tol_time
		eps_tot_old = utils.numpy2torch(self.eps_tot.vector()[:])

		while error_time > tol_time or n_step <= 2:

			# Increment time
			t += dt

			# Compute GT matrix field for viscoelastic elements
			self.compute_GT_BT_ve(dt)

			# Iterative loop settings
			tol = 1e-7
			error = 2*tol
			ite = 0
			maxiter = 40

			while error > tol and ite < maxiter:

				# Update total strain of previous iteration (eps_tot_k <-- eps_tot)
				eps_tot_k = utils.numpy2torch(self.eps_tot.vector()[:])

				# Update stress of previous iteration (stress_k <-- stress)
				self.update_stress()

				# Compute CT
				self.compute_CT(dt)

				# Compute right-hand side
				self.compute_eps_rhs(dt)

				# Build rhs
				b_rhs = do.inner(utils.dotdot(self.CT, self.eps_rhs), utils.epsilon(self.v))*self.dx
				b = do.assemble(self.b_body + sum(self.integral_neumann) + b_rhs)

				# Build lhs
				a_form = do.inner(utils.dotdot(self.CT, utils.epsilon(self.du)), utils.epsilon(self.v))*self.dx
				A = do.assemble(a_form)

				# Solve linear system
				[bc.apply(A, b) for bc in self.bcs]
				self.solver.solve(A, self.u.vector(), b)

				# Compute total strain
				self.eps_tot.assign(utils.local_projection(utils.epsilon(self.u), self.DG_3x3))
				eps_tot_torch = utils.numpy2torch(self.eps_tot.vector()[:].reshape((self.n_elems, 3, 3)))

				# Compute stress
				self.compute_stress(eps_tot_torch, dt)

				# Compute strain rates
				self.compute_eps_ve_rate(dt)

				# Compute error
				if self.theta == 1.0:
					error = 0.0
				else:
					eps_tot_k_flat = to.flatten(eps_tot_k)
					eps_tot_flat = self.eps_tot.vector()[:]
					error = np.linalg.norm(eps_tot_k_flat - eps_tot_flat) / np.linalg.norm(eps_tot_flat)

				# Increment iteration counter
				ite += 1

			# Compute strains
			self.compute_eps_ve(dt)

			# Update old viscoelastic strains
			self.update_eps_ve_old()

			# Update old viscoelastic strain rates
			self.update_eps_ve_rate_old()

			# Compute time error (to check if steady state is achieved)
			eps_tot_flat = self.eps_tot.vector()[:]
			error_time = np.linalg.norm(eps_tot_old - eps_tot_flat) / np.linalg.norm(eps_tot_flat)
			eps_tot_old = utils.numpy2torch(self.eps_tot.vector()[:])
			# self.sigma.vector()[:] = self.stress_torch.flatten()

			if save_results:
				u_vtk << (self.u, t)
				stress_vtk << (self.sigma, t)

			# Print stuff
			if verbose:
				print(n_step, t/utils.hour, ite, error, error_time)
				# print()
			n_step += 1



	def solve(self, t, dt):
		# Apply Neumann boundary condition
		self.apply_neumann_bc(t)

		# Define Dirichlet boundary conditions
		self.define_dirichlet_bc(t)

		# Compute GT and BT matrix fields for viscoelastic elements
		self.compute_GT_BT_ve(dt)

		# Iterative loop settings
		tol = 1e-7
		self.error = 2*tol
		self.ite = 0
		maxiter = 40

		while self.error > tol and self.ite < maxiter:

			# Update total strain of previous iteration (eps_tot_k <-- eps_tot)
			eps_tot_k = utils.numpy2torch(self.eps_tot.vector()[:])

			# Update stress of previous iteration (stress_k <-- stress)
			self.update_stress()

			# Compute GT and BT matrix fields for inelastic elements
			self.compute_GT_BT_ie(dt)

			# Compute CT
			self.compute_CT(dt)

			# Compute right-hand side
			self.compute_eps_rhs(dt)

			# Build rhs
			b_rhs = do.inner(utils.dotdot(self.CT, self.eps_rhs), utils.epsilon(self.v))*self.dx
			# b = do.assemble(sum(self.integral_neumann) + b_rhs)
			b = do.assemble(self.b_body + sum(self.integral_neumann) + b_rhs)

			# Build lhs
			a_form = do.inner(utils.dotdot(self.CT, utils.epsilon(self.du)), utils.epsilon(self.v))*self.dx
			A = do.assemble(a_form)

			# Solve linear system
			[bc.apply(A, b) for bc in self.bcs]
			self.solver.solve(A, self.u.vector(), b)

			# Compute total strain
			self.eps_tot.assign(utils.local_projection(utils.epsilon(self.u), self.DG_3x3))
			eps_tot_torch = utils.numpy2torch(self.eps_tot.vector()[:].reshape((self.n_elems, 3, 3)))

			# Compute stress
			self.compute_stress(eps_tot_torch, dt)

			# Increment internal variables
			self.increment_internal_variables(dt)

			# Compute strain rates
			self.compute_eps_ie_rate()
			self.compute_eps_ve_rate(dt)

			# Compute error
			if self.theta == 1.0:
				self.error = 0.0
			else:
				eps_tot_k_flat = to.flatten(eps_tot_k)
				eps_tot_flat = self.eps_tot.vector()[:]
				self.error = np.linalg.norm(eps_tot_k_flat - eps_tot_flat) / np.linalg.norm(eps_tot_flat)

			self.ite += 1

		# Update internal variables
		self.update_internal_variables()

		# Compute strains
		self.compute_eps_ie(dt)
		self.compute_eps_ve(dt)

		# Update old non-elastic strains
		self.update_eps_ie_old()
		self.update_eps_ve_old()

		# Update old non-elastic strain rates
		self.update_eps_ie_rate_old()
		self.update_eps_ve_rate_old()




	def compute_initial_hardening(self, verbose):
		# Compute initial hardening
		for elem in self.m.elems_ie:
			try:
				# Compute initial hardening parameter (alpha_0) based on initial stresses
				elem.compute_initial_hardening(self.stress_torch, Fvp_0=0.0)
				
				# Compute initial yield function values (the computed values should be equal to Fvp_0)
				I1, I2, I3, J2, J3, Sr, I1_star = elem.compute_stress_invariants(*elem.extract_stress_components(self.stress_torch))
				_ = elem.compute_Fvp(elem.alpha, I1_star, J2, Sr)

				if verbose:
					print("Fvp: ", float(max(elem.Fvp)))
					print("alpha_min: ", float(min(elem.alpha)))
					print("alpha_max: ", float(max(elem.alpha)))
					print("alpha_avg: ", float(np.average(elem.alpha)))
					print()
			except:
				pass




	def compute_GT_BT_ve(self, dt):
		"""
		Computes matrices :math:`\\mathbb{G}_T` and :math:`\\mathbf{B}_T` for the
		viscoelastic elements. That is

		.. math::

			\\mathbb{G}_T = \\sum_{i=1}^{N_{ve}} \\mathbb{G}_i
			\\quad \\text{and} \\quad
			\\mathbf{B}_T = \\sum_{i=1}^{N_{ve}} \\mathbf{B}_i

		where :math:`N_{ve}` denotes the number of viscoelastic (Kelvin-Voigt)
		elements present in the constitutive model.

		Parameters
		----------
		dt : float
			Time step size.

		Returns
		-------
		None 
		"""
		self.GT_ve_torch = to.zeros((self.n_elems, 6, 6), dtype=to.float64)
		self.BT_ve_torch = to.zeros((self.n_elems, 3, 3), dtype=to.float64)
		phi2 = dt*(1 - self.theta)
		for elem_ve in self.m.elems_ve:
			elem_ve.compute_G_B(self.stress_torch, phi2)
			self.GT_ve_torch += elem_ve.G
			self.BT_ve_torch += elem_ve.B

	def compute_GT_BT_ie(self, dt):
		"""
		Computes matrices :math:`\\mathbb{G}_{ie,T}` and :math:`\\mathbf{B}_{ie,T}` for the
		inelastic elements. That is

		.. math::

			\\mathbb{G}_{ie,T} = \\sum_{i=1}^{N_{ie}} \\mathbb{G}_i
			\\quad \\text{and} \\quad
			\\mathbf{B}_{ie,T} = \\sum_{i=1}^{N_{ie}} \\mathbf{B}_i

		where :math:`N_{ie}` denotes the number of inelastic elements present
		in the constitutive model.

		Parameters
		----------
		dt : float
			Time step size.

		Returns
		-------
		None 
		"""
		self.GT_ie_torch = to.zeros((self.n_elems, 6, 6), dtype=to.float64)
		self.BT_ie_torch = to.zeros((self.n_elems, 3, 3), dtype=to.float64)
		for elem_ie in self.m.elems_ie:
			elem_ie.compute_G_B(self.stress_torch, dt)
			self.GT_ie_torch += elem_ie.G
			self.BT_ie_torch += elem_ie.B

	def compute_CT(self, dt):
		"""
		Computes the consistent tangent matrix :math:`\\mathbb{C}_{T}`, which is given by

		.. math::

			\\mathbb{C}_{T} = \\left( 
				\\mathbb{C}_0^{-1} + \\Delta t (1 - \\theta) \\mathbb{G}_{T}
			\\right)^{-1}

		where :math:`\\mathbb{C}_0` is the stiffness matrix associated to the linear spring,
		and :math:`\\mathbb{G}_{T} = \\mathbb{G}_{ve,T} + \\mathbb{G}_{ie,T}`.

		Parameters
		----------
		dt : float
			Time step size.

		Returns
		-------
		None 
		"""
		GT = self.GT_ie_torch + self.GT_ve_torch
		self.CT_torch = to.linalg.inv(self.m.C0_inv + dt*(1-self.theta)*GT)
		self.CT.vector()[:] = to.flatten(self.CT_torch)

	def compute_eps_ve(self, dt):
		"""
		Computes the strains of all **viscoelastic** elements of the constitutive model.

		Parameters
		----------
		dt : float
			Time step size.

		Returns
		-------
			None
		"""
		self.eps_ve_torch = to.zeros((self.n_elems, 3, 3), dtype=to.float64)
		for elem_ve in self.m.elems_ve:
			elem_ve.compute_eps_ve(self.stress_torch, self.stress_k_torch, dt*(1-self.theta))
			self.eps_ve_torch += elem_ve.eps_ve

	def compute_eps_ie(self, dt):
		"""
		Computes the strains of all **inelastic** elements of the constitutive model.

		Parameters
		----------
		dt : float
			Time step size.

		Returns
		-------
			None
		"""
		self.eps_ie_torch = to.zeros((self.n_elems, 3, 3), dtype=to.float64)
		for elem_ie in self.m.elems_ie:
			elem_ie.compute_eps_ie(self.stress_torch, self.stress_k_torch, dt*(1-self.theta))
			self.eps_ie_torch += elem_ie.eps_ie

	def compute_eps_e(self):
		"""
		Computes the strains of all **elastic** elements of the constitutive model.

		Parameters
		----------
		dt : float
			Time step size.

		Returns
		-------
			None
		"""
		self.eps_e_torch = to.zeros((self.n_elems, 3, 3), dtype=to.float64)
		for elem_e in self.m.elems_e:
			elem_e.compute_eps_e(self.stress_torch)
			self.eps_e_torch += elem_e.eps_e

	def compute_eps_bar_ie(self, dt):
		"""
		Computes the eps_bar of all **inelastic** elements of the constitutive model.

		Parameters
		----------
		dt : float
			Time step size.

		Returns
		-------
			None
		"""
		self.eps_bar_ie_torch = to.zeros((self.n_elems, 3, 3), dtype=to.float64)
		for elem_ie in self.m.elems_ie:
			elem_ie.compute_eps_bar(dt*self.theta, dt*(1 - self.theta))
			self.eps_bar_ie_torch += elem_ie.eps_bar

	def compute_eps_bar_ve(self, dt):
		"""
		Computes the eps_bar of all **viscoelastic** elements of the constitutive model.

		Parameters
		----------
		dt : float
			Time step size.

		Returns
		-------
			None
		"""
		self.eps_bar_ve_torch = to.zeros((self.n_elems, 3, 3), dtype=to.float64)
		for elem_ve in self.m.elems_ve:
			elem_ve.compute_eps_bar(dt*self.theta, dt*(1 - self.theta))
			self.eps_bar_ve_torch += elem_ve.eps_bar

	def compute_eps_bar(self, dt):
		"""
		Computes the eps_bar of **all** elements of the constitutive model.

		Parameters
		----------
		dt : float
			Time step size.

		Returns
		-------
			None
		"""
		self.compute_eps_bar_ie(dt)
		self.compute_eps_bar_ve(dt)
		self.eps_bar_torch = self.eps_bar_ve_torch + self.eps_bar_ie_torch

	def compute_eps_rhs(self, dt, *args):
		"""
		Computes the term :math:`\\pmb{\\varepsilon}_{rhs}^k` of the linearized momentum equation, that is,

		.. math::

			\\pmb{\\varepsilon}_\\text{rhs}^k = \\bar{\\pmb{\\varepsilon}}_{ne}^k - \\Delta t \\left( 1 - \\theta \\right) \\left( \\mathbb{G}_{ne} : \\pmb{\\sigma}^k + \\mathbf{B}_{ne} \\right)

		Parameters
		----------
		dt : float
			Time step size.

		Returns
		-------
			None
		"""
		self.compute_eps_bar(dt)
		BT = self.BT_ie_torch + self.BT_ve_torch
		GT = self.GT_ie_torch + self.GT_ve_torch
		self.eps_rhs_torch = self.eps_bar_torch - dt*(1-self.theta)*(BT + utils.dotdot2(GT, self.stress_k_torch))
		self.eps_rhs.vector()[:] = self.eps_rhs_torch.flatten()

	def compute_stress_C0(self, eps_e):
		"""
		Compute stress tensor as

		.. math::

			\\pmb{\\sigma} = \\mathbb{C}_0 : \\pmb{\\varepsilon}_e

		.. note::
			This operation considers that :math:`\\mathbb{C}_0` is represented by Voigt notation.

		Parameters
		----------
		eps_e : torch.Tensor
			A (nelems, 3, 3) storing the elastic strain for all grid elements.
		"""
		self.stress_torch = utils.dotdot2(self.m.C0, eps_e)

	def compute_stress(self, eps_tot, dt):
		"""
		Compute stress tensor as

		.. math::

			\\pmb{\\sigma}^{k+1} = \\mathbb{C}_T :
										    \\left[
										        \\pmb{\\varepsilon}^{k+1}
										        - \\bar{\\pmb{\\varepsilon}}_{ne}^k
										        + \\Delta t (1 - \\theta)
										            \\left( 
										               \\mathbf{B}_{ne}
										               + \\mathbb{G}_{ne} : \\pmb{\\sigma}^k
										            \\right)
										    \\right]
		
		Parameters
		----------
		eps_tot : torch.Tensor
			A (nelems, 3, 3) storing the total strain for all grid elements.
		"""
		self.compute_eps_bar(dt)
		GT = self.GT_ie_torch + self.GT_ve_torch
		BT = self.BT_ie_torch + self.BT_ve_torch
		self.stress_torch = utils.dotdot2(self.CT_torch, eps_tot - self.eps_bar_torch + dt*(1-self.theta)*(BT + utils.dotdot2(GT, self.stress_k_torch)))
		self.sigma.vector()[:] = self.stress_torch.flatten()

	def update_stress(self):
		"""
		Updates stress of the previous interation, that is :math:`\\pmb{\\sigma}^k \\leftarrow \\pmb{\\sigma}^{k+1}`.
		"""
		self.stress_k_torch = self.stress_torch.clone()

	def compute_eps_ie_rate(self):
		"""
		Computes strain rates of all inelastic elements of the constitutive model.
		"""
		for elem_ie in self.m.elems_ie:
			elem_ie.compute_eps_ie_rate(self.stress_torch, return_eps_ie=False)

	def update_eps_ie_rate_old(self):
		"""
		Updates inelastic strain rates of the previous time level,
		that is :math:`\\dot{\\pmb{\\varepsilon}}^t_{ie} \\leftarrow \\dot{\\pmb{\\varepsilon}}^{t + \\Delta t}_{ie}`.
		"""
		for elem_ie in self.m.elems_ie:
			elem_ie.update_eps_ie_rate_old()

	def update_eps_ie_old(self):
		"""
		Updates inelastic strains of the previous time level,
		that is :math:`\\pmb{\\varepsilon}^t_{ie} \\leftarrow \\pmb{\\varepsilon}^{t + \\Delta t}_{ie}`.
		"""
		for elem_ie in self.m.elems_ie:
			elem_ie.update_eps_ie_old()

	def compute_eps_ve_rate(self, dt):
		"""
		Computes the strain rates of the viscoelastic elements of the constitutive model.

		Parameters
		----------
		dt : float
			Time step size.
		"""
		for elem_ve in self.m.elems_ve:
			elem_ve.compute_eps_ve_rate(self.stress_torch, dt*self.theta, return_eps_ve=False)

	def update_eps_ve_rate_old(self):
		"""
		Updates viscoelastic strain rates of the previous time level,
		that is :math:`\\dot{\\pmb{\\varepsilon}}^t_{ve} \\leftarrow \\dot{\\pmb{\\varepsilon}}^{t + \\Delta t}_{ve}`.
		"""
		for elem_ve in self.m.elems_ve:
			elem_ve.update_eps_ve_rate_old()

	def update_eps_ve_old(self):
		"""
		Updates viscoelastic strains of the previous time level,
		that is :math:`\\pmb{\\varepsilon}^t_{ve} \\leftarrow \\pmb{\\varepsilon}^{t + \\Delta t}_{ve}`.
		"""
		for elem_ve in self.m.elems_ve:
			elem_ve.update_eps_ve_old()

	def increment_internal_variables(self, dt):
		"""
		Increment internal variables.

		Parameters
		----------
		dt : float
			Time step size.
		"""
		for elem_ie in self.m.elems_ie:
			elem_ie.increment_internal_variables(self.stress_torch, self.stress_k_torch, dt)

	def update_internal_variables(self):
		"""
		Update internal variables with values of previous iteration, that is
		:math:`\\alpha_i^k \\leftarrow \\alpha_i^{k+1}.`
		"""
		for elem_ie in self.m.elems_ie:
			elem_ie.update_internal_variables()



	def apply_neumann_bc(self, t):
		"""
		It reads all Neumann boundary conditions (external loads) applied to the geometry and
		builds the right-hand side vector.

		Parameters
		----------
		t : float
			Time level.
		"""
		self.integral_neumann = []
		i = 0
		bc_neumann_list = []
		for boundary in self.input_file["boundary_conditions"]:
			if self.input_file["boundary_conditions"][boundary]["type"] == "neumann":
				bc_direction = self.input_file["boundary_conditions"][boundary]["direction"]
				bc_density = self.input_file["boundary_conditions"][boundary]["density"]
				ref_position = self.input_file["boundary_conditions"][boundary]["reference_position"]
				values = self.input_file["boundary_conditions"][boundary]["values"]
				ref_load = -np.interp(t, self.time_list, values)
				value_neumann = do.Expression(f"load_ref + rho*g*(H - x[{bc_direction}])", load_ref=ref_load, rho=bc_density, g=self.gravity, H=ref_position, degree=1)
				self.integral_neumann.append(value_neumann*self.normal*self.ds(self.grid.get_boundary_tags(boundary)))



	def define_dirichlet_bc(self, t):
		"""
		Defines a list of Dirichlet boundary conditions to be applied to the
		linear system.

		Parameters
		----------
		t : float
			Time level of the simulation.

		Returns
		-------
		bcs : list[dolfin.fem.dirichletbc.DirichletBC]
			List containing the Dirichlet boundary conditions.
		"""
		self.bcs = []
		bc_dirichlet_list = []
		for boundary in self.input_file["boundary_conditions"]:
			if self.input_file["boundary_conditions"][boundary]["type"] == "dirichlet":
				values = self.input_file["boundary_conditions"][boundary]["values"]
				value = np.interp(t, self.time_list, values)
				value_dirichlet = do.Expression("value", value=value, degree=1)
				component = self.input_file["boundary_conditions"][boundary]["component"]
				self.bcs.append(do.DirichletBC(self.CG_3x1.sub(component), value_dirichlet, self.grid.get_boundaries(), self.grid.get_boundary_tags(boundary)))


	def define_solver(self):
		"""
		Defines the solver for solving the linear system according to the specifications
		in the input_file.json.

		Returns
		-------
		solver : dolfin.cpp.la.KrylovSolver or dolfin.cpp.la.LUSolver
			The solver can be either an iterative solver (KrylovSolver) or a direct solver
			(LUSolver).
		"""
		solver_type = self.input_file["solver_settings"]["type"]
		if solver_type == "KrylovSolver":
			solver = do.KrylovSolver(
			                      	method = self.input_file["solver_settings"]["method"],
			                      	preconditioner = self.input_file["solver_settings"]["preconditioner"]
			                      )
			solver.parameters["relative_tolerance"] = self.input_file["solver_settings"]["relative_tolerance"]
		elif solver_type == "LU":
			solver = do.LUSolver(self.input_file["solver_settings"]["method"])
		else:
			raise Exception(f"Solver type {solver_type} not supported. Choose between KrylovSolver and LU.")
		return solver









