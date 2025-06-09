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
import copy
import torch as to
import dolfinx as do
from petsc4py import PETSc
import ufl
import numpy as np
import Utils as utils
from ConstitutiveModel import ConstitutiveModel
from ScreenOutput import ScreenPrinter

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
		self.x = ufl.SpatialCoordinate(self.grid.mesh)

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
		self.CG0_3x1 = do.fem.functionspace(self.grid.mesh, ("Lagrange", 1, (grid.domain_dim, )))
		self.DG0_1 = do.fem.functionspace(self.grid.mesh, ("DG", 0))
		self.DG0_3x3 = do.fem.functionspace(self.grid.mesh, ("DG", 0, (3, 3)))
		self.DG0_6x6 = do.fem.functionspace(self.grid.mesh, ("DG", 0, (6, 6)))


		# Create tensor fields
		self.C0 = do.fem.Function(self.DG0_6x6)
		self.C1 = do.fem.Function(self.DG0_6x6)
		self.CT = do.fem.Function(self.DG0_6x6)
		self.eps_tot = do.fem.Function(self.DG0_3x3)
		self.eps_rhs = do.fem.Function(self.DG0_3x3)

		self.sigma = do.fem.Function(self.DG0_3x3)
		self.sigma_0 = do.fem.Function(self.DG0_3x3)
		self.sigma.name = "Stress"

		self.sigma_v = do.fem.Function(self.DG0_1)
		self.sigma_v.name = "Mean stress"

		self.von_mises = do.fem.Function(self.DG0_1)
		self.von_mises.name = "Von Mises stress"

		# Define variational problem
		self.du = ufl.TrialFunction(self.CG0_3x1)
		self.v = ufl.TestFunction(self.CG0_3x1)
		self.ds = ufl.Measure("ds", domain=self.grid.mesh, subdomain_data=self.grid.get_boundaries())
		self.dx = ufl.Measure("dx", domain=self.grid.mesh, subdomain_data=self.grid.get_subdomains())

		# Define outward normal vectors
		self.n = ufl.FacetNormal(self.grid.mesh)
		self.normal = ufl.dot(self.n, self.v)

		# Create displacement vector
		self.u = do.fem.Function(self.CG0_3x1)
		self.u.name = "Displacement"

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
		self.density = do.fem.Function(self.DG0_1)
		self.density.x.array[:] = self.grid.get_parameter(self.input_file["body_force"]["density"])
		self.gravity = self.input_file["body_force"]["gravity"]

		self.g = [0.0, 0.0, 0.0]
		self.g[self.direction] = self.gravity
		self.body_force = self.density*do.fem.Constant(self.grid.mesh, do.default_scalar_type(tuple(self.g)))

		# Build body forces on RHS vector
		self.b_body = ufl.dot(self.body_force, self.v)*self.dx

		# Define linear solver
		self.solver = self.define_solver()

		# print("check")
		# self.test_vtk = do.io.VTKFile(self.grid.mesh.comm, os.path.join(self.operation_output_folder, "vtk", "test", "test.pvd"), "w")
		# print("check")
		# self.test_vtk.write_mesh(self.grid.mesh)
		# print("check")



	def initialize_ouput_files(self, output_folder):
		# Create output file
		self.u_vtk = do.io.XDMFFile(self.grid.mesh.comm, os.path.join(output_folder, "vtk", "displacement", "u.xdmf"), "w")
		self.stress_vtk = do.io.XDMFFile(self.grid.mesh.comm, os.path.join(output_folder, "vtk", "stress", "stress.xdmf"), "w")
		self.q_vtk = do.io.XDMFFile(self.grid.mesh.comm, os.path.join(output_folder, "vtk", "q", "q.xdmf"), "w")
		self.p_vtk = do.io.XDMFFile(self.grid.mesh.comm, os.path.join(output_folder, "vtk", "p", "p.xdmf"), "w")

		self.u_vtk.name = "Displacement"
		self.stress_vtk.name = "Stress"
		self.q_vtk.name = "Von Mises stress"
		self.p_vtk.name = "Mean stress"

		self.u_vtk.write_mesh(self.grid.mesh)
		self.stress_vtk.write_mesh(self.grid.mesh)
		self.q_vtk.write_mesh(self.grid.mesh)
		self.p_vtk.write_mesh(self.grid.mesh)


	def save_solution(self, t):

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

		self.von_mises.x.array[:] = self.grid.smoother.dot(q.numpy())
		self.q_vtk.write_function(self.von_mises, t)

		self.sigma_v.x.array[:] = self.grid.smoother.dot(p.numpy())
		self.p_vtk.write_function(self.sigma_v, t)

		self.sigma.x.array[:] = self.apply_smoother(self.stress_torch.numpy())
		self.stress_vtk.write_function(self.sigma, t)

		self.u_vtk.write_function(self.u, t)


	def apply_smoother(self, field_np):
		field_copy = field_np.copy()
		field_copy[:,0,0] = self.grid.smoother.dot(field_copy[:,0,0])
		field_copy[:,1,1] = self.grid.smoother.dot(field_copy[:,1,1])
		field_copy[:,2,2] = self.grid.smoother.dot(field_copy[:,2,2])
		field_copy[:,0,1] = field_copy[:,1,0] = self.grid.smoother.dot(field_copy[:,0,1])
		field_copy[:,0,2] = field_copy[:,2,0] = self.grid.smoother.dot(field_copy[:,0,2])
		field_copy[:,1,2] = field_copy[:,2,1] = self.grid.smoother.dot(field_copy[:,1,2])
		return field_copy.flatten()



	def initialize(self, solve_equilibrium=False, verbose=True, save_results=False, calculate_hardening=True):
		# Apply Neumann boundary condition
		self.apply_neumann_bc(self.t0)

		# Define Dirichlet boundary conditions
		self.define_dirichlet_bc(self.t0)

		# Initialize elastic stiffness matrix, C0
		self.C0.x.array[:] = to.flatten(self.m.C0)

		# Build lhs
		a = ufl.inner(utils.dotdot(self.C0, utils.epsilon(self.du)), utils.epsilon(self.v))*self.dx
		bilinear_form = do.fem.form(a)
		A = do.fem.petsc.assemble_matrix(bilinear_form, bcs=self.bcs)
		A.assemble()

		# Build rhs
		linear_form = do.fem.form(self.b_body + sum(self.integral_neumann))
		b = do.fem.petsc.assemble_vector(linear_form)
		do.fem.petsc.apply_lifting(b, [bilinear_form], [self.bcs])
		do.fem.petsc.set_bc(b, self.bcs)

		# Solve linear system
		self.solver.setOperators(A)
		self.solver.solve(b, self.u.x.petsc_vec)

		# Compute total strain
		self.eps_tot = utils.project(utils.epsilon(self.u), self.DG0_3x3)
		eps_tot_torch = utils.numpy2torch(self.eps_tot.x.array.reshape((self.n_elems, 3, 3)))

		# Compute stress
		self.compute_stress_C0(eps_tot_torch)

		if solve_equilibrium:
			self.solve_equilibrium(verbose, save_results)

		# Compute initial hardening
		if calculate_hardening:
			self.compute_initial_hardening(verbose)

		# Compute old ielastic strain rates
		self.compute_eps_ie_rate()

		# Update inelastic strain rate (Warning! Do NOT update eps_ie here, because this is wrong!)
		self.update_eps_ie_rate_old()

		# Re-initialize output files
		self.initialize_ouput_files(self.operation_output_folder)

		# Assign stress
		self.sigma.x.array[:] = self.stress_torch.flatten()



	def solve_equilibrium(self, verbose=True, save_results=False):
		# Build constitutive model for equilibrium stage
		mod_input_file = copy.deepcopy(self.input_file["constitutive_model"])
		elem_names = []
		for elem_type in mod_input_file.keys():
			for elem_name in mod_input_file[elem_type].keys():
				if mod_input_file[elem_type][elem_name]["active"] == True:
					if mod_input_file[elem_type][elem_name]["equilibrium"] == False:
						mod_input_file[elem_type][elem_name]["active"] = False
					else:
						elem_names.append(elem_name)
						mod_input_file[elem_type][elem_name]["active"] = True
		self.m = ConstitutiveModel(self.grid, mod_input_file)

		# Screen info
		screen = ScreenPrinter()
		screen.set_header_columns(["Time step", "Pseudo-time (h)", "Error time"], "center")
		screen.set_row_formats(["%s", "%.3f", "%.4e"], ["center", "center", "center"])
		screen.start_timer()
		screen.print_on_screen(screen.master_division_plus)
		screen.print_on_screen(" ")
		screen.print_on_screen(screen.master_division_plus)
		screen.print_comment(" Equilibrium Stage", "center")
		screen.print_on_screen(screen.master_division_plus)
		screen.print_comment(" ")
		screen.print_comment(" Constitutive model:")
		for elem_name in elem_names:
			screen.print_comment(f"          {elem_name}")
		screen.print_comment(" ")
		screen.print_header()
		# screen.print(screen.divider)

		self.compute_eps_ie_rate()

		# Compute old viscoelastic strain rates
		self.compute_eps_ve_rate(0)

		# Update viscoelastic strain rate (Warning! Do NOT update eps_ie here, because this is wrong!)
		self.update_eps_ve_rate_old()

		# Get maximum time step size
		dt = self.input_file["simulation_settings"]["equilibrium"]["dt_max"]
		t = self.t0

		if save_results:
			# Output folder
			equilibrium_output_folder = os.path.join(self.input_file["output"]["path"], "equilibrium")

			# Initialize output files
			self.initialize_ouput_files(equilibrium_output_folder)

			# Save output fields
			self.save_solution(t)

		# Initialize pseudo time control settings
		n_step = 1
		tol_time = self.input_file["simulation_settings"]["equilibrium"]["time_tol"]
		ite_max = self.input_file["simulation_settings"]["equilibrium"]["ite_max"]
		error_time = 2*tol_time
		eps_tot_old = utils.numpy2torch(self.eps_tot.x.array[:])

		while error_time > tol_time or n_step <= 2:

			# Increment time
			t += dt

			self.solve(0, dt)

			# Compute time error (to check if steady state is achieved)
			eps_tot_flat = self.eps_tot.x.array
			error_time = np.linalg.norm(eps_tot_old - eps_tot_flat) / np.linalg.norm(eps_tot_flat)
			eps_tot_old = utils.numpy2torch(self.eps_tot.x.array)
			self.sigma.x.array[:] = self.stress_torch.flatten()

			if save_results:
				self.save_solution(t)

			# Print stuff
			if verbose:
				screen_output_row = [f"{n_step}/{ite_max}", t/utils.hour, error_time]
				screen.print_row(screen_output_row)
			n_step += 1

			if n_step > ite_max:
				break

		# Reinitialize original constitutive model
		m_operation = ConstitutiveModel(self.grid, self.input_file["constitutive_model"])

		# Copy elements active on the equilibrium stage to 
		# the constitutive model of the operation stage
		for i, elem in enumerate(self.input_file["constitutive_model"]["elastic"].keys()):
			if self.input_file["constitutive_model"]["elastic"][elem]["active"] == True:
				if self.input_file["constitutive_model"]["elastic"][elem]["equilibrium"] == True:
					m_operation._elems_e[i] = self.m.elems_e[i]

		for i, elem in enumerate(self.input_file["constitutive_model"]["viscoelastic"].keys()):
			if self.input_file["constitutive_model"]["viscoelastic"][elem]["active"] == True:
				if self.input_file["constitutive_model"]["viscoelastic"][elem]["equilibrium"] == True:
					m_operation._elems_ve[i] = self.m.elems_ve[i]

		for i, elem in enumerate(self.input_file["constitutive_model"]["inelastic"].keys()):
			if self.input_file["constitutive_model"]["inelastic"][elem]["active"] == True:
				if self.input_file["constitutive_model"]["inelastic"][elem]["equilibrium"] == True:
					m_operation._elems_ie[i] = self.m.elems_ie[i]

		# Assign operation constitutive model to the internal object self.m
		self.m = m_operation

		screen.close()
		# screen.save_log(equilibrium_output_folder)



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
			eps_tot_k = utils.numpy2torch(self.eps_tot.x.array[:])

			# Update stress of previous iteration (stress_k <-- stress)
			self.update_stress()

			# Compute GT and BT matrix fields for inelastic elements
			self.compute_GT_BT_ie(dt)

			# Compute CT
			self.compute_CT(dt)

			# Compute right-hand side
			self.compute_eps_rhs(dt)
			

			# Build lhs
			a = ufl.inner(utils.dotdot(self.CT, utils.epsilon(self.du)), utils.epsilon(self.v))*self.dx
			bilinear_form = do.fem.form(a)
			A = do.fem.petsc.assemble_matrix(bilinear_form, bcs=self.bcs)
			A.assemble()

			# Build rhs
			b_rhs = ufl.inner(utils.dotdot(self.CT, self.eps_rhs), utils.epsilon(self.v))*self.dx
			linear_form = do.fem.form(self.b_body + sum(self.integral_neumann) + b_rhs)
			b = do.fem.petsc.assemble_vector(linear_form)
			do.fem.petsc.apply_lifting(b, [bilinear_form], [self.bcs])
			do.fem.petsc.set_bc(b, self.bcs)

			# Solve linear system
			self.solver.setOperators(A)
			self.solver.solve(b, self.u.x.petsc_vec)

			# Compute total strain
			self.eps_tot = utils.project(utils.epsilon(self.u), self.DG0_3x3)

			# self.eps_tot.assign(utils.local_projection(utils.epsilon(self.u), self.DG0_3x3))
			eps_tot_torch = utils.numpy2torch(self.eps_tot.x.array.reshape((self.n_elems, 3, 3)))

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
				eps_tot_flat = self.eps_tot.x.array
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

				if float(min(elem.alpha)) < 0:
					print("Warning! Negative hardening parameter for Desai's model.")
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
		self.CT.x.array[:] = to.flatten(self.CT_torch)

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
		self.eps_rhs.x.array[:] = self.eps_rhs_torch.flatten()

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
		self.sigma.x.array[:] = self.stress_torch.flatten()

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
		bc_neumann_list = []
		for boundary in self.input_file["boundary_conditions"]:
			if self.input_file["boundary_conditions"][boundary]["type"] == "neumann":
				i = self.input_file["boundary_conditions"][boundary]["direction"]
				rho = self.input_file["boundary_conditions"][boundary]["density"]
				H = self.input_file["boundary_conditions"][boundary]["reference_position"]
				values = self.input_file["boundary_conditions"][boundary]["values"]
				p = -np.interp(t, self.time_list, values)
				value_neumann = p + rho*self.gravity*(H - self.x[i])
				self.integral_neumann.append(value_neumann*self.normal*self.ds(self.grid.get_boundary_tag(boundary)))



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
				component = self.input_file["boundary_conditions"][boundary]["component"]
				values = self.input_file["boundary_conditions"][boundary]["values"]
				value = np.interp(t, self.time_list, values)
				dofs = do.fem.locate_dofs_topological(
					self.CG0_3x1.sub(component),
					self.grid.boundary_dim,
					self.grid.get_boundary_tags(boundary)
				)
				self.bcs.append(
					do.fem.dirichletbc(
						do.default_scalar_type(value),
						dofs,
						self.CG0_3x1.sub(component)
					)
				)

	def define_solver(self):
		"""
		Defines the solver for the linear system according to the specifications
		in the input_file.json.

		Returns
		-------
		solver : 
		"""
		solver = PETSc.KSP().create(self.grid.mesh.comm)
		solver.setType(self.input_file["solver_settings"]["solver_type"])
		solver.getPC().setType(self.input_file["solver_settings"]["solver_PC"])
		solver.setTolerances(rtol=self.input_file["solver_settings"]["rtol"], 
							 max_it=self.input_file["solver_settings"]["maxite"])
		return solver









