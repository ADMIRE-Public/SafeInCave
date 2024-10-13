"""
The class implements the iterative process to solve the non-linear equilibrium equations.
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

from ConstitutiveModel import ConstitutiveModel
from Equations import LinearMomentum
from Grid import GridHandlerGMSH
import Utils as utils
import dolfin as do
import torch as to
import numpy as np
import os
import copy

class Simulator(object):
	"""
	This class is responsible to carry out the simulation according to the
	input_file.json. It loads the mesh, builds the constitutive model, builds
	the linear momenum equation, defines the solver, defines the weak formulation
	the run the simulation.

	Parameters
	----------
	input_file : dict
		Dictionary extracted from the JSON file.
	"""
	def __init__(self, input_file):
		self.input_file = input_file

		# This is input_file to be saved
		self.input_file_to_be_saved = copy.deepcopy(input_file)

		# Output folder
		self.output_folder = input_file["output"]["path"]

		# Transient settings
		self.time_list = input_file["time_settings"]["time_list"]
		theta = input_file["time_settings"]["theta"]
	
		# Create mesh
		self.grid = GridHandlerGMSH(input_file["grid"]["name"], input_file["grid"]["path"])

		# Linear momentum balance equation
		self.eq_mom = LinearMomentum(self.grid, theta, self.input_file)

		# Save input file
		filename = os.path.join(os.path.join(input_file["output"]["path"], "input_file.json"))
		os.makedirs(os.path.dirname(filename), exist_ok=True)
		utils.save_json(self.input_file_to_be_saved, filename)


	def run(self):
		"""
		Runs transient simulation.
		"""
		if self.input_file["simulation_settings"]["equilibrium"]["active"] == True:
			self.run_equilibrium()
			utils.save_json(self.input_file_to_be_saved, os.path.join(self.output_folder, "equilibrium", "input_file.json"))
			self.run_operation()
		else:
			self.run_simulation()

		# Save input file
		utils.save_json(self.input_file_to_be_saved, os.path.join(self.output_folder, "operation", "input_file.json"))


	def run_simulation(self, verbose=True):
		"""
		Runs simulation **without** solving the equilibrium condition.
		"""
		# Pseudo time
		t = self.time_list[0]
		t_final = self.time_list[-1]

		# Get maximum time step size
		dt = self.input_file["simulation_settings"]["operation"]["dt_max"]

		# Read number of time steps to skip before saving results
		n_skip = self.input_file["simulation_settings"]["operation"]["n_skip"]

		# Perform initial computations
		self.eq_mom.initialize(verbose)

		# Save initial solution
		self.eq_mom.save_solution(t)

		# Transient simulation
		n_step = 1
		while t < t_final:

			# Increment time
			t += dt

			# Solve
			self.eq_mom.solve(t, dt)

			# Update internal variables
			self.eq_mom.update_internal_variables()

			# Compute strains
			self.eq_mom.compute_eps_ie(dt)
			self.eq_mom.compute_eps_ve(dt)

			# Update old non-elastic strains
			self.eq_mom.update_eps_ie_old()
			self.eq_mom.update_eps_ve_old()

			# Update old non-elastic strain rates
			self.eq_mom.update_eps_ie_rate_old()
			self.eq_mom.update_eps_ve_rate_old()

			# Save displacement field
			if n_step % n_skip == 0 or n_step == 1:
				self.eq_mom.save_solution(t)
				if verbose:
					print("Save step %i"%n_step)

			# Print stuff
			if verbose:
				print(n_step, f"{t_final/utils.hour}", t/utils.hour, self.eq_mom.ite, self.eq_mom.error)
				try:
					print("Fvp: ", float(max(self.eq_mom.m.elems_ie[0].Fvp)))
					print("alpha: ", float(max(self.eq_mom.m.elems_ie[0].alpha)))
				except:
					pass
				print()

			n_step += 1


	def run_equilibrium(self):
		pass

	def run_equilibrium2(self):
		"""
		Runs simulation to determine the equilibrium condition. The equilibrium condition
		only considers the **elastic** and **viscoelastic** (if present) elements of the 
		constitutive model.
		"""
		# Pseudo time
		t_0 = self.time_list[0]
		t = t_0

		# Get maximum time step size
		dt = self.input_file["simulation_settings"]["equilibrium"]["dt_max"]

		# Output folder
		equilibrium_output_folder = os.path.join(self.input_file["output"]["path"], "equilibrium")

		# Define Dirichlet boundary conditions
		bcs = self.define_dirichlet_bc(t_0)
		
		# Apply Neumann boundary conditions
		b_outer = self.apply_neumann_bc(t_0)
		
		# Build RHS vector
		b = do.assemble(self.b_body + b_outer)

		# Initialize elastic stiffness matrix, C0
		self.C0.vector()[:] = to.flatten(self.eq.m.C0)

		# Build stiffness matrix
		a_form = do.inner(utils.dotdot(self.C0, utils.epsilon(self.du)), utils.epsilon(self.v))*self.dx
		A = do.assemble(a_form)

		# Solve linear system
		[bc.apply(A, b) for bc in bcs]
		self.solver.solve(A, self.u.vector(), b)

		# Compute total strain
		self.eps_tot.assign(utils.local_projection(utils.epsilon(self.u), self.DG_3x3))

		# Compute stress
		eps_tot_torch = utils.numpy2torch(self.eps_tot.vector()[:].reshape((self.n_elems, 3, 3)))
		self.eq.compute_stress_C0(eps_tot_torch)

		# Compute old viscoelastic strain rates
		self.eq.compute_eps_ve_rate(0)

		# Update viscoelastic strain rate (Warning! Do NOT update eps_ie here, because this is wrong!)
		self.eq.update_eps_ve_rate_old()

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
			self.eq.compute_GT_BT_ve(dt)

			# Iterative loop settings
			tol = 1e-7
			error = 2*tol
			ite = 0
			maxiter = 40

			while error > tol and ite < maxiter:

				# Update total strain of previous iteration (eps_tot_k <-- eps_tot)
				eps_tot_k = utils.numpy2torch(self.eps_tot.vector()[:])

				# Update stress of previous iteration (stress_k <-- stress)
				self.eq.update_stress()

				# Compute CT
				self.eq.compute_CT(dt)
				self.CT.vector()[:] = to.flatten(self.eq.CT)

				# Compute right-hand side
				self.eq.compute_eps_rhs(dt)

				# Assign eps_rhs
				self.eps_rhs.vector()[:] = self.eq.eps_rhs.flatten()

				# Build rhs
				b_rhs = do.inner(utils.dotdot(self.CT, self.eps_rhs), utils.epsilon(self.v))*self.dx
				b = do.assemble(self.b_body + b_outer + b_rhs)

				# Build lhs
				a_form = do.inner(utils.dotdot(self.CT, utils.epsilon(self.du)), utils.epsilon(self.v))*self.dx
				A = do.assemble(a_form)

				# Solve linear system
				[bc.apply(A, b) for bc in bcs]
				self.solver.solve(A, self.u.vector(), b)

				# Compute total strain
				self.eps_tot.assign(utils.local_projection(utils.epsilon(self.u), self.DG_3x3))
				eps_tot_torch = utils.numpy2torch(self.eps_tot.vector()[:].reshape((self.n_elems, 3, 3)))

				# Compute stress
				self.eq.compute_stress(eps_tot_torch, dt)

				# Compute strain rates
				self.eq.compute_eps_ve_rate(dt)

				# Compute error
				if self.eq.theta == 1.0:
					error = 0.0
				else:
					eps_tot_k_flat = to.flatten(eps_tot_k)
					eps_tot_flat = self.eps_tot.vector()[:]
					error = np.linalg.norm(eps_tot_k_flat - eps_tot_flat) / np.linalg.norm(eps_tot_flat)

				# Increment iteration counter
				ite += 1

			# Compute strains
			self.eq.compute_eps_ve(dt)

			# Update old viscoelastic strains
			self.eq.update_eps_ve_old()

			# Update old viscoelastic strain rates
			self.eq.update_eps_ve_rate_old()

			# Compute time error (to check if steady state is achieved)
			eps_tot_flat = self.eps_tot.vector()[:]
			error_time = np.linalg.norm(eps_tot_old - eps_tot_flat) / np.linalg.norm(eps_tot_flat)
			eps_tot_old = utils.numpy2torch(self.eps_tot.vector()[:])

			# Save displacement field
			u_vtk << (self.u, t)
			self.sigma.vector()[:] = self.eq.stress.flatten()
			stress_vtk << (self.sigma, t)

			# Print stuff
			print(n_step, t/utils.hour, ite, error, error_time)
			# print()
			n_step += 1



	def run_operation(self):
		"""
		Runs transient simulation **after** the equilibrium condition. 
		"""
		# Pseudo time
		t = self.time_list[0]
		t_final = self.time_list[-1]

		# Get maximum time step size
		dt = self.input_file["simulation_settings"]["operation"]["dt_max"]

		# Read number of time steps to skip before saving results
		n_skip = self.input_file["simulation_settings"]["operation"]["n_skip"]

		# Output folder
		operation_output_folder = os.path.join(self.input_file["output"]["path"], "operation")

		# Define Dirichlet boundary conditions
		bcs = self.define_dirichlet_bc(t)

		# Apply Neumann boundary conditions
		b_outer = self.apply_neumann_bc(t)

		# Build RHS vector
		b = do.assemble(self.b_body + b_outer)

		# Compute initial hardening
		for elem in self.eq.m.elems_ie:
			try:
				# Compute initial hardening parameter (alpha_0) based on initial stresses
				elem.compute_initial_hardening(self.eq.stress, Fvp_0=0.0)
				
				# Compute initial yield function values
				I1, I2, I3, J2, J3, Sr, I1_star = elem.compute_stress_invariants(*elem.extract_stress_components(self.eq.stress))
				_ = elem.compute_Fvp(elem.alpha, I1_star, J2, Sr)
				# Fvp.vector()[:] = elem.compute_Fvp(elem.alpha, I1_star, J2, Sr)

				print("Fvp: ", float(max(elem.Fvp)))
				print("alpha_min: ", float(min(elem.alpha)))
				print("alpha_max: ", float(max(elem.alpha)))
				print("alpha_avg: ", float(np.average(elem.alpha)))
				print()
			except:
				pass

		# Compute total strain
		self.eps_tot.assign(utils.local_projection(utils.epsilon(self.u), self.DG_3x3))

		# # Compute stress
		# eps_tot_torch = utils.numpy2torch(self.eps_tot.vector()[:].reshape((self.n_elems, 3, 3)))
		# self.eq.compute_stress_C0(eps_tot_torch)

		# Compute old ielastic strain rates
		self.eq.compute_eps_ie_rate()
		# m.compute_eps_ve_rate(0)

		# Update inelastic strain rate (Warning! Do NOT update eps_ie here, because this is wrong!)
		self.eq.update_eps_ie_rate_old()
		# m.update_eps_ve_rate_old()

		# Assign stress
		self.sigma.vector()[:] = self.eq.stress.flatten()

		# Create output file
		u_vtk = do.File(os.path.join(operation_output_folder, "vtk", "displacement", "displacement.pvd"))
		stress_vtk = do.File(os.path.join(operation_output_folder, "vtk", "stress", "stress.pvd"))

		# Save output fields
		u_vtk << (self.u, t)
		stress_vtk << (self.sigma, t)

		# Transient simulation
		n_step = 1
		while t < t_final:

			# Increment time
			t += dt

			# Define Dirichlet boundary conditions
			bcs = self.define_dirichlet_bc(t)

			# Apply Neumann boundary conditions
			b_outer = self.apply_neumann_bc(t)

			# Compute GT and BT matrix fields for viscoelastic elements
			self.eq.compute_GT_BT_ve(dt)

			# Iterative loop settings
			tol = 1e-7
			error = 2*tol
			ite = 0
			maxiter = 40

			while error > tol and ite < maxiter:

				# Update total strain of previous iteration (eps_tot_k <-- eps_tot)
				eps_tot_k = utils.numpy2torch(self.eps_tot.vector()[:])

				# Update stress of previous iteration (stress_k <-- stress)
				self.eq.update_stress()

				# Compute GT and BT matrix fields for inelastic elements
				self.eq.compute_GT_BT_ie(dt)

				# Compute CT
				self.eq.compute_CT(dt)
				self.CT.vector()[:] = to.flatten(self.eq.CT)

				# Compute right-hand side
				self.eq.compute_eps_rhs(dt)

				# Assign eps_rhs
				self.eps_rhs.vector()[:] = self.eq.eps_rhs.flatten()

				# Build rhs
				b_rhs = do.inner(utils.dotdot(self.CT, self.eps_rhs), utils.epsilon(self.v))*self.dx
				b = do.assemble(self.b_body + b_outer + b_rhs)

				# Build lhs
				a_form = do.inner(utils.dotdot(self.CT, utils.epsilon(self.du)), utils.epsilon(self.v))*self.dx
				A = do.assemble(a_form)

				# Solve linear system
				[bc.apply(A, b) for bc in bcs]
				self.solver.solve(A, self.u.vector(), b)

				# Compute total strain
				self.eps_tot.assign(utils.local_projection(utils.epsilon(self.u), self.DG_3x3))
				eps_tot_torch = utils.numpy2torch(self.eps_tot.vector()[:].reshape((self.n_elems, 3, 3)))

				# Compute stress
				self.eq.compute_stress(eps_tot_torch, dt)

				# Increment internal variables
				self.eq.increment_internal_variables(dt)

				# Compute strain rates
				self.eq.compute_eps_ie_rate()
				self.eq.compute_eps_ve_rate(dt)

				# Compute error
				if self.eq.theta == 1.0:
					error = 0.0
				else:
					eps_tot_k_flat = to.flatten(eps_tot_k)
					eps_tot_flat = self.eps_tot.vector()[:]
					error = np.linalg.norm(eps_tot_k_flat - eps_tot_flat) / np.linalg.norm(eps_tot_flat)

				# Increment iteration counter
				ite += 1

			# Update internal variables
			self.eq.update_internal_variables()

			# Compute strains
			self.eq.compute_eps_ie(dt)
			self.eq.compute_eps_ve(dt)

			# Update old non-elastic strains
			self.eq.update_eps_ie_old()
			self.eq.update_eps_ve_old()

			# Update old non-elastic strain rates
			self.eq.update_eps_ie_rate_old()
			self.eq.update_eps_ve_rate_old()

			# Save displacement field
			if n_step % n_skip == 0 or n_step == 1:
				u_vtk << (self.u, t)
				self.sigma.vector()[:] = self.eq.stress.flatten()
				stress_vtk << (self.sigma, t)
				print("Save step %i"%n_step)

			# Print stuff
			print(n_step, f"{t_final/utils.hour}", t/utils.hour, ite, error)
			try:
				print("Fvp: ", float(max(self.eq.m.elems_ie[0].Fvp)))
				print("alpha: ", float(max(self.eq.m.elems_ie[0].alpha)))
			except:
				pass
			print()
			n_step += 1

