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
import torch as to
import numpy as np
from mpi4py import MPI
from Utils import numpy2torch

class Simulator(ABC):
	@abstractmethod
	def run(self):
		pass


class Simulator_TM(Simulator):
	def __init__(self, eq_mom, eq_heat, t_control, outputs, compute_elastic_response=True):
		self.eq_mom = eq_mom
		self.eq_heat = eq_heat
		self.t_control = t_control
		self.outputs = outputs
		self.compute_elastic_response = compute_elastic_response

	def run(self):
		# Output field
		for output in self.outputs:
			output.initialize()

		# Set initial temperature
		T_elems = self.eq_heat.get_T_elems()
		self.eq_mom.set_T0(T_elems)

		# Update boundary conditions
		self.eq_mom.bc.update_dirichlet(self.t_control.t)
		self.eq_mom.bc.update_neumann(self.t_control.t)

		if self.compute_elastic_response:
			# Solve elasticity
			self.eq_mom.solve_elastic_response()

			# Calculate total (elastic) strain
			eps_tot_to = self.eq_mom.compute_total_strain()

			# Compute stress
			stress_to = self.eq_mom.compute_elastic_stress(eps_tot_to)

		else:
			# Calculate total strain
			eps_tot_to = self.eq_mom.compute_total_strain()

			# Retrieve stress
			stress_to = numpy2torch(self.eq_mom.sig.x.array.reshape((self.eq_mom.n_elems, 3, 3)))

		# Set new temperature to momentum equation
		T_elems = self.eq_heat.get_T_elems()
		self.eq_mom.set_T(T_elems)
		self.eq_mom.set_T0(T_elems)

		# Calculate and eps_ie_rate_old
		self.eq_mom.compute_eps_ne_rate(stress_to, self.t_control.t)
		self.eq_mom.update_eps_ne_rate_old()


		# self.eq_heat.solve(0, self.t_control.dt)

		# Save fields
		self.eq_mom.compute_p_elems()
		self.eq_mom.compute_q_elems()
		self.eq_mom.compute_p_nodes()
		self.eq_mom.compute_q_nodes()
		output.save_fields(0)

		# Time loop
		while self.t_control.keep_looping():

			# Advance time
			self.t_control.advance_time()
			t = self.t_control.t
			dt = self.t_control.dt

			# Update boundary conditions
			self.eq_mom.bc.update_dirichlet(t)
			self.eq_mom.bc.update_neumann(t)
			self.eq_heat.bc.update_dirichlet(t)
			self.eq_heat.bc.update_neumann(t)

			# Solve heat
			self.eq_heat.solve(t, dt)

			# Set new temperature to momentum equation
			T_elems = self.eq_heat.get_T_elems()
			self.eq_mom.set_T(T_elems)

			# Iterative loop settings
			tol = 1e-6
			error = 2*tol
			ite = 0
			maxiter = 20

			while error > tol and ite < maxiter:

				# Update total strain of previous iteration (eps_tot_k <-- eps_tot)
				eps_tot_k_to = eps_tot_to.clone()

				# Update stress
				stress_k_to = stress_to.clone()

				# Build bi-linear form
				self.eq_mom.solve(stress_k_to, t, dt)

				# Compute total strain
				eps_tot_to = self.eq_mom.compute_total_strain()

				# Compute stress
				stress_to = self.eq_mom.compute_stress(eps_tot_to)

				# Increment internal variables
				self.eq_mom.increment_internal_variables(stress_to, stress_k_to, dt)

				# Compute inelastic strain rates
				self.eq_mom.compute_eps_ne_rate(stress_to, dt)

				# Compute error
				if self.eq_mom.theta == 1.0:
					error = 0.0
				else:
					eps_tot_k_flat = to.flatten(eps_tot_k_to)
					eps_tot_flat = to.flatten(eps_tot_to)
					local_error =  np.linalg.norm(eps_tot_k_flat - eps_tot_flat) / np.linalg.norm(eps_tot_flat)
					error = self.eq_mom.grid.mesh.comm.allreduce(local_error, op=MPI.SUM)

				ite += 1

			# Update internal variables
			self.eq_mom.update_internal_variables()

			# Update strain rates
			self.eq_mom.update_eps_ne_rate_old()

			# Update strain
			self.eq_mom.update_eps_ne_old(stress_to, stress_k_to, dt)

			# Print stuff
			if self.eq_mom.grid.mesh.comm.rank == 0:
				print(t/self.t_control.time_unit, ite, error)

			# Save fields
			self.eq_mom.compute_p_elems()
			self.eq_mom.compute_q_elems()
			self.eq_mom.compute_p_nodes()
			self.eq_mom.compute_q_nodes()
			for output in self.outputs:
				output.save_fields(t)

		for output in self.outputs:
			output.save_mesh()


class Simulator_M(Simulator):
	def __init__(self, eq_mom, t_control, outputs, compute_elastic_response=True):
		self.eq_mom = eq_mom
		self.t_control = t_control
		self.outputs = outputs
		self.compute_elastic_response = compute_elastic_response

	def run(self):
		# Output field
		for output in self.outputs:
			output.initialize()

		# Update boundary conditions
		self.eq_mom.bc.update_dirichlet(self.t_control.t)
		self.eq_mom.bc.update_neumann(self.t_control.t)

		if self.compute_elastic_response:
			# Solve elasticity
			self.eq_mom.solve_elastic_response()

			# Calculate total (elastic) strain
			eps_tot_to = self.eq_mom.compute_total_strain()

			# Compute stress
			stress_to = self.eq_mom.compute_elastic_stress(eps_tot_to)

		else:
			# Calculate total strain
			eps_tot_to = self.eq_mom.compute_total_strain()

			# Retrieve stress
			stress_to = numpy2torch(self.eq_mom.sig.x.array.reshape((self.eq_mom.n_elems, 3, 3)))

		# Calculate and eps_ie_rate_old
		self.eq_mom.compute_eps_ne_rate(stress_to, self.t_control.t)
		self.eq_mom.update_eps_ne_rate_old()

		# Save fields
		self.eq_mom.compute_p_elems()
		self.eq_mom.compute_q_elems()
		self.eq_mom.compute_p_nodes()
		self.eq_mom.compute_q_nodes()
		output.save_fields(0)

		# Time loop
		while self.t_control.keep_looping():

			# Advance time
			self.t_control.advance_time()
			t = self.t_control.t
			dt = self.t_control.dt

			# Update boundary conditions
			self.eq_mom.bc.update_dirichlet(t)
			self.eq_mom.bc.update_neumann(t)

			# Iterative loop settings
			tol = 1e-8
			error = 2*tol
			ite = 0
			maxiter = 40

			while error > tol and ite < maxiter:

				# Update total strain of previous iteration (eps_tot_k <-- eps_tot)
				eps_tot_k_to = eps_tot_to.clone()

				# Update stress
				stress_k_to = stress_to.clone()

				# Build bi-linear form
				self.eq_mom.solve(stress_k_to, t, dt)

				# Compute total strain
				eps_tot_to = self.eq_mom.compute_total_strain()

				# Compute stress
				stress_to = self.eq_mom.compute_stress(eps_tot_to)

				# Increment internal variables
				self.eq_mom.increment_internal_variables(stress_to, stress_k_to, dt)

				# Compute inelastic strain rates
				self.eq_mom.compute_eps_ne_rate(stress_to, dt)

				# Compute error
				if self.eq_mom.theta == 1.0:
					error = 0.0
				else:
					eps_tot_k_flat = to.flatten(eps_tot_k_to)
					eps_tot_flat = to.flatten(eps_tot_to)
					local_error =  np.linalg.norm(eps_tot_k_flat - eps_tot_flat) / np.linalg.norm(eps_tot_flat)
					error = self.eq_mom.grid.mesh.comm.allreduce(local_error, op=MPI.SUM)

				ite += 1

			# Update internal variables
			self.eq_mom.update_internal_variables()

			# Update strain rates
			self.eq_mom.update_eps_ne_rate_old()

			# Update strain
			self.eq_mom.update_eps_ne_old(stress_to, stress_k_to, dt)

			# Save fields
			self.eq_mom.compute_p_elems()
			self.eq_mom.compute_q_elems()
			self.eq_mom.compute_p_nodes()
			self.eq_mom.compute_q_nodes()
			for output in self.outputs:
				output.save_fields(t)

			# Print stuff
			if self.eq_mom.grid.mesh.comm.rank == 0:
				print(t/self.t_control.time_unit, ite, error)
				try:
					print(float(self.eq_mom.mat.elems_ne[-1].Fvp.max()))
				except:
					pass

		for output in self.outputs:
			output.save_mesh()


class Simulator_T(Simulator):
	def __init__(self, eq_heat, t_control, outputs, compute_elastic_response=True):
		self.eq_heat = eq_heat
		self.t_control = t_control
		self.outputs = outputs

	def run(self):
		# Output field
		for output in self.outputs:
			output.initialize()

		# Solve initial T field
		self.eq_heat.solve(0, self.t_control.dt)

		# Save fields
		output.save_fields(0)

		# Time loop
		while self.t_control.keep_looping():

			# Advance time
			self.t_control.advance_time()
			t = self.t_control.t
			dt = self.t_control.dt

			# Update boundary conditions
			self.eq_heat.bc.update_dirichlet(t)
			self.eq_heat.bc.update_neumann(t)

			# Solve heat
			self.eq_heat.solve(t, dt)

			# Print stuff
			if self.eq_heat.grid.mesh.comm.rank == 0:
				print(t/self.t_control.time_unit)

			# Save fields
			for output in self.outputs:
				output.save_fields(t)

		for output in self.outputs:
			output.save_mesh()