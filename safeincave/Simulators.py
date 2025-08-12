"""
Implementation of the simulators
"""
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
from abc import ABC, abstractmethod
import torch as to
import numpy as np
from mpi4py import MPI
from Utils import numpy2torch
from HeatEquation import HeatDiffusion
from MomentumEquation import LinearMomentum
from TimeHandler import TimeControllerBase
from OutputHandler import SaveFields

class Simulator(ABC):
	"""
    Abstract simulation driver interface.

    Subclasses implement a concrete `run()` method that advances one or more
    coupled PDE solvers in time, handles I/O, and updates material/internal
    variables as needed.
    """
	@abstractmethod
	def run(self):
		"""
        Execute the simulation.

        Returns
        -------
        None
        """
		pass


class Simulator_TM(Simulator):
	"""
	Run the coupled thermo–mechanical simulation.

	Workflow
	--------
	1. Initialize outputs.
	2. Initialize momentum temperature from the heat solution and update BCs.
	3. Optionally solve a purely elastic response.
	4. Initialize non-elastic rates.
	5. For each time step:

	   - Advance time and update boundary conditions for both equations.
	   - Solve the heat equation for ``(t, dt)`` and set temperatures in momentum.
	   - Iterate the momentum step (assemble/solve, update internal variables and rates).
	   - Save requested fields.

	Returns
	-------
	None
	"""
	def __init__(self, eq_mom: LinearMomentum, 
					   eq_heat: HeatDiffusion, 
					   t_control: TimeControllerBase, 
					   outputs: list[SaveFields],
					   compute_elastic_response: bool=True):
		self.eq_mom = eq_mom
		self.eq_heat = eq_heat
		self.t_control = t_control
		self.outputs = outputs
		self.compute_elastic_response = compute_elastic_response

	def run(self) -> None:
		"""
		Run the coupled thermo–mechanical simulation.

		Workflow
		--------
		1. Initialize outputs.
		2. Initialize momentum temperature history from the heat solution.
		3. Update BCs and optionally solve a purely elastic step.
		4. Initialize non-elastic rates.
		5. Time loop:

		   - Advance time, update BCs, solve heat step.
		   - Fixed-point (or single-pass) iterate the momentum step with the
		     current temperature, updating internal variables and rates.
		   - Save requested fields.

		Convergence
		-----------
		Uses a relative change in total strain between iterations as error.
		If ``theta == 1.0`` (backward Euler) or there are no non-elastic
		elements, iteration terminates immediately.

		Returns
		-------
		None

		Notes
		-----
		- Calls ``SaveFields.initialize()`` once and ``save_fields(t)`` at each
		  saved time, followed by ``save_mesh()`` after the loop.
		- Printing of progress occurs on rank 0 only.
		- The first ``output.save_fields(0)`` call targets the last ``output``
		  from the preceding loop variable; ensure all outputs are saved by
		  iterating over ``self.outputs``.
		"""
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
				elif len(self.eq_mom.mat.elems_ne) == 0:
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
	"""
    Mechanical-only simulator (linear momentum).

    Solves the momentum equation with possible non-elastic behavior using a
    θ-method loop and fixed-point iterations per step. No thermal coupling.

    Parameters
    ----------
    eq_mom : LinearMomentum
        Configured momentum equation (materials, BCs, solver set).
    t_control : TimeControllerBase
        Time controller providing `t`, `dt`, and loop control.
    outputs : list of SaveFields
        Output writers to initialize and use at each saved time.
    compute_elastic_response : bool, default=True
        If True, starts with a purely elastic solve to initialize fields.

    Attributes
    ----------
    eq_mom : LinearMomentum
    t_control : TimeControllerBase
    outputs : list[SaveFields]
    compute_elastic_response : bool
    """
	def __init__(self, eq_mom: LinearMomentum, 
					   t_control: TimeControllerBase,
					   outputs: list[SaveFields],
					   compute_elastic_response: bool=True):
		self.eq_mom = eq_mom
		self.t_control = t_control
		self.outputs = outputs
		self.compute_elastic_response = compute_elastic_response

	def run(self) -> None:
		"""
        Run the mechanical simulation.

        Workflow
        --------
        1. Initialize outputs and boundary conditions.
        2. Optionally solve a purely elastic step.
        3. Initialize non-elastic rates.
        4. For each time step: assemble/solve, update internal variables and
           rates, compute relevant quantities, and save fields.

        Convergence
        -----------
        Uses a relative change in total strain between iterations as error.
        If `theta == 1.0` or no non-elastic elements exist, iteration ends
        immediately.

        Returns
        -------
        None

        Notes
        -----
        - Printing occurs on rank 0 only.
        - The first `output.save_fields(0)` call uses the last `output`
          from the preceding loop variable.
        """
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
				elif len(self.eq_mom.mat.elems_ne) == 0:
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
	"""
    Thermal-only simulator (heat diffusion).

    Advances the heat equation with fully-implicit time loop and writes fields.

    Parameters
    ----------
    eq_heat : HeatDiffusion
        Configured heat equation (materials, BCs, solver set).
    t_control : TimeControllerBase
        Time controller providing `t`, `dt`, and loop control.
    outputs : list of SaveFields
        Output writers to initialize and use at each saved time.
    compute_elastic_response : bool, default=True
        Unused placeholder kept for interface parity.

    Attributes
    ----------
    eq_heat : HeatDiffusion
    t_control : TimeControllerBase
    outputs : list[SaveFields]
    """
	def __init__(self, eq_heat: HeatDiffusion,
					   t_control: TimeControllerBase,
					   outputs: list[SaveFields],
					   compute_elastic_response: bool=True):
		self.eq_heat = eq_heat
		self.t_control = t_control
		self.outputs = outputs

	def run(self) -> None:
		"""
        Run the thermal simulation.

        Workflow
        --------
        1. Initialize outputs.
        2. (Optionally) solve an initial step.
        3. Time loop: update BCs, solve heat equation for `(t, dt)`, and save.

        Returns
        -------
        None

        Notes
        -----
        Printing of progress occurs on rank 0 only.
        """
		# Output field
		for output in self.outputs:
			output.initialize()

		# # Solve initial T field
		# self.eq_heat.solve(0, self.t_control.dt)

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