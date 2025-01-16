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

from Equations import LinearMomentum
from Grid import GridHandlerGMSH
import Utils as utils
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

		# # Save input file
		# filename = os.path.join(os.path.join(input_file["output"]["path"], "input_file.json"))
		# os.makedirs(os.path.dirname(filename), exist_ok=True)
		# utils.save_json(self.input_file_to_be_saved, filename)

	def __save_input_file(self, filename):
		os.makedirs(os.path.dirname(filename), exist_ok=True)
		utils.save_json(self.input_file_to_be_saved, filename)


	def run(self):
		"""
		Runs transient simulation.
		"""
		solve_equilibrium = self.input_file["simulation_settings"]["equilibrium"]["active"]

		# Save input file
		self.__save_input_file(os.path.join(self.output_folder, "operation", "input_file.json"))
		if solve_equilibrium:
			self.__save_input_file(os.path.join(self.output_folder, "equilibrium", "input_file.json"))

		# Run simulation
		self.run_simulation(solve_equilibrium=solve_equilibrium, verbose=True)


	def run_simulation(self, solve_equilibrium=False, verbose=True):
		"""
		Runs simulation **without** solving the equilibrium condition.

		Parameters
		----------
		solve_equilibrium : bool
			If **True**, it calculates the equilibrium conditions considering the elastic and viscoelastic (if present) elements only. 
			If **False**, then it skips the equilibrium condition.

		verbose : bool
			Shows real time simulation info on screen.
		"""
		# Pseudo time
		t = self.time_list[0]
		t_final = self.time_list[-1]

		# Get maximum time step size
		dt = self.input_file["simulation_settings"]["operation"]["dt_max"]

		# Read number of time steps to skip before saving results
		n_skip = self.input_file["simulation_settings"]["operation"]["n_skip"]

		# Perform initial computations
		self.eq_mom.initialize(solve_equilibrium=solve_equilibrium, verbose=verbose, save_results=True)

		# Save initial solution
		self.eq_mom.save_solution(t)

		# Transient simulation
		n_step = 1
		while t < t_final:

			# Increment time
			t += dt

			# Solve
			self.eq_mom.solve(t, dt)

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


