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
from Gridx import GridHandlerGMSH
import Utils as utils
from ScreenOutput import ScreenPrinter
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
		self.input_file = copy.deepcopy(input_file)

		# This is input_file to be saved
		self.input_file_to_be_saved = copy.deepcopy(input_file)

		# Output folder
		self.output_folder = input_file["output"]["path"]

		# Transient settings
		self.time_list = input_file["time_settings"]["time_list"]
		theta = input_file["time_settings"]["theta"]
	
		# Create mesh
		grid_path = input_file["grid"]["path"]
		self.grid = GridHandlerGMSH(input_file["grid"]["name"], grid_path)

		# Linear momentum balance equation
		self.eq_mom = LinearMomentum(self.grid, theta, self.input_file)

		# # Save input file
		# filename = os.path.join(os.path.join(input_file["output"]["path"], "input_file.json"))
		# os.makedirs(os.path.dirname(filename), exist_ok=True)
		# utils.save_json(self.input_file_to_be_saved, filename)

		# Screen info
		ScreenPrinter.reset_instance()
		self.screen = ScreenPrinter()
		self.screen.print_welcome()
		self.screen.print_comment(" ")
		self.screen.print_comment(" Finite element (FE) simulator.")

		self.screen.print_comment(" ")
		self.screen.print_comment(" Results folder:")
		self.screen.print_comment(f"          {self.output_folder}")

		self.screen.print_comment(" ")
		self.screen.print_comment(" Mesh info:")
		self.screen.print_comment(f"          Location: {grid_path}")
		self.screen.print_comment(f"          Number of elements: {self.grid.n_elems}")
		self.screen.print_comment(f"          Number of nodes: {self.grid.n_nodes}")

		solver_type = input_file["solver_settings"]["solver_type"]
		solver_method = input_file["solver_settings"]["solver_PC"]
		rtol = input_file["solver_settings"]["rtol"]

		self.screen.print_comment(" ")
		self.screen.print_comment(" Solver info:")
		self.screen.print_comment(f"          Type: {solver_type}")
		self.screen.print_comment(f"          PC: {solver_method}")
		self.screen.print_comment(f"          Tol: {rtol}")
		# if solver_type == "KrylovSolver":
		# 	solver_prec = input_file["solver_settings"]["preconditioner"]
		# 	solver_tol = input_file["solver_settings"]["relative_tolerance"]
		# 	self.screen.print_comment(f"          Preconditioner: {solver_prec}")
		# 	self.screen.print_comment(f"          Tolerance: {solver_tol}")
		self.screen.print_comment(" ")

	def __save_input_file(self, filename):
		os.makedirs(os.path.dirname(filename), exist_ok=True)
		utils.save_json(self.input_file_to_be_saved, filename)


	def run(self):
		"""
		Runs transient simulation.
		"""
		solve_equilibrium = self.input_file["simulation_settings"]["equilibrium"]["active"]
		solve_operation = self.input_file["simulation_settings"]["operation"]["active"]

		# Save input file
		self.__save_input_file(os.path.join(self.output_folder, "operation", "input_file.json"))
		if solve_equilibrium:
			self.__save_input_file(os.path.join(self.output_folder, "equilibrium", "input_file.json"))

		# Run simulation
		self.run_simulation(solve_operation, solve_equilibrium, verbose=True)


	def run_simulation(self, solve_operation=True, solve_equilibrium=False, verbose=True):
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

		# Shoud calculate initial hardening such that Fvp=0 everywhere?
		hardening = self.input_file["simulation_settings"]["operation"]["hardening"]

		# Perform initial computations
		self.eq_mom.initialize(solve_equilibrium=solve_equilibrium, verbose=verbose, save_results=True, calculate_hardening=hardening)

		print("check")
		
		if solve_operation:

			# header_columns = ["Time step", "Final time (h)", "Current time (h)", "# of iters", "Non-linear error", "Save solution"],
			# header_align = "center",
			# row_formats = ["%s", "%.3f", "%.3f", "%.i", "%.4e", "%s"],
			# row_align = ["center", "center", "center", "center", "center", "center"],

			self.screen.start_timer()
			self.screen.set_header_columns(["Time step", "Final time (h)", "Current time (h)", "# of iters", "Non-linear error", "Save solution"], "center")
			self.screen.set_row_formats(["%s", "%.3f", "%.3f", "%.i", "%.4e", "%s"], ["center" for i in range(6)])

			elem_names = []
			for elem_type in self.input_file["constitutive_model"].keys():
				for elem_name in self.input_file["constitutive_model"][elem_type].keys():
					if self.input_file["constitutive_model"][elem_type][elem_name]["active"] == True:
						elem_names.append(elem_name)

			self.screen.print_on_screen(" ")
			self.screen.print_on_screen(self.screen.master_division_plus)
			self.screen.print_comment(" Operation Stage", "center")
			self.screen.print_on_screen(self.screen.master_division_plus)
			self.screen.print_comment(" ")
			self.screen.print_comment(" Constitutive model:")
			for elem_name in elem_names:
				self.screen.print_comment(f"          {elem_name}")
			self.screen.print_comment(" ")
			self.screen.print_header()

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
				save_solution = False
				if n_step % n_skip == 0 or n_step == 1:
					self.eq_mom.save_solution(t)
					save_solution = True

				# Print stuff
				if verbose:
					if save_solution:
						screen_output_row = [str(n_step), t_final/utils.hour, t/utils.hour, self.eq_mom.ite, self.eq_mom.error, "Save"]
					else:
						screen_output_row = [str(n_step), t_final/utils.hour, t/utils.hour, self.eq_mom.ite, self.eq_mom.error, "|"]
					self.screen.print_row(screen_output_row)

				n_step += 1

			self.screen.close()
			self.screen.save_log(self.output_folder)


