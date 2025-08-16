"""
Class responsible for printing data on screen during simulation
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

import time
import os
import sys
from mpi4py import MPI
from .MaterialProps import Material
from .OutputHandler import SaveFields
from .Grid import GridHandlerGMSH
from petsc4py import PETSc

def singleton(cls):
	instances = {}
	def get_instance(*args, **kwargs):
		if cls not in instances:
			instances[cls] = cls(*args, **kwargs)
		return instances[cls]
	def reset_instance():
		instances.pop(cls, None)
	get_instance.reset_instance = reset_instance
	return get_instance

@singleton
class ScreenPrinter():
	def __init__(self, grid: GridHandlerGMSH, solver: PETSc.KSP, material: Material, outputs: list[SaveFields], time_unit: str="hour"):
		self.master_division_plus = "+-----------------------------------------------------------------------------------------------+"
		self.master_division = "-----------------------------------------------------------------------------------------------"
		self.log = ""
		self.grid = grid
		self.solver = solver
		self.mat = material
		self.outputs = outputs
		self.time_unit = time_unit

		self.set_welcome()
		self.print_welcome()
		self.print_mesh_info()
		self.print_partition_info()
		self.print_solver_info()
		self.print_constitutive_model()
		self.print_output_info()
		self.begin()

	def begin(self):
		self.start_timer()
		self.set_header_columns([	
									"Step counter",
									f"dt ({self.time_unit})",
									f"t / t_final ({self.time_unit})",
									"# of iters",
									"Non-linear error"
								], "center")
		self.set_row_formats([
								"%i",
								"%.3f",
								"%s",
								"%.i",
								"%.4e",
							], ["center" for i in range(5)])
		self.print_header()

	def print_solver_info(self):
		ksp_type = self.solver.getType()
		pc_type = self.solver.getPC().getType()
		rtol, atol, divtol, max_it = self.solver.getTolerances()
		self.print_comment(" Solver info:")
		self.set_header_columns(["KSP_type", "PC_type", "  rtol  ", "max_it"], "center")
		self.print_header()
		self.set_row_formats(["%s", "%s", "%.1e", "%.i"], ["center" for i in range(4)])
		self.print_row([ksp_type, pc_type, rtol, max_it])
		self.print_on_screen(self.divider)
		self.print_comment(" ")


	def print_mesh_info(self):
		size = len(self.grid.grid_folder) - len("Location")
		self.print_comment(" Mesh info:")
		self.set_header_columns(["# of elements", "# of nodes", "Location"+size*" "], "left")
		self.print_header()
		self.set_row_formats(["%.i", "%.i", "%s"], ["left", "left", "left"])
		n_elems = MPI.COMM_WORLD.allreduce(self.grid.mesh.topology.index_map(self.grid.domain_dim).size_local, op=MPI.SUM)
		n_nodes = MPI.COMM_WORLD.allreduce(self.grid.mesh.topology.index_map(0).size_local, op=MPI.SUM)
		self.print_row([n_elems, n_nodes, self.grid.grid_folder])
		self.print_on_screen(self.divider)
		self.print_comment(" ")

	def print_partition_info(self):
		size = len(self.grid.grid_folder) - len("Location")
		self.print_comment(" Partition(s) info:")
		self.set_header_columns(["Partition #", "# of elements", "# of nodes"], "left")
		self.print_header()
		self.set_row_formats(["%.i", "%.i", "%.i"], ["center", "center", "center"])
		comm = MPI.COMM_WORLD
		for rank in range(1, comm.Get_size()):
			if comm.rank == rank:
				n_elems = self.grid.mesh.topology.index_map(self.grid.domain_dim).size_local
				n_nodes = self.grid.mesh.topology.index_map(0).size_local
				comm.send([n_elems, n_nodes], dest=0, tag=10*comm.rank)
		if comm.rank == 0:
			n_elems = self.grid.mesh.topology.index_map(self.grid.domain_dim).size_local
			n_nodes = self.grid.mesh.topology.index_map(0).size_local
			self.print_row([1, n_elems, n_nodes])
			for rank in range(1, comm.Get_size()):
				values = comm.recv(source=rank, tag=10*rank)
				self.print_row([rank+1, values[0], values[1]])
		self.print_on_screen(self.divider)
		self.print_comment(" ")

	def print_constitutive_model(self):
		elems_e_list = ""
		for elem_e in self.mat.elems_e:
			elems_e_list += elem_e.name if len(elems_e_list) == 0 else ", " + elem_e.name
		elems_ne_list = ""
		for elem_ne in self.mat.elems_ne:
			elems_ne_list += elem_ne.name if len(elems_ne_list) == 0 else ", " + elem_ne.name
		elems_th_list = ""
		for elem_th in self.mat.elems_th:
			elems_th_list += elem_th.name if len(elems_th_list) == 0 else ", " + elem_th.name

		n_elems = len(self.mat.elems_e) + len(self.mat.elems_ne) + len(self.mat.elems_th)
		if n_elems > 0:
			size = max([len(elems_e_list), len(elems_ne_list), len(elems_th_list)]) - len("List of elements")
			self.print_comment(" Constitutive model:")
			self.set_header_columns(["Element type ", "List of elements" + " "*size], "left")
			self.print_header()
			self.set_row_formats(["%s", "%s"], ["left", "left"])
			self.print_row(["elastic", elems_e_list])
			self.print_row(["non-elastic", elems_ne_list])
			self.print_row(["thermoelastic", elems_th_list])

			self.print_on_screen(self.divider)
			self.print_comment(" ")

	def print_output_info(self):
		self.print_comment(" Output info:")
		output_folder = self.outputs[0].output_folder
		size = len(output_folder)
		self.set_header_columns(["Location"+10*" ", "Field name      ", "Label name             "], "center")
		self.print_header()
		self.set_row_formats(["%s", "%s", "%s"], ["left", "left", "left"])
		self.output_folders = []
		for output in self.outputs:
			self.output_folders.append(output.output_folder)
			for field_data in output.fields_data:
				self.print_row([output.output_folder, field_data["field_name"], field_data["label_name"]])
		self.print_on_screen(self.divider)
		self.print_comment(" ")




	def set_welcome(self):
		# Generated at https://www.asciiart.eu/text-to-ascii-art with Standard font
		self.max_width = len(self.master_division_plus)
		self.welcome_text =  "+===============================================================================================+\n"
		self.welcome_text += "|   ____    _    _____ _____   ___ _   _    ____    ___     _______         ____    ___   ___   |\n"
		self.welcome_text += "|  / ___|  / \  |  ___| ____| |_ _| \ | |  / ___|  / \ \   / / ____| __   _|___ \  / _ \ / _ \  |\n"
		self.welcome_text += "|  \___ \ / _ \ | |_  |  _|    | ||  \| | | |     / _ \ \ / /|  _|   \ \ / / __) || | | | | | | |\n"
		self.welcome_text += "|   ___) / ___ \|  _| | |___   | || |\  | | |___ / ___ \ V / | |___   \ V / / __/ | |_| | |_| | |\n"
		self.welcome_text += "|  |____/_/   \_\_|   |_____| |___|_| \_|  \____/_/   \_\_/  |_____|   \_/ |_____(_)___(_)___/  |\n"
		self.welcome_text += "|                                                                                               |\n"
		self.welcome_text += "+===============================================================================================+"

	def set_row_formats(self, row_formats, row_align):
		self.row_formats = row_formats
		self.row_align = row_align

	def add_to_log(self, message):
		self.log += "\n" + message

	def print_welcome(self):
		self.print_on_screen(self.welcome_text)
		self.print_comment(" ")

	def start_timer(self):
		comm = MPI.COMM_WORLD
		comm.Barrier()
		if MPI.COMM_WORLD.rank == 0:
		    self.start = MPI.Wtime()

	def close(self):
		self.print_on_screen(self.divider)
		if MPI.COMM_WORLD.rank == 0:
			self.final = MPI.Wtime()
			cpu_time = self.final - self.start
			formatted_time = time.strftime("%H:%M:%S", time.gmtime(cpu_time))
			full_width = len(self.master_division_plus)
			message = self.format_cell(f"Total time: {formatted_time} ({cpu_time} seconds)", full_width, "right")
			self.print_on_screen(message)
			for output_folder in self.output_folders:
				self.save_log(output_folder)

	def save_log(self, output_folder):
		with open(os.path.join(output_folder, "log.txt"), 'w') as output:
			output.write(self.log)

	def print_comment(self, comment, align="left"):
		if comment != None:
			full_width = len(self.master_division_plus)
			message = "|" + self.format_cell(comment, full_width-2, align) + "|"
			self.print_on_screen(message)

	def print_on_screen(self, raw_comment):
		if MPI.COMM_WORLD.rank == 0:
			print(raw_comment)
			sys.stdout.flush()
			self.add_to_log(raw_comment)

	def set_header_columns(self, header_columns, align):
		self.header_columns = header_columns
		self.header_align = [align for i in range(len(header_columns))]
		self.widths = []
		for header_column in self.header_columns:
			self.widths.append(len(header_column))
		self.divider = self.make_divider(self.widths, "+")

	def print_header(self):
		"""
		Print the top divider, a header row (using alignments), and a divider beneath.
		"""
		# Print header row
		self.print_on_screen(self.divider)
		header_line = "| " + " | ".join(
		    self.format_cell(col, w, align)
		    for col, w, align in zip(self.header_columns, self.widths, self.header_align)
		) #+ " |"
		if self.max_width - len(header_line) - 1 > 1:
			header_line += " |"
			header_line += " " * (self.max_width - len(header_line) - 1)
			header_line += "|"
		else:
			header_line += " " * (self.max_width - len(header_line) - 1)
			header_line += "|"
		self.print_on_screen(header_line)
		self.print_on_screen(self.divider)

	def print_row(self, values):
		"""
		Print a single row of values, using the provided alignment rules.
		"""
		row_line = "| " + " | ".join(
			self.format_cell(val, w, align, text_format) 
			for val, w, align, text_format in zip(values, self.widths, self.row_align, self.row_formats)
		) #+ " |"
		if self.max_width - len(row_line) - 1 > 1:
			row_line += " |"
			row_line += " " * (self.max_width - len(row_line) - 1)
			row_line += "|"
		else:
			row_line += " " * (self.max_width - len(row_line) - 1)
			row_line += "|"
		self.print_on_screen(row_line)

	def make_divider(self, widths, middle="+"):
		"""
		Return a string for the horizontal divider line.
		Example: +--------+--------+
		"""
		segments = [ "-" * (w + 2) for w in widths ]
		divider = "+" + middle.join(segments) + "+"
		if self.max_width - len(divider) - 1 > -1:
			divider += "-" * (self.max_width - len(divider) - 1)
			divider += "+"
		else:
			divider += "-" * (self.max_width - len(divider) - 2)
		return divider


	def format_cell(self, text, width, alignment, text_format=None):
		"""
		Format a cell according to the width and alignment.
		alignment can be 'left' or 'center' (extendable to 'right').
		"""
		if alignment == 'left':
		    if text_format != None: 
		        text = text_format%text
		    return f"{text:<{width}}"
		elif alignment == 'center':
		    if text_format != None: 
		        text = text_format%text
		    return f"{text:^{width}}"
		else:
		    # Fallback (could be 'right' or any other future alignment)
		    if text_format != None: 
		        text = text_format%text
		    return f"{text:>{width}}"
