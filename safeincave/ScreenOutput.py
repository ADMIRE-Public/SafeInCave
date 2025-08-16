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
	"""
	Pretty, MPI-aware console logger for simulation metadata and progress.

	Prints a formatted header (mesh, partitions, solver settings, material
	model, and outputs) and a live table of time-stepping information.
	Only **rank 0** writes to stdout; all printed lines are also appended to
	an internal log buffer and persisted to ``log.txt`` in each output folder
	on :meth:`close`.

	Parameters
	----------
	grid : GridHandlerGMSH
	    Mesh/grid handler used to query number of nodes/elements and partition
	    information.
	solver : petsc4py.PETSc.KSP
	    Configured linear solver; its type and tolerances are reported.
	material : Material
	    Material container; its elastic, non-elastic, and thermoelastic
	    elements are listed.
	outputs : list[SaveFields]
	    Output handlers to extract output directories and registered field
	    names/labels.
	time_unit : {"second", "minute", "hour", "day", "year"}, default="hour"
	    Unit label used in the progress table headings (purely cosmetic).

	Attributes
	----------
	log : str
	    Accumulated text printed to screen.
	header_columns : list[str]
	    Column headers for the live progress table.
	row_formats : list[str]
	    C-printf format strings for each column value (e.g., ``"%.3f"``).
	row_align : list[str]
	    Alignment for each column: ``"left"`` or ``"center"`` (``"right"``
	    is accepted by :meth:`format_cell`).
	max_width : int
	    Width used to pad and cap the table rows (derived from ASCII banner).
	output_folders : list[str]
	    Paths collected from the provided :class:`SaveFields` objects.

	Notes
	-----
	- Construction performs side effects: prints welcome screen, mesh/partition
	  info, solver and material summaries, output info, and begins the table.
	- Use :meth:`close` when the run ends to finalize timing and write logs.
	"""
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

	def begin(self) -> None:
		"""
		Initialize the progress table and start timing.

		Sets default column headers and formats for the step counter, time
		increment, normalized time, nonlinear iterations, and residual error,
		then prints the table header.

		Returns
		-------
		None
		"""
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

	def print_solver_info(self) -> None:
		"""
		Print PETSc KSP/PC configuration and tolerances.

		Extracts KSP type, PC type, relative tolerance, and maximum
		iterations, then renders a small table.

		Returns
		-------
		None
		"""
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


	def print_mesh_info(self) -> None:
		"""
		Print global mesh size (elements, nodes) and location on disk.

		Gathers counts across MPI ranks and displays a one-row table with the
		mesh path.

		Returns
		-------
		None
		"""
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

	def print_partition_info(self) -> None:
		"""
		Print per-partition (MPI rank) local element/node counts.

		Rank 0 collects local sizes from all ranks and renders a table
		listing each partition with its local element and node counts.

		Returns
		-------
		None
		"""
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

	def print_constitutive_model(self) -> None:
		"""
		Print a summary of elastic, non-elastic, and thermoelastic elements.

		Aggregates element names from the material container and displays
		them grouped by type.

		Returns
		-------
		None
		"""
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

	def print_output_info(self) -> None:
		"""
		Print output destinations and the registered output fields.

		Shows, for each :class:`SaveFields` object, its output directory and
		the field name/label pairs that will be written.

		Returns
		-------
		None
		"""
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




	def set_welcome(self) -> None:
		"""
		Build the ASCII-art banner and capture width for table formatting.

		The banner determines :attr:`max_width`, which is used to pad header
		and row strings uniformly.

		Returns
		-------
		None
		"""
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

	def set_row_formats(self, row_formats: list[str], row_align: list[str]) -> None:
		"""
		Set value format strings and alignment for subsequent table rows.

		Parameters
		----------
		row_formats : list[str]
		    Format strings (printf-like) for each column value.
		row_align : list[str]
		    Alignments per column (``"left"``, ``"center"``, or ``"right"``).

		Returns
		-------
		None
		"""
		self.row_formats = row_formats
		self.row_align = row_align

	def add_to_log(self, message: str) -> None:
		"""
		Append a line to the internal log buffer.

		Parameters
		----------
		message : str
		    Line to append. A newline is prepended automatically.

		Returns
		-------
		None
		"""
		self.log += "\n" + message

	def print_welcome(self) -> None:
		"""
		Print the welcome banner and a blank spacer line.

		Returns
		-------
		None
		"""
		self.print_on_screen(self.welcome_text)
		self.print_comment(" ")

	def start_timer(self) -> None:
		"""
		Synchronize MPI ranks and start the wall-clock timer.

		On rank 0, stores :attr:`start` using :func:`MPI.Wtime`.

		Returns
		-------
		None
		"""
		comm = MPI.COMM_WORLD
		comm.Barrier()
		if MPI.COMM_WORLD.rank == 0:
		    self.start = MPI.Wtime()

	def close(self) -> None:
		"""
		Finalize logging: print total runtime and save logs to disk.

		Prints a closing divider and, on rank 0, computes wall time since
		:meth:`start_timer`, formats it as ``HH:MM:SS`` and seconds, and
		writes the accumulated :attr:`log` to ``log.txt`` in each output
		folder recorded during :meth:`print_output_info`.

		Returns
		-------
		None
		"""
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

	def save_log(self, output_folder: str) -> None:
		"""
		Write the accumulated log to ``log.txt`` inside ``output_folder``.

		Parameters
		----------
		output_folder : str
		    Destination directory. Created beforehand elsewhere.

		Returns
		-------
		None
		"""
		with open(os.path.join(output_folder, "log.txt"), 'w') as output:
			output.write(self.log)

	def print_comment(self, comment: str, align: str="left") -> None:
		"""
		Print a single comment line boxed to the full table width.

		Parameters
		----------
		comment : str or None
		    Text to print. If ``None``, nothing is printed.
		align : {"left", "center", "right"}, default="left"
		    Horizontal alignment of the text within the box.

		Returns
		-------
		None
		"""
		if comment != None:
			full_width = len(self.master_division_plus)
			message = "|" + self.format_cell(comment, full_width-2, align) + "|"
			self.print_on_screen(message)

	def print_on_screen(self, raw_comment: str) -> None:
		"""
		Print a raw string to stdout on rank 0 and add it to the log.

		Parameters
		----------
		raw_comment : str
		    A full line to print (including any box characters).

		Returns
		-------
		None
		"""
		if MPI.COMM_WORLD.rank == 0:
			print(raw_comment)
			sys.stdout.flush()
			self.add_to_log(raw_comment)

	def set_header_columns(self, header_columns: str, align: str) -> None:
		"""
		Configure table header columns and compute the divider line.

		Parameters
		----------
		header_columns : list[str]
		    Column titles to display in the header row.
		align : {"left", "center", "right"}
		    Alignment to apply to **all** header columns.

		Returns
		-------
		None
		"""
		self.header_columns = header_columns
		self.header_align = [align for i in range(len(header_columns))]
		self.widths = []
		for header_column in self.header_columns:
			self.widths.append(len(header_column))
		self.divider = self.make_divider(self.widths, "+")

	def print_header(self) -> None:
		"""
		Print the top divider, a formatted header row, and a divider below.

		Returns
		-------
		None
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

	def print_row(self, values: str) -> None:
		"""
		Print a single row of values under the current header.

		Values are formatted with :attr:`row_formats` and aligned by
		:attr:`row_align`.

		Parameters
		----------
		values : list[object]
		    One value per column.

		Returns
		-------
		None
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

	def make_divider(self, widths: list[int], middle: str="+") -> None:
		"""
		Build a horizontal divider string matching the current column widths.

		Example
		-------
		``+--------+--------+``

		Parameters
		----------
		widths : list[int]
		    Content widths for each column (without padding).
		middle : str, default="+"
		    String used between column segments.

		Returns
		-------
		str
		    Divider line capped/padded to :attr:`max_width`.
		"""
		segments = [ "-" * (w + 2) for w in widths ]
		divider = "+" + middle.join(segments) + "+"
		if self.max_width - len(divider) - 1 > -1:
			divider += "-" * (self.max_width - len(divider) - 1)
			divider += "+"
		else:
			divider += "-" * (self.max_width - len(divider) - 2)
		return divider


	def format_cell(self, text, width: int, alignment: str, text_format: str=None) -> None:
		"""
		Format a cell string with width, alignment, and optional printf-format.

		Parameters
		----------
		text : object
		    Value to render. Converted to string (or via ``text_format``).
		width : int
		    Visible width reserved for the cell (without border chars).
		alignment : {"left", "center", "right"}
		    Horizontal alignment.
		text_format : str or None, default=None
		    If provided, a printf-like format string applied to ``text``
		    before alignment (e.g., ``\"%.3f\"``).

		Returns
		-------
		str
		    Formatted cell string of length at most ``width``.
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
