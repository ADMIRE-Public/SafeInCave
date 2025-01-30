"""
Useful to read vtk files and post-process results.
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

import time
import os

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
	def __init__(self):
		self.master_division_plus = "+-----------------------------------------------------------------------------------------------+"
		self.master_division = "-----------------------------------------------------------------------------------------------"

		self.log = ""

		self.set_welcome()

	def set_welcome(self):
		self.max_width = len(self.master_division_plus)
		self.welcome_text =  "+===============================================================================================+\n"
		self.welcome_text += "|   ____    _    _____ _____   ___ _   _    ____    ___     _______         _   ____    ___     |\n"
		self.welcome_text += "|  / ___|  / \  |  ___| ____| |_ _| \ | |  / ___|  / \ \   / / ____| __   _/ | |___ \  / _ \    |\n"
		self.welcome_text += "|  \___ \ / _ \ | |_  |  _|    | ||  \| | | |     / _ \ \ / /|  _|   \ \ / / |   __) || | | |   |\n"
		self.welcome_text += "|   ___) / ___ \|  _| | |___   | || |\  | | |___ / ___ \ V / | |___   \ V /| |_ / __/ | |_| |   |\n"
		self.welcome_text += "|  |____/_/   \_\_|   |_____| |___|_| \_|  \____/_/   \_\_/  |_____|   \_/ |_(_)_____(_)___/    |\n"
		self.welcome_text += "|                                                                                               |\n"
		self.welcome_text += "+===============================================================================================+\n"
		self.welcome_text += "|                                                                                               |\n"
		self.welcome_text += "| Herminio Tasinafo Honorio (main developer)                                                    |\n"
		self.welcome_text += "| Hadi Hajibeygi (PI)                                                                           |\n"
		self.welcome_text += "|                                                                                               |\n"
		self.welcome_text += "+===============================================================================================+"

	def set_header_columns(self, header_columns, align):
		self.header_columns = header_columns
		self.header_align = [align for i in range(len(header_columns))]
		self.widths = []
		for header_column in self.header_columns:
			self.widths.append(len(header_column))
		self.divider = self.make_divider(self.widths, "+")

	def set_row_formats(self, row_formats, row_align):
		self.row_formats = row_formats
		self.row_align = row_align

	def start_timer(self):
		self.start = time.time()

	def add_to_log(self, message):
		self.log += "\n" + message

	def print_welcome(self):
		self.print_on_screen(self.welcome_text)

	def close(self):
		self.print_on_screen(self.divider)
		self.final = time.time()
		cpu_time = self.final - self.start
		formatted_time = time.strftime("%H:%M:%S", time.gmtime(cpu_time))
		full_width = len(self.master_division_plus)
		message = self.format_cell(f"Total time: {formatted_time} ({cpu_time} seconds)", full_width, "right")
		self.print_on_screen(message)

	def save_log(self, output_folder):
		with open(os.path.join(output_folder, "log.txt"), 'w') as output:
			output.write(self.log)

	def print_comment(self, comment, align="left"):
		if comment != None:
			full_width = len(self.master_division_plus)
			message = "|" + self.format_cell(comment, full_width-2, align) + "|"
			self.print_on_screen(message)

	def print_on_screen(self, raw_comment):
		print(raw_comment)
		self.add_to_log(raw_comment)

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


if __name__ == '__main__':
	screen = ScreenPrinter(
			header_columns = ["Time step", "Final time (h)", "Current time (h)", "# of iters", "Non-linear error", "Save solution"],
			header_align = "center",
			row_formats = ["%i", "%.3f", "%.3f", "%i", "%.4e", "%s"],
			row_align = ["center", "center", "center", "center", "center", "center"],
			comment = "Running equilibrium stage"
	)

	screen.print_header()
	screen.print_row([1,2,3,4,5,"Yes"])
	screen.print_row([1,2,3,4,5,"Yes"])
	screen.print_row([1,2,3,4,5,"Yes"])
	screen.print_row([1,2,3,4,5,"Yes"])
	screen.close()