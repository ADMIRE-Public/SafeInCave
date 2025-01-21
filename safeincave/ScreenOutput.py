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

class ScreenPrinter():
	def __init__(self, header_columns, header_align, row_formats, row_align, comment=None):
		self.start = time.time()
		self.header_columns = header_columns
		self.header_align = [header_align for i in range(len(header_columns))]
		self.row_formats = row_formats
		self.row_align = row_align
		self.comment = comment
		self.master_division = "+-----------------------------------------------------------------------------------------------+"

		self.set_welcome()

		self.widths = []
		for header_column in self.header_columns:
			self.widths.append(len(header_column))

		self.divider = self.make_divider(self.widths, "+")


	def set_welcome(self):
		self.max_width = len(self.master_division)
		self.welcome_text =  "+-----------------------------------------------------------------------------------------------+\n"
		self.welcome_text += "|   ____    _    _____ _____   ___ _   _    ____    ___     _______         _   ____    ___     |\n"
		self.welcome_text += "|  / ___|  / \  |  ___| ____| |_ _| \ | |  / ___|  / \ \   / / ____| __   _/ | |___ \  / _ \    |\n"
		self.welcome_text += "|  \___ \ / _ \ | |_  |  _|    | ||  \| | | |     / _ \ \ / /|  _|   \ \ / / |   __) || | | |   |\n"
		self.welcome_text += "|   ___) / ___ \|  _| | |___   | || |\  | | |___ / ___ \ V / | |___   \ V /| |_ / __/ | |_| |   |\n"
		self.welcome_text += "|  |____/_/   \_\_|   |_____| |___|_| \_|  \____/_/   \_\_/  |_____|   \_/ |_(_)_____(_)___/    |\n"
		self.welcome_text += "|                                                                                               |\n"
		self.welcome_text += "+-----------------------------------------------------------------------------------------------+\n"

	def print_welcome(self):
		print(self.welcome_text)


	def close(self):
		print(self.divider)
		self.final = time.time()
		cpu_time = self.final - self.start
		formatted_time = time.strftime("%H:%M:%S", time.gmtime(cpu_time))
		full_width = len(self.divider)
		# print(f"Total time: {formatted_time} ({cpu_time} seconds)")
		message = f"Total time: {formatted_time} ({cpu_time} seconds)"
		print("|" + self.format_cell(message, full_width-2, "right") + "|")
		print(self.master_division)
		print()


	def print_comment(self):
		if self.comment != None:
			full_width = len(self.divider)
			print(self.make_divider(self.widths, "-"))
			print("|" + self.format_cell(self.comment, full_width-2, "left") + "|")
			# print(self.comment)


	def print_header(self):
		"""
		Print the top divider, a header row (using alignments), and a divider beneath.
		"""

		self.print_comment()
		print(self.divider)

		# Print header row
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
		print(header_line)
		print(self.divider)


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
		print(row_line)

	def make_divider(self, widths, middle="+"):
		"""
		Return a string for the horizontal divider line.
		Example: +--------+--------+
		"""
		segments = [ "-" * (w + 2) for w in widths ]
		# segments = [ "-" * 95 ]
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