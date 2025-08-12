"""
Class responsible for saving results
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

import dolfinx as do
import shutil
import os

class SaveFields():
	def __init__(self, eq):
		self.eq = eq
		self.fields_data = []
		self.output_fields = []

	def set_output_folder(self, output_folder):
		self.output_folder = output_folder

	def add_output_field(self, field_name : str, label_name : str):
		data = {
			"field_name": field_name,
			"label_name": label_name,
		}
		self.fields_data.append(data)

	def initialize(self):
		for field_data in self.fields_data:
			field_name = field_data["field_name"]
			output_field = do.io.XDMFFile(self.eq.grid.mesh.comm, os.path.join(self.output_folder, field_name, f"{field_name}.xdmf"), "w")
			output_field.write_mesh(self.eq.grid.mesh)
			self.output_fields.append(output_field)

	def save_fields(self, t : float):
		for i, field_data in enumerate(self.fields_data):
			field = getattr(self.eq, field_data["field_name"])
			field.name = field_data["label_name"]
			self.output_fields[i].write_function(field, t)

	def save_mesh(self):
		mesh_origin_file = os.path.join(self.eq.grid.grid_folder, f"{self.eq.grid.geometry_name}.msh")
		mesh_destination_folder = os.path.join(self.output_folder, "mesh")
		if not os.path.exists(mesh_destination_folder):
			os.makedirs(mesh_destination_folder, exist_ok=True)
		shutil.copy(mesh_origin_file, mesh_destination_folder)




