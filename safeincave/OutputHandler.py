import dolfinx as do
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


