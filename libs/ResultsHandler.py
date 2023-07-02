import meshio
import os
import numpy as np
import pandas as pd
from abc import ABC, abstractmethod

# class Saver():
# 	@abstractmethod
# 	def record(self):
# 		pass

# 	@abstractmethod
# 	def save(self):
# 		pass

# class AverageSaver(Saver):
# 	def __init__(self, dx, field_name, field, time_handler, output_folder):
# 		self.saver = TensorSaver(field_name, dx)
# 		self.field = field
# 		self.time_handler = time_handler
# 		self.output_folder = output_folder

# 	def record(self):
# 		self.saver.record_average(self.field, self.time_handler.time)

# 	def save(self):
# 		self.saver.save(os.path.join(self.output_folder, "avg"))

# class VtkSaver(Saver):
# 	def __init__(self, field_name, field, time_handler, output_folder):
# 		from fenics import File
# 		self.field = field
# 		self.time_handler = time_handler
# 		self.output_folder = output_folder
# 		self.vtk = File(os.path.join(output_folder, f"{field_name}.pvd"))
# 		# self.vtk = File(os.path.join(output_folder, "vtk", f"{field_name}.pvd"))

# 	def record(self):
# 		self.vtk << (self.field, self.time_handler.time)

# 	def save(self):
# 		pass

class TensorSaver():
	def __init__(self, name, dx):
		import fenics as fe
		self.fe = fe
		self.name = name
		self.dx = dx
		self.vol = fe.assemble(1*self.dx)
		self.tensor_data = {
			"Time": [],
			"00": [],
			"01": [],
			"02": [],
			"10": [],
			"11": [],
			"12": [],
			"20": [],
			"21": [],
			"22": []
		}

	def get_average(self, sigma_ij):
		return self.fe.assemble(sigma_ij*self.dx)/self.vol

	def record_average(self, tensor, t):
		self.tensor_data["Time"].append(t)
		for i in range(3):
			for j in range(3):
				sigma_ij = tensor[i,j]
				avg_value = self.get_average(sigma_ij)
				self.tensor_data[f"{i}{j}"].append(avg_value)

	def save(self, output_folder):
		if not os.path.exists(output_folder):
			os.makedirs(output_folder)
		self.df = pd.DataFrame(self.tensor_data)
		self.df.to_excel(os.path.join(output_folder, f"{self.name}.xlsx"))
		# df.to_csv(os.path.join(output_folder, f"{self.name}.csv"))

class ScalarSaver():
	def __init__(self, name, dx):
		import fenics as fe
		self.fe = fe
		self.name = name
		self.dx = dx
		self.vol = fe.assemble(1*self.dx)
		self.scalar_data = {
			"Time": [],
			"Scalar": []
		}

	def get_average(self, scalar):
		return self.fe.assemble(scalar*self.dx)/self.vol

	def record_average(self, scalar, t):
		self.scalar_data["Time"].append(t)
		avg_value = self.get_average(scalar)
		self.scalar_data["Scalar"].append(avg_value)

	def save(self, output_folder):
		if not os.path.exists(output_folder):
			os.makedirs(output_folder)
		self.df = pd.DataFrame(self.scalar_data)
		self.df.to_excel(os.path.join(output_folder, f"{self.name}.xlsx"))
		# df.to_csv(os.path.join(output_folder, f"{self.name}.csv"))

class ResultsHandler(object):
	def __init__(self, output_folder, vtu_file_name):
		self.output_folder = output_folder
		self.vtu_file_name = vtu_file_name

		self.read_gmsh()
		self.build_tags()
		self.read_vtu()

	def read_gmsh(self):
		file_list = os.listdir()
		gmsh_files = [f for f in os.listdir(self.output_folder) if f.endswith(".msh")]
		if len(gmsh_files) > 1:
			raise Exception("Too many .msh files. Opening %s."%gmsh_files[0])
		self.gmsh_file = meshio.read(os.path.join(self.output_folder, gmsh_files[0]))

	def build_tags(self):
		self.gmsh_tags = {1:{}, 2:{}, 3:{}}
		for key, value in self.gmsh_file.field_data.items():
			self.gmsh_tags[value[1]][key] = value[0]

	def read_vtu(self):
		self.vtu_file = meshio.read(os.path.join(self.output_folder, self.vtu_file_name))

	def get_group_indices(self, group_name):
		bound_id = self.gmsh_tags[1][group_name]
		indices = np.where(self.gmsh_file.cell_data["gmsh:physical"][0] == bound_id)[0]
		edges_on_boundary = self.gmsh_file.cells[0].data[indices].flatten()
		return list(set(edges_on_boundary))


	def get_results(self, group_name, field_name):
		group_idx = self.get_group_indices(group_name)
		coords = self.vtu_file.points[group_idx]
		field = self.vtu_file.point_data[field_name][group_idx]
		x = coords[:,0]
		y = coords[:,1]
		z = coords[:,2]
		return x, y, z, field

class ResultsReader(object):
	def __init__(self, output_folder, pvd_file=None, pvd_file_index=0, gmsh_file_index=0):
		self.output_folder = output_folder
		self.pvd_file_index = pvd_file_index
		self.pvd_file = pvd_file
		self.gmsh_file_index = gmsh_file_index

		self.load_inital_files()
		self.build_tags()
		# self.read_vtu()

	def load_inital_files(self):
		folder_files = os.listdir(self.output_folder)

		vtu_files = self.get_file_names(folder_files, ".vtu")
		pvd_files = self.get_file_names(folder_files, ".pvd")
		gmsh_files = self.get_file_names(folder_files, ".msh")

		# pvd_file = self.check_file_names(file_list=pvd_files, file_index=0)
		# pvd_file = self.check_file_names(file_list=pvd_files, file_index=self.pvd_file_index)
		gmsh_file = self.check_file_names(file_list=gmsh_files, file_index=self.gmsh_file_index)

		# print()
		# print(self.pvd_file)
		# print()

		self.read_gmsh(gmsh_file)
		self.read_pvd(self.pvd_file)

	def get_file_names(self, folder_files, extension):
		return [f for f in folder_files if f.endswith(extension)]

	def check_file_names(self, file_list, file_index):
		if len(file_list) > 1:
			ext = file_list[0].split(".")[1]
			print("WARNING: Too many .%s files. Opening %s."%(ext, file_list[file_index]))
		return file_list[file_index]

	def read_gmsh(self, gmsh_file):
		self.gmsh_file = meshio.read(os.path.join(self.output_folder, gmsh_file))

	def read_pvd(self, pvd_file_name):
		self.timestep_dict = {}
		for line in open(os.path.join(self.output_folder, pvd_file_name)).readlines()[3:-2]:
			line_split = line.split(" ")
			timestep = float(line_split[5][line_split[5].find("=")+2:-1])
			vtu_file_name = line_split[7][line_split[7].find("=")+2:-1]
			self.timestep_dict[timestep] = vtu_file_name

	def build_tags(self):
		self.gmsh_tags = {1:{}, 2:{}, 3:{}}
		for key, value in self.gmsh_file.field_data.items():
			self.gmsh_tags[value[1]][key] = value[0]

	def read_vtu(self):
		self.vtu_file = meshio.read(os.path.join(self.output_folder, self.vtu_file_name))

	def load_points_of_interest(self, fun):
		first_vtu_file_name = list(self.timestep_dict.values())[0]
		vtu_file = meshio.read(os.path.join(self.output_folder, first_vtu_file_name))
		self.points_of_interest = []
		for i, point in enumerate(vtu_file.points):
			is_of_interest = fun(point[0], point[1], point[2])
			if is_of_interest:
				self.points_of_interest.append(i)

	def get_results_over_all_times(self, field_name):
		field = []
		time = []
		for timestep, vtu_file_name in self.timestep_dict.items():
			try:
				vtu_file = meshio.read(os.path.join(self.output_folder, vtu_file_name))
				field.append(vtu_file.point_data[field_name][self.points_of_interest])
				time.append(timestep)
			except:
				print(f"Couldn't read file {vtu_file_name}.")
		coords = vtu_file.points[self.points_of_interest]
		x = coords[:,0]
		y = coords[:,1]
		z = coords[:,2]
		return x, y, z, time, np.array(field)

	def get_results_at_timelist(self, field_name, time_list):
		field = []
		time = []
		for timestep in time_list:
			vtu_file_name = self.timestep_dict[timestep]
			vtu_file = meshio.read(os.path.join(self.output_folder, vtu_file_name))
			field.append(vtu_file.point_data[field_name][self.points_of_interest])
		coords = vtu_file.points[self.points_of_interest]
		x = coords[:,0]
		y = coords[:,1]
		z = coords[:,2]
		return x, y, z, time, np.array(field)

	def get_group_indices(self, group_name):
		bound_id = self.gmsh_tags[1][group_name]
		indices = np.where(self.gmsh_file.cell_data["gmsh:physical"][0] == bound_id)[0]
		edges_on_boundary = self.gmsh_file.cells[0].data[indices].flatten()
		return list(set(edges_on_boundary))

	def get_results(self, group_name, field_name):
		group_idx = self.get_group_indices(group_name)
		coords = self.vtu_file.points[group_idx]
		field = self.vtu_file.point_data[field_name][group_idx]
		x = coords[:,0]
		y = coords[:,1]
		z = coords[:,2]
		return x, y, z, field