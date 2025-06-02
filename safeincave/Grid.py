"""
It reads gmsh grids and gives access to relevant mesh information.
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

from fenics import Mesh, MeshFunction, SubDomain, near
import numpy as np
import dolfin as do
import torch as to
from scipy.sparse import csr_matrix
import meshio
import os

class GridHandlerGMSH(object):
	"""
	This class is responsible to read a Gmsh grid and build the internal
	data structure in a convenient way.

	Parameters
	----------
	geometry_name : str
		Name of the .geo file (usually *geom.geo*).
	grid_folder : str
		Path to where the geometry and grid are located.
	"""
	def __init__(self, geometry_name, grid_folder):
		self.grid_folder = grid_folder
		self.geometry_name = geometry_name

		self.load_mesh()
		self.build_tags()
		self.read_grid_dimensions()
		self.load_subdomains()
		self.load_boundaries()
		self.build_box_dimensions()
		self.__extract_grid_data()
		self.build_smoother()

	def __tetrahedron_volume(self, x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4):
	    volume = abs((1/6) * ((x2 - x1) * ((y3 - y1)*(z4 - z1) - (z3 - z1)*(y4 - y1)) + 
	                 (y2 - y1) * ((z3 - z1)*(x4 - x1) - (x3 - x1)*(z4 - z1)) + 
	                 (z2 - z1) * ((x3 - x1)*(y4 - y1) - (y3 - y1)*(x4 - x1))))
	    return volume

	def __compute_volumes(self):
		conn = self.mesh.cells()
		coord = self.mesh.coordinates()
		self.volumes = np.zeros(self.mesh.num_cells())
		for i in range(self.mesh.num_cells()):
			nodes = conn[i]
			x1, y1, z1 = coord[nodes[0], 0], coord[nodes[0], 1], coord[nodes[0], 2]
			x2, y2, z2 = coord[nodes[1], 0], coord[nodes[1], 1], coord[nodes[1], 2]
			x3, y3, z3 = coord[nodes[2], 0], coord[nodes[2], 1], coord[nodes[2], 2]
			x4, y4, z4 = coord[nodes[3], 0], coord[nodes[3], 1], coord[nodes[3], 2]
			self.volumes[i] = self.__tetrahedron_volume(x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4)

	def __build_node_elem_stencil(self):
		conn = self.mesh.cells()
		coord = self.mesh.coordinates()
		n_nodes = coord.shape[0]
		self.stencil = [[] for i in range(n_nodes)]
		for elem, elem_conn in enumerate(conn):
			for node in elem_conn:
				if elem not in self.stencil[node]:
					self.stencil[node].append(elem)

	def build_smoother(self):
		self.__compute_volumes()
		self.__build_node_elem_stencil()
		n_nodes = self.mesh.num_vertices()
		n_elems = self.mesh.num_cells()
		A_row, A_col, A_data = [], [], []
		for node in range(n_nodes):
			vol = self.volumes[self.stencil[node]].sum()
			for elem in self.stencil[node]:
				A_row.append(node)
				A_col.append(elem)
				A_data.append(self.volumes[elem]/vol)
		self.A_csr = csr_matrix((A_data, (A_row, A_col)), shape=(n_nodes, n_elems))
		B_row, B_col, B_data = [], [], []
		for elem, nodes in enumerate(self.mesh.cells()):
			for node in nodes:
				B_row.append(elem)
				B_col.append(node)
				B_data.append(1/len(nodes))
		self.B_csr = csr_matrix((B_data, (B_row, B_col)), shape=(n_elems, n_nodes))
		self.smoother = self.B_csr.dot(self.A_csr)

	def load_mesh(self):
		"""
		Reads the .xml file containing the mesh.
		"""
		file_name_xml = os.path.join(self.grid_folder, self.geometry_name+".xml")
		self.mesh = Mesh(file_name_xml)

	def build_tags(self):
		"""
		Reads the mesh tags for lines (1D), surfaces(2D) and volumes (3D) and
		creates a dictionary associating the entity name to an unique integer.
		"""
		file_name_msh = os.path.join(self.grid_folder, self.geometry_name+".msh")
		grid = meshio.read(file_name_msh)
		self.tags = {1:{}, 2:{}, 3:{}}
		for key, value in grid.field_data.items():
			self.tags[value[1]][key] = value[0]
		self.dolfin_tags = self.tags

	def read_grid_dimensions(self):
		"""
		Reads the grid dimensions.
		"""
		self.domain_dim = self.mesh.topology().dim()
		self.boundary_dim = self.domain_dim - 1

	def load_subdomains(self):
		"""
		Renumbers subdomain (3D) tag numbers.
		"""
		file_name_physical_region_xml = os.path.join(self.grid_folder, self.geometry_name+"_physical_region.xml")
		subdomains0 = MeshFunction("size_t", self.mesh, file_name_physical_region_xml)
		self.subdomains = MeshFunction("size_t", self.mesh, self.domain_dim)
		self.subdomains.set_all(0)
		for i, value in enumerate(self.tags[self.domain_dim].items()):
			self.subdomains.array()[subdomains0.array() == value[1]] = i + 1
			self.dolfin_tags[self.domain_dim][value[0]] = i + 1
		# print("Hello")
		# print(self.subdomains.array())


	def load_boundaries(self):
		"""
		Renumbers boundaries (2D) tag numbers.
		"""
		file_name_facet_region_xml = os.path.join(self.grid_folder, self.geometry_name+"_facet_region.xml")
		boundaries0 = MeshFunction("size_t", self.mesh, file_name_facet_region_xml)
		self.boundaries = MeshFunction("size_t", self.mesh, self.mesh.topology().dim() - 1)
		self.boundaries.set_all(0)
		for i, value in enumerate(self.tags[self.boundary_dim].items()):
			self.boundaries.array()[boundaries0.array() == value[1]] = i + 1
			self.dolfin_tags[self.boundary_dim][value[0]] = i + 1

	def build_box_dimensions(self):
		"""
		Reads box dimensions, that is, maximum *x*, *y* and *z* coordinates.
		"""
		self.Lx = self.mesh.coordinates()[:,0].max() - self.mesh.coordinates()[:,0].min()
		self.Ly = self.mesh.coordinates()[:,1].max() - self.mesh.coordinates()[:,1].min()
		self.Lz = self.mesh.coordinates()[:,2].max() - self.mesh.coordinates()[:,2].min()

	def get_boundaries(self):
		"""
		Get mesh boundaries. It is used when applying Dirichlet boundary conditions.

		Returns
		-------
		boundaries : dolfin.cpp.mesh.MeshFunctionSizet
			Mesh boundaries.
		"""
		return self.boundaries

	def get_boundary_tags(self, BOUNDARY_NAME):
		"""
		Get boundary tag (integer) corresponding to *BOUNDARY_NAME*.

		Parameters
		----------
		BOUNDARY_NAME : str
			Name of the boundary.

		Returns
		-------
		tag_number : int
			Integer representing BOUNDARY_NAME
		"""
		if BOUNDARY_NAME == None:
			return None
		else:
			tag_number = self.dolfin_tags[self.boundary_dim][BOUNDARY_NAME]
			return tag_number

	def get_boundary_names(self):
		"""
		Provides the names of all the boundaries.

		Returns
		-------
		boundary_names : dict_keys
			List of strings containing the boundary names.
		"""
		boundary_names = self.dolfin_tags[self.boundary_dim].keys()
		return boundary_names

	def get_subdomain_tags(self, DOMAIN_NAME):
		"""
		Get subdomain tag (integer) corresponding to *DOMAIN_NAME*.

		Parameters
		----------
		DOMAIN_NAME : str
			Name of the subdomain.

		Returns
		-------
		tag_number : int
			Integer representing DOMAIN_NAME
		"""
		tag_number = self.dolfin_tags[self.domain_dim][DOMAIN_NAME]
		return tag_number

	def get_subdomains(self):
		"""
		Get mesh subdomains. It can be used for solving different models in different subdomains.

		Returns
		-------
		subdomains : dolfin.cpp.mesh.MeshFunctionSizet
			Mesh subdomains.
		"""
		return self.subdomains

	def get_subdomain_names(self):
		"""
		Provides the names of all the subdomains.

		Returns
		-------
		subdomain_names : dict_keys
			List of strings containing the subdomain names.
		"""
		subdomain_names = self.dolfin_tags[self.domain_dim].keys()
		return subdomain_names

	def __extract_grid_data(self):
		self.region_names = list(self.get_subdomain_names())
		self.n_regions = len(self.region_names)
		self.n_elems = self.mesh.num_cells()
		self.n_nodes = self.mesh.num_vertices()
		self.region_indices = {}
		self.tags_dict = {}
		for i in range(len(self.region_names)):
			self.region_indices[self.region_names[i]] = []
			tag = self.get_subdomain_tags(self.region_names[i])
			self.tags_dict[tag] = self.region_names[i]

		for cell in do.cells(self.mesh):
			region_marker = self.subdomains[cell]
			self.region_indices[self.tags_dict[region_marker]].append(cell.index())

	def get_parameter(self, param):
		if type(param) == int or type(param) == float:
			return to.tensor([param for i in range(self.n_elems)])
		elif len(param) == self.n_regions:
			param_to = to.zeros(self.n_elems)
			for i, region in enumerate(self.region_indices.keys()):
				param_to[self.region_indices[region]] = param[i]
			return param_to
		elif len(param) == self.n_elems:
			if type(param) == to.Tensor:
				return param
			else:
				return to.tensor(param)
		else:
			raise Exception("Size of parameter list does not match neither # of elements nor # of regions.")



# class GridHandlerFEniCS(object):
# 	"""
# 	This class is responsible to read a FEniCS mesh of a brick shape (not 
# 	necessarily a cube)	and build the internal data structure in a convenient way.

# 	Parameters
# 	----------
# 	geometry_name : str
# 		Name of the .geo file (usually *geom.geo*).
# 	grid_folder : str
# 		Path to where the geometry and grid are located.
# 	"""
# 	def __init__(self, mesh):
# 		self.mesh = mesh
# 		self.domain_dim = self.mesh.topology().dim()
# 		self.dolfin_tags = {1:{}, 2:{}, 3:{}}

# 		self.build_box_dimensions()
# 		self.build_grid_dimensions()
# 		self.build_boundaries()
# 		self.build_dolfin_tags()
# 		self.build_subdomains()
# 		self.__extract_grid_data()
# 		self.build_smoother()

# 	def __tetrahedron_volume(self, x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4):
# 	    volume = abs((1/6) * ((x2 - x1) * ((y3 - y1)*(z4 - z1) - (z3 - z1)*(y4 - y1)) + 
# 	                 (y2 - y1) * ((z3 - z1)*(x4 - x1) - (x3 - x1)*(z4 - z1)) + 
# 	                 (z2 - z1) * ((x3 - x1)*(y4 - y1) - (y3 - y1)*(x4 - x1))))
# 	    return volume

# 	def __compute_volumes(self):
# 		conn = self.mesh.cells()
# 		coord = self.mesh.coordinates()
# 		self.volumes = np.zeros(self.mesh.num_cells())
# 		for i in range(self.mesh.num_cells()):
# 			nodes = conn[i]
# 			x1, y1, z1 = coord[nodes[0], 0], coord[nodes[0], 1], coord[nodes[0], 2]
# 			x2, y2, z2 = coord[nodes[1], 0], coord[nodes[1], 1], coord[nodes[1], 2]
# 			x3, y3, z3 = coord[nodes[2], 0], coord[nodes[2], 1], coord[nodes[2], 2]
# 			x4, y4, z4 = coord[nodes[3], 0], coord[nodes[3], 1], coord[nodes[3], 2]
# 			self.volumes[i] = self.__tetrahedron_volume(x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4)

# 	def __build_node_elem_stencil(self):
# 		conn = self.mesh.cells()
# 		coord = self.mesh.coordinates()
# 		n_nodes = coord.shape[0]
# 		self.stencil = [[] for i in range(n_nodes)]
# 		for elem, elem_conn in enumerate(conn):
# 			for node in elem_conn:
# 				if elem not in self.stencil[node]:
# 					self.stencil[node].append(elem)

# 	def build_smoother(self):
# 		self.__compute_volumes()
# 		self.__build_node_elem_stencil()
# 		n_nodes = self.mesh.num_vertices()
# 		n_elems = self.mesh.num_cells()
# 		A_row, A_col, A_data = [], [], []
# 		for node in range(n_nodes):
# 			vol = self.volumes[self.stencil[node]].sum()
# 			for elem in self.stencil[node]:
# 				A_row.append(node)
# 				A_col.append(elem)
# 				A_data.append(self.volumes[elem]/vol)
# 		self.A_csr = csr_matrix((A_data, (A_row, A_col)), shape=(n_nodes, n_elems))
# 		B_row, B_col, B_data = [], [], []
# 		for elem, nodes in enumerate(self.mesh.cells()):
# 			for node in nodes:
# 				B_row.append(elem)
# 				B_col.append(node)
# 				B_data.append(1/len(nodes))
# 		self.B_csr = csr_matrix((B_data, (B_row, B_col)), shape=(n_elems, n_nodes))
# 		self.smoother = self.B_csr.dot(self.A_csr)

# 	def build_grid_dimensions(self):
# 		"""
# 		Reads the grid dimensions.
# 		"""
# 		self.domain_dim = self.mesh.topology().dim()
# 		self.boundary_dim = self.domain_dim - 1

# 	def build_box_dimensions(self):
# 		"""
# 		Reads box dimensions, that is, maximum *x*, *y* and *z* coordinates.
# 		"""
# 		self.Lx = self.mesh.coordinates()[:,0].max() - self.mesh.coordinates()[:,0].min()
# 		self.Ly = self.mesh.coordinates()[:,1].max() - self.mesh.coordinates()[:,1].min()
# 		self.Lz = self.mesh.coordinates()[:,2].max() - self.mesh.coordinates()[:,2].min()

# 	def build_boundaries(self):
# 		"""
# 		Builds list of boundaries for faces EAST, WEST, NORTH, SOUTH, BOTTOM and TOP.
# 		"""
# 		TOL = 1E-14
# 		Lx = self.Lx
# 		Ly = self.Ly
# 		Lz = self.Lz
# 		class boundary_facet_WEST(SubDomain):
# 			def inside(self, x, on_boundary):
# 				return on_boundary and near(x[0], 0.0, TOL)
# 		class boundary_facet_EAST(SubDomain):
# 			def inside(self, x, on_boundary):
# 				return on_boundary and near(x[0], Lx, TOL)
# 		class boundary_facet_SOUTH(SubDomain):
# 			def inside(self, x, on_boundary):
# 				return on_boundary and near(x[1], 0, TOL)
# 		class boundary_facet_NORTH(SubDomain):
# 			def inside(self, x, on_boundary):
# 				return on_boundary and near(x[1], Ly, TOL)
# 		class boundary_facet_BOTTOM(SubDomain):
# 			def inside(self, x, on_boundary):
# 				return on_boundary and near(x[2], 0.0, TOL)
# 		class boundary_facet_TOP(SubDomain):
# 			def inside(self, x, on_boundary):
# 				return on_boundary and near(x[2], Lz, TOL)
# 		self.boundaries = MeshFunction("size_t", self.mesh, self.mesh.topology().dim()-1)
# 		self.boundaries.set_all(0)
# 		boundary_facet_WEST().mark(self.boundaries, 1)
# 		boundary_facet_EAST().mark(self.boundaries, 2)
# 		boundary_facet_SOUTH().mark(self.boundaries, 3)
# 		boundary_facet_NORTH().mark(self.boundaries, 4)
# 		boundary_facet_BOTTOM().mark(self.boundaries, 5)
# 		boundary_facet_TOP().mark(self.boundaries, 6)

# 	def build_subdomains(self):
# 		"""
# 		Builds subdomains, which in this case is just one.
# 		"""
# 		self.subdomains = MeshFunction("size_t", self.mesh, self.domain_dim)
# 		self.subdomains.set_all(0)
# 		class Omega(SubDomain):
# 			def inside(self, x, on_boundary):
# 				return True
# 		Omega().mark(self.subdomains, 0)

# 	def build_dolfin_tags(self):
# 		"""
# 		Defines tags for boundaries (2D) and subdomains (3D).
# 		"""
# 		self.dolfin_tags[2]["WEST"] = 1
# 		self.dolfin_tags[2]["EAST"] = 2
# 		self.dolfin_tags[2]["SOUTH"] = 3
# 		self.dolfin_tags[2]["NORTH"] = 4
# 		self.dolfin_tags[2]["BOTTOM"] = 5
# 		self.dolfin_tags[2]["TOP"] = 6
# 		self.dolfin_tags[3]["BODY"] = 1


# 	def get_boundaries(self):
# 		"""
# 		Get mesh boundaries. It is used when applying Dirichlet boundary conditions.

# 		Returns
# 		-------
# 		boundaries : dolfin.cpp.mesh.MeshFunctionSizet
# 			Mesh boundaries.
# 		"""
# 		return self.boundaries

# 	def get_boundary_tags(self, BOUNDARY_NAME):
# 		"""
# 		Get boundary tag (integer) corresponding to *BOUNDARY_NAME*.

# 		Parameters
# 		----------
# 		BOUNDARY_NAME : str
# 			Name of the boundary.

# 		Returns
# 		-------
# 		tag_number : int
# 			Integer representing BOUNDARY_NAME
# 		"""
# 		if BOUNDARY_NAME == None:
# 			return None
# 		else:
# 			return self.dolfin_tags[self.boundary_dim][BOUNDARY_NAME]

# 	def get_boundary_names(self):
# 		"""
# 		Provides the names of all the boundaries.

# 		Returns
# 		-------
# 		boundary_names : dict_keys
# 			List of strings containing the boundary names.
# 		"""
# 		boundary_names = self.dolfin_tags[self.boundary_dim].keys()
# 		return boundary_names

# 	def get_subdomain_tags(self, DOMAIN_NAME):
# 		"""
# 		Get subdomain tag (integer) corresponding to *DOMAIN_NAME*.

# 		Parameters
# 		----------
# 		DOMAIN_NAME : str
# 			Name of the subdomain.

# 		Returns
# 		-------
# 		tag_number : int
# 			Integer representing DOMAIN_NAME
# 		"""
# 		return self.dolfin_tags[self.domain_dim][DOMAIN_NAME]

# 	def get_subdomain_names(self):
# 		"""
# 		Provides the names of all the subdomains.

# 		Returns
# 		-------
# 		subdomain_names : dict_keys
# 			List of strings containing the subdomain names.
# 		"""
# 		subdomain_names = self.dolfin_tags[self.domain_dim].keys()
# 		return subdomain_names

# 	def get_subdomains(self):
# 		"""
# 		Get mesh subdomains. It can be used for solving different models in different subdomains.

# 		Returns
# 		-------
# 		subdomains : dolfin.cpp.mesh.MeshFunctionSizet
# 			Mesh subdomains.
# 		"""
# 		return self.subdomains

# 	def __extract_grid_data(self):
# 		self.region_names = list(self.get_subdomain_names())
# 		self.n_regions = len(self.region_names)
# 		self.n_elems = self.mesh.num_cells()
# 		self.n_nodes = self.mesh.num_vertices()
# 		self.region_indices = {}
# 		self.tags_dict = {}
# 		for i in range(len(self.region_names)):
# 			self.region_indices[self.region_names[i]] = []
# 			tag = self.get_subdomain_tags(self.region_names[i])
# 			self.tags_dict[tag] = self.region_names[i]

# 		for cell in do.cells(self.mesh):
# 			region_marker = self.subdomains[cell]
# 			self.region_indices[self.tags_dict[region_marker]].append(cell.index())

# 	def get_parameter(self, param):
# 		if type(param) == int or type(param) == float:
# 			return to.tensor([param for i in range(self.n_elems)])
# 		elif len(param) == self.n_regions:
# 			param_to = to.zeros(self.n_elems)
# 			for i, region in enumerate(self.region_indices.keys()):
# 				param_to[self.region_indices[region]] = param[i]
# 			return param_to
# 		elif len(param) == self.n_elems:
# 			return to.tensor(param)
# 		else:
# 			raise Exception("Size of parameter list does not match neither # of elements nor # of regions.")
