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

from mpi4py import MPI
from dolfinx.io import gmshio
from dolfinx.mesh import meshtags
from dolfinx import mesh
import numpy as np
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
		self.comm = MPI.COMM_WORLD
		self.rank = self.comm.rank

		self.load_mesh()
		self.build_tags()
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
		conn = self.mesh.topology.connectivity(3, 0).array.reshape(-1, 4)
		coord = self.mesh.geometry.x
		self.volumes = np.zeros(self.n_elems)
		for i in range(self.n_elems):
			nodes = conn[i]
			x1, y1, z1 = coord[nodes[0], 0], coord[nodes[0], 1], coord[nodes[0], 2]
			x2, y2, z2 = coord[nodes[1], 0], coord[nodes[1], 1], coord[nodes[1], 2]
			x3, y3, z3 = coord[nodes[2], 0], coord[nodes[2], 1], coord[nodes[2], 2]
			x4, y4, z4 = coord[nodes[3], 0], coord[nodes[3], 1], coord[nodes[3], 2]
			self.volumes[i] = self.__tetrahedron_volume(x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4)

	def __build_node_elem_stencil(self):
		conn = self.mesh.topology.connectivity(3, 0).array.reshape(-1, 4)
		coord = self.mesh.geometry.x
		self.stencil = [[] for i in range(self.n_nodes)]
		for elem, elem_conn in enumerate(conn):
			for node in elem_conn:
				if elem not in self.stencil[node]:
					self.stencil[node].append(elem)

	def build_smoother(self):
		self.__compute_volumes()
		self.__build_node_elem_stencil()
		A_row, A_col, A_data = [], [], []
		for node in range(self.n_nodes):
			vol = self.volumes[self.stencil[node]].sum()
			for elem in self.stencil[node]:
				A_row.append(node)
				A_col.append(elem)
				A_data.append(self.volumes[elem]/vol)
		self.A_csr = csr_matrix((A_data, (A_row, A_col)), shape=(self.n_nodes, self.n_elems))
		conn = self.mesh.topology.connectivity(3, 0).array.reshape(-1, 4)
		B_row, B_col, B_data = [], [], []
		for elem, nodes in enumerate(conn):
			for node in nodes:
				B_row.append(elem)
				B_col.append(node)
				B_data.append(1/len(nodes))
		self.B_csr = csr_matrix((B_data, (B_row, B_col)), shape=(self.n_elems, self.n_nodes))
		self.smoother = self.B_csr.dot(self.A_csr)

	def load_mesh(self):
		self.mesh, self.subdomains, self.boundaries = gmshio.read_from_msh(
													    os.path.join(self.grid_folder, f"{self.geometry_name}.msh"),
													    self.comm,
													    rank=0
													)
		self.domain_dim = self.mesh.topology.dim
		self.boundary_dim = self.domain_dim - 1
		self.n_elems = self.mesh.topology.index_map(self.domain_dim).size_local + len(self.mesh.topology.index_map(self.domain_dim).ghosts)
		self.n_nodes = self.mesh.topology.index_map(0).size_local + len(self.mesh.topology.index_map(0).ghosts)

	def build_tags(self):
		grid = meshio.read(os.path.join(self.grid_folder, self.geometry_name+".msh"))
		self.tags = {1:{}, 2:{}, 3:{}}
		for key, value in grid.field_data.items():
			self.tags[value[1]][key] = value[0]
		self.dolfin_tags = self.tags

	def load_subdomains(self):
		self.subdomain_tags = {}
		for subdomain_name in self.get_subdomain_names():
			self.subdomain_tags[subdomain_name] = []


	def load_boundaries(self):
		self.boundary_tags = {}

		for boundary_name in self.get_boundary_names():
			self.boundary_tags[boundary_name] = []
		
		tag_to_name = {fd: name for name, fd in self.dolfin_tags[2].items()}
		boundary_facets = mesh.exterior_facet_indices(self.mesh.topology)
		for i, facet in zip(boundary_facets, self.boundaries.values):
			boundary_name = tag_to_name[facet]
			self.boundary_tags[boundary_name].append(i)


	def build_box_dimensions(self):
		self.Lx = self.mesh.geometry.x[:,0].max() - self.mesh.geometry.x[:,0].min()
		self.Ly = self.mesh.geometry.x[:,1].max() - self.mesh.geometry.x[:,1].min()
		self.Lz = self.mesh.geometry.x[:,2].max() - self.mesh.geometry.x[:,2].min()

	def get_boundaries(self):
		return self.boundaries

	def get_boundary_tags(self, BOUNDARY_NAME):
		if BOUNDARY_NAME == None:
			return None
		else:
			return self.boundary_tags[BOUNDARY_NAME]

	def get_boundary_tag(self, BOUNDARY_NAME):
		if BOUNDARY_NAME == None:
			return None
		else:
			tag_number = self.dolfin_tags[self.boundary_dim][BOUNDARY_NAME]
			return tag_number

	def get_boundary_names(self):
		boundary_names = list(self.dolfin_tags[self.boundary_dim].keys())
		return boundary_names

	def get_subdomain_tag(self, DOMAIN_NAME):
		tag_number = self.dolfin_tags[self.domain_dim][DOMAIN_NAME]
		return tag_number

	def get_subdomains(self):
		return self.subdomains

	def get_subdomain_names(self):
		subdomain_names = list(self.dolfin_tags[self.domain_dim].keys())
		return subdomain_names

	def __extract_grid_data(self):
		self.region_names = self.get_subdomain_names()
		self.n_regions = len(self.region_names)
		self.region_indices = {}
		self.tags_dict = {}

		for i in range(len(self.region_names)):
			self.region_indices[self.region_names[i]] = []
			tag = self.get_subdomain_tag(self.region_names[i])
			self.tags_dict[tag] = self.region_names[i]

		for cell in range(self.n_elems):
			region_marker = self.subdomains.values[cell]
			self.region_indices[self.tags_dict[region_marker]].append(cell)

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



class GridHandlerFEniCS(object):
	def __init__(self, mesh):
		self.mesh = mesh
		self.domain_dim = self.mesh.geometry.dim
		self.boundary_dim = self.domain_dim - 1
		self.n_elems = self.mesh.topology.index_map(self.domain_dim).size_local
		self.n_nodes = self.mesh.topology.index_map(0).size_global
		self.dolfin_tags = {1:{}, 2:{}, 3:{}}

		self.build_box_dimensions()
		self.build_boundaries()
		self.build_dolfin_tags()
		self.load_boundaries()
		self.build_subdomains()
		self.__extract_grid_data()
		# self.build_smoother()

	def build_box_dimensions(self):
		self.Lx = self.mesh.geometry.x[:,0].max() - self.mesh.geometry.x[:,0].min()
		self.Ly = self.mesh.geometry.x[:,1].max() - self.mesh.geometry.x[:,1].min()
		self.Lz = self.mesh.geometry.x[:,2].max() - self.mesh.geometry.x[:,2].min()

	def build_boundaries(self):
		TOL = 1E-10
		Lx = self.Lx
		Ly = self.Ly
		Lz = self.Lz

		boundaries = [	(1, lambda x: np.isclose(x[0], 0., rtol=TOL)),
						(2, lambda x: np.isclose(x[0], self.Lx, rtol=TOL)),
						(3, lambda x: np.isclose(x[1], 0., rtol=TOL)),
						(4, lambda x: np.isclose(x[1], self.Ly, rtol=TOL)),
						(5, lambda x: np.isclose(x[2], 0., rtol=TOL)),
						(6, lambda x: np.isclose(x[2], self.Lz, rtol=TOL))]

		facet_indices, facet_markers = [], []
		for (marker, locator) in boundaries:
			facets = mesh.locate_entities(self.mesh, self.boundary_dim, locator)
			facet_indices.append(facets)
			facet_markers.append(np.full_like(facets, marker))
		facet_indices = np.hstack(facet_indices).astype(np.int32)
		facet_markers = np.hstack(facet_markers).astype(np.int32)
		sorted_facets = np.argsort(facet_indices)
		self.boundaries = mesh.meshtags(self.mesh, self.boundary_dim, facet_indices[sorted_facets], facet_markers[sorted_facets])


	def load_boundaries(self):
		self.boundary_tags = {}

		for boundary_name in self.get_boundary_names():
			self.boundary_tags[boundary_name] = []
		
		tag_to_name = {fd: name for name, fd in self.dolfin_tags[2].items()}
		self.mesh.topology.create_connectivity(self.mesh.topology.dim-1, self.mesh.topology.dim)
		boundary_facets = mesh.exterior_facet_indices(self.mesh.topology)
		for i, facet in zip(boundary_facets, self.boundaries.values):
			boundary_name = tag_to_name[facet]
			self.boundary_tags[boundary_name].append(i)


	def build_subdomains(self):
		cell_indices = [i for i in range(self.n_elems)]
		cell_markers = [0]*self.n_elems
		self.subdomains = mesh.meshtags(self.mesh, self.domain_dim, cell_indices, cell_markers)

	def build_dolfin_tags(self):
		self.dolfin_tags[2]["WEST"] = 1
		self.dolfin_tags[2]["EAST"] = 2
		self.dolfin_tags[2]["SOUTH"] = 3
		self.dolfin_tags[2]["NORTH"] = 4
		self.dolfin_tags[2]["BOTTOM"] = 5
		self.dolfin_tags[2]["TOP"] = 6
		self.dolfin_tags[3]["BODY"] = 1


	def get_boundaries(self):
		return self.boundaries

	def get_boundary_tag(self, BOUNDARY_NAME):
		if BOUNDARY_NAME == None:
			return None
		else:
			return self.dolfin_tags[self.boundary_dim][BOUNDARY_NAME]

	def get_boundary_tags(self, BOUNDARY_NAME):
		if BOUNDARY_NAME == None:
			return None
		else:
			return self.boundary_tags[BOUNDARY_NAME]

	def get_boundary_names(self):
		boundary_names = list(self.dolfin_tags[self.boundary_dim].keys())
		return boundary_names

	def get_subdomain_tag(self, DOMAIN_NAME):
		return self.dolfin_tags[self.domain_dim][DOMAIN_NAME]

	def get_subdomain_names(self):
		subdomain_names = list(self.dolfin_tags[self.domain_dim].keys())
		return subdomain_names

	def get_subdomains(self):
		return self.subdomains

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
			return to.tensor(param)
		else:
			raise Exception("Size of parameter list does not match neither # of elements nor # of regions.")
