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
	Handler for reading a Gmsh mesh into DOLFINx and exposing convenient
	grid-related utilities (tags, regions, volumes, smoothers).

	The constructor loads the mesh, builds tag maps, extracts subdomain and
	boundary metadata, computes bounding-box dimensions, derives region-wise
	element indices, and constructs node–element smoothing operators.

	Parameters
	----------
	geometry_name : str
	    Base name (without extension) of the `.msh` file to read.
	grid_folder : str
	    Directory where the mesh file `{geometry_name}.msh` resides.

	Attributes
	----------
	grid_folder : str
	    Mesh directory provided at construction.
	geometry_name : str
	    Mesh base name provided at construction.
	comm : MPI.Comm
	    MPI communicator (defaults to `MPI.COMM_WORLD`).
	rank : int
	    Rank of the current process.
	mesh : dolfinx.mesh.Mesh
	    Loaded DOLFINx mesh.
	subdomains : dolfinx.mesh.MeshTags
	    Cell (volume) tags read from the mesh.
	boundaries : dolfinx.mesh.MeshTags
	    Facet (surface) tags read from the mesh.
	domain_dim : int
	    Topological dimension of the mesh cells.
	boundary_dim : int
	    Dimension of boundary entities (`domain_dim - 1`).
	n_elems : int
	    Number of local cells including ghosts.
	n_nodes : int
	    Number of local vertices including ghosts.
	tags : dict[int, dict[str, int]]
	    Mapping `dimension -> {name -> tag_id}` parsed from the Gmsh file.
	dolfin_tags : dict
	    Alias of `tags` for convenience.
	subdomain_tags : dict[str, list[int]]
	    Placeholder for subdomain-wise cell indices (filled later).
	boundary_tags : dict[str, list[int]]
	    Mapping from boundary name to list of exterior facet indices.
	Lx, Ly, Lz : float
	    Extents of the mesh bounding box in x, y, z.
	region_names : list[str]
	    List of subdomain (region) names.
	n_regions : int
	    Number of regions.
	region_indices : dict[str, list[int]]
	    Mapping from region name to list of cell indices in that region.
	tags_dict : dict[int, str]
	    Reverse mapping `{tag_id -> region_name}`.
	volumes : numpy.ndarray
	    Per-cell volumes, shape `(n_elems,)`. Created by `build_smoother()`.
	stencil : list[list[int]]
	    Node-to-element adjacency (per-node list of incident cell indices).
	A_csr : scipy.sparse.csr_matrix
	    Node-to-element averaging weights, shape `(n_nodes, n_elems)`.
	B_csr : scipy.sparse.csr_matrix
	    Element-to-node averaging weights, shape `(n_elems, n_nodes)`.
	smoother : scipy.sparse.csr_matrix
	    Element-wise smoother `B_csr @ A_csr`, shape `(n_elems, n_elems)`.

	Notes
	-----
	- Assumes a tetrahedral volume mesh for volume and centroid calculations.
	- Counts (`n_elems`, `n_nodes`) include local ghosts for parallel runs.
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
		"""
		Compute the absolute volume of a tetrahedron from its 4 vertices.

		Parameters
		----------
		x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4 : float
		    Coordinates of the four tetrahedron vertices.

		Returns
		-------
		float
		    Absolute volume of the tetrahedron.

		Notes
		-----
		Uses the scalar triple product formula:
		:math:`V = |((a-b) \\cdot ((c-b) \\times (d-b)))/6|`.
		"""
		volume = abs((1/6) * ((x2 - x1) * ((y3 - y1)*(z4 - z1) - (z3 - z1)*(y4 - y1)) + 
		             (y2 - y1) * ((z3 - z1)*(x4 - x1) - (x3 - x1)*(z4 - z1)) + 
		             (z2 - z1) * ((x3 - x1)*(y4 - y1) - (y3 - y1)*(x4 - x1))))
		return volume

	def __compute_volumes(self):
		"""
        Compute per-cell volumes for all tetrahedral elements.

        Parameters
        ----------
        None

        Returns
        -------
        None

        Side Effects
        ------------
        volumes : numpy.ndarray
            Sets `self.volumes` with shape `(n_elems,)`.

        Notes
        -----
        Extracts connectivity (3→0) and coordinates from the DOLFINx mesh and
        applies `__tetrahedron_volume` to each cell.
        """
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
		"""
        Build node-to-element adjacency lists.

        Parameters
        ----------
        None

        Returns
        -------
        None

        Side Effects
        ------------
        stencil : list[list[int]]
            Sets `self.stencil` such that `self.stencil[i]` contains the
            indices of elements shared by node `i`.
        """
		conn = self.mesh.topology.connectivity(3, 0).array.reshape(-1, 4)
		coord = self.mesh.geometry.x
		self.stencil = [[] for i in range(self.n_nodes)]
		for elem, elem_conn in enumerate(conn):
			for node in elem_conn:
				if elem not in self.stencil[node]:
					self.stencil[node].append(elem)

	def build_smoother(self):
		"""
        Construct element–node smoothing operators and their product.

        Parameters
        ----------
        None

        Returns
        -------
        None

        Side Effects
        ------------
        A_csr : scipy.sparse.csr_matrix
            Node-to-element weights, shape `(n_nodes, n_elems)`.
        B_csr : scipy.sparse.csr_matrix
            Element-to-node weights, shape `(n_elems, n_nodes)`.
        smoother : scipy.sparse.csr_matrix
            Product `B_csr @ A_csr`, shape `(n_elems, n_elems)`.

        Notes
        -----
        - `A_csr[i, e] = vol_e / sum(vol_e' for e' in stencil[i])`.
        - `B_csr[e, i] = 1/4` for tetrahedra (uniform average over the 4 nodes).
        """
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
		"""
        Load mesh and tag metadata from a Gmsh `.msh` file.

        Parameters
        ----------
        None

        Returns
        -------
        None

        Side Effects
        ------------
        mesh : dolfinx.mesh.Mesh
            Sets `self.mesh`.
        subdomains : dolfinx.mesh.MeshTags
            Sets `self.subdomains` (cell tags).
        boundaries : dolfinx.mesh.MeshTags
            Sets `self.boundaries` (facet tags).
        domain_dim : int
            Sets `self.domain_dim` from mesh topology.
        boundary_dim : int
            Sets `self.boundary_dim = domain_dim - 1`.
        n_elems, n_nodes : int
            Sets local counts including ghosts.

        Notes
        -----
        Uses `gmshio.read_from_msh(os.path.join(grid_folder, f"{geometry_name}.msh"), comm, rank=0)`.
        """
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
		"""
        Parse Gmsh field data into a dimension→name→tag mapping.

        Parameters
        ----------
        None

        Returns
        -------
        None

        Side Effects
        ------------
        dolfin_tags : dict[int, dict[str, int]]
            Populates `self.dolfin_tags` with entries for dims 1, 2, 3.

        Notes
        -----
        Reads the `.msh` via `meshio.read` to access `field_data`.
        """
		grid = meshio.read(os.path.join(self.grid_folder, self.geometry_name+".msh"))
		# self.tags = {1:{}, 2:{}, 3:{}}
		# for key, value in grid.field_data.items():
		# 	self.tags[value[1]][key] = value[0]
		# self.dolfin_tags = self.tags
		self.dolfin_tags = {1:{}, 2:{}, 3:{}}
		for key, value in grid.field_data.items():
			self.dolfin_tags[value[1]][key] = value[0]

	def load_subdomains(self):
		"""
        Initialize container for subdomain cell indices per region name.

        Parameters
        ----------
        None

        Returns
        -------
        None

        Side Effects
        ------------
        subdomain_tags : dict[str, list[int]]
            Initializes empty lists keyed by subdomain names.
        """
		self.subdomain_tags = {}
		for subdomain_name in self.get_subdomain_names():
			self.subdomain_tags[subdomain_name] = []


	def load_boundaries(self):
		"""
        Build a map from boundary name to exterior facet indices.

        Parameters
        ----------
        None

        Returns
        -------
        None

        Side Effects
        ------------
        boundary_tags : dict[str, list[int]]
            Populates with facet indices for each named boundary.

        Notes
        -----
        Uses `mesh.exterior_facet_indices(self.mesh.topology)` and the facet
        tag values in `self.boundaries` to assign names via `self.dolfin_tags[2]`.
        """
		self.boundary_tags = {}

		for boundary_name in self.get_boundary_names():
			self.boundary_tags[boundary_name] = []
		
		tag_to_name = {fd: name for name, fd in self.dolfin_tags[2].items()}
		boundary_facets = mesh.exterior_facet_indices(self.mesh.topology)
		for i, facet in zip(boundary_facets, self.boundaries.values):
			boundary_name = tag_to_name[facet]
			self.boundary_tags[boundary_name].append(i)


	def build_box_dimensions(self):
		"""
        Compute axis-aligned bounding-box extents (Lx, Ly, Lz).

        Parameters
        ----------
        None

        Returns
        -------
        None

        Side Effects
        ------------
        Lx, Ly, Lz : float
            Sets the extents along x, y, z.
        """
		self.Lx = self.mesh.geometry.x[:,0].max() - self.mesh.geometry.x[:,0].min()
		self.Ly = self.mesh.geometry.x[:,1].max() - self.mesh.geometry.x[:,1].min()
		self.Lz = self.mesh.geometry.x[:,2].max() - self.mesh.geometry.x[:,2].min()

	def get_boundaries(self):
		"""
        Return the facet `MeshTags` object.

        Returns
        -------
        dolfinx.mesh.MeshTags
            The boundary tags read from the mesh.
        """
		return self.boundaries

	def get_boundary_tags(self, BOUNDARY_NAME):
		"""
        Get list of exterior facet indices for a named boundary.

        Parameters
        ----------
        BOUNDARY_NAME : str or None
            Boundary name as defined in the Gmsh field data. If `None`,
            returns `None`.

        Returns
        -------
        list[int] or None
            Facet indices on the exterior boundary corresponding to the name,
            or `None` if `BOUNDARY_NAME` is `None`.
        """
		if BOUNDARY_NAME == None:
			return None
		else:
			return self.boundary_tags[BOUNDARY_NAME]

	def get_boundary_tag(self, BOUNDARY_NAME):
		"""
        Get the integer tag ID for a named boundary.

        Parameters
        ----------
        BOUNDARY_NAME : str or None
            Boundary name. If `None`, returns `None`.

        Returns
        -------
        int or None
            Integer tag in `self.dolfin_tags[self.boundary_dim]`, or `None`.
        """
		if BOUNDARY_NAME == None:
			return None
		else:
			tag_number = self.dolfin_tags[self.boundary_dim][BOUNDARY_NAME]
			return tag_number

	def get_boundary_names(self):
		"""
        List all boundary names present in the mesh.

        Returns
        -------
        list[str]
            Boundary names from `self.dolfin_tags[self.boundary_dim]`.
        """
		boundary_names = list(self.dolfin_tags[self.boundary_dim].keys())
		return boundary_names

	def get_subdomain_tag(self, DOMAIN_NAME):
		"""
        Get the integer tag ID for a named subdomain (cell region).

        Parameters
        ----------
        DOMAIN_NAME : str
            Subdomain name.

        Returns
        -------
        int
            Integer tag for the subdomain in `self.dolfin_tags[self.domain_dim]`.
        """
		tag_number = self.dolfin_tags[self.domain_dim][DOMAIN_NAME]
		return tag_number

	def get_subdomains(self):
		"""
        Return the cell `MeshTags` object.

        Returns
        -------
        dolfinx.mesh.MeshTags
            The subdomain (cell) tags read from the mesh.
        """
		return self.subdomains

	def get_subdomain_names(self):
		"""
        List all subdomain (region) names present in the mesh.

        Returns
        -------
        list[str]
            Subdomain names from `self.dolfin_tags[self.domain_dim]`.
        """
		subdomain_names = list(self.dolfin_tags[self.domain_dim].keys())
		return subdomain_names

	def __extract_grid_data(self):
		"""
        Build region indices for elements based on cell tags.

        Parameters
        ----------
        None

        Returns
        -------
        None

        Side Effects
        ------------
        region_names : list[str]
            Sets from `get_subdomain_names()`.
        n_regions : int
            Total number of regions.
        region_indices : dict[str, list[int]]
            Maps each region name to a list of cell indices.
        tags_dict : dict[int, str]
            Reverse mapping from integer tag to region name.

        Notes
        -----
        Iterates over local elements; if running in parallel, these are
        local (including ghosts).
        """
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
		"""
        Expand a parameter specification to per-element values.

        Parameters
        ----------
        param : int or float or sequence or torch.Tensor
            - Scalar (`int`/`float`): broadcast to all elements.
            - Sequence of length `n_regions`: values per region in the
              order of `self.region_indices.keys()`.
            - Sequence or tensor of length `n_elems`: per-element values.

        Returns
        -------
        torch.Tensor
            1D tensor of length `n_elems` with parameter values.

        Raises
        ------
        Exception
            If `param` length is neither `n_regions` nor `n_elems`.

        Notes
        -----
        Converts Python sequences to `torch.Tensor` when needed. For the
        region-wise case, elements are filled according to
        `self.region_indices[region]` for each region name.
        """
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
