# Copyright (c) 2026, The SafeInCave Developers
#
# SPDX-License-Identifier: BSD-3-Clause

from mpi4py import MPI
from dolfinx.io import gmshio
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
        volume = abs(
            (1 / 6)
            * (
                (x2 - x1) * ((y3 - y1) * (z4 - z1) - (z3 - z1) * (y4 - y1))
                + (y2 - y1) * ((z3 - z1) * (x4 - x1) - (x3 - x1) * (z4 - z1))
                + (z2 - z1) * ((x3 - x1) * (y4 - y1) - (y3 - y1) * (x4 - x1))
            )
        )
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
            self.volumes[i] = self.__tetrahedron_volume(
                x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4
            )

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
                A_data.append(self.volumes[elem] / vol)
        self.A_csr = csr_matrix(
            (A_data, (A_row, A_col)), shape=(self.n_nodes, self.n_elems)
        )
        conn = self.mesh.topology.connectivity(3, 0).array.reshape(-1, 4)
        B_row, B_col, B_data = [], [], []
        for elem, nodes in enumerate(conn):
            for node in nodes:
                B_row.append(elem)
                B_col.append(node)
                B_data.append(1 / len(nodes))
        self.B_csr = csr_matrix(
            (B_data, (B_row, B_col)), shape=(self.n_elems, self.n_nodes)
        )
        self.smoother = self.B_csr.dot(self.A_csr)
        self._cg1_space_cache = {}

    def smooth_dg0_to_cg1(self, dg0_func):
        """
        Average a DG0 field (scalar/vector/tensor) to CG1 nodes via A_csr component-wise.

        Parameters
        ----------
        dg0_func : dolfinx.fem.Function
            A function in a DG0 scalar, vector, or tensor space.

        Returns
        -------
        dolfinx.fem.Function
            Function in an equivalent CG1 space with smoothed nodal values.

        Notes
        -----
        Uses the pre-built volume-weighted averaging matrix `self.A_csr`
        (shape (n_nodes, n_elems)) applied independently to each component.
        CG1 function spaces are cached to avoid repeated construction.
        """
        import dolfinx as do

        ufl_elem = dg0_func.function_space.ufl_element()
        shape = ufl_elem.reference_value_shape

        if shape not in self._cg1_space_cache:
            if shape == ():
                cg1_space = do.fem.functionspace(self.mesh, ("Lagrange", 1))
            else:
                cg1_space = do.fem.functionspace(self.mesh, ("Lagrange", 1, shape))
            self._cg1_space_cache[shape] = cg1_space
        else:
            cg1_space = self._cg1_space_cache[shape]

        arr = dg0_func.x.array.reshape(self.n_elems, -1)
        smoothed = np.stack([self.A_csr.dot(arr[:, k]) for k in range(arr.shape[1])], axis=-1)

        cg1_func = do.fem.Function(cg1_space)
        cg1_func.x.array[:] = smoothed.flatten()
        return cg1_func

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
            rank=0,
        )
        self.domain_dim = self.mesh.topology.dim
        self.boundary_dim = self.domain_dim - 1
        self.n_elems = self.mesh.topology.index_map(self.domain_dim).size_local + len(
            self.mesh.topology.index_map(self.domain_dim).ghosts
        )
        self.n_nodes = self.mesh.topology.index_map(0).size_local + len(
            self.mesh.topology.index_map(0).ghosts
        )

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
        grid = meshio.read(os.path.join(self.grid_folder, self.geometry_name + ".msh"))
        # self.tags = {1:{}, 2:{}, 3:{}}
        # for key, value in grid.field_data.items():
        # 	self.tags[value[1]][key] = value[0]
        # self.dolfin_tags = self.tags
        self.dolfin_tags = {1: {}, 2: {}, 3: {}}
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
        self.Lx = self.mesh.geometry.x[:, 0].max() - self.mesh.geometry.x[:, 0].min()
        self.Ly = self.mesh.geometry.x[:, 1].max() - self.mesh.geometry.x[:, 1].min()
        self.Lz = self.mesh.geometry.x[:, 2].max() - self.mesh.geometry.x[:, 2].min()

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
        if BOUNDARY_NAME is None:
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
        if BOUNDARY_NAME is None:
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
        if type(param) is int or type(param) is float:
            return to.tensor([param for i in range(self.n_elems)])
        elif len(param) == self.n_regions:
            param_to = to.zeros(self.n_elems)
            for i, region in enumerate(self.region_indices.keys()):
                param_to[self.region_indices[region]] = param[i]
            return param_to
        elif len(param) == self.n_elems:
            if type(param) is to.Tensor:
                return param
            else:
                return to.tensor(param)
        else:
            raise Exception(
                "Size of parameter list does not match neither # of elements nor # of regions."
            )


# ============================================================================
# GridHandlerPythonScript: Generate mesh via external Python scripts
# ============================================================================


class GridHandlerPythonScript(GridHandlerGMSH):
    """
    Handler for generating mesh via external Python scripts and reading into DOLFINx.

    This handler loads and executes a Python script containing a mesh generation
    function, passes user-provided parameters to that function, and uses the
    returned DOLFINx mesh objects to initialize the handler (compatible with
    GridHandlerGMSH interface).

    The script must define a `main(**kwargs)` function that returns a tuple of
    (mesh, cell_tags, facet_tags), which are DOLFINx Mesh and MeshTags objects.

    Parameters
    ----------
    script_path : str
        Path to the Python script containing the mesh generation function.
    parameters : dict, optional
        Dictionary of keyword arguments to pass to main(). Default is empty dict {}.
    function_name : str, optional
        Deprecated; ignored. For compatibility only. The script must define main().

    Attributes
    ----------
    script_path : str
        Path to the mesh generation script.
    parameters : dict
        Parameters passed to main().
    mesh : dolfinx.mesh.Mesh
        Loaded DOLFINx mesh from the generated mesh.
    subdomains : dolfinx.mesh.MeshTags
        Cell (volume) tags from the generated mesh.
    boundaries : dolfinx.mesh.MeshTags
        Facet (surface) tags from the generated mesh.

    Notes
    -----
    - Inherits from GridHandlerGMSH for compatibility with downstream code.
    - No .msh file is read from disk; mesh is generated in-memory by the script.
    - Script execution is isolated to avoid polluting the caller's namespace.
    - Error handling provides detailed messages if script loading or main() execution fails.

    Examples
    --------
    >>> handler = GridHandlerPythonScript(
    ...     script_path="path/to/cavern2D_mesh_generation.py",
    ...     parameters={"cavern_radius": 5.0}
    ... )
    >>> print(handler.mesh)  # Use like any GridHandlerGMSH
    """

    def __init__(self, script_path, parameters=None, function_name=None):
        """Initialize GridHandlerPythonScript.
        
        Args:
            script_path: Path to the mesh generation script.
            parameters: Dictionary of parameters to pass to main().
            function_name: Deprecated; ignored for backwards compatibility.
        """
        self.script_path = script_path
        self.parameters = parameters or {}
        self.comm = MPI.COMM_WORLD
        self.rank = self.comm.rank
        
        # Set grid_folder for compatibility with downstream code
        from pathlib import Path
        script_path_obj = Path(script_path)
        self.grid_folder = str(script_path_obj.parent)
        self.geometry_name = script_path_obj.stem

        # Generate mesh by executing the script
        mesh, cell_tags, facet_tags = self._generate_mesh_from_script()

        # Store DOLFINx objects (mimic GridHandlerGMSH interface)
        self.mesh = mesh
        self.subdomains = cell_tags
        self.boundaries = facet_tags

        # Extract metadata from mesh (compatible with GridHandlerGMSH)
        self.domain_dim = self.mesh.topology.dim
        self.boundary_dim = self.domain_dim - 1
        self.n_elems = self.mesh.topology.index_map(self.domain_dim).size_local + len(
            self.mesh.topology.index_map(self.domain_dim).ghosts
        )
        self.n_nodes = self.mesh.topology.index_map(0).size_local + len(
            self.mesh.topology.index_map(0).ghosts
        )

        # Build metadata structures (reuse GridHandlerGMSH logic)
        # Skip load_mesh() since we already have the mesh
        self._build_tags_from_mesh()
        self._load_boundaries_from_mesh()
        self._build_dolfin_tags()
        self.build_box_dimensions()
        # Override __extract_grid_data to use our already-extracted tags
        self._extract_grid_data_from_tags()
        self.build_smoother()

    def _build_tags_from_mesh(self):
        """
        Extract tags from already-loaded mesh objects with proper name mapping.
        Uses boundary group metadata from the mesh generation script if available.
        """
        # Try to get boundary group names from the script module that was just executed
        tag_value_to_name = {}
        try:
            import importlib
            import sys
            from pathlib import Path
            
            # Get the module name from the script path
            script_module_name = Path(self.script_path).stem
            
            # Check if the module is in sys.modules (it should be after execution)
            for module_name, module in sys.modules.items():
                if script_module_name in module_name or module_name.endswith(script_module_name):
                    if hasattr(module, "_MESH_METADATA"):
                        tag_value_to_name = module._MESH_METADATA
                        break
        except Exception as e:
            pass  # Fall back to generic names
        
        # Initialize tags dictionary from DOLFINx MeshTags
        self.tags = {}
        self.dolfin_tags = self.tags
        self.subdomain_tags = {}
        self.tags_dict = {}
        
        # Extract cell tags (3D)
        if self.subdomains is not None:
            cell_values = self.subdomains.values
            cell_indices = self.subdomains.indices
            
            for i, tag_value in enumerate(sorted(set(cell_values))):
                # Use boundary group name if available, otherwise generic name
                if tag_value in tag_value_to_name:
                    tag_name = tag_value_to_name[tag_value]
                else:
                    tag_name = f"region_{tag_value}"
                self.subdomain_tags[tag_name] = list(cell_indices[cell_values == tag_value])
                self.tags_dict[tag_value] = tag_name
        
        # Set up basic tag structure for compatibility
        if 3 not in self.tags:
            self.tags[3] = {}
        for tag_value, tag_name in self.tags_dict.items():
            self.tags[3][tag_name] = tag_value

    def _load_boundaries_from_mesh(self):
        """
        Extract boundaries from already-loaded mesh objects with proper name mapping.
        """
        # Try to get boundary group names from the script module
        tag_value_to_name = {}
        try:
            import sys
            from pathlib import Path
            
            script_module_name = Path(self.script_path).stem
            for module_name, module in sys.modules.items():
                if script_module_name in module_name or module_name.endswith(script_module_name):
                    if hasattr(module, "_MESH_METADATA"):
                        tag_value_to_name = module._MESH_METADATA
                        break
        except Exception:
            pass
        
        self.boundary_tags = {}
        
        # Extract facet tags (2D for 3D mesh)
        if self.boundaries is not None:
            facet_values = self.boundaries.values
            facet_indices = self.boundaries.indices
            
            for tag_value in set(facet_values):
                # Use boundary group name if available, otherwise generic name
                if tag_value in tag_value_to_name:
                    tag_name = tag_value_to_name[tag_value]
                else:
                    tag_name = f"boundary_{tag_value}"
                self.boundary_tags[tag_name] = list(facet_indices[facet_values == tag_value])
        
        # Set up region_indices and tags_dict for compatibility
        self.region_names = list(self.subdomain_tags.keys())
        self.n_regions = len(self.region_names)
        self.region_indices = self.subdomain_tags

    def _extract_grid_data_from_tags(self):
        """
        Build region indices for elements based on cell tags.
        This is a simplified version that uses tags we've already extracted,
        avoiding the need to call get_subdomain_tag() which requires dolfin_tags
        to be fully populated.
        """
        # Use the tags and subdomain_tags we've already extracted
        # region_indices is already set to subdomain_tags in _load_boundaries_from_mesh
        # Just ensure tags_dict is properly mapped
        pass  # Everything is already done in _build_tags_from_mesh and _load_boundaries_from_mesh

    def _build_dolfin_tags(self):
        """
        Build dolfin_tags dictionary from boundary_tags and subdomain_tags.
        This mimics the behavior of GridHandlerGMSH.get_tags_from_file().
        """
        self.dolfin_tags = {1: {}, 2: {}, 3: {}}
        
        # Add subdomain (volume) tags (dimension domain_dim, e.g., 3 for 3D)
        # Create reverse mapping: tag_name -> tag_value
        for tag_name, tag_value in self.tags_dict.items():
            # tags_dict maps: tag_value -> tag_name, so we need to reverse it
            self.dolfin_tags[self.domain_dim][tag_name] = tag_value
        
        # Add boundary tags (dimension boundary_dim, e.g., 2 for 3D boundaries)
        if self.boundaries is not None:
            facet_values = self.boundaries.values
            facet_indices = self.boundaries.indices
            
            for tag_name, facet_list in self.boundary_tags.items():
                # Find the tag value from the boundaries MeshTags
                for face_val in set(facet_values):
                    face_indices = facet_indices[facet_values == face_val]
                    # Check if any of our facets match this tag value
                    if len(facet_list) > 0 and face_val in facet_values:
                        # Check if facets with this value are in our list
                        if any(idx in facet_list for idx in face_indices[:min(10, len(face_indices))]):
                            self.dolfin_tags[self.boundary_dim][tag_name] = face_val
                            break
        self.tags_dict = {i: name for i, name in enumerate(self.region_names)}

    def _generate_mesh_from_script(self):
        """
        Load and execute the mesh generation script, returning DOLFINx objects.

        The script must define a `main()` function that returns (mesh, cell_tags, facet_tags).

        Returns
        -------
        tuple[mesh, cell_tags, facet_tags]
            DOLFINx Mesh, cell MeshTags, and facet MeshTags.

        Raises
        ------
        FileNotFoundError
            If the script file does not exist.
        AttributeError
            If the script does not define a main() function.
        RuntimeError
            If main() execution fails or returns invalid types.
        """
        import importlib.util
        from pathlib import Path

        script_path_obj = Path(self.script_path).resolve()

        # Check file exists
        if not script_path_obj.exists():
            raise FileNotFoundError(
                f"Mesh generation script not found: {script_path_obj}"
            )

        # Load the script as a module
        spec = importlib.util.spec_from_file_location(
            script_path_obj.stem, script_path_obj
        )
        if spec is None or spec.loader is None:
            raise RuntimeError(
                f"Failed to load script module from {script_path_obj}"
            )
        module = importlib.util.module_from_spec(spec)
        
        # Register module in sys.modules before executing so it can be found by the MeshGenerator class
        import sys
        sys.modules[spec.name] = module
        
        spec.loader.exec_module(module)

        # Get the main() function
        if not hasattr(module, "main"):
            available = [
                name
                for name in dir(module)
                if callable(getattr(module, name)) and not name.startswith("_")
            ]
            raise AttributeError(
                f"Script must define a main() function. Available callables: "
                f"{', '.join(available)}"
            )

        func = getattr(module, "main")

        # Execute main() with provided parameters
        # YAML-invoked scripts should not save mesh files to disk by default
        try:
            if self.rank == 0:
                print(
                    f"[GridHandlerPythonScript] Executing main({self.parameters})"
                )
            result = func(save_to_disk=False, **self.parameters)
        except Exception as e:
            raise RuntimeError(
                f"main() execution failed with error: {e}"
            ) from e

        # Validate return type
        if not isinstance(result, tuple) or len(result) != 3:
            raise RuntimeError(
                f"Expected main() to return (mesh, cell_tags, facet_tags) tuple, "
                f"got {type(result)}"
            )

        mesh, cell_tags, facet_tags = result

        # Validate types
        from dolfinx.mesh import Mesh, MeshTags

        if not isinstance(mesh, Mesh):
            raise RuntimeError(
                f"First return value must be dolfinx.mesh.Mesh, got {type(mesh)}"
            )
        if not isinstance(cell_tags, MeshTags):
            raise RuntimeError(
                f"Second return value must be dolfinx.mesh.MeshTags, got {type(cell_tags)}"
            )
        if not isinstance(facet_tags, MeshTags):
            raise RuntimeError(
                f"Third return value must be dolfinx.mesh.MeshTags, got {type(facet_tags)}"
            )

        if self.rank == 0:
            domain_dim = mesh.topology.dim
            n_cells = mesh.topology.index_map(domain_dim).size_local + len(mesh.topology.index_map(domain_dim).ghosts)
            n_vertices = mesh.topology.index_map(0).size_local + len(mesh.topology.index_map(0).ghosts)
            print(
                f"[GridHandlerPythonScript] Mesh generated: {n_cells} cells, "
                f"{n_vertices} vertices"
            )

        return mesh, cell_tags, facet_tags
