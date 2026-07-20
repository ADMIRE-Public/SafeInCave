# Copyright (c) 2026, The SafeInCave Developers
#
# SPDX-License-Identifier: BSD-3-Clause

import ufl
from dolfinx import fem
from mpi4py import MPI
import numpy as np


class CavernVolumeComputer:
    """
    Compute the volume of a cavern from a tagged boundary surface.

    The cavern wall is extracted from a 3D tetrahedral mesh, oriented with
    outward-facing triangle normals, and integrated as a closed triangulated
    surface. Optional displacement fields can be supplied to compute the
    deformed volume.

    Parameters
    ----------
    grid : GridHandlerGMSH-like
        Grid object exposing the mesh, boundary tags, and facet markers.
    boundary_name : str
        Name/label of the cavern-wall boundary in the mesh tags.
    reference_point : list of float, optional
        Point used as the origin for tetrahedral volume integration. If
        ``None``, the area-weighted surface centroid is used. Required when
        ``sym_scale`` is greater than 1.
    sym_scale : int, default=1
        Symmetry multiplier for meshes representing part of a cavern. Valid
        values are ``1`` (full cavern), ``2`` (half cavern), and ``4``
        (quarter cavern).

    Attributes
    ----------
    grid : GridHandlerGMSH-like
        Stored grid reference.
    boundary_name : str
        Boundary label used to extract the cavern wall.
    sym_scale : int
        Validated symmetry multiplier.
    reference_point : list of float or ndarray
        Point used for volume integration.
    coords_wall : ndarray
        Coordinates of the wall vertices.
    conn_wall : ndarray
        Triangle connectivity indexing into ``coords_wall``.
    ids_wall : ndarray
        Global vertex ids of the wall vertices.
    wall_cells : ndarray
        Incident cell id used to evaluate displacement at each wall vertex.
    """

    def __init__(
        self,
        grid,
        boundary_name: str,
        reference_point: list[float, float, float] = None,
        sym_scale: int = 1,
    ):
        self.grid = grid
        self.boundary_name = boundary_name
        self.sym_scale = self.validate_sym_scale(sym_scale)
        self.reference_point = reference_point
        self.verify_reference_point()
        cavern_data = self.__extract_cavern_surface_from_grid(self.boundary_name)
        self.coords_wall, self.conn_wall, self.ids_wall = cavern_data
        if reference_point is None:
            self.reference_point = self.__surface_centroid(
                self.coords_wall, self.conn_wall
            )
        self.gather_cells_for_wall_vertices()

    def verify_reference_point(self):
        if self.sym_scale != 1 and self.reference_point is None:
            raise ValueError(f"""
Reference point should be provided when using symmetry (sym_scale = {self.sym_scale}).
Additionally, it must be on the intersection of all symmetry planes.
                             """)

    def validate_sym_scale(self, sym_scale):
        """
        Validate that sym_scale is a positive integer and that the mesh represents a fraction 1/sym_scale of the full cavern.
        It can be either 1 (full cavern), 2 (half cavern), or 4 (quarter cavern).
        """
        if sym_scale not in [1, 2, 4]:
            raise ValueError(f"sym_scale must be 1, 2, or 4. Got {sym_scale}.")
        return sym_scale

    def calculate_normals(self):
        normals = np.zeros((len(self.conn_wall), 3))
        for i, tri in enumerate(self.conn_wall):
            p0, p1, p2 = self.coords_wall[tri]
            normal = np.cross(p1 - p0, p2 - p0)
            normals[i] = normal / np.linalg.norm(normal)
        return normals

    def calculate_cavern_midpoint(self, direction: int, u=None):
        if u is None:
            coords = self.coords_wall
        else:
            disp_wall = u.eval(self.coords_wall, self.wall_cells)
            coords = self.coords_wall + disp_wall
        min_coord = np.min(coords[:, direction])
        max_coord = np.max(coords[:, direction])
        return (min_coord + max_coord) / 2.0

    def gather_cells_for_wall_vertices(self):
        tdim = self.grid.mesh.topology.dim
        self.grid.mesh.topology.create_connectivity(0, tdim)
        v2c = self.grid.mesh.topology.connectivity(0, tdim)
        self.wall_cells = np.empty(len(self.ids_wall), dtype=np.int32)
        for k, v in enumerate(self.ids_wall):
            links = v2c.links(v)
            if len(links) == 0:
                raise RuntimeError(f"Vertex {v} has no incident cell.")
            self.wall_cells[k] = links[0]

    def compute(self, u=None):
        if u is None:
            coords = self.coords_wall
        else:
            disp_wall = u.eval(self.coords_wall, self.wall_cells)
            coords = self.coords_wall + disp_wall
        volume = 0.0
        for tri in self.conn_wall:
            v0 = coords[tri][0] - self.reference_point
            v1 = coords[tri][1] - self.reference_point
            v2 = coords[tri][2] - self.reference_point
            volume += np.dot(v0, np.cross(v1, v2)) / 6.0
        return abs(volume) * self.sym_scale

    def __extract_cavern_surface_from_grid(self, boundary_name: str):
        """
        Return (coords_wall, tris_local, wall_ids) for the named boundary.

        The returned triangles are oriented so that their normals point outward
        from the 3D domain (i.e. away from the adjacent tetrahedral cell).

        Parameters
        ----------
        grid : GridHandlerGMSH-like
            Object exposing:
            - grid.mesh
            - grid.get_boundary_tag(boundary_name)
            - grid.boundaries or grid.facet_tags
        boundary_name : str
            Name of the boundary to extract.

        Returns
        -------
        coords_wall : (n_wall, 3) ndarray
            Coordinates of wall vertices, in the same order as wall_ids.
        tris_local : (n_tris, 3) ndarray
            Triangle connectivity indexing into coords_wall, with consistent
            outward orientation.
        wall_ids : (n_wall,) ndarray
            Global vertex ids of the wall vertices.
        """
        mesh = self.grid.mesh
        tdim = mesh.topology.dim
        fdim = tdim - 1

        if tdim != 3:
            raise RuntimeError("This function currently assumes a 3D tetrahedral mesh.")

        tag = self.grid.get_boundary_tag(boundary_name)

        mt = getattr(self.grid, "boundaries", None) or getattr(
            self.grid, "facet_tags", None
        )
        if mt is None:
            raise RuntimeError(
                "Grid does not expose facet tags as 'boundaries' or 'facet_tags'."
            )

        # Facets carrying the requested boundary tag
        facets = mt.indices[mt.values == tag]
        if len(facets) == 0:
            raise RuntimeError(
                f"No facets found for boundary '{boundary_name}' (tag={tag})."
            )

        # Required connectivities
        mesh.topology.create_connectivity(fdim, 0)  # facet -> vertices
        mesh.topology.create_connectivity(fdim, tdim)  # facet -> cell(s)
        mesh.topology.create_connectivity(tdim, 0)  # cell  -> vertices

        f2v = mesh.topology.connectivity(fdim, 0)
        f2c = mesh.topology.connectivity(fdim, tdim)
        c2v = mesh.topology.connectivity(tdim, 0)

        coords_all = mesh.geometry.x

        tri_global = []
        wall_set = set()

        for f in facets:
            facet_vertices = np.array(f2v.links(f), dtype=np.int64)
            if len(facet_vertices) != 3:
                # Expected for tetrahedral boundary facets
                continue

            attached_cells = f2c.links(f)
            if len(attached_cells) != 1:
                raise RuntimeError(
                    f"Boundary facet {f} should have exactly one adjacent cell, got {len(attached_cells)}."
                )

            cell = attached_cells[0]
            cell_vertices = np.array(c2v.links(cell), dtype=np.int64)

            # The tetra vertex opposite to this facet
            opposite = np.setdiff1d(cell_vertices, facet_vertices)
            if len(opposite) != 1:
                raise RuntimeError(
                    f"Could not identify unique opposite vertex for facet {f} in cell {cell}."
                )
            opposite_vertex = opposite[0]

            v0, v1, v2 = facet_vertices
            p0 = coords_all[v0]
            p1 = coords_all[v1]
            p2 = coords_all[v2]
            p_opp = coords_all[opposite_vertex]

            # Normal from current ordering
            n = np.cross(p1 - p0, p2 - p0)

            # Vector from facet toward the cell interior (toward the opposite tetra vertex)
            to_interior = p_opp - p0

            # If normal points toward the interior, flip triangle orientation
            if np.dot(n, to_interior) > 0.0:
                facet_vertices = np.array([v0, v2, v1], dtype=np.int64)

            tri_global.append(facet_vertices)
            wall_set.update(facet_vertices.tolist())

        if not tri_global:
            raise RuntimeError(
                f"No triangular facets found for boundary '{boundary_name}' (tag={tag})."
            )

        wall_ids = np.array(sorted(wall_set), dtype=np.int64)
        gid2lid = {gid: i for i, gid in enumerate(wall_ids)}

        tris_local = np.array(
            [[gid2lid[v] for v in tri] for tri in tri_global], dtype=np.int32
        )
        coords_wall = coords_all[wall_ids]
        return coords_wall, tris_local, wall_ids

    def __surface_centroid(self, coordinates, triangles):
        """Calculate the centroid of a surface defined by triangles."""
        total_area = 0
        weighted_sum = np.zeros(3)
        for tri in triangles:
            p0, p1, p2 = coordinates[tri]
            center = (p0 + p1 + p2) / 3.0
            area = np.linalg.norm(np.cross(p1 - p0, p2 - p0)) / 2.0
            weighted_sum += center * area
            total_area += area
        return weighted_sum / total_area


class HeatFluxComputer:
    """
    Compute integrated heat flux through a named cavern boundary.

    Parameters
    ----------
    grid : GridHandlerGMSH-like
        Grid object exposing the mesh, boundary tags, and facet markers.
    boundary_name : str
        Name/label of the boundary where heat flux is integrated.

    Attributes
    ----------
    grid : GridHandlerGMSH-like
        Stored grid reference.
    boundary_name : str
        Boundary label used to query mesh tags.
    ds : ufl.Measure
        Exterior-facet measure with grid boundary markers.
    n : ufl.FacetNormal
        Outward unit normal on the mesh facets.
    """

    def __init__(self, grid, boundary_name: str):
        self.grid = grid
        self.boundary_name = boundary_name

        self.ds = ufl.Measure(
            "ds", domain=self.grid.mesh, subdomain_data=self.grid.get_boundaries()
        )
        self.n = ufl.FacetNormal(self.grid.mesh)

    def compute(self, dt: float, T: any = None, kappa: any = None) -> float:
        if T is None or kappa is None:
            print(T, kappa)
            return 0.0
        else:
            Q_form = fem.form(
                ufl.dot(kappa * ufl.grad(T), self.n)
                * self.ds(self.grid.get_boundary_tag(self.boundary_name))
            )
            Q = fem.assemble_scalar(Q_form)
            Q = -self.grid.mesh.comm.allreduce(Q, op=MPI.SUM)
            return dt * Q
