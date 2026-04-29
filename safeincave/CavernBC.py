# Copyright 2025 The safeincave community.
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
from abc import ABC, abstractmethod
import CoolProp.CoolProp as CP
from .Thermodynamics import CavernThermodynamics
import ufl
from dolfinx import fem
from mpi4py import MPI
from .Utils import save_json
import numpy as np
import os


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


class Cavern(ABC):
    """
    Abstract base class for cavern boundary-condition models.

    Parameters
    ----------
    cavern_name : str
        Name/label of the cavern boundary in the mesh tags.
    fluid : str, optional
        CoolProp fluid name. If provided, it is validated against CoolProp's
        ``HEOS`` backend.
    h_conv : float, optional
        Convective heat-transfer coefficient used by the heat boundary handler.

    Attributes
    ----------
    cavern_name : str
        Boundary label associated with the cavern.
    fluid : str or None
        CoolProp fluid identifier.
    h_conv : float or None
        Convective heat-transfer coefficient.
    type : str or None
        Cavern type identifier set by subclasses.
    """

    def __init__(self, cavern_name: str, fluid: str = None, h_conv: float = None):
        self.cavern_name = cavern_name
        self.check_fluid(fluid)
        self.fluid = fluid
        self.h_conv = h_conv
        self.type = None

    def check_fluid(self, fluid: str) -> None:
        if fluid is not None:
            try:
                CP.AbstractState("HEOS", fluid)
            except ValueError:
                raise ValueError(f"Fluid '{fluid}' not recognized by CoolProp.")

    @abstractmethod
    def update_cavern(self, t: float) -> None:
        pass

    @abstractmethod
    def record_data(self, t: float) -> None:
        pass

    def calculate_initial_condition(self) -> None:
        pass

    def calculate_heat(self, T: any = None, kappa: any = None) -> None:
        pass


class CavernHandler:
    """
    Container and dispatcher for cavern boundary-condition models.

    This class stores caverns by model type and forwards volume, heat,
    time-update, recording, and export operations to the registered caverns.

    Attributes
    ----------
    caverns_T : list[Cavern_T]
        Caverns with prescribed time-dependent temperature.
    caverns_PT : list[Cavern_PT]
        Caverns with prescribed time-dependent pressure and temperature.
    caverns_MFlux : list[Cavern_MassFlux]
        Caverns advanced from a prescribed mass-flux history.
    output_folder : str or None
        Folder where cavern history data are written by :meth:`save_caverns_data`.
    """

    def __init__(self):
        self.caverns_T = []
        self.caverns_PT = []
        self.caverns_MFlux = []
        self.output_folder = None

    def add_cavern(self, cavern: Cavern) -> None:
        if cavern.type == "Cavern_T":
            self.caverns_T.append(cavern)
        elif cavern.type == "Cavern_PT":
            self.caverns_PT.append(cavern)
        elif cavern.type == "Cavern_MFlux":
            self.caverns_MFlux.append(cavern)
        else:
            raise ValueError(f"Cavern type {cavern.type} not supported")

    def set_output_folder(self, output_folder: str) -> None:
        self.output_folder = output_folder

    def calculate_volumes(self, u: any = None) -> None:
        for cavern in self.caverns_MFlux + self.caverns_PT:
            cavern.calculate_volume(u)

    def calculate_total_heat(self, dt: float, T: any = None, kappa: any = None) -> None:
        for cavern in self.caverns_MFlux + self.caverns_PT:
            cavern.calculate_heat(dt, T, kappa)

    def update_caverns(
        self, t: float, dt: float = None, u: any = None, T: any = None
    ) -> None:
        for cavern in self.caverns_T:
            cavern.update_cavern(t, dt)

        for cavern in self.caverns_PT:
            cavern.update_cavern(t, dt)

        for cavern in self.caverns_MFlux:
            cavern.update_cavern(t, dt)

    def calculate_initial_conditions(self) -> None:
        for cavern in self.caverns_MFlux:
            cavern.calculate_initial_condition()

    def record_cavern_data(self, t: float) -> None:
        for cavern in self.caverns_T:
            cavern.record_data(t)

        for cavern in self.caverns_PT:
            cavern.record_data(t)

        for cavern in self.caverns_MFlux:
            cavern.record_data(t)

    def save_caverns_data(self):
        cavern_data = {}
        for cavern in self.caverns_T + self.caverns_PT + self.caverns_MFlux:
            cavern_data[cavern.cavern_name] = {}
            cavern_data[cavern.cavern_name]["Type"] = cavern.type
            cavern_data[cavern.cavern_name]["Fluid"] = cavern.fluid
            cavern_data[cavern.cavern_name]["Data"] = cavern.create_data()
        if self.output_folder is not None:
            save_json(cavern_data, os.path.join(self.output_folder, "cavern_data.json"))


class Cavern_T(Cavern):
    """
    Cavern model with prescribed time-dependent temperature.

    Parameters
    ----------
    cavern_name : str
        Name/label of the cavern boundary in the mesh tags.
    T_values : list of float
        Prescribed cavern temperatures over time.
    time_values : list of float
        Times corresponding to ``T_values``.
    h_conv : float, optional
        Convective heat-transfer coefficient used by the heat boundary handler.

    Attributes
    ----------
    type : str
        Always ``"Cavern_T"``.
    T_values : list of float
        Time-sampled prescribed temperatures.
    time_values : list of float
        Sample times corresponding to ``T_values``.
    T : float
        Current cavern temperature.
    T_hist : list of float
        Recorded temperature history.
    t_hist : list of float
        Recorded time history.
    """

    def __init__(
        self, cavern_name: str, T_values: list, time_values: list, h_conv: float = None
    ):
        super().__init__(cavern_name, None, h_conv)
        self.type = "Cavern_T"
        self.T_values = T_values
        self.time_values = time_values
        self.T = self.T_values[0]

        # Initialize histories
        self.T_hist = [self.T]
        self.t_hist = [self.time_values[0]]

    def update_cavern(self, t: float, dt: float) -> None:
        self.T = np.interp(t, self.time_values, self.T_values)
        if self.T <= 0.0:
            raise ValueError(f"T must be > 0, got {self.T}")

    def record_data(self, t: float) -> None:
        self.T_hist.append(self.T)
        self.t_hist.append(t)

    def create_data(self, output_dir: str):
        pass


class Cavern_PT(Cavern):
    """
    Cavern model with prescribed time-dependent pressure and temperature.

    Pressure and temperature are interpolated from user-provided histories.
    The current density is evaluated with CoolProp using the absolute pressure
    ``P + P_atm`` and temperature ``T``.

    Parameters
    ----------
    grid : GridHandlerGMSH-like
        Grid object used for heat-flux and volume computations.
    cavern_name : str
        Name/label of the cavern boundary in the mesh tags.
    fluid : str
        CoolProp fluid name.
    sym_scale : int
        Symmetry multiplier for volume computation. Valid values are ``1``,
        ``2``, and ``4``.
    reference_point : list of float, optional
        Point used by :class:`CavernVolumeComputer` for volume integration.
    P_values : list of float
        Prescribed gauge pressure values over time, in Pa.
    T_values : list of float
        Prescribed temperature values over time, in K.
    time_values : list of float
        Times corresponding to ``P_values`` and ``T_values``.
    ref_pos : float, default=0.0
        Reference coordinate used by coupled cavern models.
    direction : int, default=0
        Coordinate direction associated with ``ref_pos`` and gravity.
    g : float, default=-9.81
        Gravitational acceleration.
    h_conv : float, optional
        Convective heat-transfer coefficient used by the heat boundary handler.
    P_atm : float, default=101325.0
        Atmospheric pressure added to gauge pressure before CoolProp calls.

    Attributes
    ----------
    type : str
        Always ``"Cavern_PT"``.
    P : float
        Current gauge pressure.
    T : float
        Current temperature.
    density : float
        Current fluid density.
    Q : float
        Integrated heat over the current time step.
    heat : HeatFluxComputer
        Helper used to integrate heat flux on the cavern boundary.
    cvc : CavernVolumeComputer
        Helper used to compute cavern volume.
    V_hist, M_hist, P_hist, T_hist, Q_hist, density_hist, t_hist : list
        Recorded histories for volume, mass, pressure, temperature, heat,
        density, and time.
    """

    def __init__(
        self,
        *,
        grid: any,
        cavern_name: str,
        fluid: str,
        sym_scale: int,
        reference_point: list[float, float, float] = None,
        P_values: list,  # Gauge pressure values (Pa)
        T_values: list,  # Temperature values (K)
        time_values: list,
        ref_pos: float = 0.0,
        direction: int = 0,
        g: float = -9.81,
        h_conv: float = None,
        P_atm: float = 101325.0,  # Atmospheric pressure in Pa
    ):
        super().__init__(cavern_name, fluid, h_conv)
        self.type = "Cavern_PT"
        self.P_values = P_values
        self.T_values = T_values
        self.time_values = time_values
        self.P_atm = P_atm
        self.ref_pos = ref_pos
        self.direction = direction
        self.gravity = g
        self.Q = 0.0

        self.AS = CP.AbstractState("HEOS", self.fluid)
        self.P = self.P_values[0]
        self.T = self.T_values[0]
        self.AS.update(CP.PT_INPUTS, self.P + self.P_atm, self.T)
        self.density = self.AS.rhomass()

        # Initialize heat flux computer
        self.heat = HeatFluxComputer(grid=grid, boundary_name=self.cavern_name)

        # Initialize cavern volume computer
        self.cvc = CavernVolumeComputer(
            grid=grid,
            boundary_name=self.cavern_name,
            reference_point=reference_point,
            sym_scale=sym_scale,
        )

        # Initialize histories
        self.V_hist = []
        self.M_hist = []
        self.P_hist = []
        self.T_hist = []
        self.Q_hist = []
        self.density_hist = []
        self.t_hist = []

    def calculate_volume(self, u: any = None) -> None:
        if u is None:
            self.V = self.cvc.compute()
        else:
            self.V = self.cvc.compute(u)

    def calculate_heat(self, dt: float, T: any = None, kappa: any = None) -> None:
        self.Q = self.heat.compute(dt, T, kappa)

    def update_cavern(self, t: float, dt: float) -> None:
        self.P = np.interp(t, self.time_values, self.P_values)
        self.T = np.interp(t, self.time_values, self.T_values)
        if self.P + self.P_atm <= 0.0:
            raise ValueError(
                f"Absolute pressure must be > 0, got {self.P + self.P_atm}"
            )
        if self.T <= 0.0:
            raise ValueError(f"Temperature must be > 0, got {self.T}")
        self.AS.update(CP.PT_INPUTS, self.P + self.P_atm, self.T)
        self.density = self.AS.rhomass()

    def record_data(self, t: float) -> None:
        self.density_hist.append(self.density)
        self.V_hist.append(self.V)
        self.M_hist.append(self.density * self.V)
        self.P_hist.append(self.P)
        self.T_hist.append(self.T)
        self.Q_hist.append(self.Q)
        self.t_hist.append(t)

    def create_data(self):
        data = {}
        data["Time"] = self.t_hist
        data["Pressure"] = self.P_hist
        data["Temperature"] = self.T_hist
        data["Density"] = self.density_hist
        data["Volume"] = self.V_hist
        data["Mass"] = self.M_hist
        data["Heat"] = self.Q_hist
        return data


class Cavern_MassFlux(Cavern):
    """
    Cavern model advanced from a prescribed mass-flux history.

    The cavern mass is updated from the interpolated mass flux and time-step
    size. Pressure, temperature, and density are then solved with
    :class:`CavernThermodynamics`, including heat exchange and volume change.

    Parameters
    ----------
    grid : GridHandlerGMSH-like
        Grid object used for heat-flux and volume computations.
    cavern_name : str
        Name/label of the cavern boundary in the mesh tags.
    sym_scale : int
        Symmetry multiplier for volume computation. Valid values are ``1``,
        ``2``, and ``4``.
    reference_point : list of float, optional
        Point used by :class:`CavernVolumeComputer` for volume integration.
    fluid : str
        CoolProp fluid name.
    P_init : float
        Initial gauge pressure, in Pa.
    T_init : float
        Initial fluid temperature, in K.
    T_in : float
        Temperature of injected fluid, in K.
    Mflux_values : list of float
        Prescribed mass-flux values over time.
    time_values : list of float
        Times corresponding to ``Mflux_values``.
    direction : int, default=2
        Coordinate direction used to calculate the cavern midpoint.
    g : float, default=-9.81
        Gravitational acceleration.
    h_conv : float, optional
        Convective heat-transfer coefficient used by the heat boundary handler.
    P_atm : float, default=101325.0
        Atmospheric pressure used to convert gauge pressure to absolute pressure.

    Attributes
    ----------
    type : str
        Always ``"Cavern_MFlux"``.
    Mflux : float
        Current mass-flux value.
    P : float
        Current gauge pressure.
    T : float
        Current temperature.
    density : float
        Current fluid density.
    Q : float
        Integrated heat over the current time step.
    heat : HeatFluxComputer
        Helper used to integrate heat flux on the cavern boundary.
    cvc : CavernVolumeComputer
        Helper used to compute cavern volume.
    model : CavernThermodynamics
        Thermodynamic model used to update pressure, temperature, and density.
    V_hist, M_hist, P_hist, T_hist, Q_hist, density_hist, t_hist : list
        Recorded histories for volume, mass, pressure, temperature, heat,
        density, and time.
    """

    def __init__(
        self,
        *,
        grid: any,
        cavern_name: str,
        sym_scale: int,
        reference_point: list[float, float, float] = None,
        fluid: str,
        P_init: float,  # Initial gauge pressure (Pa)
        T_init: float,  # Initial fluid temperature (K)
        T_in: float,  # Temperature of injected fluid (K)
        Mflux_values: list,
        time_values: list,
        direction: int = 2,
        g: float = -9.81,
        h_conv: float = None,
        P_atm: float = 101325.0,  # Atmospheric pressure in Pa
    ):
        super().__init__(cavern_name, fluid, h_conv)
        self.type = "Cavern_MFlux"
        self.Mflux_values = Mflux_values
        self.time_values = time_values
        self.Mflux = self.Mflux_values[0]
        self.direction = direction
        self.gravity = g
        self.P_atm = P_atm
        self.P_init = P_init
        self.T_init = T_init
        self.T_in = T_in
        self.Q = 0.0

        # Initial gauge pressure and temperature
        self.P = P_init
        self.T = T_init
        self.P0 = P_init
        self.T0 = T_init

        # Initialize heat flux computer
        self.heat = HeatFluxComputer(grid=grid, boundary_name=self.cavern_name)

        # Initialize cavern volume computer
        self.cvc = CavernVolumeComputer(
            grid=grid,
            boundary_name=self.cavern_name,
            reference_point=reference_point,
            sym_scale=sym_scale,
        )
        self.ref_pos = self.cvc.calculate_cavern_midpoint(direction=self.direction)

        # Initialize thermodynamic model
        self.model = CavernThermodynamics(self.fluid)

        # Calculate initial density
        AS = CP.AbstractState("HEOS", self.fluid)
        AS.update(CP.PT_INPUTS, self.P + self.P_atm, self.T)
        self.density = AS.rhomass()

        # Initialize histories
        self.V_hist = []
        self.M_hist = []
        self.P_hist = []
        self.T_hist = []
        self.Q_hist = []
        self.density_hist = []
        self.t_hist = []

    def calculate_volume(self, u: any = None) -> None:
        if u is None:
            self.V = self.cvc.compute()
            self.ref_pos = self.cvc.calculate_cavern_midpoint(direction=self.direction)
        else:
            self.V = self.cvc.compute(u)
            self.ref_pos = self.cvc.calculate_cavern_midpoint(
                direction=self.direction, u=u
            )

    def calculate_heat(self, dt: float, T: any = None, kappa: any = None) -> None:
        self.Q = self.heat.compute(dt, T, kappa)
        # self.Q = 0.0

    def update_cavern(self, t: float, dt: float) -> None:
        Mflux = np.interp(t, self.time_values, self.Mflux_values)
        self.M = self.M0 + Mflux * dt
        dm = self.M - self.M0
        self.P, self.T, self.density = self.model.solve(
            dm=dm,
            Q_in=self.Q,
            T_in=self.T_in,
            P0=self.P0 + self.P_atm,  # Convert to absolute pressure
            T0=self.T0,
            V0=self.V0,
            V1=self.V,
        )
        if self.P <= 0.0:
            raise ValueError(f"P must be > 0, got {self.P}")
        if self.T <= 0.0:
            raise ValueError(f"T must be > 0, got {self.T}")

        # Convert to gauge pressure
        self.P -= self.P_atm
        # print(self.P/1e6, self.T, self.density, self.M0, dm, self.V)

    def record_data(self, t: float) -> None:
        self.density_hist.append(self.density)
        self.V_hist.append(self.V)
        self.M_hist.append(self.V * self.density)
        self.P_hist.append(self.P)
        self.T_hist.append(self.T)
        self.Q_hist.append(self.Q)
        self.t_hist.append(t)

    def calculate_initial_condition(self):
        self.V0 = self.V
        self.P0 = self.P
        self.T0 = self.T
        self.M0 = self.V * self.density

    def create_data(self):
        data = {}
        data["Time"] = self.t_hist
        data["Pressure"] = self.P_hist
        data["Temperature"] = self.T_hist
        data["Density"] = self.density_hist
        data["Volume"] = self.V_hist
        data["Mass"] = self.M_hist
        data["Heat"] = self.Q_hist
        return data
