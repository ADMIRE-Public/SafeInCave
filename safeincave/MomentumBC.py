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
from __future__ import annotations
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from .MomentumEquation import LinearMomentum
    
from abc import ABC
import numpy as np
import dolfinx as do
import ufl

class GeneralBC(ABC):
	"""
    Base container for time-dependent boundary-condition data.

    Subclasses (e.g., :class:`DirichletBC`, :class:`NeumannBC`) set the
    concrete fields and `type` identifier.

    Attributes
    ----------
    boundary_name : str or None
        Name/label of the boundary as defined by mesh tags.
    type : str or None
        Boundary condition type identifier (e.g., ``"dirichlet"``, ``"neumann"``).
    values : list[float] or None
        Boundary values sampled at the times in ``time_values``.
    time_values : list[float] or None
        Monotonically increasing times associated with ``values``.
    """
	def __init__(self):
		self.boundary_name = None
		self.type = None
		self.values = None
		self.time_values = None


class DirichletBC(GeneralBC):
	"""
    Time-dependent Dirichlet (essential) boundary condition for one component.

    Parameters
    ----------
    boundary_name : str
        Named boundary in the mesh tags (handled by the grid object).
    component : int
        Component index of the vector field to constrain (e.g., 0=x, 1=y, 2=z).
    values : list[float]
        Prescribed values over time (interpolated with :func:`numpy.interp`).
    time_values : list[float]
        Times corresponding to ``values`` (must be monotone).

    Attributes
    ----------
    type : str
        Always ``"dirichlet"``.
    component : int
        Constrained component index.
    values, time_values : list[float]
        Stored time history for interpolation.
    boundary_name : str
        Stored boundary label.
    """
	def __init__(self, boundary_name: str, component: int, values: list, time_values: list):
		self.boundary_name = boundary_name
		self.type = "dirichlet"
		self.values = values
		self.time_values = time_values
		self.component = component

class NeumannBC(GeneralBC):
	"""
    Time-dependent Neumann (traction/pressure) boundary condition with hydrostatics.

    The applied boundary value is modeled as
    ``p(t) + ρ g (H - x[i])``, where ``p(t)`` is a time-dependent pressure,
    ``ρ`` is density, ``g`` is gravity (default negative for downward),
    ``H`` is a reference elevation, and ``x[i]`` is the spatial coordinate
    in the chosen direction.

    Parameters
    ----------
    boundary_name : str
        Named boundary in the mesh tags.
    direction : int
        Coordinate index used for hydrostatic variation (0=x, 1=y, 2=z).
    density : float
        Fluid/solid density ``ρ``.
    ref_pos : float
        Reference elevation ``H``.
    values : list[float]
        Time samples for the base pressure ``p(t)`` (before sign).
    time_values : list[float]
        Times corresponding to ``values`` (must be monotone).
    g : float, default=-9.81
        Gravitational acceleration (sign included).

    Attributes
    ----------
    type : str
        Always ``"neumann"``.
    direction : int
        Axis index for hydrostatic term.
    density : float
        Stored density.
    ref_pos : float
        Stored reference elevation.
    gravity : float
        Stored gravitational acceleration.
    boundary_name, values, time_values : as given
        Stored metadata and time history.
    """
	def __init__(self, boundary_name: str, direction: int, density: float, ref_pos: float, values: list, time_values: list, g=-9.81):
		self.boundary_name = boundary_name
		self.type = "neumann"
		self.values = values
		self.time_values = time_values
		self.direction = direction
		self.density = density
		self.ref_pos = ref_pos
		self.gravity = g


class BcHandler():
	"""
    Boundary-condition handler for a linear momentum problem.

    Stores user-defined BC objects, organizes them by type, and converts them
    into DOLFINx/UFL objects at a given time ``t`` for assembly.

    Parameters
    ----------
    equation : LinearMomentum
        Momentum equation object providing:
        - ``grid`` with mesh, boundary dimension, and tag queries,
        - ``get_uV()`` to access the vector function space (for Dirichlet DOFs),
        - ``normal`` (outward unit normal vector for boundary integrals),
        - ``ds`` boundary measure.

    Attributes
    ----------
    eq : LinearMomentum
        Stored equation reference.
    dirichlet_boundaries : list[DirichletBC]
        Registered Dirichlet BCs.
    neumann_boundaries : list[NeumannBC]
        Registered Neumann BCs.
    dirichlet_bcs : list[dolfinx.fem.DirichletBC] (set by :meth:`update_dirichlet`)
        DOLFINx Dirichlet BC objects for the current time.
    neumann_bcs : list[ufl.Form] (set by :meth:`update_neumann`)
        UFL surface terms to add to the right-hand side.
    x : ufl.core.expr.Expr
        Spatial coordinate vector ``x = SpatialCoordinate(mesh)`` used for hydrostatics.
    """
	def __init__(self, equation: LinearMomentum):
		self.eq = equation
		self.dirichlet_boundaries = []
		self.neumann_boundaries = []
		self.x = ufl.SpatialCoordinate(self.eq.grid.mesh)

	def reset_boundary_conditions(self) -> None:
		"""
        Clear all registered boundary conditions.

        Returns
        -------
        None
        """
		self.dirichlet_boundaries = []
		self.neumann_boundaries = []

	def add_boundary_condition(self, bc : GeneralBC) -> None:
		"""
        Register a boundary condition instance.

        Parameters
        ----------
        bc : GeneralBC
            One of :class:`DirichletBC` or :class:`NeumannBC`.

        Returns
        -------
        None

        Raises
        ------
        Exception
            If the boundary condition type is not supported.
        """
		if bc.type == "dirichlet":
			self.dirichlet_boundaries.append(bc)
		elif bc.type == "neumann":
			self.neumann_boundaries.append(bc)
		else:
			raise Exception(f"Boundary type {bc.type} not supported.")

	def update_dirichlet(self, t: float) -> None:
		"""
        Build Dirichlet BC objects at time ``t``.

        Parameters
        ----------
        t : float
            Current simulation time.

        Returns
        -------
        None

        Side Effects
        ------------
        Populates :attr:`dirichlet_bcs` with
        :class:`dolfinx.fem.DirichletBC` constructed by:
        - locating DOFs on each boundary for the target component, and
        - interpolating the prescribed value via :func:`numpy.interp`.
        """
		self.dirichlet_bcs = []
		for bc in self.dirichlet_boundaries:
			value = np.interp(t, bc.time_values, bc.values)
			dofs = do.fem.locate_dofs_topological(
				self.eq.get_uV().sub(bc.component),
				self.eq.grid.boundary_dim,
				self.eq.grid.get_boundary_tags(bc.boundary_name)
			)
			self.dirichlet_bcs.append(
				do.fem.dirichletbc(
					do.default_scalar_type(value),
					dofs,
					self.eq.get_uV().sub(bc.component)
				)
			)

	def update_neumann(self, t: float) -> None:
		"""
        Build Neumann contributions (tractions/pressures) at time ``t``.

        For each :class:`NeumannBC`, the boundary term is constructed as
        ``(p(t) + ρ g (H - x[i])) * n * ds(tag)``, where ``n`` is the outward
        unit normal, and ``ds(tag)`` is the boundary measure for the target
        boundary.

        Parameters
        ----------
        t : float
            Current simulation time.

        Returns
        -------
        None

        Side Effects
        ------------
        Populates :attr:`neumann_bcs` with UFL surface integrals to be added
        to the right-hand side form.
        """
		self.neumann_bcs = []
		for bc in self.neumann_boundaries:
			i = bc.direction
			rho = bc.density
			H = bc.ref_pos
			p = -np.interp(t, bc.time_values, bc.values)
			value_neumann = p + rho*bc.gravity*(H - self.x[i])
			self.neumann_bcs.append(value_neumann*self.eq.normal*self.eq.ds(self.eq.grid.get_boundary_tag(bc.boundary_name)))

