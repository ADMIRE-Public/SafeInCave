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

from safeincave.CavernBC import CavernHandler

if TYPE_CHECKING:
    from .MassEquation import MassDeformablePorousMedia
    
from abc import ABC
import numpy as np
import dolfinx as do
import ufl


class GeneralBC(ABC):
	"""
    Base container for boundary-condition data (time-dependent).

    Parameters
    ----------
    boundary_name : str
        Name/label of the boundary as defined in the mesh tags handled by
        the equation's grid (e.g., ``eq.grid.get_boundary_tags(boundary_name)``).
    values : list of float
        Values of the boundary condition at the times given by ``time_values``.
        Interpreted piecewise-linearly via :func:`numpy.interp`.
    time_values : list of float
        Monotonically increasing times associated with ``values``.

    Attributes
    ----------
    boundary_name : str
        Boundary label used to query mesh tags.
    values : list of float
        Time-sampled boundary values.
    time_values : list of float
        Sample times corresponding to ``values``.
    type : str or None
        Boundary condition type identifier (set by subclasses).
    """
	def __init__(self, boundary_name: str, values: list, time_values: list):
		self.boundary_name = boundary_name
		self.values = values
		self.time_values = time_values
		self.type = None


class DirichletBC(GeneralBC):
	"""
    Time-dependent Dirichlet boundary condition (essential BC).

    Parameters
    ----------
    boundary_name : str
        Named boundary in the mesh tags.
    values : list of float
        Prescribed values over time.
    time_values : list of float
        Times corresponding to ``values``.

    Attributes
    ----------
    type : str
        Always ``"dirichlet"``.
    """
	def __init__(self, boundary_name: str, values: list, time_values: list):
		super().__init__(boundary_name, values, time_values)
		self.type = "dirichlet"

class NeumannBC(GeneralBC):
	"""
    Time-dependent Neumann boundary condition (natural BC / flux).

    Parameters
    ----------
    boundary_name : str
        Named boundary in the mesh tags.
    values : list of float
        Flux/intensity values over time.
    time_values : list of float
        Times corresponding to ``values``.

    Attributes
    ----------
    type : str
        Always ``"neumann"``.
    """
	def __init__(self, boundary_name: str, values: list, time_values: list):
		super().__init__(boundary_name, values, time_values)
		self.type = "neumann"



class BcHandler():
    """
    Boundary-condition handler for mass balance equation.

    This class stores user-defined boundary conditions, organizes them by
    type (Dirichlet, Neumann), and converts them into DOLFINx/UFL
    objects at a given time ``t`` for assembly.

    Parameters
    ----------
    equation : MassDeformablePorousMedia
        Equation object providing:

        - ``V``: function space of the primary field.
        - ``grid``: mesh/grid handler with methods
          :meth:`get_boundary_tags` and :meth:`get_boundary_tag`.
        - Variational symbols, e.g., ``T_`` (test function), ``dT`` (trial
          increment), and the measure ``ds``.

    Attributes
    ----------
    eq : MassDeformablePorousMedia
        Stored equation reference.
    dirichlet_boundaries : list[DirichletBC]
        Registered Dirichlet BCs.
    neumann_boundaries : list[NeumannBC]
        Registered Neumann BCs.
    dirichlet_bcs : list[dolfinx.fem.DirichletBC]
        DOLFINx Dirichlet BC objects at the current time (set by :meth:`update_dirichlet`).
    neumann_bcs : list[ufl.Form]
        Neumann contributions to the linear form (set by :meth:`update_neumann`).

    Notes
    -----
    Time dependence is handled via piecewise-linear interpolation using
    :func:`numpy.interp` between ``values`` and ``time_values`` stored in each
    boundary-condition object.
    """
    def __init__(self, equation: MassDeformablePorousMedia):
        self.eq = equation
        self.dirichlet_boundaries = []
        self.neumann_boundaries = []

    def reset_boundary_conditions(self) -> None:
        """
        Clear all registered boundary conditions of all types.

        Returns
        -------
        None
        """
        self.dirichlet_boundaries = []
        self.neumann_boundaries = []

    def add_boundary_condition(self, bc : GeneralBC) -> None:
        """
        Register a boundary condition by its type.

        Parameters
        ----------
        bc : GeneralBC
            One of :class:`DirichletBC`, :class:`NeumannBC`, or :class:`RobinBC`.

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

    def update_bcs(self, t: float) -> None:
        """
        Update all boundary-condition objects for a given time.

        This builds the runtime lists:
        ``dirichlet_bcs``, ``neumann_bcs``.

        Parameters
        ----------
        t : float
            Current simulation time.

        Returns
        -------
        None
        """
        self.update_dirichlet(t)
        self.update_neumann(t)

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
        :class:`dolfinx.fem.DirichletBC` objects constructed from:
        - located dofs on the boundary,
        - interpolated Dirichlet value at time ``t``.
        """
        self.dirichlet_bcs = []
        for bc in self.dirichlet_boundaries:
            value = np.interp(t, bc.time_values, bc.values)
            dofs = do.fem.locate_dofs_topological(
                self.eq.V,
                self.eq.grid.boundary_dim,
                self.eq.grid.get_boundary_tags(bc.boundary_name)
            )
            self.dirichlet_bcs.append(
                do.fem.dirichletbc(
                    do.default_scalar_type(value),
                    dofs,
                    self.eq.V
                )
            )


    def update_neumann(self, t: float) -> None:
        """
        Build Neumann contributions to the linear form at time ``t``.

        Parameters
        ----------
        t : float
            Current simulation time.

        Returns
        -------
        None

        Side Effects
        ------------
        Populates :attr:`neumann_bcs` with UFL terms of the form
        ``value * eq.T_ * ds(tag)`` to be added to the right-hand side.
        """
        self.neumann_bcs = []
        for bc in self.neumann_boundaries:
            value = np.interp(t, bc.time_values, bc.values)
            self.neumann_bcs.append(value*self.eq.T_*self.eq.ds(self.eq.grid.get_boundary_tag(bc.boundary_name)))
	

    def update_cavern_bcs(self, cavern_handler: CavernHandler):
        self.cavern_bcs = []
        for cavern in cavern_handler.caverns_PT + cavern_handler.caverns_MFlux:
            value = cavern.P
            dofs = do.fem.locate_dofs_topological(
                self.eq.V,
                self.eq.grid.boundary_dim,
                self.eq.grid.get_boundary_tags(cavern.cavern_name)
            )
            self.cavern_bcs.append(
                do.fem.dirichletbc(
                    do.default_scalar_type(value),
                    dofs,
                    self.eq.V
                )
            )

