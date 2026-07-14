# Copyright (c) 2026, The SafeInCave Developers
#
# SPDX-License-Identifier: BSD-3-Clause

from __future__ import annotations
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from .CavernBC import CavernHandler

from .BC.base import GeneralBC, DirichletBC, NeumannBC
from .BC.cavern_utils import compute_cavern_pressure_load
import numpy as np
import dolfinx as do
import ufl


class BcHandler:
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

    def __init__(self):
        self.dirichlet_boundaries = []
        self.neumann_boundaries = []

    def set_uV(self, uV):
        self.uV = uV

    def set_boundary_dim(self, boundary_dim):
        self.boundary_dim = boundary_dim

    def set_boudary_tags(self, boundary_tags):
        self.boundary_tags = boundary_tags

    def set_dolfin_tags(self, dolfin_tags):
        self.dolfin_tags = dolfin_tags

    def set_normal(self, normal):
        self.normal = normal

    def set_ds(self, ds):
        self.ds = ds

    def set_spatial_coordinates(self, mesh):
        self.x = ufl.SpatialCoordinate(mesh)

    def reset_boundary_conditions(self) -> None:
        """
        Clear all registered boundary conditions.

        Returns
        -------
        None
        """
        self.dirichlet_boundaries = []
        self.neumann_boundaries = []

    def add_boundary_condition(self, bc: GeneralBC) -> None:
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
                self.uV.sub(bc.component),
                self.boundary_dim,
                self.boundary_tags[bc.boundary_name],
            )
            self.dirichlet_bcs.append(
                do.fem.dirichletbc(
                    do.default_scalar_type(value), dofs, self.uV.sub(bc.component)
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
            value_neumann = p + rho * bc.gravity * (H - self.x[i])
            self.neumann_bcs.append(
                value_neumann
                * self.normal
                * self.ds(self.dolfin_tags[self.boundary_dim][bc.boundary_name])
            )

    def update_cavern_bcs(self, cavern_handler: CavernHandler):
        self.cavern_bcs = []
        for cavern in cavern_handler.caverns_PT + cavern_handler.caverns_MFlux:
            load = compute_cavern_pressure_load(cavern, self.x)
            self.cavern_bcs.append(
                load
                * self.normal
                * self.ds(self.dolfin_tags[self.boundary_dim][cavern.cavern_name])
            )
