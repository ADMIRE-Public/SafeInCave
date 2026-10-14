# Copyright (c) 2026, The SafeInCave Developers
#
# SPDX-License-Identifier: BSD-3-Clause

from __future__ import annotations
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from ...Cavern import CavernHandler

from ..base import GeneralBC
from ..cavern_utils import compute_cavern_pressure_load
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
        self.mesh = mesh
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

    # The BC objects and UFL terms below are built once and their time-varying
    # magnitudes live in dolfinx Constants updated in place. Rebuilding them
    # with the value baked in as a Python float changed every form's
    # signature each step, so each time step re-ran the FFCx JIT (a fresh
    # compilation on a cold cache) and solvers caching compiled forms on term
    # identity (LinearMomentumNewton._bc_forms_key) rebuilt them every step.
    # They are rebuilt only when the registered boundaries/caverns change.

    def update_dirichlet(self, t: float) -> None:
        """
        Set Dirichlet BC values at time ``t``.

        Side Effects
        ------------
        :attr:`dirichlet_bcs` holds one :class:`dolfinx.fem.DirichletBC` per
        registered Dirichlet boundary (DOFs located once per boundary set),
        whose value is the prescribed history interpolated at ``t``.
        """
        key = tuple(map(id, self.dirichlet_boundaries))
        if getattr(self, "_dirichlet_key", None) != key:
            self._dirichlet_values = []
            self.dirichlet_bcs = []
            for bc in self.dirichlet_boundaries:
                V = self.uV.sub(bc.component)
                dofs = do.fem.locate_dofs_topological(
                    V, self.boundary_dim, self.boundary_tags[bc.boundary_name]
                )
                value = do.fem.Constant(self.mesh, do.default_scalar_type(0.0))
                self._dirichlet_values.append(value)
                self.dirichlet_bcs.append(do.fem.dirichletbc(value, dofs, V))
            self._dirichlet_key = key
        for bc, value in zip(self.dirichlet_boundaries, self._dirichlet_values):
            value.value = np.interp(t, bc.time_values, bc.values)

    def update_neumann(self, t: float) -> None:
        """
        Set Neumann contributions (tractions/pressures) at time ``t``.

        Each :class:`NeumannBC` contributes ``(p(t) + ρ g (H - x[i])) * n *
        ds(tag)``, with ``p(t)`` the (sign-flipped) interpolated history held
        in a Constant.

        Side Effects
        ------------
        :attr:`neumann_bcs` holds the UFL surface integrals to add to the
        right-hand side form.
        """
        key = tuple(map(id, self.neumann_boundaries))
        if getattr(self, "_neumann_key", None) != key:
            self._neumann_values = []
            self.neumann_bcs = []
            for bc in self.neumann_boundaries:
                p = do.fem.Constant(self.mesh, do.default_scalar_type(0.0))
                value_neumann = p + bc.density * bc.gravity * (bc.ref_pos - self.x[bc.direction])
                self._neumann_values.append(p)
                self.neumann_bcs.append(
                    value_neumann
                    * self.normal
                    * self.ds(self.dolfin_tags[self.boundary_dim][bc.boundary_name])
                )
            self._neumann_key = key
        for bc, p in zip(self.neumann_boundaries, self._neumann_values):
            p.value = -np.interp(t, bc.time_values, bc.values)

    def update_cavern_bcs(self, cavern_handler: CavernHandler):
        caverns = cavern_handler.caverns_PT + cavern_handler.caverns_MFlux
        key = tuple(map(id, caverns))
        if getattr(self, "_cavern_key", None) != key:
            self._cavern_pressures = []
            self.cavern_bcs = []
            for cavern in caverns:
                P = do.fem.Constant(self.mesh, do.default_scalar_type(0.0))
                load = compute_cavern_pressure_load(cavern, self.x, pressure=P)
                self._cavern_pressures.append(P)
                self.cavern_bcs.append(
                    load
                    * self.normal
                    * self.ds(self.dolfin_tags[self.boundary_dim][cavern.cavern_name])
                )
            self._cavern_key = key
        for cavern, P in zip(caverns, self._cavern_pressures):
            P.value = cavern.P
