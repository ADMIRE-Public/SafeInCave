# Copyright (c) 2026, The SafeInCave Developers
#
# SPDX-License-Identifier: BSD-3-Clause

from __future__ import annotations
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from ...Cavern import CavernHandler

from ..base import GeneralBC
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
        # Constant-backed BC caches. The DOLFINx BC objects / UFL terms are
        # built ONCE with a persistent fem.Constant per time-varying load;
        # per-step updates only assign Constant.value. This keeps the
        # dirichlet_bcs/neumann_bcs/cavern_bcs list objects identity-stable,
        # which is the hard prerequisite for caching compiled forms that
        # reference them (a form compiled against a load baked in as a float
        # would freeze that load at its compile-time value).
        self._dirichlet_consts = None
        self._neumann_consts = None
        self._cavern_cache_key = None
        self._cavern_consts = None

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
        self._dirichlet_consts = None
        self._neumann_consts = None

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
            self._dirichlet_consts = None
        elif bc.type == "neumann":
            self.neumann_boundaries.append(bc)
            self._neumann_consts = None
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
        On first call (or after the registered BC set changes), builds
        :attr:`dirichlet_bcs` once: DOFs are located per boundary and each
        :class:`dolfinx.fem.DirichletBC` wraps a persistent
        :class:`dolfinx.fem.Constant`. Subsequent calls only assign
        ``Constant.value`` from :func:`numpy.interp`, keeping the BC objects
        (and the list) identity-stable for cached compiled forms.
        """
        if self._dirichlet_consts is None:
            mesh = self.uV.mesh
            self._dirichlet_consts = []
            self.dirichlet_bcs = []
            for bc in self.dirichlet_boundaries:
                const = do.fem.Constant(mesh, do.default_scalar_type(0.0))
                dofs = do.fem.locate_dofs_topological(
                    self.uV.sub(bc.component),
                    self.boundary_dim,
                    self.boundary_tags[bc.boundary_name],
                )
                self.dirichlet_bcs.append(
                    do.fem.dirichletbc(const, dofs, self.uV.sub(bc.component))
                )
                self._dirichlet_consts.append(const)
        for bc, const in zip(self.dirichlet_boundaries, self._dirichlet_consts):
            const.value = np.interp(t, bc.time_values, bc.values)

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
        On first call (or after the registered BC set changes), builds
        :attr:`neumann_bcs` once as UFL surface integrals whose time-varying
        pressure lives in a persistent :class:`dolfinx.fem.Constant`.
        Subsequent calls only assign ``Constant.value``, keeping the UFL
        terms (and the list) identity-stable for cached compiled forms.
        """
        if self._neumann_consts is None:
            mesh = self.uV.mesh
            self._neumann_consts = []
            self.neumann_bcs = []
            for bc in self.neumann_boundaries:
                i = bc.direction
                rho = bc.density
                H = bc.ref_pos
                p_const = do.fem.Constant(mesh, do.default_scalar_type(0.0))
                value_neumann = p_const + rho * bc.gravity * (H - self.x[i])
                self.neumann_bcs.append(
                    value_neumann
                    * self.normal
                    * self.ds(self.dolfin_tags[self.boundary_dim][bc.boundary_name])
                )
                self._neumann_consts.append(p_const)
        for bc, p_const in zip(self.neumann_boundaries, self._neumann_consts):
            p_const.value = -np.interp(t, bc.time_values, bc.values)

    def update_cavern_bcs(self, cavern_handler: CavernHandler):
        """
        Build/update cavern pressure loads. The UFL terms are built once per
        cavern set with persistent :class:`dolfinx.fem.Constant` objects for
        the two time-varying scalars (gauge pressure ``P`` and fluid density);
        per-call updates only assign ``Constant.value``. A change in the
        cavern set rebuilds the terms (invalidating any cached forms that
        reference :attr:`cavern_bcs`).
        """
        caverns = cavern_handler.caverns_PT + cavern_handler.caverns_MFlux
        cache_key = tuple(id(cavern) for cavern in caverns)
        if self._cavern_cache_key != cache_key:
            mesh = self.uV.mesh
            self._cavern_cache_key = cache_key
            self._cavern_consts = []
            self.cavern_bcs = []
            for cavern in caverns:
                p_const = do.fem.Constant(mesh, do.default_scalar_type(0.0))
                rho_const = do.fem.Constant(mesh, do.default_scalar_type(0.0))
                load = p_const + rho_const * cavern.gravity * (
                    cavern.ref_pos - self.x[cavern.direction]
                )
                self.cavern_bcs.append(
                    load
                    * self.normal
                    * self.ds(self.dolfin_tags[self.boundary_dim][cavern.cavern_name])
                )
                self._cavern_consts.append((p_const, rho_const))
        for cavern, (p_const, rho_const) in zip(caverns, self._cavern_consts):
            p_const.value = -cavern.P
            rho_const.value = cavern.density
