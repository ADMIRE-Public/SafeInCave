# Copyright (c) 2026, The SafeInCave Developers
#
# SPDX-License-Identifier: BSD-3-Clause

"""
Base class for equation solvers.

This module provides common functionality shared between different equation
types (heat diffusion, momentum, etc.). Both HeatDiffusion and LinearMomentum
import from here to reduce code duplication.
"""

from abc import ABC, abstractmethod
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from ..MaterialProps import Material
    from ..Grid import GridHandlerGMSH

import dolfinx as do
from petsc4py import PETSc


class EquationBase(ABC):
    """
    Abstract base class for PDE equation solvers.

    Provides common functionality for setting up function spaces, measures,
    solvers, and material properties. Subclasses implement equation-specific
    variational forms and solve routines.

    Parameters
    ----------
    grid : GridHandlerGMSH
        Grid/mesh handler providing the DOLFINx mesh and mesh tags.
    solver_name : str, default="cg"
        PETSc solver type (e.g., "cg", "gmres").
    preconditioner : str, default="asm"
        PETSc preconditioner type.
    rtol : float, default=1e-12
        Relative tolerance for the linear solver.
    max_it : int, default=100
        Maximum number of linear solver iterations.

    Attributes
    ----------
    grid : GridHandlerGMSH
        Input grid handler.
    solver_name : str
        PETSc solver type.
    preconditioner : str
        PETSc preconditioner type.
    rtol : float
        Relative tolerance for solver.
    max_it : int
        Maximum solver iterations.
    solver : petsc4py.PETSc.KSP
        Linear solver instance.
    mat : Material
        Material properties container.
    """

    def __init__(
        self,
        grid: "GridHandlerGMSH",
        solver_name: str = "cg",
        preconditioner: str = "asm",
        rtol: float = 1e-12,
        max_it: int = 100,
    ):
        self.grid = grid
        self.solver_name = solver_name
        self.preconditioner = preconditioner
        self.rtol = rtol
        self.max_it = max_it

        self.build_solver()
        self.create_function_spaces()
        self.create_ds_dx()

    def build_solver(self) -> None:
        """
        Create and configure the PETSc linear solver.

        Returns
        -------
        None

        Side Effects
        ------------
        Creates :attr:`solver` as a configured PETSc.KSP instance.
        """
        self.solver = PETSc.KSP().create(self.grid.mesh.comm)
        self.solver.setType(self.solver_name)
        self.solver.getPC().setType(self.preconditioner)
        self.solver.setTolerances(rtol=self.rtol, max_it=self.max_it)

    @abstractmethod
    def create_function_spaces(self) -> None:
        """
        Create finite element function spaces.

        This is equation-specific and must be implemented by subclasses.

        Returns
        -------
        None

        Side Effects
        ------------
        Defines equation-specific function space attributes.
        """
        pass

    def create_ds_dx(self) -> None:
        """
        Create UFL measures with subdomain data for integration.

        Returns
        -------
        None

        Side Effects
        ------------
        Defines :attr:`ds` for boundary integrals and :attr:`dx` for domain
        integrals, using mesh tags from :attr:`grid`.
        """
        import ufl
        self.ds = ufl.Measure(
            "ds", domain=self.grid.mesh, subdomain_data=self.grid.get_boundaries()
        )
        self.dx = ufl.Measure(
            "dx", domain=self.grid.mesh, subdomain_data=self.grid.get_subdomains()
        )

    def set_material(self, material: "Material") -> None:
        """
        Attach material properties and initialize FE fields.

        Parameters
        ----------
        material : Material
            Container with per-element material properties.

        Returns
        -------
        None

        Side Effects
        ------------
        Calls :meth:`initialize` to copy material arrays into FE functions.
        """
        self.mat = material
        self.initialize()

    @abstractmethod
    def initialize(self) -> None:
        """
        Copy material arrays into FE functions.

        This is equation-specific and must be implemented by subclasses.

        Returns
        -------
        None

        Side Effects
        ------------
        Writes material property arrays into FE function storage.
        """
        pass

    def set_solver(self, solver: PETSc.KSP) -> None:
        """
        Set the PETSc linear solver.

        Parameters
        ----------
        solver : petsc4py.PETSc.KSP
            Preconfigured Krylov solver.

        Returns
        -------
        None
        """
        self.solver = solver

    @abstractmethod
    def set_boundary_conditions(self, bc) -> None:
        """
        Set the boundary-condition handler.

        This is equation-specific and must be implemented by subclasses.

        Parameters
        ----------
        bc : BcHandler
            Boundary-condition manager.

        Returns
        -------
        None
        """
        pass
