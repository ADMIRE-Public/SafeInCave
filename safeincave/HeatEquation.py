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
    from .HeatBC import BcHandler
    # try:
    #     from .HeatBC import BcHandler
    # except:
    #     from HeatBC import BcHandler

from abc import ABC, abstractmethod
import dolfinx as do
import ufl
from petsc4py import PETSc
import torch as to
# from .MaterialProps import Material
# from .Grid import GridHandlerGMSH
# from .Utils import numpy2torch, project
try:
    from .MaterialProps import Material
    from .Grid import GridHandlerGMSH
    from .Utils import numpy2torch, project
except:
    from MaterialProps import Material
    from Grid import GridHandlerGMSH
    from Utils import numpy2torch, project




class HeatDiffusion():
    """
    Transient heat-diffusion solver on a DOLFINx mesh.

    This class sets up the finite-element spaces, fields, measures, and linear
    system for the implicit θ-method (backward Euler when θ=1). It assembles
    the bilinear and linear forms, applies time-dependent boundary conditions,
    and solves for the temperature increment/field at each time step.

    Parameters
    ----------
    grid : GridHandlerGMSH
        Grid/mesh handler providing the DOLFINx mesh and mesh tags
        (subdomains/boundaries).

    Attributes
    ----------
    grid : GridHandlerGMSH
        Input grid handler.
    DG0_1 : dolfinx.fem.FunctionSpace
        Discontinuous Galerkin (degree 0) space for cell-wise material fields.
    V : dolfinx.fem.FunctionSpace
        Lagrange (degree 1) space for temperature.
    n_elems : int
        Number of (local + ghost) cells in the mesh.
    n_nodes : int
        Number of (local + ghost) nodes in the mesh.
    dt : dolfinx.fem.Constant
        Time-step parameter (as a scalar Constant on the mesh).
    dT : ufl.trialfunctions.TrialFunction
        Trial function in space `V` (temperature increment).
    T_ : ufl.testfunctions.TestFunction
        Test function in space `V`.
    ds, dx : ufl.Measure
        Boundary and volume measures with subdomain data.
    k, rho, cp : dolfinx.fem.Function
        Cell-wise material properties (DG0_1): conductivity, density, heat capacity.
    T_old, T, X : dolfinx.fem.Function
        Temperature at previous step, current temperature, and solver vector
        (all in `V`).
    mat : Material
        Material properties container (set via :meth:`set_material`).
    solver : petsc4py.PETSc.KSP
        Linear solver (set via :meth:`set_solver`).
    bc : BcHandler
        Boundary-condition handler (set via :meth:`set_boundary_conditions`).

    Notes
    -----
    - Voigt/tensor conventions are irrelevant here; only scalar temperature is solved.
    - The class expects that a solver (:class:`PETSc.KSP`), a material
      (:class:`Material`), and a BC handler (:class:`BcHandler`) are set
      before calling :meth:`solve`.
    """
    def __init__(self, grid: GridHandlerGMSH):
        self.grid = grid

        self.create_function_spaces()

        self.n_elems = self.DG0_1.dofmap.index_map.size_local + len(self.DG0_1.dofmap.index_map.ghosts)
        self.n_nodes = self.V.dofmap.index_map.size_local + len(self.V.dofmap.index_map.ghosts)
        self.dt = do.fem.Constant(self.grid.mesh, 1.0)

        self.create_trial_test_functions()
        self.create_ds_dx()
        self.create_fenicsx_fields()
        # self.create_pytorch_fields()

    def set_material(self, material: Material) -> None:
        """
        Attach material properties and initialize FE fields.

        Parameters
        ----------
        material : Material
            Container with per-element tensors/vectors (`k`, `density`, `cp`).

        Returns
        -------
        None

        Side Effects
        ------------
        Calls :meth:`initialize` to copy material arrays into FE functions.
        """
        self.mat = material
        self.initialize()

    def set_solver(self, solver: PETSc.KSP) -> None:
        """
        Set the PETSc linear solver.

        Parameters
        ----------
        solver : petsc4py.PETSc.KSP
            Configured Krylov solver (e.g., with preconditioner and tolerances).

        Returns
        -------
        None
        """
        self.solver = solver

    def set_boundary_conditions(self, bc: BcHandler) -> None:
        """
        Set the boundary-condition handler.

        Parameters
        ----------
        bc : BcHandler
            Boundary-condition manager providing Dirichlet, Neumann, and Robin
            contributions and update routines.

        Returns
        -------
        None
        """
        self.bc = bc

    def create_trial_test_functions(self) -> None:
        """
        Create trial and test functions in `V`.

        Returns
        -------
        None

        Side Effects
        ------------
        Defines :attr:`dT` and :attr:`T_` as UFL trial/test functions.
        """
        self.dT = ufl.TrialFunction(self.V)
        self.T_ = ufl.TestFunction(self.V)

    def create_function_spaces(self) -> None:
        """
        Build function spaces for materials (DG0) and temperature (P1).

        Returns
        -------
        None

        Side Effects
        ------------
        Defines :attr:`DG0_1` and :attr:`V`.
        """
        self.DG0_1 = do.fem.functionspace(self.grid.mesh, ("DG", 0))
        self.V = do.fem.functionspace(self.grid.mesh, ("Lagrange", 1))

    def create_ds_dx(self) -> None:
        """
        Create UFL measures with subdomain data for integration.

        Returns
        -------
        None

        Side Effects
        ------------
        Defines :attr:`ds` for boundary integrals and :attr:`dx` for domain integrals,
        using mesh tags from :attr:`grid`.
        """
        self.ds = ufl.Measure("ds", domain=self.grid.mesh, subdomain_data=self.grid.get_boundaries())
        self.dx = ufl.Measure("dx", domain=self.grid.mesh, subdomain_data=self.grid.get_subdomains())

    def create_fenicsx_fields(self) -> None:
        """
        Allocate DOLFINx Functions for material and temperature fields.

        Returns
        -------
        None

        Side Effects
        ------------
        Creates :attr:`k`, :attr:`rho`, :attr:`cp` in `DG0_1` and
        :attr:`T_old`, :attr:`T`, :attr:`X` in `V`.
        """
        self.k = do.fem.Function(self.DG0_1)
        self.rho = do.fem.Function(self.DG0_1)
        self.cp = do.fem.Function(self.DG0_1)
        self.T_old = do.fem.Function(self.V)
        self.T = do.fem.Function(self.V)
        self.X = do.fem.Function(self.V)

    def initialize(self) -> None:
        """
        Copy material arrays into FE functions (`k`, `rho`, `cp`).

        Returns
        -------
        None

        Side Effects
        ------------
        Writes into `self.k.x.array`, `self.rho.x.array`, and `self.cp.x.array`
        from the arrays stored in :attr:`mat`.
        """
        self.k.x.array[:] = self.mat.k
        self.rho.x.array[:] = self.mat.density
        self.cp.x.array[:] = self.mat.cp

    def split_solution(self) -> None:
        """
        Assign the solver vector `X` to the temperature field `T`.

        Returns
        -------
        None

        Notes
        -----
        This sets `self.T = self.X` (rebinds the reference), rather than copying
        underlying arrays. If you need a deep copy, copy `self.X.x.array` into
        `self.T.x.array`.
        """
        self.T = self.X

    def update_T_old(self) -> None:
        """
        Copy the current temperature into the `T_old` storage.

        Returns
        -------
        None

        Side Effects
        ------------
        Overwrites `self.T_old.x.array` with `self.T.x.array`.
        """
        self.T_old.x.array[:] = self.T.x.array

    def set_initial_T(self, T_field: to.Tensor) -> None:
        """
        Initialize both current and previous temperature fields.

        Parameters
        ----------
        T_field : torch.Tensor
            Temperature values for all nodal DOFs (shape `(n_nodes,)` or compatible).

        Returns
        -------
        None

        Side Effects
        ------------
        Writes into `self.T_old.x.array` and `self.T.x.array`.
        """
        self.T_old.x.array[:] = T_field
        self.T.x.array[:] = T_field

    def get_T_elems(self) -> None:
        """
        Project the nodal temperature to DG0 and return as a torch tensor.

        Returns
        -------
        torch.Tensor
            Cell-wise constant temperature values of length `n_elems`.

        Notes
        -----
        Uses :func:`Utils.project` to project `T` onto `DG0_1` and
        :func:`Utils.numpy2torch` to convert to a Torch tensor.
        """
        T_elems = project(self.T, self.DG0_1)
        return numpy2torch(T_elems.x.array)


    def solve(self, t: float, dt: float) -> None:
        """
        Assemble and solve one implicit time step.

        The method:
        1. Update boundary conditions at time ``t`` via :class:`BcHandler`.
        2. Set the time-step constant ``dt``.
        3. Assemble the bilinear form :math:`(\\rho c_p / \\Delta t) (dT, T_) + (k \\nabla dT, \\nabla T_)` plus Robin terms.
        4. Assemble the right-hand side :math:`(\\rho c_p / \\Delta t) (T_{old}, T_)` plus Neumann and Robin terms.
        5. Apply Dirichlet BCs, solve the linear system for :attr:`X`, and update :attr:`T` and :attr:`T_old`.

        Parameters
        ----------
        t : float
            Current time at which boundary conditions are evaluated.
        dt : float
            Time-step size to advance.

        Returns
        -------
        None

        Raises
        ------
        AttributeError
            If solver, BC handler, or material fields are not set.

        Notes
        -----
        - The linear system is solved via the user-provided :class:`PETSc.KSP`.
        - Dirichlet conditions are enforced strongly; Neumann/Robin are added
          to the forms via the handler lists.
        """
        # Update boundary conditions
        self.bc.update_bcs(t)

        self.dt.value = dt

        # Build bilinear form
        a = (self.rho*self.cp*self.dT*self.T_/self.dt + self.k*ufl.dot(ufl.grad(self.dT), ufl.grad(self.T_)))*self.dx
        a += sum(self.bc.robin_bcs_a)
        bilinear_form = do.fem.form(a)
        A = do.fem.petsc.assemble_matrix(bilinear_form, bcs=self.bc.dirichlet_bcs)
        A.assemble()

        # Build linear form
        L = (self.rho*self.cp*self.T_old*self.T_/self.dt)*self.dx + sum(self.bc.neumann_bcs) + sum(self.bc.robin_bcs_b)
        linear_form = do.fem.form(L)
        b = do.fem.petsc.assemble_vector(linear_form)
        do.fem.petsc.apply_lifting(b, [bilinear_form], [self.bc.dirichlet_bcs])
        b.ghostUpdate(addv=PETSc.InsertMode.ADD_VALUES, mode=PETSc.ScatterMode.REVERSE)
        do.fem.petsc.set_bc(b, self.bc.dirichlet_bcs)
        b.ghostUpdate(addv=PETSc.InsertMode.INSERT, mode=PETSc.ScatterMode.FORWARD)

        # Solve linear system
        self.solver.setOperators(A)
        self.solver.solve(b, self.X.x.petsc_vec)
        self.X.x.scatter_forward()
        self.split_solution()

        # Update old temperature field
        self.update_T_old()


