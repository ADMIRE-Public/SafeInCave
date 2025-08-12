"""
Discretization of the momentum balance equations
"""
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
    from MomentumBC import BcHandler

from abc import ABC, abstractmethod
import dolfinx as do
from dolfinx.fem import petsc as fem_petsc
import basix
import ufl
from petsc4py import PETSc
import torch as to
from MaterialProps import Material
from Grid import GridHandlerGMSH
import Utils as ut

class LinearMomentumBase(ABC):
    """
    Abstract base for a thermo-(visco)elastic linear momentum solver on a DOLFINx mesh.

    Sets up common FE spaces, measures, and fields; provides utilities for
    assembling body forces, computing invariants, and coordinating inelastic
    elements via the `Material` container. Concrete subclasses supply the
    variational forms and solve routines.

    Parameters
    ----------
    grid : GridHandlerGMSH
        Mesh/grid handler that provides the DOLFINx mesh and meshtags.
    theta : float
        Time integration parameter: 0 for fully implicit, 0.5 for Crank-Nicolson, 1 for explicit.

    Attributes
    ----------
    grid : GridHandlerGMSH
        Input grid handler.
    theta : float
        Time integration parameter.
    DG0_1, CG1_1 : dolfinx.fem.FunctionSpace
        Scalar DG0 (per element) and CG1 (per node) spaces.
    CG1_3x1 : dolfinx.fem.FunctionSpace
        Vector CG1 space of size equal to the spatial dimension.
    DG0_3x3, DG0_6x6 : dolfinx.fem.FunctionSpace
        Tensor DG0 spaces for 3×3 and 6×6 (Voigt) fields.
    n_elems : int
        Number of local+ghost elements.
    n_nodes : int
        Number of local+ghost nodes.
    u : dolfinx.fem.Function
        Displacement field (vector).
    sig, eps_tot : dolfinx.fem.Function
        Stress and total strain (DG0 3×3 tensors).
    q_nodes, q_elems, p_nodes, p_elems : dolfinx.fem.Function
        Von Mises magnitude and pressure in node/element spaces.
    Temp, T0 : torch.Tensor
        Current and reference temperatures per element, shape ``(n_elems,)``.
    normal : ufl.core.expr.Expr
        Test-function-weighted outward normal used for Neumann terms.
    ds, dx : ufl.Measure
        Boundary and domain measures with subdomain data.
    X : dolfinx.fem.Function
        Solution vector (same space as :meth:`get_uV()`).
    mat : Material
        Material container (set via :meth:`set_material`).
    solver : petsc4py.PETSc.KSP
        PETSc linear solver (set via :meth:`set_solver`).
    bc : BcHandler
        Boundary-condition handler (set via :meth:`set_boundary_conditions`).
    """
    def __init__(self, grid: GridHandlerGMSH, theta: float):
        self.grid = grid
        self.theta = theta

        self.create_function_spaces()
        self.create_ds_dx()

        self.n_elems = self.DG0_1.dofmap.index_map.size_local + len(self.DG0_1.dofmap.index_map.ghosts)
        self.n_nodes = self.CG1_1.dofmap.index_map.size_local + len(self.CG1_1.dofmap.index_map.ghosts)

        self.commom_fields()
        self.create_fenicsx_fields()
        self.create_pytorch_fields()

    def commom_fields(self) -> None:
        """
        Allocate common storage for temperature, stresses, and strains.

        Returns
        -------
        None

        Side Effects
        ------------
        Initializes tensors/functions:
        `T0`, `Temp`, `sig`, `eps_tot`, `u`, `q_elems`, `q_nodes`, `p_elems`, `p_nodes`.
        """
        self.T0 = to.zeros(self.n_elems, dtype=to.float64)
        self.Temp = to.zeros(self.n_elems, dtype=to.float64)
        self.sig = do.fem.Function(self.DG0_3x3)
        self.eps_tot = do.fem.Function(self.DG0_3x3)
        self.u = do.fem.Function(self.CG1_3x1)
        self.q_elems = do.fem.Function(self.DG0_1)
        self.q_nodes = do.fem.Function(self.CG1_1)
        self.p_elems = do.fem.Function(self.DG0_1)
        self.p_nodes = do.fem.Function(self.CG1_1)

    def set_material(self, material: Material) -> None:
        """
        Attach a material model and initialize FE fields from it.

        Parameters
        ----------
        material : Material
            Material container with elastic and non-elastic elements.

        Returns
        -------
        None

        Side Effects
        ------------
        Calls :meth:`initialize`.
        """
        self.mat = material
        self.initialize()

    def set_T(self, T: to.Tensor) -> None:
        """
        Set the current element-wise temperature.

        Parameters
        ----------
        T : torch.Tensor
            Temperature per element, shape ``(n_elems,)``.

        Returns
        -------
        None
        """
        self.Temp = T

    def set_T0(self, T0: to.Tensor) -> None:
        """
        Set the reference element-wise temperature.

        Parameters
        ----------
        T0 : torch.Tensor
            Reference temperature per element, shape ``(n_elems,)``.

        Returns
        -------
        None
        """
        self.T0 = T0

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

    def set_boundary_conditions(self, bc: BcHandler) -> None:
        """
        Set the boundary-condition handler.

        Parameters
        ----------
        bc : BcHandler
            Handler providing Dirichlet/Neumann terms and updates.

        Returns
        -------
        None
        """
        self.bc = bc

    def create_function_spaces(self) -> None:
        """
        Create function spaces used by the formulation.

        Returns
        -------
        None

        Side Effects
        ------------
        Defines :attr:`CG1_3x1`, :attr:`DG0_1`, :attr:`CG1_1`,
        :attr:`DG0_3x3`, and :attr:`DG0_6x6`.
        """
        self.CG1_3x1 = do.fem.functionspace(self.grid.mesh, ("Lagrange", 1, (self.grid.domain_dim, )))
        self.DG0_1 = do.fem.functionspace(self.grid.mesh, ("DG", 0))
        self.CG1_1 = do.fem.functionspace(self.grid.mesh, ("Lagrange", 1))
        self.DG0_3x3 = do.fem.functionspace(self.grid.mesh, ("DG", 0, (3, 3)))
        self.DG0_6x6 = do.fem.functionspace(self.grid.mesh, ("DG", 0, (6, 6)))

    def create_ds_dx(self) -> None:
        """
        Create boundary and domain measures with subdomain data.

        Returns
        -------
        None

        Side Effects
        ------------
        Defines :attr:`ds` and :attr:`dx` from grid meshtags.
        """
        self.ds = ufl.Measure("ds", domain=self.grid.mesh, subdomain_data=self.grid.get_boundaries())
        self.dx = ufl.Measure("dx", domain=self.grid.mesh, subdomain_data=self.grid.get_subdomains())

    def create_normal(self) -> None:
        """
        Create a test-function-weighted outward normal for surface terms.

        Returns
        -------
        None

        Side Effects
        ------------
        Defines :attr:`normal` as ``dot(FacetNormal(mesh), self.u_)``.
        """
        n = ufl.FacetNormal(self.grid.mesh)
        self.normal = ufl.dot(n, self.u_)

    def build_body_force(self, g: list) -> None:
        """
        Build the body-force linear form ``∫ ρ g · u_ dx``.

        Parameters
        ----------
        g : list of float
            Gravity/body acceleration vector components.

        Returns
        -------
        None

        Side Effects
        ------------
        Defines :attr:`b_body` as a UFL form for the right-hand side.
        """
        density = do.fem.Function(self.DG0_1)
        density.x.array[:] = self.mat.density
        body_force = density*do.fem.Constant(self.grid.mesh, do.default_scalar_type(tuple(g)))
        self.b_body = ufl.dot(body_force, self.u_)*self.dx

	# def compute_q_nodes(self) -> do.fem.Function:
	# 	dev = self.sig - (1/3)*ufl.tr(self.sig)*ufl.Identity(3)
	# 	q_form = ufl.sqrt((3/2)*ufl.inner(dev, dev))
	# 	self.q_nodes = ut.project(q_form, self.CG1_1)

	# def compute_q_elems(self) -> do.fem.Function:
	# 	dev = self.sig - (1/3)*ufl.tr(self.sig)*ufl.Identity(3)
	# 	q_form = ufl.sqrt((3/2)*ufl.inner(dev, dev))
	# 	self.q_elems = ut.project(q_form, self.DG0_1)

    def compute_q_nodes(self) -> None:
        """
        Compute von Mises equivalent stress and smooth it to nodes.

        Returns
        -------
        None

        Side Effects
        ------------
        Sets :attr:`q_nodes` by applying a node-element averaging matrix
        (:attr:`grid.A_csr`) to the element-wise von Mises values.
        """
        stress = ut.numpy2torch(self.sig.x.array.reshape((self.n_elems, 3, 3)))
        I1 = stress[:,0,0] + stress[:,1,1] + stress[:,2,2]
        I2 = stress[:,0,0]*stress[:,1,1] + stress[:,1,1]*stress[:,2,2] + stress[:,0,0]*stress[:,2,2] - stress[:,0,1]**2 - stress[:,0,2]**2 - stress[:,1,2]**2
        J2 = (1/3)*I1**2 - I2
        q_to = to.sqrt(3*J2)
        self.q_nodes.x.array[:] = self.grid.A_csr.dot(q_to.numpy())

    def compute_q_elems(self) -> None:
        """
        Compute von Mises equivalent stress and smooth it to elements.

        Returns
        -------
        None

        Side Effects
        ------------
        Sets :attr:`q_elems` by applying :attr:`grid.smoother` to nodal values.
        """
        stress = ut.numpy2torch(self.sig.x.array.reshape((self.n_elems, 3, 3)))
        I1 = stress[:,0,0] + stress[:,1,1] + stress[:,2,2]
        I2 = stress[:,0,0]*stress[:,1,1] + stress[:,1,1]*stress[:,2,2] + stress[:,0,0]*stress[:,2,2] - stress[:,0,1]**2 - stress[:,0,2]**2 - stress[:,1,2]**2
        J2 = (1/3)*I1**2 - I2
        q_to = to.sqrt(3*J2)
        self.q_elems.x.array[:] = self.grid.smoother.dot(q_to.numpy())

    def compute_total_strain(self) -> to.Tensor:
        """
        Project total small-strain tensor to DG0 and return as torch.

        Returns
        -------
        torch.Tensor
            Total strain per element, shape ``(n_elems, 3, 3)``.

        Notes
        -----
        Uses :func:`ut.project` on ``ε(u)``.
        """
        self.eps_tot = ut.project(ut.epsilon(self.u), self.DG0_3x3)
        eps_to = ut.numpy2torch(self.eps_tot.x.array.reshape((self.n_elems, 3, 3)))
        return eps_to

    def compute_eps_th(self) -> to.Tensor:
        """
        Compute element-wise thermal strain by aggregating thermoelastic elements.

        Returns
        -------
        torch.Tensor
            Thermal strain per element, shape ``(n_elems, 3, 3)``.
        """
        eps_th = to.zeros((self.n_elems, 3, 3), dtype=to.float64)
        deltaT = self.Temp - self.T0
        for elem_th in self.mat.elems_th:
        	elem_th.compute_eps_th(deltaT)
        	eps_th += elem_th.eps_th
        return eps_th

    def compute_eps_ne_k(self, dt: float) -> to.Tensor:
        """
        Compute predictor of non-elastic strain at the previous iteration k.

        Parameters
        ----------
        dt : float
            Time-step size.

        Returns
        -------
        torch.Tensor
            Predicted non-elastic strain per element, shape ``(n_elems, 3, 3)``.
        """
        eps_ne_k = to.zeros((self.n_elems, 3, 3), dtype=to.float64)
        for elem_ne in self.mat.elems_ne:
        	elem_ne.compute_eps_ne_k(dt*self.theta, dt*(1 - self.theta))
        	eps_ne_k += elem_ne.eps_ne_k
        return eps_ne_k

    def compute_eps_ne_rate(self, stress: to.Tensor, dt: float) -> None:
        """
        Update non-elastic strain rate for all non-elastic elements.

        Parameters
        ----------
        stress : torch.Tensor
            Stress per element, shape ``(n_elems, 3, 3)``.
        dt : float
            Time-step size.

        Returns
        -------
        None
        """
        for elem_ne in self.mat.elems_ne:
        	elem_ne.compute_eps_ne_rate(stress, dt*self.theta, self.Temp, return_eps_ne=False)

    def update_eps_ne_rate_old(self) -> None:
        """
        Update non-elastic strain rate from the previous time step “old”.

        Returns
        -------
        None
        """
        for elem_ne in self.mat.elems_ne:
        	elem_ne.update_eps_ne_rate_old()

    def update_eps_ne_old(self, stress: to.Tensor, stress_k: to.Tensor, dt: float) -> None:
        """
        Update non-elastic strain tensor from the previous time step “old”.

        Parameters
        ----------
        stress : torch.Tensor
            Stress at current iteration k+1, shape ``(n_elems, 3, 3)``.
        stress_k : torch.Tensor
            Stress from previous iteration k, shape ``(n_elems, 3, 3)``.
        dt : float
            Time-step size.

        Returns
        -------
        None
        """
        for elem_ne in self.mat.elems_ne:
        	elem_ne.update_eps_ne_old(stress, stress_k, dt*(1-self.theta))

    def increment_internal_variables(self, stress: to.Tensor, stress_k: to.Tensor, dt: float) -> None:
        """
        Increment material internal variables (e.g., hardening).

        Parameters
        ----------
        stress : torch.Tensor
        stress_k : torch.Tensor
        dt : float

        Returns
        -------
        None
        """
        for elem_ne in self.mat.elems_ne:
        	elem_ne.increment_internal_variables(stress, stress_k, dt)

    def update_internal_variables(self) -> None:
        """
        Commit internal variables at the end of a time step.

        Returns
        -------
        None
        """
        for elem_ne in self.mat.elems_ne:
        	elem_ne.update_internal_variables()

    def create_solution_vector(self) -> None:
        """
        Allocate the solution function `X` in the primary space.

        Returns
        -------
        None

        Side Effects
        ------------
        Sets :attr:`X` as a `dolfinx.fem.Function(self.V)`.
        """
        self.X = do.fem.Function(self.V)

    def run_after_solve(self) -> None:
        """
        Optional hook called after each linear solve.

        Returns
        -------
        None
        """
        pass

    @abstractmethod
    def compute_CT(self, dt: float, stress_k: to.Tensor):
        """
        Build the consistent tangent operator (per element).

        Parameters
        ----------
        dt : float
            Time-step size.
        stress_k : torch.Tensor
            Stress from previous iteration k, shape ``(n_elems, 3, 3)``.

        Returns
        -------
        None
        """
        pass

    @abstractmethod
    def compute_eps_rhs(self, dt: float, stress_k: to.Tensor, eps_k: to.Tensor):
        """
        Compute the right-hand-side strain term used in the linear form.

        Parameters
        ----------
        dt : float
            Time-step size.
        stress_k : torch.Tensor
            Intermediate stress, shape ``(n_elems, 3, 3)``.
        eps_k : torch.Tensor
            Optional strain input for schemes that need it.

        Returns
        -------
        None
        """
        pass

    @abstractmethod
    def compute_elastic_stress(self, eps_e: to.Tensor):
        """
        Compute elastic stress from elastic strain using `C`.

        Parameters
        ----------
        eps_e : torch.Tensor
            Elastic strain, shape ``(n_elems, 3, 3)``.

        Returns
        -------
        torch.Tensor
            Stress, shape ``(n_elems, 3, 3)``.
        """
        pass

    @abstractmethod
    def compute_stress(self, eps_tot: to.Tensor, eps_rhs: to.Tensor, p: to.Tensor):
        """
        Compute stress from total strain and RHS strain (and optionally pressure).

        Parameters
        ----------
        eps_tot : torch.Tensor
            Total strain, shape ``(n_elems, 3, 3)``.
        eps_rhs : torch.Tensor
            RHS strain term, shape ``(n_elems, 3, 3)``.
        p : torch.Tensor
            Optional pressure term.

        Returns
        -------
        torch.Tensor
            Stress, shape ``(n_elems, 3, 3)``.
        """
        pass

    @abstractmethod
    def create_fenicsx_fields(self) -> None:
        """
        Create FE functions specific to the concrete formulation.

        Returns
        -------
        None
        """
        pass

    @abstractmethod
    def create_pytorch_fields(self) -> None:
        """
        Create torch tensors specific to the concrete formulation.

        Returns
        -------
        None
        """
        pass

    @abstractmethod
    def create_trial_test_functions(self) -> None:
        """
        Create UFL trial and test functions in the primary space.

        Returns
        -------
        None
        """
        pass

    @abstractmethod
    def get_uV(self) -> do.fem.FunctionSpace:
        """
        Return the primary function space for displacements.

        Returns
        -------
        dolfinx.fem.FunctionSpace
            Vector function space used for `u`.
        """
        pass

    @abstractmethod
    def initialize(self) -> None:
        """
        Initialize FE fields from the material container.

        Returns
        -------
        None
        """
        pass

    @abstractmethod
    def split_solution(self) -> None:
        """
        Assign the computed solution `X` to the primary field (e.g., `u`).

        Returns
        -------
        None
        """
        pass

    @abstractmethod
    def split_solution(self) -> None:
        """
        Duplicate abstract declaration (kept for API compatibility).
        """
        pass

    @abstractmethod
    def compute_p_nodes(self) -> None:
        """
        Compute nodal pressure (mean stress) from the stress field.

        Returns
        -------
        None
        """
        pass

    @abstractmethod
    def solve_elastic_response(self) -> None:
        """
        Solve the purely elastic problem (e.g., for initialization).

        Returns
        -------
        None
        """
        pass

    @abstractmethod
    def solve(self) -> None:
        """
        Assemble and solve one step of the inelastic problem.

        Returns
        -------
        None
        """
        pass





class LinearMomentum(LinearMomentumBase):
    """
    Linear momentum formulation with thermo-(visco)elastic tangent.

    Implements the concrete FE fields, consistent tangent assembly, right-hand
    side strain, and linear solves for both elastic and inelastic steps.
    """
    def __init__(self, grid: GridHandlerGMSH, theta: float):
        """
        Initialize spaces, measures, fields, and solution vector.

        Parameters
        ----------
        grid : GridHandlerGMSH
        theta : float
        """
        super().__init__(grid, theta)
        self.V = self.CG1_3x1
        self.create_trial_test_functions()
        self.create_normal()
        self.create_solution_vector()

    def create_fenicsx_fields(self) -> None:
        """
        Allocate FE functions specific to this formulation.

        Returns
        -------
        None

        Side Effects
        ------------
        Creates :attr:`C`, :attr:`CT` (DG0 6×6 tangents), and :attr:`eps_rhs`
        (DG0 3×3 RHS strain).
        """
        self.C = do.fem.Function(self.DG0_6x6)
        self.CT = do.fem.Function(self.DG0_6x6)
        self.eps_rhs = do.fem.Function(self.DG0_3x3)

    def create_pytorch_fields(self) -> None:
        """
        Allocate torch tensors specific to this formulation.

        Returns
        -------
        None

        Side Effects
        ------------
        Creates :attr:`eps_rhs_to` with shape ``(n_elems, 3, 3)``.
        """
        self.eps_rhs_to = to.zeros((self.n_elems, 3, 3))

    def create_trial_test_functions(self) -> None:
        """
        Create UFL trial/test functions for displacement.

        Returns
        -------
        None

        Side Effects
        ------------
        Defines :attr:`du` (trial) and :attr:`u_` (test) in :attr:`V`.
        """
        self.du = ufl.TrialFunction(self.V)
        self.u_ = ufl.TestFunction(self.V)

    def get_uV(self) -> do.fem.FunctionSpace:
        """
        Return the primary displacement function space.

        Returns
        -------
        dolfinx.fem.FunctionSpace
        """
        return self.V

    def initialize(self) -> None:
        """
        Initialize elastic tangent from the material container.

        Returns
        -------
        None

        Side Effects
        ------------
        Flattens and copies `mat.C` into :attr:`C`.
        """
        self.C.x.array[:] = to.flatten(self.mat.C)

    def compute_CT(self, stress_k: to.Tensor, dt: float) -> None:
        """
        Assemble consistent tangent operator `CT` for the current step.

        Parameters
        ----------
        stress_k : torch.Tensor
            Stress from previous iteration k, shape ``(n_elems, 3, 3)``.
        dt : float
            Time-step size.

        Returns
        -------
        None

        Side Effects
        ------------
        Updates material operators (`G`, `B`, `CT`) and copies to FE field :attr:`CT`.
        """
        self.mat.compute_G_B(stress_k, dt, self.theta, self.Temp)
        self.mat.compute_CT(dt, self.theta)
        self.CT.x.array[:] = to.flatten(self.mat.CT)

    def compute_elastic_stress(self, eps_e: to.Tensor) -> to.Tensor:
        """
        Compute elastic Cauchy stress using the elastic stiffness `C`.

        Parameters
        ----------
        eps_e : torch.Tensor
            Elastic strain, shape ``(n_elems, 3, 3)``.

        Returns
        -------
        torch.Tensor
            Stress, shape ``(n_elems, 3, 3)``.

        Side Effects
        ------------
        Copies the stress into :attr:`sig`.
        """
        stress_to = ut.dotdot_torch(self.mat.C, eps_e)
        self.sig.x.array[:] = to.flatten(stress_to)
        return stress_to

    def compute_stress(self, eps_tot_to: to.Tensor, *_) -> to.Tensor:
        """
        Compute stress using the consistent tangent and RHS strain.

        Parameters
        ----------
        eps_tot_to : torch.Tensor
            Total strain, shape ``(n_elems, 3, 3)``.
        *_
            Unused extra arguments (kept for signature compatibility).

        Returns
        -------
        torch.Tensor
            Stress, shape ``(n_elems, 3, 3)``.

        Side Effects
        ------------
        Copies the stress into :attr:`sig`.
        """
        stress_to = ut.dotdot_torch(self.mat.CT, eps_tot_to - self.eps_rhs_to)
        self.sig.x.array[:] = to.flatten(stress_to)
        return stress_to

    def compute_eps_rhs(self, dt: float, stress_k: to.Tensor) -> None:
        """
        Compute the right-hand-side strain term for the variational form.

        Parameters
        ----------
        dt : float
            Time-step size.
        stress_k : torch.Tensor
            Intermediate stress, shape ``(n_elems, 3, 3)``.

        Returns
        -------
        None

        Side Effects
        ------------
        Sets :attr:`eps_rhs_to` (torch) and :attr:`eps_rhs` (FE field).
        """
        eps_ne_k = self.compute_eps_ne_k(dt)
        eps_th = self.compute_eps_th()
        self.eps_rhs_to = eps_ne_k + eps_th - dt*(1 - self.theta)*(self.mat.B + ut.dotdot_torch(self.mat.G, stress_k))
        self.eps_rhs.x.array[:] = to.flatten(self.eps_rhs_to)

    def solve_elastic_response(self) -> None:
        """
        Solve the purely elastic boundary-value problem.

        Returns
        -------
        None

        Side Effects
        ------------
        - Assembles and solves the linear system with :math:`C`.
        - Updates :attr:`X` and calls :meth:`split_solution`.
        """
        # Build bilinear form
        a = ufl.inner(ut.dotdot_ufl(self.C, ut.epsilon(self.du)), ut.epsilon(self.u_))*self.dx
        bilinear_form = do.fem.form(a)
        A = do.fem.petsc.assemble_matrix(bilinear_form, bcs=self.bc.dirichlet_bcs)
        A.assemble()

        # Build linear form
        linear_form = do.fem.form(self.b_body + sum(self.bc.neumann_bcs))
        b = fem_petsc.assemble_vector(linear_form)
        fem_petsc.apply_lifting(b, [bilinear_form], [self.bc.dirichlet_bcs])
        b.ghostUpdate(addv=PETSc.InsertMode.ADD_VALUES, mode=PETSc.ScatterMode.REVERSE)
        fem_petsc.set_bc(b, self.bc.dirichlet_bcs)
        b.ghostUpdate(addv=PETSc.InsertMode.INSERT, mode=PETSc.ScatterMode.FORWARD)

        # Solve linear system
        self.solver.setOperators(A)
        self.solver.solve(b, self.X.x.petsc_vec)
        self.X.x.scatter_forward()
        self.split_solution()

    def split_solution(self) -> None:
        """
        Assign displacement solution `X` to the primary field `u`.

        Returns
        -------
        None
        """
        self.u = self.X

	# def compute_p_nodes(self) -> do.fem.Function:
	# 	self.p_nodes = ut.project(ufl.tr(self.sig)/3, self.CG1_1)

	# def compute_p_elems(self) -> do.fem.Function:
	# 	# self.p_elems = ut.project(ufl.tr(self.sig)/3, self.DG0_1)
	# 	stress_to = ut.numpy2torch(self.sig.x.array.reshape((self.n_elems, 3, 3)))
	# 	p_to = to.einsum("kii->k", stress_to)
	# 	self.p_elems.x.array[:] = to.flatten(p_to)

    def compute_p_nodes(self) -> None:
        """
        Compute nodal pressure ``p = tr(σ)/3`` via node-element averaging.

        Returns
        -------
        None

        Side Effects
        ------------
        Writes to :attr:`p_nodes`.
        """
        stress = ut.numpy2torch(self.sig.x.array.reshape((self.n_elems, 3, 3)))
        I1 = stress[:,0,0] + stress[:,1,1] + stress[:,2,2]
        p_to = I1/3
        self.p_nodes.x.array[:] = self.grid.A_csr.dot(p_to)

    def compute_p_elems(self) -> None:
        """
        Compute element pressure by smoothing the nodal trace of stress.

        Returns
        -------
        None

        Side Effects
        ------------
        Writes to :attr:`p_elems`.
        """
        stress_to = ut.numpy2torch(self.sig.x.array.reshape((self.n_elems, 3, 3)))
        p_to = to.einsum("kii->k", stress_to)
        p_to = self.grid.smoother.dot(p_to.numpy())
        self.p_elems.x.array[:] = p_to

    def solve(self, stress_k_to: to.Tensor, t: float, dt: float) -> None:
        """
        Assemble and solve one implicit time step for the inelastic problem.

        Parameters
        ----------
        stress_k_to : torch.Tensor
            Stress at previous iteration k, shape ``(n_elems, 3, 3)``.
        t : float
            Current time (used by BC handler externally).
        dt : float
            Time-step size.

        Returns
        -------
        None

        Side Effects
        ------------
        - Builds `CT` and `eps_rhs`, assembles and solves the linear system.
        - Updates :attr:`X`, calls :meth:`split_solution`, then :meth:`run_after_solve`.
        """

        # Compute consistent tangent matrix
        self.compute_CT(stress_k_to, dt)

        # Compute right-hand side epsilon
        self.compute_eps_rhs(dt, stress_k_to)

        # Build bilinear form
        a = ufl.inner(ut.dotdot_ufl(self.CT, ut.epsilon(self.du)), ut.epsilon(self.u_))*self.dx
        bilinear_form = do.fem.form(a)
        A = fem_petsc.assemble_matrix(bilinear_form, bcs=self.bc.dirichlet_bcs)
        A.assemble()

        # Build linear form
        b_rhs = ufl.inner(ut.dotdot_ufl(self.CT, self.eps_rhs), ut.epsilon(self.u_))*self.dx
        linear_form = do.fem.form(self.b_body + sum(self.bc.neumann_bcs) + b_rhs)
        b = fem_petsc.assemble_vector(linear_form)
        fem_petsc.apply_lifting(b, [bilinear_form], [self.bc.dirichlet_bcs])
        b.ghostUpdate(addv=PETSc.InsertMode.ADD_VALUES, mode=PETSc.ScatterMode.REVERSE)
        fem_petsc.set_bc(b, self.bc.dirichlet_bcs)
        b.ghostUpdate(addv=PETSc.InsertMode.INSERT, mode=PETSc.ScatterMode.FORWARD)

        # Solve linear system
        self.solver.setOperators(A)
        self.solver.solve(b, self.X.x.petsc_vec)
        self.X.x.scatter_forward()
        self.split_solution()

        self.run_after_solve()



