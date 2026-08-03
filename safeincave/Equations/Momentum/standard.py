# Copyright (c) 2026, The SafeInCave Developers
#
# SPDX-License-Identifier: BSD-3-Clause

from __future__ import annotations
from typing import TYPE_CHECKING
import warnings

if TYPE_CHECKING:
    from ...Mesh.Grid import GridHandlerGMSH

import dolfinx as do
import ufl
from petsc4py import PETSc
import torch as to
from dolfinx.fem import petsc as fem_petsc
from .base import LinearMomentumBase
from ...Utils import dotdot_ufl, dotdot_torch, epsilon

class LinearMomentum(LinearMomentumBase):
    """
    Linear momentum formulation with thermo-(visco)elastic tangent.

    Implements the concrete FE fields, consistent tangent assembly, right-hand
    side strain, and linear solves for both elastic and inelastic steps.
    """

    def __init__(
        self,
        grid: GridHandlerGMSH,
        theta: float,
        solver_name: str = "gmres",
        preconditioner: str = "asm",
        rtol: float = 1e-12,
        max_it: int = 100,
    ):
        """
        Initialize spaces, measures, fields, and solution vector.

        Parameters
        ----------
        grid : GridHandlerGMSH
        theta : float
        """
        super().__init__(grid, theta, solver_name, preconditioner, rtol, max_it)
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
        self.eps_0 = do.fem.Function(self.DG0_3x3)

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
        self.eps0_to = to.zeros((self.n_elems, 3, 3))

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

    def compute_CT(self, dt: float, stress_k: to.Tensor) -> None:
        """
        Assemble consistent tangent operator `CT` for the current step.

        Parameters
        ----------
        dt : float
            Time-step size.
        stress_k : torch.Tensor
            Stress from previous iteration k, shape ``(n_elems, 3, 3)``.


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
        stress_to = dotdot_torch(self.mat.C, eps_e)
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
        stress_to = dotdot_torch(
            self.mat.CT, self.eps0_to + eps_tot_to - self.eps_rhs_to
        )
        # stress_to = dotdot_torch(self.mat.CT, eps_tot_to - self.eps_rhs_to)
        self.sig.x.array[:] = to.flatten(stress_to)
        return stress_to

    def apply_initial_stress(self, sig0: to.Tensor) -> None:
        """
        Set the initial element-wise stress.

        Parameters
        ----------
        sig0 : torch.Tensor
            Initial stress per element, shape ``(n_elems, 3, 3)``.

        Returns
        -------
        None

        Side Effects
        ------------
        Sets :attr:`sig` field.
        """
        self.eps0_to = dotdot_torch(self.mat.C_inv, sig0)
        self.eps_0.x.array[:] = to.flatten(self.eps0_to)
        # self.eps_tot.x.array[:] = to.flatten(self.eps0_to)
        self.sig.x.array[:] = to.flatten(sig0)
        # Kept for the exact-trial reconstruction in compute_eps_ne_rate.
        self.sig0_to = sig0.clone()
        self.sig0.x.array[:] = to.flatten(sig0)

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
        self.eps_rhs_to = (
            eps_ne_k
            + eps_th
            - dt * (1 - self.theta) * (self.mat.B + dotdot_torch(self.mat.G, stress_k))
        )
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
        # Compiled forms and PETSc objects are cached; only values change.
        key = self._bc_forms_key()
        if getattr(self, "_elastic_forms_key", None) != key:
            a = (
                ufl.inner(dotdot_ufl(self.C, epsilon(self.du)), epsilon(self.u_))
                * self.dx
            )
            bilinear_form = do.fem.form(a)
            b_rhs = ufl.inner(self.sig0, epsilon(self.u_)) * self.dx
            linear_form = do.fem.form(
                self.b_body
                + sum(self.bc.neumann_bcs)
                + sum(self.bc.cavern_bcs)
                + b_rhs
            )
            A = fem_petsc.create_matrix(bilinear_form)
            b = fem_petsc.create_vector(linear_form)
            self._elastic_forms = (bilinear_form, linear_form, A, b)
            self._elastic_forms_key = key
        bilinear_form, linear_form, A, b = self._elastic_forms

        A.zeroEntries()
        fem_petsc.assemble_matrix(A, bilinear_form, bcs=self.bc.dirichlet_bcs)
        A.assemble()

        with b.localForm() as b_local:
            b_local.set(0.0)
        fem_petsc.assemble_vector(b, linear_form)
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
        # Cached fem.Expression interpolated into the persistent p_elems
        # Function (self.sig is persistent, so the Expression stays valid) —
        # avoids the per-call Function/Expression allocation of Utils.project.
        if getattr(self, "_p_elems_expr", None) is None:
            self._p_elems_expr = do.fem.Expression(
                ufl.tr(self.sig) / 3, self.DG0_1.element.interpolation_points()
            )
        self.p_elems.interpolate(self._p_elems_expr)

    def _bc_forms_key(self) -> tuple:
        """
        Identity key of the BC terms referenced by cached compiled forms.

        BC Constants are updated in place (see BcHandler), so cached forms
        stay valid as long as the term objects themselves are unchanged; a
        change in the registered BC or cavern set produces new term objects
        and must invalidate the cached forms.
        """
        return (
            tuple(map(id, self.bc.dirichlet_bcs)),
            tuple(map(id, self.bc.neumann_bcs)),
            tuple(map(id, getattr(self.bc, "cavern_bcs", []))),
        )

    def _ensure_solve_forms(self) -> None:
        """
        Compile the bilinear/linear forms of :meth:`solve` once and allocate
        the persistent PETSc matrix/vector. All time- and iteration-varying
        data (CT, eps_rhs, eps_0, BC loads) live in Function/Constant objects
        updated in place, so re-assembly with the same compiled forms is
        exact.
        """
        key = self._bc_forms_key()
        if getattr(self, "_solve_forms_key", None) == key:
            return
        a = ufl.inner(dotdot_ufl(self.CT, epsilon(self.du)), epsilon(self.u_)) * self.dx
        bilinear_form = do.fem.form(a)
        b_rhs = (
            ufl.inner(dotdot_ufl(self.CT, self.eps_rhs - self.eps_0), epsilon(self.u_))
            * self.dx
        )
        linear_form = do.fem.form(
            self.b_body + sum(self.bc.neumann_bcs) + sum(self.bc.cavern_bcs) + b_rhs
        )
        A = fem_petsc.create_matrix(bilinear_form)
        b = fem_petsc.create_vector(linear_form)
        self._solve_forms = (bilinear_form, linear_form, A, b)
        self._solve_forms_key = key

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
        self.compute_CT(dt, stress_k_to)

        # Compute right-hand side epsilon
        self.compute_eps_rhs(dt, stress_k_to)

        # Compiled forms and PETSc objects are cached; only values change.
        self._ensure_solve_forms()
        bilinear_form, linear_form, A, b = self._solve_forms

        A.zeroEntries()
        fem_petsc.assemble_matrix(A, bilinear_form, bcs=self.bc.dirichlet_bcs)
        A.assemble()

        with b.localForm() as b_local:
            b_local.set(0.0)
        fem_petsc.assemble_vector(b, linear_form)
        fem_petsc.apply_lifting(b, [bilinear_form], [self.bc.dirichlet_bcs])
        b.ghostUpdate(addv=PETSc.InsertMode.ADD_VALUES, mode=PETSc.ScatterMode.REVERSE)
        fem_petsc.set_bc(b, self.bc.dirichlet_bcs)
        b.ghostUpdate(addv=PETSc.InsertMode.INSERT, mode=PETSc.ScatterMode.FORWARD)

        # Solve linear system
        self.solver.setOperators(A)
        self.solver.solve(b, self.X.x.petsc_vec)

        # A silently failed linear solve feeds garbage displacements into the
        # nonlinear loop and can masquerade as material-model divergence.
        reason = self.solver.getConvergedReason()
        if reason < 0:
            warnings.warn(
                f"Linear solver did not converge (KSP reason={reason}, "
                f"iterations={self.solver.getIterationNumber()}, "
                f"residual={self.solver.getResidualNorm():.3e}). "
                "Consider a direct solver (preonly/lu) or higher max_it.",
                RuntimeWarning,
            )

        self.X.x.scatter_forward()
        self.split_solution()

        self.run_after_solve()
