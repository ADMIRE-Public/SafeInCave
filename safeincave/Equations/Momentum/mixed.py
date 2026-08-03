# Copyright (c) 2026, The SafeInCave Developers
#
# SPDX-License-Identifier: BSD-3-Clause

from __future__ import annotations
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    pass

import basix
import dolfinx as do
import ufl
from petsc4py import PETSc
import torch as to
from dolfinx.fem import petsc as fem_petsc
from .base import LinearMomentumBase
from ...Utils import numpy2torch, dotdot_ufl, dotdot_torch, epsilon, project, ModelML

class LinearMomentumMixed(LinearMomentumBase):
    def __init__(
        self,
        grid,
        theta,
        stab_scaling: float = 1.0,
        solver_name: str = "gmres",
        preconditioner: str = "asm",
        rtol: float = 1e-12,
        max_it: int = 100,
    ):
        super().__init__(grid, theta, solver_name, preconditioner, rtol, max_it)
        Vue = basix.ufl.element(
            "CG", self.grid.mesh.basix_cell(), 1, shape=(3,)
        )  # displacement finite element
        Vpe = basix.ufl.element(
            "CG", self.grid.mesh.basix_cell(), 1
        )  # mean stress finite element
        el_mixed = basix.ufl.mixed_element([Vue, Vpe])
        self.V = do.fem.functionspace(self.grid.mesh, el_mixed)
        self.create_trial_test_functions()
        self.create_normal()
        self.create_solution_vector()

        assert stab_scaling >= 0.0, "stab_scaling must be >= 0.0"
        if stab_scaling == 0.0:
            pass
        else:
            self.calculate_h(stab_scaling)

    def calculate_h(self, stab_scaling: float = 0.0) -> None:
        model_h = ModelML()
        h = stab_scaling * model_h.compute_mesh_h(self.grid.mesh)
        self.h_cell_2.x.array[:] = h**2

    def create_fenicsx_fields(self):
        self.CT_tilde = do.fem.Function(self.DG0_6x6)
        self.C_tilde = do.fem.Function(self.DG0_6x6)
        self.C_tilde_inv = do.fem.Function(self.DG0_6x6)
        self.T_vol = do.fem.Function(self.DG0_1)
        self.B_vol = do.fem.Function(self.DG0_1)
        self.eps_rhs_tilde = do.fem.Function(self.DG0_3x3)
        self.eps_ne_vol = do.fem.Function(self.DG0_1)
        self.eps_th_vol = do.fem.Function(self.DG0_1)
        self.eps_0_tilde = do.fem.Function(self.DG0_3x3)
        self.K = do.fem.Function(self.DG0_1)
        self.E = do.fem.Function(self.DG0_1)
        self.E_star = do.fem.Function(self.DG0_1)
        self.h_cell_2 = do.fem.Function(self.DG0_1)
        self.p_k = do.fem.Function(self.CG1_1)

    def create_pytorch_fields(self):
        self.eps_rhs_to = to.zeros((self.n_elems, 3, 3))
        self.eps_rhs_tilde_to = to.zeros((self.n_elems, 3, 3))
        self.eps_0_tilde_to = to.zeros((self.n_elems, 3, 3))

    def create_trial_test_functions(self):
        self.du, self.dp = ufl.TrialFunctions(self.V)
        self.u_, self.p_ = ufl.TestFunctions(self.V)

    def get_uV(self):
        return self.V.sub(0)

    def initialize(self) -> None:
        self.C_tilde.x.array[:] = self.mat.C_tilde.flatten()
        self.C_tilde_inv.x.array[:] = self.mat.C_tilde_inv.flatten()
        self.K.x.array[:] = self.mat.K
        self.E.x.array[:] = self.mat.E
        self.E_star.x.array[:] = self.mat.E

    def compute_CT(self, dt, stress_k):
        self.mat.compute_G_B(stress_k, dt, self.theta, self.Temp)
        self.mat.compute_T_IT()
        self.mat.compute_Bvol_Tvol(stress_k, dt)
        self.mat.compute_Gtilde_Btilde(stress_k, dt)
        self.mat.compute_CT_tilde(dt, self.theta)
        self.CT_tilde.x.array[:] = to.flatten(self.mat.CT_tilde)
        self.T_vol.x.array[:] = self.mat.T_vol
        self.B_vol.x.array[:] = self.mat.B_vol

    def compute_elastic_stress(self, eps_e: to.Tensor) -> to.Tensor:
        I_3x3 = to.eye(3).expand(self.n_elems, -1, -1)
        eps_e_tilde = self.compute_eps_tilde(eps_e)
        stress_to = (
            dotdot_torch(self.mat.C_tilde, eps_e_tilde)
            + self.p_to[:, None, None] * I_3x3
        )
        self.sig.x.array[:] = to.flatten(stress_to)
        return stress_to

    def compute_stress(self, eps_tot):
        eps_tilde = self.compute_eps_tilde(eps_tot)
        I_3x3 = to.eye(3).expand(self.n_elems, -1, -1)
        pI = self.p_to[:, None, None] * I_3x3
        eps_aux = eps_tilde - self.eps_rhs_tilde_to + self.eps_0_tilde_to
        eps_aux += dotdot_torch(self.mat.C_tilde_inv, pI)
        stress_to = dotdot_torch(self.mat.CT_tilde, eps_aux)
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
        # self.eps0_to = dotdot_torch(self.mat.C_inv, sig0)
        one_over_2G = self.mat.C_tilde_inv[:,0,0]
        self.eps_0_tilde_to = sig0 * one_over_2G[:, None, None]
        self.eps_0_tilde.x.array[:] = to.flatten(self.eps_0_tilde_to)
        self.sig.x.array[:] = to.flatten(sig0)
        # Kept for the exact-trial reconstruction in compute_eps_ne_rate.
        self.sig0_to = sig0.clone()


    def compute_eps_k_ne_vol(self, eps_k_ne):
        eps_k_ne_vol = to.einsum("bii->b", eps_k_ne)
        return eps_k_ne_vol

    def compute_eps_k_tilde(self, eps_k_ne, eps_k_ne_vol):
        I_3x3 = to.eye(3).expand(self.n_elems, -1, -1)
        eps_k_ne_tilde = eps_k_ne - (1 / 3) * eps_k_ne_vol[:, None, None] * I_3x3
        return eps_k_ne_tilde

    def compute_eps_tilde(self, eps):
        I_3x3 = to.eye(3).expand(self.n_elems, -1, -1)
        eps_vol = to.einsum("bii->b", eps)[:, None, None]
        return eps - (1 / 3) * eps_vol * I_3x3

    def compute_eps_rhs(self, dt, stress_k, eps_k_tilde):
        self.eps_rhs_tilde_to = eps_k_tilde - dt * (1 - self.theta) * (
            self.mat.B_tilde + dotdot_torch(self.mat.G_tilde, stress_k)
        )
        self.eps_rhs_tilde.x.array[:] = to.flatten(self.eps_rhs_tilde_to)

    def split_solution(self):
        self.u = self.X.sub(0).collapse()
        self.p_nodes = self.X.sub(1).collapse()

        self.u.x.scatter_forward()
        self.p_nodes.x.scatter_forward()

        # Project mean stress to Gauss points
        p_elems = project(self.p_nodes, self.DG0_1)
        self.p_to = numpy2torch(p_elems.x.array)

    def compute_p_nodes(self) -> do.fem.Function:
        return self.p_nodes

    def compute_p_elems(self) -> do.fem.Function:
        self.p_elems = project(ufl.tr(self.sig) / 3, self.DG0_1)
        # self.p_elems = project(self.p_nodes, self.DG0_1)
        return self.p_elems

    def compute_moduli(self, stress_to):
        """
        Update the stabilization modulus E_star from the current principal
        stress/strain ratio.

        Principal values are computed with the closed-form symmetric-3x3
        formula (robust for repeated eigenvalues, e.g. the isotropic
        geostatic state). epsil_1 is clamped away from zero to avoid
        dividing by a vanishing strain at the very first step.
        """
        strain_to = self.compute_total_strain()

        def _sym3x3_eigvals(A: to.Tensor) -> to.Tensor:
            """
            Closed-form (trigonometric) eigenvalues of symmetric 3x3 tensors,
            ascending, shape (N, 3). Replaces torch.linalg.eigvalsh, which is
            both orders of magnitude slower for large batches of tiny matrices
            and prone to non-convergence on (near-)isotropic tensors.
            """
            q = to.diagonal(A, dim1=-2, dim2=-1).sum(-1) / 3.0
            p1 = A[:, 0, 1] ** 2 + A[:, 0, 2] ** 2 + A[:, 1, 2] ** 2
            d0 = A[:, 0, 0] - q
            d1 = A[:, 1, 1] - q
            d2 = A[:, 2, 2] - q
            p2 = d0**2 + d1**2 + d2**2 + 2.0 * p1
            p = to.sqrt(to.clamp(p2 / 6.0, min=0.0))
            p_safe = to.clamp(p, min=1e-30)
            # det(B) with B = (A - q I)/p, via explicit cofactor expansion
            detB = (
                d0 * (d1 * d2 - A[:, 1, 2] ** 2)
                - A[:, 0, 1] * (A[:, 0, 1] * d2 - A[:, 1, 2] * A[:, 0, 2])
                + A[:, 0, 2] * (A[:, 0, 1] * A[:, 1, 2] - d1 * A[:, 0, 2])
            ) / p_safe**3
            r = to.clamp(detB / 2.0, min=-1.0, max=1.0)
            phi = to.acos(r) / 3.0
            e_hi = q + 2.0 * p * to.cos(phi)
            e_lo = q + 2.0 * p * to.cos(phi + 2.0 * to.pi / 3.0)
            e_mid = 3.0 * q - e_hi - e_lo
            return to.stack([e_lo, e_mid, e_hi], dim=1)

        principal_stresses = _sym3x3_eigvals(stress_to)
        principal_strains = _sym3x3_eigvals(strain_to)
        sigma_1 = principal_stresses[:, 0]
        sigma_2 = principal_stresses[:, 1]
        sigma_3 = principal_stresses[:, 2]
        epsil_1 = principal_strains[:, 0]
        epsil_1_safe = to.where(
            epsil_1.abs() < 1e-12, to.full_like(epsil_1, 1e-12), epsil_1
        )
        nu = self.mat.elems_e[0].nu
        E_star_1 = (sigma_1 - nu * (sigma_2 + sigma_3)) / epsil_1_safe
        # Undeformed elements (e.g. the very first iteration) give E_star = 0,
        # which blows up the K/E_star stabilization term — fall back to the
        # elastic modulus there.
        E_elastic = self.mat.elems_e[0].E
        E_star_1 = to.where(
            E_star_1.abs() < 1e-6 * E_elastic, E_elastic, E_star_1
        )
        self.E_star.x.array[:] = E_star_1

    def solve_elastic_response(self):
        # Build bilinear form
        eps_u = epsilon(self.du)
        eps_tilde = eps_u - (1 / 3) * ufl.tr(eps_u) * ufl.Identity(3)
        a = (
            ufl.inner(
                dotdot_ufl(
                    self.C_tilde,
                    eps_tilde + dotdot_ufl(self.C_tilde_inv, self.dp * ufl.Identity(3)),
                ),
                epsilon(self.u_),
            )
            * self.dx
        )
        a += (self.dp * self.p_ - self.K * ufl.tr(epsilon(self.du)) * self.p_) * self.dx
        a += (
            3
            * ufl.dot(
                self.h_cell_2 * (self.K / self.E) * ufl.grad(self.dp), ufl.grad(self.p_)
            )
        ) * self.dx
        bilinear_form = do.fem.form(a)
        A = do.fem.petsc.assemble_matrix(bilinear_form, bcs=self.bc.dirichlet_bcs)
        A.assemble()

        # Build linear form
        linear_form = do.fem.form(
            self.b_body + 
            sum(self.bc.neumann_bcs) + 
            sum(self.bc.cavern_bcs) +
            ufl.inner(self.sig0, epsilon(self.u_)) * self.dx
        )
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

    def solve(self, stress_k_to, t, dt):
        # Calculate stabilization moduli
        self.compute_moduli(stress_k_to)

        # Compute consistent tangent matrix
        self.compute_CT(dt, stress_k_to)

        # Compute epsilons
        eps_ne_k = self.compute_eps_ne_k(dt)
        eps_ne_k_vol = self.compute_eps_k_ne_vol(eps_ne_k)
        eps_k_tilde_to = self.compute_eps_k_tilde(eps_ne_k, eps_ne_k_vol)

        # Compute right-hand side epsilon
        self.compute_eps_rhs(dt, stress_k_to, eps_k_tilde_to)

        # Compute non-elastic volumetric strain
        self.eps_ne_vol.x.array[:] = eps_ne_k_vol

        # Compute thermoelastic strains
        eps_th_to = self.compute_eps_th()
        eps_th_vol_to = to.einsum("bii->b", eps_th_to)
        self.eps_th_vol.x.array[:] = eps_th_vol_to

        # phi2 enters the UFL forms as a fem.Constant so changing dt does NOT
        # require recompiling them: only the Constant's value is updated.
        phi2 = dt * (1 - self.theta)
        if not hasattr(self, "_phi2_const"):
            self._phi2_const = do.fem.Constant(self.grid.mesh, float(phi2))
        else:
            self._phi2_const.value = float(phi2)

        # Build (and cache) the bi-linear form. All spatially varying
        # coefficients are fem.Functions updated in place, so the compiled
        # form stays valid across iterations and time steps.
        if not hasattr(self, "_bilinear_form"):
            eps_u = epsilon(self.du)
            eps_tilde = eps_u - (1 / 3) * ufl.tr(eps_u) * ufl.Identity(3)
            a = (
                ufl.inner(
                    dotdot_ufl(
                        self.CT_tilde,
                        eps_tilde
                        + dotdot_ufl(self.C_tilde_inv, self.dp * ufl.Identity(3)),
                    ),
                    epsilon(self.u_),
                )
                * self.dx
            )
            a += (
                (1 + self._phi2_const * self.K * self.T_vol) * self.dp * self.p_
                - self.K * ufl.tr(epsilon(self.du)) * self.p_
            ) * self.dx
            a += (
                3
                * ufl.dot(
                    self.h_cell_2 * (self.K / self.E_star) * ufl.grad(self.dp),
                    ufl.grad(self.p_),
                )
            ) * self.dx
            self._bilinear_form = do.fem.form(a)
            self._A_mat = do.fem.petsc.create_matrix(self._bilinear_form)
        bilinear_form = self._bilinear_form
        A = self._A_mat
        A.zeroEntries()
        do.fem.petsc.assemble_matrix(A, bilinear_form, bcs=self.bc.dirichlet_bcs)
        A.assemble()

        # Build (and cache) the linear form. BC loads live in in-place
        # updated fem.Constants (see BcHandler), so the compiled form stays
        # valid while the BC/cavern term objects are unchanged; the identity
        # key rebuilds it if the BC or cavern set itself changes.
        bc_key = (
            tuple(map(id, self.bc.neumann_bcs)),
            tuple(map(id, getattr(self.bc, "cavern_bcs", []))),
        )
        if getattr(self, "_linear_form_key", None) != bc_key:
            b_u = (
                ufl.inner(
                    dotdot_ufl(self.CT_tilde, self.eps_rhs_tilde - self.eps_0_tilde),
                    epsilon(self.u_),
                )
                * self.dx
            )
            b_p = (
                self.K
                * (
                    self._phi2_const * (self.T_vol * self.p_k + self.B_vol)
                    - self.eps_ne_vol
                    - self.eps_th_vol
                )
                * self.p_
                * self.dx
            )
            self._linear_form = do.fem.form(
                self.b_body
                + sum(self.bc.neumann_bcs)
                + sum(self.bc.cavern_bcs)
                + b_u
                + b_p
            )
            self._b_vec = fem_petsc.create_vector(self._linear_form)
            self._linear_form_key = bc_key
        linear_form = self._linear_form
        b = self._b_vec
        with b.localForm() as b_local:
            b_local.set(0.0)
        fem_petsc.assemble_vector(b, linear_form)
        do.fem.petsc.apply_lifting(b, [bilinear_form], [self.bc.dirichlet_bcs])
        b.ghostUpdate(addv=PETSc.InsertMode.ADD_VALUES, mode=PETSc.ScatterMode.REVERSE)
        do.fem.petsc.set_bc(b, self.bc.dirichlet_bcs)
        b.ghostUpdate(addv=PETSc.InsertMode.INSERT, mode=PETSc.ScatterMode.FORWARD)

        # Solve linear system
        self.solver.setOperators(A)
        self.solver.solve(b, self.X.x.petsc_vec)
        self.X.x.scatter_forward()
        self.split_solution()

        # Update p_k
        self.p_k.x.array[:] = self.p_nodes.x.array
        self.p_k.x.scatter_forward()

        self.run_after_solve()
