# Copyright (c) 2026, The SafeInCave Developers
#
# SPDX-License-Identifier: BSD-3-Clause

from __future__ import annotations
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from ...MomentumBC import BcHandler

from ..base import EquationBase
from abc import ABC, abstractmethod
import warnings
import dolfinx as do
from dolfinx.fem import petsc as fem_petsc
import basix
import ufl
from petsc4py import PETSc
import torch as to
from ...Materials.Material import Material
from ...Mesh.Grid import GridHandlerGMSH
from ...Utils import numpy2torch, project, epsilon, dotdot_torch, dotdot_ufl
from ...Mesh.MeshParameter import ModelML

class LinearMomentumBase(EquationBase, ABC):
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
    q_elems, p_elems : dolfinx.fem.Function
        Von Mises magnitude and pressure (DG0 scalar fields, smoothable to CG1 via smooth_output flag).
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

    def __init__(
        self,
        grid: GridHandlerGMSH,
        theta: float,
        solver_name: str = "gmres",
        preconditioner: str = "asm",
        rtol: float = 1e-12,
        max_it: int = 100,
    ):
        super().__init__(grid, solver_name, preconditioner, rtol, max_it)
        self.theta = theta

        self.n_elems = self.DG0_1.dofmap.index_map.size_local + len(
            self.DG0_1.dofmap.index_map.ghosts
        )
        self.n_nodes = self.CG1_1.dofmap.index_map.size_local + len(
            self.CG1_1.dofmap.index_map.ghosts
        )

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
        `T0`, `Temp`, `sig`, `eps_tot`, `u`, `q_elems`, `p_elems`.
        """
        self.T0 = to.zeros(self.n_elems, dtype=to.float64)
        self.Temp = to.zeros(self.n_elems, dtype=to.float64)
        self.sig = do.fem.Function(self.DG0_3x3)
        self.sig0 = do.fem.Function(self.DG0_3x3)
        self.eps_tot = do.fem.Function(self.DG0_3x3)
        self.eps_0 = do.fem.Function(self.DG0_3x3)
        self.u = do.fem.Function(self.CG1_3x1)
        self.q_elems = do.fem.Function(self.DG0_1)
        self.p_elems = do.fem.Function(self.DG0_1)
        self.principal_stresses = do.fem.Function(self.DG0_3x1)
        self.principal_stress_directions = do.fem.Function(self.DG0_3x3)
        self.principal_strains = do.fem.Function(self.DG0_3x1)
        self.principal_strain_directions = do.fem.Function(self.DG0_3x3)

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
        bc.set_uV(self.get_uV())
        bc.set_boundary_dim(self.grid.boundary_dim)
        bc.set_boudary_tags(self.grid.boundary_tags)
        bc.set_dolfin_tags(self.grid.dolfin_tags)
        bc.set_normal(self.normal)
        bc.set_ds(self.ds)
        bc.set_spatial_coordinates(self.grid.mesh)
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
        self.CG1_3x1 = do.fem.functionspace(
            self.grid.mesh, ("Lagrange", 1, (self.grid.domain_dim,))
        )
        self.DG0_1 = do.fem.functionspace(self.grid.mesh, ("DG", 0))
        self.CG1_1 = do.fem.functionspace(self.grid.mesh, ("Lagrange", 1))
        self.DG0_3x1 = do.fem.functionspace(self.grid.mesh, ("DG", 0, (3,)))
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
        self.ds = ufl.Measure(
            "ds", domain=self.grid.mesh, subdomain_data=self.grid.get_boundaries()
        )
        self.dx = ufl.Measure(
            "dx", domain=self.grid.mesh, subdomain_data=self.grid.get_subdomains()
        )

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
        body_force = density * do.fem.Constant(
            self.grid.mesh, do.default_scalar_type(tuple(g))
        )
        self.b_body = ufl.dot(body_force, self.u_) * self.dx

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
        stress = numpy2torch(self.sig.x.array.reshape((self.n_elems, 3, 3)))
        I1 = stress[:, 0, 0] + stress[:, 1, 1] + stress[:, 2, 2]
        I2 = (
            stress[:, 0, 0] * stress[:, 1, 1]
            + stress[:, 1, 1] * stress[:, 2, 2]
            + stress[:, 0, 0] * stress[:, 2, 2]
            - stress[:, 0, 1] ** 2
            - stress[:, 0, 2] ** 2
            - stress[:, 1, 2] ** 2
        )
        J2 = (1 / 3) * I1**2 - I2
        q_to = to.sqrt(3 * J2)
        self.q_elems.x.array[:] = self.grid.smoother.dot(q_to.numpy())

    def compute_principal_stresses(self) -> None:
        """
        Compute principal stress magnitudes and directions from the full
        Cauchy stress tensor using closed-form eigenvalue decomposition
        (robust for repeated eigenvalues, e.g. isotropic stress states).

        Principal stresses are sorted descending: σ₁ ≥ σ₂ ≥ σ₃
        (most tensile / least compressive first).

        Principal directions are stored as a 3×3 tensor where each column
        is the eigenvector corresponding to (σ₁, σ₂, σ₃).

        Returns
        -------
        None

        Side Effects
        ------------
        Sets :attr:`principal_stresses` (DG0_3x1) and
        :attr:`principal_stress_directions` (DG0_3x3).
        """
        stress = numpy2torch(self.sig.x.array.reshape((self.n_elems, 3, 3)))

        # Closed-form eigenvalues and eigenvectors for symmetric 3x3 tensors.
        # Uses the trigonometric (Cardano) formula for eigenvalues, then
        # solves (A - λI) v = 0 via cross-product of two non-parallel rows.
        # This is both faster and more robust than iterative methods for
        # large batches of small matrices.

        # --- Step 1: Eigenvalues via trigonometric formula ---
        q = to.diagonal(stress, dim1=-2, dim2=-1).sum(-1) / 3.0  # trace/3
        p1 = stress[:, 0, 1] ** 2 + stress[:, 0, 2] ** 2 + stress[:, 1, 2] ** 2
        d0 = stress[:, 0, 0] - q
        d1 = stress[:, 1, 1] - q
        d2 = stress[:, 2, 2] - q
        p2 = d0**2 + d1**2 + d2**2 + 2.0 * p1
        p = to.sqrt(to.clamp(p2 / 6.0, min=0.0))
        p_safe = to.clamp(p, min=1e-30)

        # det(B) with B = (A - qI) / p, via explicit cofactor expansion
        detB = (
            d0 * (d1 * d2 - stress[:, 1, 2] ** 2)
            - stress[:, 0, 1] * (stress[:, 0, 1] * d2 - stress[:, 1, 2] * stress[:, 0, 2])
            + stress[:, 0, 2] * (stress[:, 0, 1] * stress[:, 1, 2] - d1 * stress[:, 0, 2])
        ) / p_safe**3
        r = to.clamp(detB / 2.0, min=-1.0, max=1.0)
        phi = to.acos(r) / 3.0

        eig_hi = q + 2.0 * p * to.cos(phi)               # σ₁ (largest)
        eig_lo = q + 2.0 * p * to.cos(phi + 2.0 * to.pi / 3.0)  # σ₃ (smallest)
        eig_mid = 3.0 * q - eig_hi - eig_lo               # σ₂ (middle)

        # Stack descending: σ₁ ≥ σ₂ ≥ σ₃
        eigenvalues = to.stack([eig_hi, eig_mid, eig_lo], dim=1)  # (n_elems, 3)

        # --- Step 2: Eigenvectors ---
        # For each eigenvalue λ, solve (A - λI) v = 0.
        # Use the cross-product of two non-parallel rows of (A - λI).
        # This avoids singular matrix issues and is robust for repeated roots.
        n_elems = self.n_elems
        eye_3 = to.eye(3, device=stress.device, dtype=stress.dtype)
        eigenvectors = to.zeros((n_elems, 3, 3), dtype=stress.dtype)

        for i in range(3):
            lam = eigenvalues[:, i]  # (n_elems,)
            # (A - λI)
            A_shifted = stress - lam[:, None, None] * eye_3[None, :, :]  # (n_elems, 3, 3)

            # Find two non-parallel rows by taking cross products
            # Row0 × Row1, Row1 × Row2, Row2 × Row0 — pick the one with largest norm
            r0 = A_shifted[:, 0, :]  # (n_elems, 3)
            r1 = A_shifted[:, 1, :]  # (n_elems, 3)
            r2 = A_shifted[:, 2, :]  # (n_elems, 3)

            c01 = to.linalg.cross(r0, r1)  # (n_elems, 3)
            c12 = to.linalg.cross(r1, r2)  # (n_elems, 3)
            c20 = to.linalg.cross(r2, r0)  # (n_elems, 3)

            norms01 = to.linalg.vector_norm(c01, dim=1)
            norms12 = to.linalg.vector_norm(c12, dim=1)
            norms20 = to.linalg.vector_norm(c20, dim=1)

            # Stack and pick the cross product with largest norm per element
            candidates = to.stack([c01, c12, c20], dim=1)  # (n_elems, 3, 3)
            candidate_norms = to.stack([norms01, norms12, norms20], dim=1)  # (n_elems, 3)
            best_idx = to.argmax(candidate_norms, dim=1)  # (n_elems,)

            # Gather the best candidate per element
            batch_indices = to.arange(n_elems, device=stress.device)
            vec = candidates[batch_indices, best_idx]  # (n_elems, 3)

            # Normalize
            vec_norm = to.linalg.vector_norm(vec, dim=1, keepdim=True)
            vec_norm_safe = to.clamp(vec_norm, min=1e-30)
            vec = vec / vec_norm_safe

            eigenvectors[:, :, i] = vec

        # --- Step 3: Store to FE fields ---
        # Principal stress magnitudes as a 3-component vector (σ₁, σ₂, σ₃)
        self.principal_stresses.x.array[:] = to.flatten(eigenvalues)

        # Principal stress directions as a 3x3 tensor (columns = eigenvectors)
        self.principal_stress_directions.x.array[:] = to.flatten(eigenvectors)

    def compute_principal_strains(self) -> None:
        """
        Compute principal strain magnitudes and directions from the full
        small-strain tensor using closed-form eigenvalue decomposition
        (robust for repeated eigenvalues, e.g. isotropic strain states).

        Principal strains are sorted descending: ε₁ ≥ ε₂ ≥ ε₃.

        Principal directions are stored as a 3×3 tensor where each column
        is the eigenvector corresponding to (ε₁, ε₂, ε₃).

        Returns
        -------
        None

        Side Effects
        ------------
        Sets :attr:`principal_strains` (DG0_3x1) and
        :attr:`principal_strain_directions` (DG0_3x3).
        """
        strain = numpy2torch(self.eps_tot.x.array.reshape((self.n_elems, 3, 3)))

        # Closed-form eigenvalues and eigenvectors for symmetric 3x3 tensors.
        # Uses the trigonometric (Cardano) formula for eigenvalues, then
        # solves (A - λI) v = 0 via cross-product of two non-parallel rows.
        # This is both faster and more robust than iterative methods for
        # large batches of small matrices.

        # --- Step 1: Eigenvalues via trigonometric formula ---
        q = to.diagonal(strain, dim1=-2, dim2=-1).sum(-1) / 3.0  # trace/3
        p1 = strain[:, 0, 1] ** 2 + strain[:, 0, 2] ** 2 + strain[:, 1, 2] ** 2
        d0 = strain[:, 0, 0] - q
        d1 = strain[:, 1, 1] - q
        d2 = strain[:, 2, 2] - q
        p2 = d0**2 + d1**2 + d2**2 + 2.0 * p1
        p = to.sqrt(to.clamp(p2 / 6.0, min=0.0))
        p_safe = to.clamp(p, min=1e-30)

        # det(B) with B = (A - qI) / p, via explicit cofactor expansion
        detB = (
            d0 * (d1 * d2 - strain[:, 1, 2] ** 2)
            - strain[:, 0, 1] * (strain[:, 0, 1] * d2 - strain[:, 1, 2] * strain[:, 0, 2])
            + strain[:, 0, 2] * (strain[:, 0, 1] * strain[:, 1, 2] - d1 * strain[:, 0, 2])
        ) / p_safe**3
        r = to.clamp(detB / 2.0, min=-1.0, max=1.0)
        phi = to.acos(r) / 3.0

        eig_hi = q + 2.0 * p * to.cos(phi)               # ε₁ (largest)
        eig_lo = q + 2.0 * p * to.cos(phi + 2.0 * to.pi / 3.0)  # ε₃ (smallest)
        eig_mid = 3.0 * q - eig_hi - eig_lo               # ε₂ (middle)

        # Stack descending: ε₁ ≥ ε₂ ≥ ε₃
        eigenvalues = to.stack([eig_hi, eig_mid, eig_lo], dim=1)  # (n_elems, 3)

        # --- Step 2: Eigenvectors ---
        # For each eigenvalue λ, solve (A - λI) v = 0.
        # Use the cross-product of two non-parallel rows of (A - λI).
        # This avoids singular matrix issues and is robust for repeated roots.
        n_elems = self.n_elems
        eye_3 = to.eye(3, device=strain.device, dtype=strain.dtype)
        eigenvectors = to.zeros((n_elems, 3, 3), dtype=strain.dtype)

        for i in range(3):
            lam = eigenvalues[:, i]  # (n_elems,)
            # (A - λI)
            A_shifted = strain - lam[:, None, None] * eye_3[None, :, :]  # (n_elems, 3, 3)

            # Find two non-parallel rows by taking cross products
            # Row0 × Row1, Row1 × Row2, Row2 × Row0 — pick the one with largest norm
            r0 = A_shifted[:, 0, :]  # (n_elems, 3)
            r1 = A_shifted[:, 1, :]  # (n_elems, 3)
            r2 = A_shifted[:, 2, :]  # (n_elems, 3)

            c01 = to.linalg.cross(r0, r1)  # (n_elems, 3)
            c12 = to.linalg.cross(r1, r2)  # (n_elems, 3)
            c20 = to.linalg.cross(r2, r0)  # (n_elems, 3)

            norms01 = to.linalg.vector_norm(c01, dim=1)
            norms12 = to.linalg.vector_norm(c12, dim=1)
            norms20 = to.linalg.vector_norm(c20, dim=1)

            # Stack and pick the cross product with largest norm per element
            candidates = to.stack([c01, c12, c20], dim=1)  # (n_elems, 3, 3)
            candidate_norms = to.stack([norms01, norms12, norms20], dim=1)  # (n_elems, 3)
            best_idx = to.argmax(candidate_norms, dim=1)  # (n_elems,)

            # Gather the best candidate per element
            batch_indices = to.arange(n_elems, device=strain.device)
            vec = candidates[batch_indices, best_idx]  # (n_elems, 3)

            # Normalize
            vec_norm = to.linalg.vector_norm(vec, dim=1, keepdim=True)
            vec_norm_safe = to.clamp(vec_norm, min=1e-30)
            vec = vec / vec_norm_safe

            eigenvectors[:, :, i] = vec

        # --- Step 3: Store to FE fields ---
        # Principal strain magnitudes as a 3-component vector (ε₁, ε₂, ε₃)
        self.principal_strains.x.array[:] = to.flatten(eigenvalues)

        # Principal directions as a 3x3 tensor (columns = eigenvectors)
        self.principal_strain_directions.x.array[:] = to.flatten(eigenvectors)

    def compute_total_strain(self) -> to.Tensor:
        """
        Project total small-strain tensor to DG0 and return as torch.

        Returns
        -------
        torch.Tensor
            Total strain per element, shape ``(n_elems, 3, 3)``.

        Notes
        -----
        Uses :func:`project` on ``ε(u)``.
        """
        self.eps_tot = project(epsilon(self.u), self.DG0_3x3)
        eps_to = numpy2torch(self.eps_tot.x.array.reshape((self.n_elems, 3, 3)))
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
            elem_ne.compute_eps_ne_k(dt * self.theta, dt * (1 - self.theta))
            eps_ne_k += elem_ne.eps_ne_k
        return eps_ne_k

    def compute_eps_ne_rate(
        self, stress: to.Tensor, dt: float, eps_tot: to.Tensor | None = None
    ) -> None:
        """
        Update non-elastic strain rate for all non-elastic elements.

        Parameters
        ----------
        stress : torch.Tensor
            Stress per element, shape ``(n_elems, 3, 3)``.
        dt : float
            Time-step size.
        eps_tot : torch.Tensor, optional
            Current-iteration total strain, shape ``(n_elems, 3, 3)``. Passed
            through (with the matching thermal strain) only to elements that
            opt in via ``uses_exact_trial = True``, letting them reconstruct
            the elastic trial stress directly from this iteration's strain
            and the start-of-step committed inelastic strain -- avoiding the
            one-iteration staleness of reconstructing it from a stored rate.

        Returns
        -------
        None
        """
        eps_th = None
        for elem_ne in self.mat.elems_ne:
            kwargs = {}
            if eps_tot is not None and getattr(elem_ne, "uses_exact_trial", False):
                if eps_th is None:
                    eps_th = self.compute_eps_th()
                kwargs["eps_tot"] = eps_tot
                kwargs["eps_th"] = eps_th
                # Initial stress (apply_initial_stress) enters the
                # formulation as an eigenstrain, so it is NOT contained in
                # eps_tot; the exact trial reconstruction must superpose it.
                sig0 = getattr(self, "sig0_to", None)
                if sig0 is not None:
                    kwargs["sig0"] = sig0
            elem_ne.compute_eps_ne_rate(
                stress, dt * self.theta, self.Temp, return_eps_ne=False, **kwargs
            )

    def update_eps_ne_rate_old(self) -> None:
        """
        Update non-elastic strain rate from the previous time step “old”.

        Returns
        -------
        None
        """
        for elem_ne in self.mat.elems_ne:
            elem_ne.update_eps_ne_rate_old()

    def update_eps_ne_old(
        self, stress: to.Tensor, stress_k: to.Tensor, dt: float
    ) -> None:
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
            elem_ne.update_eps_ne_old(stress, stress_k, dt * (1 - self.theta))

    def increment_internal_variables(
        self, stress: to.Tensor, stress_k: to.Tensor, dt: float
    ) -> None:
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
    def compute_p_elems(self) -> None:
        """
        Compute element pressure (mean stress) from the stress field.

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

    @abstractmethod
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
        pass

