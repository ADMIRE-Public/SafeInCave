# Copyright (c) 2026, The SafeInCave Developers
#
# SPDX-License-Identifier: BSD-3-Clause

"""Incremental Newton-Raphson linear momentum equation."""

from __future__ import annotations
from typing import TYPE_CHECKING
import warnings

if TYPE_CHECKING:
    from ...Mesh.Grid import GridHandlerGMSH

import numpy as np
import basix
import dolfinx as do
import ufl
from petsc4py import PETSc
import torch as to
from dolfinx.fem import petsc as fem_petsc
from .standard import LinearMomentum
from ...Utils import dotdot_ufl, dotdot_torch, epsilon, numpy2torch


class LinearMomentumNewton(LinearMomentum):
    """
    Momentum equation solved with incremental full Newton-Raphson.

    Commercial-code style solve path (Abaqus/COMSOL): per time step, the
    residual ``R(u) = f_int(σ(ε(u))) − f_ext`` is driven to zero by Newton
    corrections ``K_T δu = −R`` with the consistent algorithmic tangent
    ``D_ep`` assembled DIRECTLY into the Jacobian — bypassing the staggered
    path's ``(C⁻¹+φ₂G)⁻¹`` operator whose structure caps plastic softening
    at half the consistent tangent (and with it, the convergence rate at
    linear). Quadratic convergence: typically 3–6 iterations per step.

    Constitutive state placement (decisive for the convergence rate):

    - ``order=1`` (P1 displacement): strain is constant per cell, so stress
      and tangent live in DG0 — one state point per cell, identical to the
      staggered layout, and the Jacobian is the exact residual derivative.
    - ``order=2`` (P2, what COMSOL uses on the cavern benchmarks — avoids
      volumetric locking): strain varies within a cell, so the state lives
      at the cells' QUADRATURE POINTS (basix quadrature elements, degree 2)
      and all Newton forms integrate with the matching quadrature rule.
      With DG0 midpoint sampling the Jacobian would be inconsistent and
      Newton degrades to linear convergence (verified on cavern2D).

    In quadrature mode the constitutive batch dimension is
    ``n_state_points = n_cells × n_qp`` — size ``Material`` and its elements
    with ``eq.n_state_points`` (equal to ``n_elems`` in DG0 mode). Cell-level
    DG0 mirrors of stress/strain are maintained for outputs and logging
    (cell averages); per-cell input fields (initial stress, temperature) can
    be expanded to state points with :meth:`expand_cell_field`.

    Constitutive requirements: every non-elastic element must expose the
    incremental-Newton interface (``supports_incremental_newton``,
    ``compute_stress_and_tangent``, ``commit_increment``); at most one
    non-secondary (primary) element is supported until the creep+plasticity
    composition lands. Elastic-only materials converge in exactly one
    iteration.

    The class inherits the staggered LinearMomentum for field bookkeeping
    and output helpers; its nonlinear path is the ``constitutive_update`` /
    ``assemble_residual`` / ``solve_increment`` triple driven by
    ``Simulator_MNewton``. The Newton path REQUIRES a direct solve
    ("preonly"+"lu"; MUMPS when available) — gmres+asm stalls silently
    (KSP reason −3) on the plastic benchmarks.
    """

    def __init__(
        self,
        grid: GridHandlerGMSH,
        theta: float = 0.5,
        order: int = 1,
        state_space: str = "auto",
        quadrature_degree: int = 2,
        solver_name: str = "preonly",
        preconditioner: str = "lu",
        rtol: float = 1e-12,
        max_it: int = 500,
    ):
        super().__init__(grid, theta, solver_name, preconditioner, rtol, max_it)
        self.order = int(order)
        if state_space == "auto":
            state_space = "dg0" if self.order == 1 else "quadrature"
        if state_space not in ("dg0", "quadrature"):
            raise ValueError(f"Unknown state_space: {state_space!r}")
        self.state_space = state_space
        self.quadrature_degree = int(quadrature_degree)

        if self.order != 1:
            # Higher-order displacement space; rebuild everything derived
            # from V.
            self.V = do.fem.functionspace(
                grid.mesh, ("Lagrange", self.order, (grid.domain_dim,))
            )
            self.create_trial_test_functions()
            self.create_normal()
            self.create_solution_vector()
            self.u = do.fem.Function(self.V)

        if solver_name != "preonly":
            warnings.warn(
                "LinearMomentumNewton expects a direct solve (preonly+lu); "
                "Krylov solvers are known to stall silently on the plastic "
                "benchmarks (KSP reason -3).",
                RuntimeWarning,
            )
        else:
            try:
                self.solver.getPC().setFactorSolverType("mumps")
            except Exception:
                pass  # PETSc build without MUMPS: default LU is fine

        self._setup_state_space()
        self._newton_forms_key = None
        # Convergence state consumed by NewtonResidualCriterion.
        self._newton_state = None

    # ------------------------------------------------------------------
    # State space (where stress/tangent live)
    # ------------------------------------------------------------------

    def _setup_state_space(self) -> None:
        mesh = self.grid.mesh
        if self.state_space == "dg0":
            # Same layout as the staggered path: reuse the inherited DG0
            # fields directly.
            self.n_state_points = self.n_elems
            self.sig_state = self.sig
            self.Dep = do.fem.Function(self.DG0_6x6)
            self.dxm = self.dx
            self._strain_space = self.DG0_3x3
            self._q_cell_dofs = None
            return

        qd = self.quadrature_degree
        cell = mesh.basix_cell()

        def qspace(shape):
            elem = basix.ufl.quadrature_element(
                cell, value_shape=shape, scheme="default", degree=qd
            )
            return do.fem.functionspace(mesh, elem)

        self.Q_1 = qspace(())
        self.Q_3x3 = qspace((3, 3))
        self.Q_6x6 = qspace((6, 6))
        imap = self.Q_1.dofmap.index_map
        self.n_state_points = imap.size_local + len(imap.ghosts)
        # (n_cells, n_qp) map from cells to their state-point indices, for
        # cell-average reductions and per-cell expansions.
        self._q_cell_dofs = np.asarray(self.Q_1.dofmap.list, dtype=np.int64)

        self.sig_state = do.fem.Function(self.Q_3x3)
        self.Dep = do.fem.Function(self.Q_6x6)
        self.eps_tot_state = do.fem.Function(self.Q_3x3)
        # Torch mirrors sized by state points (apply_initial_stress
        # overwrites eps0_to; this is the no-initial-stress default).
        self.eps0_to = to.zeros((self.n_state_points, 3, 3), dtype=to.float64)
        # All Newton forms must integrate on the state points' rule.
        self.dxm = ufl.Measure(
            "dx",
            domain=mesh,
            subdomain_data=self.grid.get_subdomains(),
            metadata={"quadrature_degree": qd, "quadrature_scheme": "default"},
        )
        self._strain_space = self.Q_3x3

    def cell_average(self, values: to.Tensor) -> to.Tensor:
        """Reduce a per-state-point tensor to a per-cell average."""
        if self._q_cell_dofs is None:
            return values
        return values[self._q_cell_dofs].mean(dim=1)

    def expand_cell_field(self, values: to.Tensor) -> to.Tensor:
        """
        Expand a per-cell tensor to state points (replication). Exact for
        cell-wise constant data (e.g. a uniform initial stress); fields that
        vary within cells should be sampled at the state points instead.
        """
        if self._q_cell_dofs is None:
            return values
        out_shape = (self.n_state_points,) + tuple(values.shape[1:])
        out = to.zeros(out_shape, dtype=values.dtype)
        out[to.from_numpy(self._q_cell_dofs)] = values[:, None]
        return out

    # ------------------------------------------------------------------
    # Overrides for the quadrature-state batch size
    # ------------------------------------------------------------------

    def initialize(self) -> None:
        """Fill the DG0 elastic tangent (used by the inherited elastic
        initializer) from the material, cell-averaged in quadrature mode."""
        self.C.x.array[:] = to.flatten(self.cell_average(self.mat.C))

    def build_body_force(self, g: list) -> None:
        if self.state_space == "dg0":
            super().build_body_force(g)
            return
        density = do.fem.Function(self.DG0_1)
        density.x.array[:] = self.cell_average(
            to.as_tensor(self.mat.density, dtype=to.float64)
        ).numpy()
        body_force = density * do.fem.Constant(
            self.grid.mesh, do.default_scalar_type(tuple(g))
        )
        self.b_body = ufl.dot(body_force, self.u_) * self.dx

    def apply_initial_stress(self, sig0: to.Tensor) -> None:
        """
        Set the initial stress (sized ``n_state_points``; a per-cell field
        can be expanded with :meth:`expand_cell_field`).
        """
        if sig0.shape[0] == self.n_elems and self.n_state_points != self.n_elems:
            sig0 = self.expand_cell_field(sig0)
        self.eps0_to = dotdot_torch(self.mat.C_inv, sig0)
        self.sig0_to = sig0.clone()
        self.sig_state.x.array[:] = to.flatten(sig0)
        sig0_cell = self.cell_average(sig0)
        self.eps_0.x.array[:] = to.flatten(self.cell_average(self.eps0_to))
        self.sig.x.array[:] = to.flatten(sig0_cell)
        self.sig0.x.array[:] = to.flatten(sig0_cell)

    def compute_total_strain(self) -> to.Tensor:
        """Strain at the state points (cached Expression); the inherited DG0
        ``eps_tot`` mirror is refreshed with cell averages for outputs."""
        if self.state_space == "dg0":
            return super().compute_total_strain()
        if getattr(self, "_eps_state_expr_key", None) != id(self.u):
            self._eps_state_expr = do.fem.Expression(
                epsilon(self.u), self._strain_space.element.interpolation_points()
            )
            self._eps_state_expr_key = id(self.u)
        self.eps_tot_state.interpolate(self._eps_state_expr)
        eps_to = numpy2torch(
            self.eps_tot_state.x.array.reshape((self.n_state_points, 3, 3))
        )
        self.eps_tot.x.array[:] = to.flatten(self.cell_average(eps_to))
        return eps_to

    def compute_eps_th(self) -> to.Tensor:
        """Thermal strain at the state points (thermoelastic elements are
        sized ``n_state_points`` on this path)."""
        eps_th = to.zeros((self.n_state_points, 3, 3), dtype=to.float64)
        deltaT = self.Temp - self.T0
        for elem_th in self.mat.elems_th:
            elem_th.compute_eps_th(deltaT)
            eps_th += elem_th.eps_th
        return eps_th

    def compute_elastic_stress(self, eps_e: to.Tensor) -> to.Tensor:
        stress_to = dotdot_torch(self.mat.C, eps_e)
        self.sig_state.x.array[:] = to.flatten(stress_to)
        self.sig.x.array[:] = to.flatten(self.cell_average(stress_to))
        return stress_to

    def compute_yield_mode(self) -> None:
        if self.state_space == "dg0":
            super().compute_yield_mode()
            return
        if not hasattr(self, "mat"):
            return
        combined = None
        for elem in getattr(self.mat, "elems_ne", []):
            mode = getattr(elem, "yield_mode", None)
            if mode is None:
                continue
            combined = mode if combined is None else to.maximum(combined, mode)
        if combined is None:
            return
        cell_mode = combined[to.from_numpy(self._q_cell_dofs)].amax(dim=1)
        self.yield_mode_elems.x.array[:] = cell_mode.to(to.float64).numpy()

    # ------------------------------------------------------------------
    # Newton machinery
    # ------------------------------------------------------------------

    def init_newton(self) -> None:
        """
        Compile the Jacobian and residual forms once and allocate persistent
        PETSc objects. All time/iteration-varying data (Dep, sig, BC loads)
        live in in-place-updated Function/Constant objects; the identity key
        rebuilds only if the BC/cavern term set itself changes.
        """
        key = self._bc_forms_key()
        if self._newton_forms_key == key:
            return
        self.du_sol = do.fem.Function(self.V)
        J = (
            ufl.inner(dotdot_ufl(self.Dep, epsilon(self.du)), epsilon(self.u_))
            * self.dxm
        )
        self._J_form = do.fem.form(J)
        self._fint_form = do.fem.form(
            ufl.inner(self.sig_state, epsilon(self.u_)) * self.dxm
        )
        self._fext_form = do.fem.form(
            self.b_body + sum(self.bc.neumann_bcs) + sum(self.bc.cavern_bcs)
        )
        self._A_newton = fem_petsc.create_matrix(self._J_form)
        self._fint_vec = fem_petsc.create_vector(self._fint_form)
        self._fext_vec = fem_petsc.create_vector(self._fext_form)
        self._b_newton = fem_petsc.create_vector(self._fint_form)
        self.solver.setOperators(self._A_newton)
        self._newton_forms_key = key

    def _newton_elements(self) -> tuple:
        """
        Partition the constitutive elements for the Newton path.

        Returns ``(primary, creep_elems, all_elems)``:

        - ``primary``: the (at most one) non-secondary element implementing
          the incremental-Newton interface (return-mapped stress + D_ep);
        - ``creep_elems``: legacy rate-based elements (creep/viscoplastic),
          composed staggered-within-Newton via the theta scheme — their
          viscous increment at the previous stress iterate is subtracted
          from the trial and their tangent G folds into the Jacobian as
          ``(D_ep⁻¹ + φ₂·G_visc)⁻¹`` (superlinear instead of quadratic,
          costing ~1-2 extra iterations).
        """
        elems = list(getattr(self.mat, "elems_ne", []))
        capable = [
            e for e in elems if getattr(e, "supports_incremental_newton", False)
        ]
        creep_elems = [
            e for e in elems if not getattr(e, "supports_incremental_newton", False)
        ]
        primaries = [
            e for e in capable if not getattr(e, "newton_secondary", False)
        ]
        if len(primaries) > 1:
            raise NotImplementedError(
                "Multiple primary Newton-capable elements are not supported; "
                "found: "
                + ", ".join(getattr(e, "name", type(e).__name__) for e in primaries)
            )
        return (primaries[0] if primaries else None), creep_elems, elems

    def constitutive_update(self, dt: float) -> tuple:
        """
        One constitutive pass at the current displacement iterate:
        build the elastic trial strain
        ``ε_el = ε₀ + ε(u) − ε_th − Σ ε_ne_old`` (initial stress enters as
        the eigenstrain ε₀ = C⁻¹:σ₀) at the state points, call the primary
        element's return mapping, and write stress and consistent tangent
        into the FE fields used by the residual/Jacobian forms.

        Returns ``(stress, eps_tot)`` as CELL-LEVEL torch tensors (averages
        in quadrature mode) for logging/outputs.

        Creep composition (staggered-within-Newton, tangent-consistent):
        rate-based elements contribute a theta-scheme viscous increment
        evaluated at the PREVIOUS stress iterate, subtracted from the trial
        before the plastic return; their tangents fold into the Jacobian as
        ``(D_ep⁻¹ + φ₂·G_visc)⁻¹`` (assembled singular-safely).
        """
        primary, creep_elems, elems = self._newton_elements()
        eps_tot = self.compute_total_strain()
        eps_th = self.compute_eps_th()
        eps_ne = to.zeros_like(eps_tot)
        for elem in elems:
            eps_ne = eps_ne + elem.eps_ne_old
        eps_el_trial = self.eps0_to + eps_tot - eps_th - eps_ne

        phi1 = dt * self.theta
        phi2 = dt * (1.0 - self.theta)
        sigma_k = getattr(self, "_sigma_k", None)
        if creep_elems:
            if sigma_k is None:
                sigma_k = dotdot_torch(self.mat.C, eps_el_trial)
            # Viscous increment at the previous stress iterate.
            eps_visc = to.zeros_like(eps_tot)
            for elem in creep_elems:
                elem.compute_eps_ne_rate(sigma_k, phi1, self.Temp)
                eps_visc = eps_visc + phi1 * elem.eps_ne_rate_old + phi2 * elem.eps_ne_rate
            eps_el_eff = eps_el_trial - eps_visc
        else:
            eps_el_eff = eps_el_trial

        if primary is not None:
            sigma, D_ep = primary.compute_stress_and_tangent(
                eps_el_eff, self.Temp, dt
            )
        else:
            sigma = dotdot_torch(self.mat.C, eps_el_eff)
            D_ep = self.mat.C

        if creep_elems:
            # Jacobian composition (D_ep^-1 + phi2*G_visc)^-1, written as
            # (I + D_ep phi2 G_visc)^-1 D_ep so a singular D_ep (apex) is
            # handled without inversion.
            G_visc = to.zeros_like(D_ep)
            for elem in creep_elems:
                elem.compute_G_B(sigma_k, dt, self.theta, self.Temp)
                G_visc = G_visc + elem.G
            eye = to.eye(6, dtype=to.float64).expand_as(D_ep)
            D_ep = to.linalg.solve(eye + phi2 * (D_ep @ G_visc), D_ep)

        self._sigma_k = sigma.detach().clone()
        self.sig_state.x.array[:] = to.flatten(sigma)
        self.Dep.x.array[:] = to.flatten(D_ep)
        if self.state_space != "dg0":
            self.sig.x.array[:] = to.flatten(self.cell_average(sigma))
        return self.cell_average(sigma), self.cell_average(eps_tot)

    def initialize_rates(self) -> None:
        """
        Seed the creep elements' rate carryover (``eps_ne_rate_old``) at the
        current (initial) stress — matching the staggered simulator's
        startup, so the first step's trapezoidal viscous increment carries
        the correct theta-weighted initial rate.
        """
        _, creep_elems, _ = self._newton_elements()
        if not creep_elems:
            return
        sigma = numpy2torch(
            self.sig_state.x.array.reshape((self.n_state_points, 3, 3))
        )
        for elem in creep_elems:
            elem.compute_eps_ne_rate(sigma, 0.0, self.Temp)
            elem.update_eps_ne_rate_old()

    def commit_step(self, dt: float) -> None:
        """
        Commit the accepted step's inelastic increments: creep elements
        first (final rate at the converged stress, trapezoidal increment
        into ``eps_ne_old``, rate carried over), then the Newton-capable
        elements' ``commit_increment``.
        """
        primary, creep_elems, elems = self._newton_elements()
        phi1 = dt * self.theta
        phi2 = dt * (1.0 - self.theta)
        if creep_elems:
            sigma = numpy2torch(
                self.sig_state.x.array.reshape((self.n_state_points, 3, 3))
            )
            for elem in creep_elems:
                elem.compute_eps_ne_rate(sigma, phi1, self.Temp)
                elem.eps_ne_old = (
                    elem.eps_ne_old
                    + phi1 * elem.eps_ne_rate_old
                    + phi2 * elem.eps_ne_rate
                )
                elem.eps_ne_k = elem.eps_ne_old.clone()
                elem.update_internal_variables()
                elem.update_eps_ne_rate_old()
        for elem in elems:
            elem.commit_increment()

    def _dirichlet_dofs(self) -> np.ndarray:
        """Local DOF indices constrained by Dirichlet BCs."""
        indices = []
        for bc in self.bc.dirichlet_bcs:
            dofs = bc.dof_indices()
            if isinstance(dofs, tuple):
                dofs = dofs[0]
            indices.append(np.asarray(dofs, dtype=np.int64))
        if not indices:
            return np.empty(0, dtype=np.int64)
        return np.unique(np.concatenate(indices))

    def _global_norm(self, values: np.ndarray) -> float:
        """MPI-safe L2 norm of locally owned entries."""
        local_sq = float(np.dot(values, values))
        return float(self.grid.mesh.comm.allreduce(local_sq) ** 0.5)

    def assemble_residual(self) -> tuple:
        """
        Assemble the Newton residual with the canonical dolfinx incremental
        Dirichlet pattern:

            b = f_int − f_ext
            apply_lifting(b, [J], [bcs], x0=[u], alpha=−1)
            set_bc(b, bcs, u, −1)

        so that solving ``J δ = b`` and updating ``u ← u − δ`` lands u
        exactly on u_D(t) in the first iteration (later corrections are
        homogeneous at constrained DOFs).

        Returns ``(r_norm, ref_norm)``: the free-DOF residual norm and the
        reference force scale ``max(free-DOF ‖f_ext‖, ‖f_int‖ incl.
        reactions)`` (same normalization as the force_displacement
        criterion).
        """
        self.init_newton()

        with self._fint_vec.localForm() as v:
            v.set(0.0)
        fem_petsc.assemble_vector(self._fint_vec, self._fint_form)
        self._fint_vec.ghostUpdate(
            addv=PETSc.InsertMode.ADD_VALUES, mode=PETSc.ScatterMode.REVERSE
        )
        with self._fext_vec.localForm() as v:
            v.set(0.0)
        fem_petsc.assemble_vector(self._fext_vec, self._fext_form)
        self._fext_vec.ghostUpdate(
            addv=PETSc.InsertMode.ADD_VALUES, mode=PETSc.ScatterMode.REVERSE
        )

        b = self._b_newton
        self._fint_vec.copy(b)
        b.axpy(-1.0, self._fext_vec)

        fem_petsc.apply_lifting(
            b,
            [self._J_form],
            [self.bc.dirichlet_bcs],
            x0=[self.u.x.petsc_vec],
            alpha=-1.0,
        )
        b.ghostUpdate(
            addv=PETSc.InsertMode.ADD_VALUES, mode=PETSc.ScatterMode.REVERSE
        )
        fem_petsc.set_bc(b, self.bc.dirichlet_bcs, self.u.x.petsc_vec, -1.0)
        b.ghostUpdate(
            addv=PETSc.InsertMode.INSERT, mode=PETSc.ScatterMode.FORWARD
        )

        dirichlet = self._dirichlet_dofs()
        n_owned = self.V.dofmap.index_map.size_local * self.V.dofmap.index_map_bs
        r = b.array[:n_owned].copy()
        f_ext = self._fext_vec.array[:n_owned].copy()
        f_int = self._fint_vec.array[:n_owned]
        owned_dirichlet = dirichlet[dirichlet < n_owned]
        if owned_dirichlet.size:
            r[owned_dirichlet] = 0.0
            f_ext[owned_dirichlet] = 0.0
        r_norm = self._global_norm(r)
        ref_norm = max(
            self._global_norm(f_ext), self._global_norm(np.asarray(f_int)), 1e-16
        )
        return r_norm, ref_norm

    def solve_increment(self, apply: bool = True) -> float:
        """
        Assemble the Jacobian and solve one Newton correction ``J δ = b``
        into ``du_sol``; when ``apply`` is True update ``u ← u − δ``
        (dolfinx sign convention, matching assemble_residual) — pass False
        when a line search will choose the step length. A failed direct
        solve raises — on the Newton path it must trigger the driver's
        time-step cutback, not limp on.

        Returns ``‖δu‖`` (MPI-safe).
        """
        A = self._A_newton
        A.zeroEntries()
        fem_petsc.assemble_matrix(A, self._J_form, bcs=self.bc.dirichlet_bcs)
        A.assemble()

        self.solver.solve(self._b_newton, self.du_sol.x.petsc_vec)
        reason = self.solver.getConvergedReason()
        if reason < 0:
            raise RuntimeError(
                f"Newton linear solve failed (KSP reason={reason}); "
                "the step should be retried with a smaller dt."
            )
        self.du_sol.x.scatter_forward()

        if apply:
            self.u.x.array[:] -= self.du_sol.x.array
            self.u.x.scatter_forward()

        n_owned = self.V.dofmap.index_map.size_local * self.V.dofmap.index_map_bs
        return self._global_norm(self.du_sol.x.array[:n_owned].copy())
