# Copyright (c) 2026, The SafeInCave Developers
#
# SPDX-License-Identifier: BSD-3-Clause

"""
Newton-step solver for momentum equations: assembles and solves one correction
via K·δu = -R (preonly + LU factorization).

Adapted from the dolfinx-external-operator tutorial's ``solvers.py``
(MIT, A. Latyshev & J. Hale).

Generic dolfinx/PETSc plumbing with no dependency on any particular solver
backend (JAX, torch, ...) — shared by any Newton-type momentum solve.
"""

from __future__ import annotations

import ufl
from petsc4py import PETSc

from dolfinx import fem


class MomentumNewtonSolver:
    """Assemble-and-solve one Newton correction step for momentum equations.
    
    Solves the linearized system K·δu = -R where:
    - K is the Jacobian of the residual (dR)
    - R is the momentum residual
    - u is the current displacement iterate
    - δu is the correction increment
    
    Uses PETSc KSP with preonly+LU (MUMPS when available).
    """

    def __init__(self, dR: ufl.Form, R: ufl.Form, u: fem.Function, bcs=None):
        """Initialize Newton solver for one momentum correction.
        
        Parameters
        ----------
        dR : ufl.Form
            Jacobian form (derivative of residual w.r.t. displacement).
        R : ufl.Form
            Momentum residual form.
        u : fem.Function
            Current displacement iterate (used for BC lifting).
        bcs : list[DirichletBC], optional
            Boundary conditions (default: none).
        """
        self.u = u
        self.bcs = bcs or []
        self.b_form = fem.form(R)
        self.A_form = fem.form(dR)
        self.b = fem.petsc.create_vector(self.b_form)
        self.A = fem.petsc.create_matrix(self.A_form)
        self.comm = u.function_space.mesh.comm
        self.solver = PETSc.KSP().create(self.comm)
        self.solver.setType("preonly")
        self.solver.getPC().setType("lu")
        # MUMPS factorizes noticeably faster than PETSc's built-in serial LU
        # on the quadrature-space P2 systems this solver assembles (see
        # FullNewtonRaphsonSolver's LinearMomentumNewton, which sets the same
        # option). Falls back to the default factorizer on PETSc builds
        # without MUMPS.
        try:
            self.solver.getPC().setFactorSolverType("mumps")
        except Exception:
            pass  # PETSc build without MUMPS: default LU is fine
        self.solver.setOperators(self.A)

    def assemble_vector(self) -> None:
        """
        Assemble the Newton right-hand side with BC lifting relative to the
        current iterate ``u`` (the step increment): after the correction,
        BC dofs of ``u`` equal the (possibly inhomogeneous) BC values and
        subsequent corrections there vanish.
        """
        with self.b.localForm() as b_local:
            b_local.set(0.0)
        fem.petsc.assemble_vector(self.b, self.b_form)
        fem.petsc.apply_lifting(
            self.b, [self.A_form], [self.bcs], x0=[self.u.x.petsc_vec], alpha=1.0
        )
        self.b.ghostUpdate(addv=PETSc.InsertMode.ADD, mode=PETSc.ScatterMode.REVERSE)
        fem.petsc.set_bc(self.b, self.bcs, x0=self.u.x.petsc_vec, alpha=1.0)

    def assemble_matrix(self) -> None:
        """Assemble the Newton Jacobian matrix (system matrix K)."""
        self.A.zeroEntries()
        fem.petsc.assemble_matrix(self.A, self.A_form, bcs=self.bcs)
        self.A.assemble()

    def solve(self, du: fem.Function) -> None:
        """Solve K·δu = -R for the correction increment du."""
        self.solver.solve(self.b, du.x.petsc_vec)

    def __del__(self):
        self.solver.destroy()
        self.A.destroy()
        self.b.destroy()
