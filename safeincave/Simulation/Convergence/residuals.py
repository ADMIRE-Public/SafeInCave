# Copyright (c) 2026, The SafeInCave Developers
from __future__ import annotations
#
# SPDX-License-Identifier: BSD-3-Clause

"""Shared numeric helpers for convergence criteria (force residuals, internal force vectors, norms)."""

from typing import Any,  Optional, Dict, Any
import torch as to
import numpy as np
from mpi4py import MPI
from petsc4py import PETSc
import ufl
from dolfinx import fem as do_fem
from dolfinx.fem import petsc as fem_petsc
from ...Utils import epsilon

def _compute_external_load_vector_norm(momentum_eq: Any) -> float:
    """
    Compute L2 norm of external load vector for residual normalization.

    Assembles and computes the L2 norm of all external loads (body forces,
    Neumann boundary conditions, and cavern loads) to use as a reference scale
    for residual error computation in convergence checks.

    Parameters
    ----------
    momentum_eq : LinearMomentumBase
        Momentum equation solver instance (provides FEniCS fields and BCs).

    Returns
    -------
    float
        L2 norm of the assembled external load vector. Protected against zero by
        floor of 1e-16 to enable safe normalization in error metrics.

    Notes
    -----
    Collects contributions from:
    - Body forces: ∫ ρ g · u_ dx
    - Neumann boundary conditions: natural tractions
    - Cavern-specific loads

    This function is MPI-safe (uses collective norm operations).
    """
    # Build the external load vector
    linear_form = do_fem.form(
        momentum_eq.b_body + sum(momentum_eq.bc.neumann_bcs) + sum(momentum_eq.bc.cavern_bcs)
    )
    f_ext = fem_petsc.assemble_vector(linear_form)
    f_ext.ghostUpdate(addv=PETSc.InsertMode.ADD, mode=PETSc.ScatterMode.REVERSE)

    # Zero the rows at Dirichlet-constrained DOFs: loads applied directly at
    # constrained DOFs are carried by reactions and play no role in the
    # free-DOF equilibrium the residual criterion measures.
    loads = f_ext.array.copy()
    dirichlet_dofs = _dirichlet_dof_indices(momentum_eq)
    if dirichlet_dofs.size:
        loads[dirichlet_dofs] = 0.0

    norm = float(np.linalg.norm(loads))

    # Protect against zero norm
    norm = max(norm, 1e-16)

    return norm


def _dirichlet_dof_indices(momentum_eq: Any) -> np.ndarray:
    """
    Collect the (locally owned) DOF indices constrained by Dirichlet BCs.

    Used to exclude constrained DOFs from force-residual norms: the
    out-of-balance force at a constrained DOF is the reaction force, not a
    convergence defect.
    """
    indices = []
    for bc in momentum_eq.bc.dirichlet_bcs:
        dofs = bc.dof_indices()
        # dolfinx returns (indices, num_owned); accept a bare array too.
        if isinstance(dofs, tuple):
            dofs = dofs[0]
        indices.append(np.asarray(dofs, dtype=np.int64))
    if not indices:
        return np.empty(0, dtype=np.int64)
    return np.unique(np.concatenate(indices))


def _compute_internal_force_vector(
    momentum_eq: Any, stress_field: to.Tensor
) -> to.Tensor:
    """
    Assemble internal force vector from stress field.

    Assembles the internal force vector q_int = ∫ σ : ∇u dx by computing
    the divergence of the stress tensor and integrating over the domain.

    Parameters
    ----------
    momentum_eq : LinearMomentumBase
        Momentum equation solver instance (provides FEniCS fields and mesh).
    stress_field : torch.Tensor
        Stress tensor per element, shape ``(n_elems, 3, 3)``.

    Returns
    -------
    torch.Tensor
        Internal force vector, same shape/structure as displacement DoFs.

    Notes
    -----
    Uses the current stiffness matrix assembly pattern without time stepping:
    q_int = ∫ σ : ∇u_ dx

    For efficiency, stores stress in the DG0 stress field (momentum_eq.sig) before
    assembly to ensure consistent integration.
    """
    # Update stress field for assembly
    momentum_eq.sig.x.array[:] = to.flatten(stress_field)

    # Build internal force form: ∫ σ : ∇u_ dx
    # This represents the internal force contribution from current stress state
    internal_force_form = ufl.inner(momentum_eq.sig, epsilon(momentum_eq.u_)) * momentum_eq.dx

    # Assemble as a vector (right-hand side)
    assembled_form = do_fem.form(internal_force_form)
    internal_force_vec = fem_petsc.assemble_vector(assembled_form)
    internal_force_vec.ghostUpdate(addv=PETSc.InsertMode.ADD, mode=PETSc.ScatterMode.REVERSE)

    # Convert to torch tensor for consistency with other methods
    internal_force_tensor = to.from_numpy(internal_force_vec.array.copy()).double()

    return internal_force_tensor


def _compute_force_residual(
    momentum_eq: Any, stress_field: to.Tensor
) -> to.Tensor:
    """
    Compute mechanical force residual vector R = P_ext - q_int.

    Computes the out-of-balance force residual as the difference between applied
    external loads (P_ext) and internal forces (q_int) derived from the current
    stress state. Used to check mechanical equilibrium in convergence checks.

    Parameters
    ----------
    momentum_eq : LinearMomentumBase
        Momentum equation solver instance.
    stress_field : torch.Tensor
        Current stress state, shape ``(n_elems, 3, 3)``.
    Returns
    -------
    tuple[torch.Tensor, float]
        Residual vector R = P_ext - q_int (Dirichlet rows zeroed), and the
        reference force norm to normalize it by (max of free-DOF external
        load norm and full internal-force norm including reactions).

    Notes
    -----
    Residual computation:
    - P_ext: assembled from external loads (body forces, Neumann BCs, cavern loads)
    - q_int: internal forces computed from divergence of stress
    - R = P_ext - q_int (small R indicates equilibrium)

    For a converged solution, both R should be small (residual criterion)
    and the displacement increment should be small (strain criterion).
    """
    # Assemble external load vector
    external_load_form = do_fem.form(
        momentum_eq.b_body + sum(momentum_eq.bc.neumann_bcs) + sum(momentum_eq.bc.cavern_bcs)
    )
    external_loads_vec = fem_petsc.assemble_vector(external_load_form)
    external_loads_vec.ghostUpdate(addv=PETSc.InsertMode.ADD, mode=PETSc.ScatterMode.REVERSE)

    # Get internal force vector
    internal_forces_tensor = _compute_internal_force_vector(momentum_eq, stress_field)

    # Convert external loads to torch tensor for residual computation
    external_loads_tensor = to.from_numpy(external_loads_vec.array.copy()).double()

    # Compute residual: R = P_ext - q_int (out-of-balance forces)
    residual_vector = external_loads_tensor - internal_forces_tensor

    # Zero the residual at Dirichlet-constrained DOFs: the imbalance there
    # is the reaction force, which equilibrium never drives to zero and
    # which must not count against convergence.
    dirichlet_dofs = _dirichlet_dof_indices(momentum_eq)
    if dirichlet_dofs.size:
        residual_vector[dirichlet_dofs] = 0.0

    # Reference force scale for normalizing the residual. The free-DOF
    # external load alone is a bad reference: it is ~zero for
    # displacement-controlled problems AND for fully confined load cases
    # (e.g. a geostatic stage with rollers on every loaded boundary), where
    # the entire applied load is carried directly by reactions. The full
    # internal force vector (constrained rows INCLUDED, i.e. reactions)
    # always carries the physical force scale of the system, so take the
    # larger of the two.
    free_loads = external_loads_tensor.clone()
    if dirichlet_dofs.size:
        free_loads[dirichlet_dofs] = 0.0
    reference_norm = max(
        float(to.linalg.norm(free_loads)),
        float(to.linalg.norm(internal_forces_tensor)),
    )

    return residual_vector, reference_norm


def _compute_vector_norm(vector: to.Tensor) -> float:
    """
    Compute L2 norm of a tensor vector.

    Generic utility to compute the L2 norm of any displacement, residual,
    or other vector field.

    Parameters
    ----------
    vector : torch.Tensor
        Vector to compute norm for (displacement, residual, etc.).

    Returns
    -------
    float
        L2 norm of the vector. Protected against zero by floor of 1e-16
        to enable safe relative error computation (normalization).

    Notes
    -----
    This norm is used in convergence criteria:
    - error_strain = ||u_new - u_old|| / ||u_old||
    - error_displacement = ||u_correction|| / ||u_total||

    For torch tensors, computes torch.norm(u).item() and converts to float.
    """
    if isinstance(vector, to.Tensor):
        computed_norm = to.norm(vector).item()
    else:
        # Handle numpy arrays if needed
        computed_norm = np.linalg.norm(vector)

    # Protect against zero norm to prevent division by zero
    computed_norm = max(computed_norm, 1e-16)

    return float(computed_norm)


def _initialize_step_displacement(momentum_eq: Any) -> to.Tensor:
    """
    Capture and store displacement field at time step start.

    Captures the displacement state at the beginning of a time step to enable
    computation of the cumulative displacement increment throughout the step.
    Used by displacement correction criterion for robust convergence checking.

    Parameters
    ----------
    momentum_eq : LinearMomentumBase
        Momentum equation solver instance.

    Returns
    -------
    torch.Tensor
        Copy of displacement array at time step start.

    Notes
    -----
    - Must be called ONCE per time step, before nonlinear iteration loop
    - Stores displacement as a FEniCS vector array copy (numpy format)
    - Used internally by DisplacementIncrementCriterion
    - MPI-safe: each rank stores its own DOFs plus ghost elements
    """
    # Capture displacement state at step start for total increment tracking
    displacement_step_start = momentum_eq.u.x.array.copy()
    momentum_eq.u_step_start = displacement_step_start
    return displacement_step_start


def resolve_convergence_criterion(
    convergence_criterion: str | "ConvergenceCriterion" = "strain_based",
) -> "ConvergenceCriterion":
    """
    Resolve convergence criterion selector into concrete strategy.

    Supported names:
    - "strain_based"
    - "force_residual"
    - "displacement_increment"
    - "force_displacement" (combined force residual + displacement increment)

    Parameters
    ----------
    convergence_criterion : str
        Selector name.

    Returns
    -------
    ConvergenceCriterion
        Concrete criterion strategy instance.
    """
    if hasattr(convergence_criterion, "compute_error") and hasattr(
        convergence_criterion, "is_converged"
    ):
        return convergence_criterion

    criterion_key = str(convergence_criterion).strip().lower()

    if criterion_key == "strain_based":
        return StrainBasedCriterion()
    if criterion_key == "force_residual":
        return ForceResidualCriterion()
    if criterion_key == "displacement_increment":
        return DisplacementIncrementCriterion()
    if criterion_key == "force_displacement":
        return ForceDisplacementCriterion()

    raise ValueError(
        "Unknown convergence_criterion. Supported values: "
        "'strain_based', 'force_residual', "
        "'displacement_increment', 'force_displacement'."
    )


def initialize_convergence_state(
    momentum_eq: Any,
    convergence_criterion: "ConvergenceCriterion",
) -> None:
    """
    Initialize criterion state at time-step start.

    Parameters
    ----------
    momentum_eq : LinearMomentumBase
        Momentum equation solver instance.
    convergence_criterion : Any
        Resolved criterion strategy.
    """
    convergence_criterion.initialize(momentum_eq)


def compute_error_from_criterion(
    momentum_eq: Any,
    convergence_criterion: "ConvergenceCriterion",
) -> float:
    """
    Compute convergence error from the selected criterion strategy.

    Parameters
    ----------
    momentum_eq : LinearMomentumBase
        Momentum equation solver instance.
    convergence_criterion : Any
        Resolved criterion strategy.

    Returns
    -------
    float
        Raw criterion error value.
    """
    scalar_error = float(convergence_criterion.compute_error(momentum_eq))
    return scalar_error


class ConvergenceErrorHandler:
    """
    Facade for convergence strategy selection and error evaluation.

    This class centralizes all convergence-related selection/dispatch logic so
    simulator drivers only need one import from this module.
    """

    def __init__(
        self,
        convergence_criterion: str | "ConvergenceCriterion" = "strain_based",
        tol: Optional[float] = None,
        plastic_consistency_tolerance: float = 1e-4,
    ):
        self.criterion = resolve_convergence_criterion(convergence_criterion)
        self.error: float = 0.0
        self.not_converged_error: bool = True
        self.below_max_iterations: bool = True
        self.last_raw_error: float = 0.0
        # If tol is None, convergence is evaluated with criterion.is_converged(error).
        # If tol is provided, it overrides criterion tolerance for loop control.
        self.tol: Optional[float] = tol
        self.ite: int = 0
        self.maxiter: Optional[int] = None
        # Independent equilibrium gate, applied regardless of the named
        # criterion: a plastic material's yield-consistency residual
        # (max(F)+/f_c across elems_ne exposing `consistency_error`) must
        # also be small before the step is accepted. Without this, a metric
        # like strain-based can be satisfied trivially on a small load
        # increment while the stress still violates the yield surface by a
        # non-negligible margin. Tolerance calibration on the cavern2D
        # benchmark: 1e-3 lets per-step residuals compound to a 10%+ bias
        # over a long load ramp (do not loosen); 1e-4 is validated stable
        # across the full ramp (~3.4% vs COMSOL, plateaued, ~20-25
        # iters/step); 1e-6 adds ~1.5-2x more iterations with no measurable
        # accuracy gain. Kept out of individual criteria so their reported
        # error stays a truthful, single-purpose metric -- this is a safety
        # net on top, not folded into any one criterion's number.
        self.plastic_consistency_tolerance = plastic_consistency_tolerance
        self.plastic_consistency_error: float = 0.0

    def initialize_step(
        self,
        momentum_eq: Any,
        maxiter: Optional[int] = None,
    ) -> None:
        """Initialize criterion state at the start of a time step."""
        self.error = 0.0
        self.not_converged_error = True
        self.ite = 0
        self.maxiter = maxiter
        self.below_max_iterations = True if maxiter is None else (self.ite < maxiter)
        initialize_convergence_state(
            momentum_eq,
            self.criterion,
        )

    def compute_error(self, momentum_eq: Any) -> float:
        """Compute raw criterion error and store for diagnostics."""
        raw_error = compute_error_from_criterion(momentum_eq, self.criterion)
        self.error = float(raw_error)
        self.last_raw_error = self.error

        local_consistency = 0.0
        for elem_ne in getattr(momentum_eq.mat, "elems_ne", []):
            local_consistency = max(
                local_consistency, float(getattr(elem_ne, "consistency_error", 0.0))
            )
        self.plastic_consistency_error = float(
            momentum_eq.grid.mesh.comm.allreduce(local_consistency, op=MPI.MAX)
        )

        if os.environ.get("SAFEINCAVE_DEBUG_CONV"):
            print(
                f"[conv] ite={self.ite:3d} error={self.error:.3e} "
                f"consistency={self.plastic_consistency_error:.3e}",
                flush=True,
            )

        return self.error

    def evaluate(
        self,
        momentum_eq: Any,
        ite: Optional[int] = None,
        maxiter: Optional[int] = None,
    ) -> bool:
        """Compute error and refresh the convergence-state boolean."""
        self.compute_error(momentum_eq)
        if ite is not None or maxiter is not None:
            self.update_below_max_iterations(ite=ite, maxiter=maxiter)
        return self.update_not_converged_error()

    def update_not_converged_error(self) -> bool:
        """Update and return the current convergence-state boolean."""
        if not math.isfinite(self.error) or not math.isfinite(self.plastic_consistency_error):
            # NaN/Inf means the nonlinear iteration is unstable and must continue/retry.
            self.not_converged_error = True
        elif self.tol is None:
            self.not_converged_error = not self.criterion.is_converged(self.error)
        else:
            self.not_converged_error = self.error > self.tol

        if not self.not_converged_error:
            # Equilibrium gate: don't accept a step where the criterion's
            # metric (e.g. strain change) looks converged but a plastic
            # material's yield-consistency residual is still large.
            if self.plastic_consistency_error > self.plastic_consistency_tolerance:
                self.not_converged_error = True

        return self.not_converged_error

    def update_below_max_iterations(
        self,
        ite: Optional[int] = None,
        maxiter: Optional[int] = None,
    ) -> bool:
        """Update and return whether the current iteration is below max iterations."""
        if maxiter is not None:
            self.maxiter = maxiter
        if ite is not None:
            self.ite = ite

        if self.maxiter is None:
            self.below_max_iterations = True
        else:
            self.below_max_iterations = self.ite < self.maxiter

        return self.below_max_iterations

    def increment_iteration(self) -> int:
        """Increment iteration counter and refresh max-iteration state."""
        self.ite += 1
        self.update_below_max_iterations(ite=self.ite)
        return self.ite

    def not_converged(self) -> bool:
        """Return unified nonlinear-loop continuation condition."""
        return self.not_converged_error and self.below_max_iterations

    def get_tolerance(self) -> float:
        """Return active criterion tolerance used for convergence checks."""
        if self.tol is not None:
            return float(self.tol)
        return float(self.criterion.tolerance)

    def get_last_raw_error(self) -> float:
        """Return most recent raw criterion value (not zeroed by convergence logic)."""
        return self.last_raw_error

    def set_criterion(self, convergence_criterion: str) -> None:
        """Swap convergence strategy at runtime."""
        self.criterion = resolve_convergence_criterion(convergence_criterion)


def _is_momentum_solver_instance(obj: Any) -> bool:
    """Return True if object provides the momentum solver interface used here."""
