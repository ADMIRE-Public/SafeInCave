# Copyright (c) 2026, The SafeInCave Developers
#
# SPDX-License-Identifier: BSD-3-Clause

from __future__ import annotations
from typing import Any
import torch as to
import numpy as np
from petsc4py import PETSc
import ufl
from dolfinx import fem as do_fem
from dolfinx.fem import petsc as fem_petsc
from ...Utils import epsilon

"""Shared numeric helpers for convergence criteria (force residuals, internal force vectors, norms)."""

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
    f_ext = _assemble_external_load_vector(momentum_eq)

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


def _bc_terms_key(momentum_eq: Any) -> tuple:
    """
    Identity key of the BC terms referenced by the cached residual forms.

    BC loads live in in-place-updated fem.Constants (see BcHandler), so the
    compiled forms stay valid while the term objects are unchanged; a change
    of the registered BC or cavern set produces new term objects and must
    invalidate the cache.
    """
    return (
        tuple(map(id, momentum_eq.bc.neumann_bcs)),
        tuple(map(id, getattr(momentum_eq.bc, "cavern_bcs", []))),
    )


def _assemble_external_load_vector(momentum_eq: Any):
    """
    Assemble the external load vector (body + Neumann + cavern) into a
    cached PETSc vector using a form compiled once per BC-term set.
    """
    key = _bc_terms_key(momentum_eq)
    if getattr(momentum_eq, "_resid_ext_key", None) != key:
        form = do_fem.form(
            momentum_eq.b_body
            + sum(momentum_eq.bc.neumann_bcs)
            + sum(momentum_eq.bc.cavern_bcs)
        )
        momentum_eq._resid_ext_form = form
        momentum_eq._resid_ext_vec = fem_petsc.create_vector(form)
        momentum_eq._resid_ext_key = key
    f_ext = momentum_eq._resid_ext_vec
    with f_ext.localForm() as f_local:
        f_local.set(0.0)
    fem_petsc.assemble_vector(f_ext, momentum_eq._resid_ext_form)
    f_ext.ghostUpdate(addv=PETSc.InsertMode.ADD, mode=PETSc.ScatterMode.REVERSE)
    return f_ext


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

    # Internal force form ∫ σ : ε(u_) dx, compiled once (momentum_eq.sig is
    # a persistent Function updated in place above).
    if getattr(momentum_eq, "_resid_int_form", None) is None:
        internal_force_form = (
            ufl.inner(momentum_eq.sig, epsilon(momentum_eq.u_)) * momentum_eq.dx
        )
        momentum_eq._resid_int_form = do_fem.form(internal_force_form)
        momentum_eq._resid_int_vec = fem_petsc.create_vector(momentum_eq._resid_int_form)
    internal_force_vec = momentum_eq._resid_int_vec
    with internal_force_vec.localForm() as f_local:
        f_local.set(0.0)
    fem_petsc.assemble_vector(internal_force_vec, momentum_eq._resid_int_form)
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
    # Assemble external load vector (cached compiled form)
    external_loads_vec = _assemble_external_load_vector(momentum_eq)

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

