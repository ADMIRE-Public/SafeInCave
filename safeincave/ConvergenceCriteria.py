# Copyright (c) 2026, The SafeInCave Developers
#
# SPDX-License-Identifier: BSD-3-Clause

"""
Modular convergence criteria for nonlinear finite element iterations.

This module provides abstract and concrete implementations of modular convergence
criteria that can be composed together based on problem physics and requirements.

**Philosophy**:
Instead of hardcoded criterion combinations, criteria are independent and can be
composed flexibly:
- Use strain-based for elastic/viscous problems
- Add residual criterion for nonlinear equilibrium checking
- Add displacement correction for robust implicit iteration
- Use ForceDisplacementCriterion for robust dual checks

**Criteria Available**:
- `StrainBasedCriterion`: Relative strain error (classical FEM approach)
- `ForceResidualCriterion`: Force residual normalized by external loads
- `DisplacementIncrementCriterion`: Newton step size relative to total increment
- `ForceDisplacementCriterion`: Combined force residual + displacement increment criterion

**Example Usage**:
```python
from safeincave.ConvergenceCriteria import (
    ForceResidualCriterion, DisplacementIncrementCriterion, ForceDisplacementCriterion
)

# Compose robust implicit criterion from modular components
residual_check = ForceResidualCriterion(tolerance=1e-3)
displacement_check = DisplacementIncrementCriterion(tolerance=1e-2)

criterion = ForceDisplacementCriterion(
    force_tolerance=1e-3,
    displacement_tolerance=1e-2,
)

criterion.initialize(momentum_eq)

# In iteration loop:
error = criterion.compute_error(momentum_eq)
converged = criterion.is_converged(error)
```
"""

from __future__ import annotations
from typing import TYPE_CHECKING, List, Optional, Dict, Any
from abc import ABC, abstractmethod
import torch as to
import numpy as np
from mpi4py import MPI
from petsc4py import PETSc
import ufl
from dolfinx import fem as do_fem
from dolfinx.fem import petsc as fem_petsc
from .Utils import epsilon

if TYPE_CHECKING:
    from MomentumEquation import LinearMomentumBase


# ============================================================================
# HELPER FUNCTIONS - Convergence Computation Utilities
# ============================================================================
# These functions encapsulate FEniCS/PETSc operations needed for convergence
# checking. They operate on LinearMomentumBase instances and are used by all
# criteria classes.


def _compute_external_load_vector_norm(momentum_eq: LinearMomentumBase) -> float:
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

    # Apply lifting for Dirichlet BCs (project to null space of essential BCs)
    fem_petsc.apply_lifting(
        f_ext,
        [
            do_fem.form(
                ufl.inner(ufl.grad(momentum_eq.du), ufl.grad(momentum_eq.u_))
                * momentum_eq.dx
            )
        ],
        [momentum_eq.bc.dirichlet_bcs],
    )
    f_ext.ghostUpdate(addv=PETSc.InsertMode.ADD, mode=PETSc.ScatterMode.REVERSE)

    # Compute norm
    norm = f_ext.norm()

    # Protect against zero norm
    norm = max(norm, 1e-16)

    return norm


def _compute_internal_force_vector(
    momentum_eq: LinearMomentumBase, stress_field: to.Tensor
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
    momentum_eq: LinearMomentumBase, stress_field: to.Tensor
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
    torch.Tensor
        Residual vector R = P_ext - q_int, same shape as displacement DoFs.

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

    # Apply lifting for Dirichlet BCs (project to null space)
    fem_petsc.apply_lifting(
        external_loads_vec,
        [
            do_fem.form(
                ufl.inner(ufl.grad(momentum_eq.du), ufl.grad(momentum_eq.u_))
                * momentum_eq.dx
            )
        ],
        [momentum_eq.bc.dirichlet_bcs],
    )
    external_loads_vec.ghostUpdate(addv=PETSc.InsertMode.ADD, mode=PETSc.ScatterMode.REVERSE)

    # Get internal force vector
    internal_forces_tensor = _compute_internal_force_vector(momentum_eq, stress_field)

    # Convert external loads to torch tensor for residual computation
    external_loads_tensor = to.from_numpy(external_loads_vec.array.copy()).double()

    # Compute residual: R = P_ext - q_int (out-of-balance forces)
    residual_vector = external_loads_tensor - internal_forces_tensor

    return residual_vector


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


def _initialize_step_displacement(momentum_eq: LinearMomentumBase) -> to.Tensor:
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
    convergence_criterion: str = "strain_based",
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
    criterion_key = convergence_criterion.strip().lower()

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
    momentum_eq: LinearMomentumBase,
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
    momentum_eq: LinearMomentumBase,
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

    def __init__(self, convergence_criterion: str = "strain_based", tol: float = 1e-8):
        self.criterion = resolve_convergence_criterion(convergence_criterion)
        self.error: float = 0.0
        self.not_converged_error: bool = True
        self.below_max_iterations: bool = True
        self.last_raw_error: float = 0.0
        self.tol: float = tol
        self.ite: int = 0
        self.maxiter: Optional[int] = None

    def initialize_step(
        self,
        momentum_eq: LinearMomentumBase,
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

    def compute_error(self, momentum_eq: LinearMomentumBase) -> float:
        """Compute raw criterion error and store for diagnostics."""
        raw_error = compute_error_from_criterion(momentum_eq, self.criterion)
        self.error = float(raw_error)
        self.last_raw_error = self.error
        return self.error

    def evaluate(
        self,
        momentum_eq: LinearMomentumBase,
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
        self.not_converged_error = self.error > self.tol
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
        return float(self.criterion.tolerance)

    def get_last_raw_error(self) -> float:
        """Return most recent raw criterion value (not zeroed by convergence logic)."""
        return self.last_raw_error

    def set_criterion(self, convergence_criterion: str) -> None:
        """Swap convergence strategy at runtime."""
        self.criterion = resolve_convergence_criterion(convergence_criterion)


def _is_momentum_solver_instance(obj: Any) -> bool:
    """Return True if object provides the momentum solver interface used here."""
    return hasattr(obj, "compute_total_strain") and hasattr(obj, "u") and hasattr(obj, "mat")


class ConvergenceCriterion(ABC):
    """
    Abstract base class for convergence criteria.

    Each criterion computes a single error metric that can be checked against
    a tolerance to create flexible convergence strategies.

    Parameters
    ----------
    tolerance : float, optional
        Convergence tolerance for this criterion. Default: 1e-7.
    name : str, optional
        Descriptive name for diagnostics/logging. Default: class name.

    Attributes
    ----------
    tolerance : float
        Convergence tolerance.
    name : str
        Criterion name.
    history : list
        Record of error values per iteration.
    """

    def __init__(self, tolerance: float = 1e-7, name: Optional[str] = None):
        self.tolerance = tolerance
        self.name = name or self.__class__.__name__
        self.history: list = []

    @abstractmethod
    def initialize(self, momentum_eq: LinearMomentumBase) -> None:
        """
        Initialize criterion at step start.

        Called once per time step, before entering the nonlinear iteration loop.

        Parameters
        ----------
        momentum_eq : LinearMomentumBase
            Momentum equation solver instance.

        Returns
        -------
        None
        """
        pass

    @abstractmethod
    def compute_error(self, *args, **kwargs) -> float:
        """
        Compute error metric for current iteration.

        Returns
        -------
        float
            Dimensionless error value. Returns 0 if criterion not applicable
            (e.g., elastic-only problem for strain criterion).
        """
        pass

    @abstractmethod
    def is_converged(self, error: float) -> bool:
        """
        Check if error is below tolerance.

        Parameters
        ----------
        error : float
            Error from compute_error().

        Returns
        -------
        bool
            True if error ≤ tolerance.
        """
        pass

    def get_history(self) -> list:
        """Return iteration history."""
        return self.history

    def reset_history(self) -> None:
        """Reset iteration history."""
        self.history = []

    def __repr__(self) -> str:
        return f"{self.name}(tol={self.tolerance:.2e})"


class StrainBasedCriterion(ConvergenceCriterion):
    """
    Strain-based convergence criterion.

    Error metric:
    $$\\text{error} = \\frac{\\|\\varepsilon_\\text{new} - \\varepsilon_\\text{old}\\|}{\\|\\varepsilon_\\text{new}\\|}$$

    This implementation preserves prior behavior,
    including MPI SUM reduction and momentum-equation-tracked previous strain.

    Parameters
    ----------
    tolerance : float, optional
        Strain error tolerance. Default: 1e-7.
    """

    def __init__(self, tolerance: float = 1e-7, name: Optional[str] = None):
        super().__init__(tolerance=tolerance, name=name or "strain_based")

    def initialize(self, momentum_eq: LinearMomentumBase) -> None:
        """Initialize tracked strain state at step start."""
        momentum_eq._strain_previous = momentum_eq.compute_total_strain().clone()
        self.reset_history()

    def compute_error(self, *args, **kwargs) -> float:
        """
        Compute strain-based error: ||Δε|| / ||ε_current||.

        Parameters
        ----------
        Supports call styles:
        1) compute_error(momentum_eq)
        2) compute_error(momentum_eq, strain_previous, strain_current)

        Returns
        -------
        float
            Relative strain change error. Returns 0 if not applicable.
        """
        if len(args) == 1 and _is_momentum_solver_instance(args[0]):
            momentum_eq = args[0]
            strain_current = momentum_eq.compute_total_strain()
            if hasattr(momentum_eq, "_strain_previous"):
                strain_previous = momentum_eq._strain_previous
            else:
                strain_previous = strain_current
        elif len(args) == 3 and _is_momentum_solver_instance(args[0]):
            momentum_eq = args[0]
            strain_previous = args[1]
            strain_current = args[2]
        else:
            raise ValueError(
                "StrainBasedCriterion.compute_error expects (momentum_eq) or "
                "(momentum_eq, strain_previous, strain_current)."
            )

        # Early exit for explicit time integration
        if momentum_eq.theta == 1.0:
            error_value = 0.0
            self.history.append(error_value)
            return error_value

        # Early exit for purely elastic material
        if len(momentum_eq.mat.elems_ne) == 0:
            error_value = 0.0
            self.history.append(error_value)
            return error_value

        strain_previous_flat = to.flatten(strain_previous)
        strain_current_flat = to.flatten(strain_current)

        denom = max(np.linalg.norm(strain_current_flat), 1e-16)
        local_error = np.linalg.norm(strain_previous_flat - strain_current_flat) / denom
        error = momentum_eq.grid.mesh.comm.allreduce(local_error, op=MPI.SUM)

        # Update tracked state for next call.
        momentum_eq._strain_previous = strain_current.clone()

        error_value = float(error)
        self.history.append(error_value)
        return float(error_value)

    def is_converged(self, error: float) -> bool:
        """Check if strain error is below tolerance."""
        return error <= self.tolerance


class ForceResidualCriterion(ConvergenceCriterion):
    """
    Force residual convergence criterion (equilibrium check).

    Error metric (normalized out-of-balance forces):
    $$\\text{error}_\\text{residual} = \\frac{\\|R\\|}{\\|P_\\text{ext}\\|}$$

    where R = P_ext - q_int (out-of-balance forces).

    **Purpose**:
    - Checks mechanical equilibrium
    - Detects stalled iterations (large residual despite small displacement)
    - Robust for problems with varying load magnitudes

    **Advantages**:
    - Direct equilibrium measure
    - Scale-robust (normalized by load)
    - Detects residual stalling

    **Disadvantages**:
    - Requires residual assembly (expensive)
    - Must be combined with displacement criterion for full convergence check

    **Typical Tolerance**: 1e-3 to 1e-4 (0.1% to 0.01% of load magnitude)

    Parameters
    ----------
    tolerance : float, optional
        Residual tolerance (as fraction of external load). Default: 1e-3.
    """

    def __init__(self, tolerance: float = 1e-3):
        super().__init__(tolerance=tolerance, name="ForceResidual")

    def initialize(self, momentum_eq: LinearMomentumBase) -> None:
        """Initialize residual criterion."""
        self.reset_history()

    def compute_error(self, *args, **kwargs) -> float:
        """
        Compute residual-based error: ||R|| / ||P_ext||.

        Parameters
        ----------
        Supports two call styles:
        1) compute_error(momentum_eq, stress)
        2) compute_error(momentum_eq)  # stress inferred from momentum_eq.sig

        Returns
        -------
        float
            Residual error. Dimensionless, typically in [0, ∞).
        """
        if len(args) == 2 and _is_momentum_solver_instance(args[0]):
            momentum_eq = args[0]
            stress = args[1]
        elif len(args) == 1 and _is_momentum_solver_instance(args[0]):
            momentum_eq = args[0]
            stress = to.as_tensor(
                momentum_eq.sig.x.array.reshape((momentum_eq.n_elems, 3, 3))
            )
        else:
            raise ValueError(
                "ForceResidualCriterion.compute_error expects (momentum_eq) or "
                "(momentum_eq, stress)."
            )

        # Compute residual R = P_ext - q_int
        residual_vector = _compute_force_residual(momentum_eq, stress)
        residual_norm = _compute_vector_norm(residual_vector)

        # Normalize by external load norm for scale-robustness
        external_load_norm = _compute_external_load_vector_norm(momentum_eq)
        external_load_norm_safe = max(external_load_norm, 1e-16)

        error_value = residual_norm / external_load_norm_safe
        self.history.append(error_value)
        return float(error_value)

    def is_converged(self, error: float) -> bool:
        """Check if residual error is below tolerance."""
        return error <= self.tolerance


class DisplacementIncrementCriterion(ConvergenceCriterion):
    """
    Displacement increment convergence criterion (Newton step measure).

    Error metric (Newton step relative to total increment):
    $$\\text{error}_\\text{displacement} = \\frac{\\|\\Delta u_\\text{correction}\\|}{\\|\\Delta u_\\text{total}\\|}$$

    where:
    - Δu_correction = ||u^(k+1) - u^(k)|| (Newton iteration step in current iteration)
    - Δu_total = ||u^(k+1) - u_step_start|| (cumulative displacement since step start)

    **Purpose**:
    - Measures iteration progress relative to step size
    - Scale-independent (correction vs. total, not vs. old)
    - Detects convergence stalling (very small corrections)

    **Advantages**:
    - Scale-robust (independent of problem size)
    - Captures iteration difficulty relative to step magnitude
    - Robust for large deformations

    **Disadvantages**:
    - Requires tracking u_step_start (requires initialize_step_displacements())
    - Can be strict for small steps (ratio may remain large)
    - Should be combined with residual criterion

    **Typical Tolerance**: 1e-2 (1% of total step increment)

    Parameters
    ----------
    tolerance : float, optional
        Displacement correction tolerance. Default: 1e-2.
    """

    def __init__(self, tolerance: float = 1e-2):
        super().__init__(tolerance=tolerance, name="DisplacementIncrement")
        self.u_step_start: Optional[to.Tensor] = None
        self.u_previous: Optional[to.Tensor] = None

    def initialize(self, momentum_eq: LinearMomentumBase) -> None:
        """
        Initialize displacement criterion and capture starting displacement.

        Called once per time step to store the initial displacement state
        for tracking the cumulative increment throughout the step.

        Parameters
        ----------
        momentum_eq : LinearMomentumBase
            Momentum equation solver instance.

        Returns
        -------
        None

        Side Effects
        -----------
        Stores initial displacement via _initialize_step_displacement().
        """
        self.u_step_start = _initialize_step_displacement(momentum_eq)
        self.u_previous = to.as_tensor(momentum_eq.u.x.array.copy())
        self.reset_history()

    def compute_error(self, *args, **kwargs) -> float:
        """
        Compute displacement correction error.

        Parameters
        ----------
        Supports two call styles:
        1) compute_error(momentum_eq, u_new, u_old)
        2) compute_error(momentum_eq)  # u_old tracked internally

        Returns
        -------
        float
            Displacement correction error. Dimensionless, in [0, ∞).
        """
        if len(args) == 3 and _is_momentum_solver_instance(args[0]):
            momentum_eq = args[0]
            u_new = args[1]
            u_old = args[2]
            self.u_previous = u_new.clone()
        elif len(args) == 1 and _is_momentum_solver_instance(args[0]):
            momentum_eq = args[0]
            u_new = to.as_tensor(momentum_eq.u.x.array.copy())
            if self.u_previous is None:
                u_old = u_new
            else:
                u_old = self.u_previous
            self.u_previous = u_new.clone()
        else:
            raise ValueError(
                "DisplacementIncrementCriterion.compute_error expects "
                "(momentum_eq) or (momentum_eq, u_new, u_old)."
            )

        # Compute Newton iteration correction (single step)
        displacement_correction = u_new - u_old
        correction_norm = _compute_vector_norm(displacement_correction)

        # Compute cumulative displacement increment from step start
        displacement_total = u_new - self.u_step_start
        total_increment_norm = _compute_vector_norm(displacement_total)

        # Ratio: single step vs. cumulative (scale-independent)
        total_increment_norm_safe = max(total_increment_norm, 1e-16)
        error_value = correction_norm / total_increment_norm_safe

        self.history.append(error_value)
        return float(error_value)

    def is_converged(self, error: float) -> bool:
        """Check if displacement correction error is below tolerance."""
        return error <= self.tolerance


class ForceDisplacementCriterion(ConvergenceCriterion):
    """
    Explicit combined criterion: force residual + displacement increment.

    The combined scalar error is defined as:
    max(error_force / tol_force, error_displacement / tol_displacement)

    Convergence is achieved when the combined error is <= 1.0, which is
    equivalent to satisfying both component tolerances.

    Parameters
    ----------
    force_tolerance : float, optional
        Tolerance for force residual criterion. Default: 1e-3.
    displacement_tolerance : float, optional
        Tolerance for displacement increment criterion. Default: 1e-2.
    name : str, optional
        Name for diagnostics. Default: "force_displacement".
    """

    def __init__(
        self,
        force_tolerance: float = 1e-3,
        displacement_tolerance: float = 1e-2,
        name: Optional[str] = None,
    ):
        criterion_name = name or "force_displacement"
        # tolerance=1 means both criteria satisfy their own tolerances
        super().__init__(tolerance=1.0, name=criterion_name)
        self.force_criterion = ForceResidualCriterion(tolerance=force_tolerance)
        self.displacement_criterion = DisplacementIncrementCriterion(
            tolerance=displacement_tolerance
        )
        self.component_history: List[Dict[str, float]] = []

    def initialize(self, momentum_eq: LinearMomentumBase) -> None:
        """
        Initialize component criteria.

        Parameters
        ----------
        momentum_eq : LinearMomentumBase
            Momentum equation solver.

        Returns
        -------
        None
        """
        self.force_criterion.initialize(momentum_eq)
        self.displacement_criterion.initialize(momentum_eq)
        self.reset_history()
        self.component_history = []

    def compute_error(self, momentum_eq: LinearMomentumBase) -> float:
        """
        Compute combined normalized error from both component criteria.

        Parameters
        ----------
        momentum_eq : LinearMomentumBase
            Momentum equation solver.

        Returns
        -------
        float
            Combined normalized error.
        """
        force_error = self.force_criterion.compute_error(momentum_eq)
        displacement_error = self.displacement_criterion.compute_error(momentum_eq)

        force_ratio = force_error / max(self.force_criterion.tolerance, 1e-16)
        displacement_ratio = displacement_error / max(
            self.displacement_criterion.tolerance, 1e-16
        )
        combined_error = max(force_ratio, displacement_ratio)

        self.component_history.append(
            {
                "force_error": float(force_error),
                "displacement_error": float(displacement_error),
                "force_ratio": float(force_ratio),
                "displacement_ratio": float(displacement_ratio),
            }
        )
        self.history.append(float(combined_error))
        return float(combined_error)

    def is_converged(self, error: float) -> bool:
        """
        Check composite convergence based on logical combination rule.

        Parameters
        ----------
        error : float
            Combined normalized error.

        Returns
        -------
        bool
            True when error <= 1.0.
        """
        return error <= self.tolerance

    def get_convergence_info(self) -> Dict[str, Any]:
        """
        Return detailed convergence diagnostics.

        Returns
        -------
        dict
            Composite convergence information with per-criterion history.
        """
        return {
            "n_iterations": len(self.history),
            "criterion": self.name,
            "combined_tolerance": self.tolerance,
            "force_tolerance": self.force_criterion.tolerance,
            "displacement_tolerance": self.displacement_criterion.tolerance,
            "combined_history": self.history,
            "component_history": self.component_history,
        }

