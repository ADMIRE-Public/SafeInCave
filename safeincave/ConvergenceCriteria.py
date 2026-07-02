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
- Combine multiple criteria with AND/OR logic via CompositeCriterion

**Criteria Available**:
- `StrainBasedCriterion`: Relative strain error (classical FEM approach)
- `ForceResidualCriterion`: Force residual normalized by external loads
- `DisplacementIncrementCriterion`: Newton step size relative to total increment
- `CompositeCriterion`: Combines multiple criteria with AND/OR logic

**Example Usage**:
```python
from safeincave.ConvergenceCriteria import (
    ForceResidualCriterion, DisplacementIncrementCriterion, CompositeCriterion
)

# Compose robust implicit criterion from modular components
residual_check = ForceResidualCriterion(tolerance=1e-3)
displacement_check = DisplacementIncrementCriterion(tolerance=1e-2)

criterion = CompositeCriterion(
    [residual_check, displacement_check],
    combine_logic='AND'  # Both must converge
)

criterion.initialize(momentum_eq)

# In iteration loop:
errors = criterion.compute_error(momentum_eq, stress, u_new, u_old)
converged = criterion.is_converged(*errors)
```
"""

from __future__ import annotations
from typing import TYPE_CHECKING, Tuple, List, Optional, Dict, Any
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
    momentum_eq: LinearMomentumBase, stress_field: to.Tensor, displacement: to.Tensor = None
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
    displacement : torch.Tensor, optional
        Displacement field (reserved for future extensions like energy norm).
        Default is None.

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


# ============================================================================
# LEGACY WRAPPER FUNCTIONS - Backward-Compatible Hardcoded Logic
# ============================================================================
# These functions encapsulate the exact hardcoded convergence logic currently
# used in Simulators.py to enable refactoring without functional changes.
# These are for backward compatibility only and replicate exact behavior.


def compute_strain_error_legacy(
    momentum_eq: LinearMomentumBase,
    strain_previous: to.Tensor,
    strain_current: to.Tensor
) -> float:
    """
    Compute strain-based error using legacy hardcoded logic.
    
    Exactly replicates the convergence check from Simulators.py.
    Returns 0 if explicit time integration or no inelastic elements.
    Uses MPI collective to synchronize across ranks.
    
    Parameters
    ----------
    momentum_eq : LinearMomentumBase
        Momentum equation solver instance.
    strain_previous : torch.Tensor
        Total strain from previous iteration, shape (n_elems, 3, 3).
    strain_current : torch.Tensor
        Total strain from current iteration, shape (n_elems, 3, 3).
    
    Returns
    -------
    float
        Relative strain error across all MPI ranks (uses MPI.SUM).
    
    Notes
    -----
    This function is provided for backward compatibility during refactoring.
    The hardcoded logic it replaces:
    - Returns 0 if theta=1.0 (explicit time integration)
    - Returns 0 if no inelastic elements
    - Otherwise computes ||ε_k - ε|| / ||ε|| with MPI synchronization
    """
    # Early exit for explicit time integration
    if momentum_eq.theta == 1.0:
        return 0.0
    
    # Early exit for purely elastic material
    if len(momentum_eq.mat.elems_ne) == 0:
        return 0.0
    
    # Compute local error
    eps_tot_k_flat = to.flatten(strain_previous)
    eps_tot_flat = to.flatten(strain_current)
    local_error = np.linalg.norm(
        eps_tot_k_flat - eps_tot_flat
    ) / np.linalg.norm(eps_tot_flat)
    
    # Synchronize across MPI ranks (uses SUM to match original behavior)
    error = momentum_eq.grid.mesh.comm.allreduce(
        local_error, op=MPI.SUM
    )
    
    return float(error)


class ConvergenceCriterion(ABC):
    """
    Abstract base class for convergence criteria.

    Each criterion computes a single error metric that can be checked against
    a tolerance. Criteria can be composed together via CompositeCriterion
    to create flexible convergence strategies.

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
    Classical strain-based convergence criterion.

    Error metric:
    $$\\text{error} = \\frac{\\|\\varepsilon_\\text{new} - \\varepsilon_\\text{old}\\|}{\\|\\varepsilon_\\text{new}\\|}$$

    **Advantages**:
    - Simple, interpretable (% change in strain)
    - Fast (no residual assembly)
    - Works for elastic and viscous problems

    **Disadvantages**:
    - Doesn't check equilibrium (residual)
    - Scale-dependent (strain magnitude affects convergence)

    **When to use**:
    - Elastic-only simulations
    - Viscous material response
    - Coupled with residual criterion for robustness

    Parameters
    ----------
    tolerance : float, optional
        Strain error tolerance. Default: 1e-7.
    """

    def __init__(self, tolerance: float = 1e-7):
        super().__init__(tolerance=tolerance, name="StrainBased")
        self.theta: float = 1.0  # Time integration parameter (0=implicit, 1=explicit)
        self.inelastic_element_count: int = 0  # Number of inelastic elements

    def initialize(self, momentum_eq: LinearMomentumBase) -> None:
        """Initialize strain criterion from momentum equation state."""
        self.theta = momentum_eq.theta
        self.inelastic_element_count = len(momentum_eq.mat.elems_ne)
        self.reset_history()

    def compute_error(
        self, strain_previous: to.Tensor, strain_current: to.Tensor
    ) -> float:
        """
        Compute strain-based error: ||Δε|| / ||ε_current||.

        Computes the relative change in total strain between consecutive iterations.
        Returns 0 if explicit time integration or purely elastic problem (no iteration).

        Parameters
        ----------
        strain_previous : torch.Tensor
            Total strain from previous iteration k, shape (n_elems, 3, 3).
        strain_current : torch.Tensor
            Total strain from current iteration k+1, shape (n_elems, 3, 3).

        Returns
        -------
        float
            Relative strain change error. Returns 0 if not applicable.
        """
        # Skip criterion if explicit time integration or purely elastic material
        if self.theta == 1.0 or self.inelastic_element_count == 0:
            error_value = 0.0
            self.history.append(error_value)
            return error_value

        strain_previous_flat = to.flatten(strain_previous)
        strain_current_flat = to.flatten(strain_current)

        strain_change_norm = to.norm(strain_current_flat - strain_previous_flat).item()
        strain_current_norm = max(to.norm(strain_current_flat).item(), 1e-16)

        error_value = strain_change_norm / strain_current_norm
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

    def compute_error(
        self, momentum_eq: LinearMomentumBase, stress: to.Tensor
    ) -> float:
        """
        Compute residual-based error: ||R|| / ||P_ext||.

        Parameters
        ----------
        momentum_eq : LinearMomentumBase
            Momentum equation solver (for assembly).
        stress : torch.Tensor
            Current stress field, shape (n_elems, 3, 3).

        Returns
        -------
        float
            Residual error. Dimensionless, typically in [0, ∞).
        """
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
        self.reset_history()

    def compute_error(
        self, momentum_eq: LinearMomentumBase, u_new: to.Tensor, u_old: to.Tensor
    ) -> float:
        """
        Compute displacement correction error.

        Parameters
        ----------
        momentum_eq : LinearMomentumBase
            Momentum equation solver (for norm computation).
        u_new : torch.Tensor
            Displacement at iteration k+1, shape (n_dofs,).
        u_old : torch.Tensor
            Displacement at iteration k, shape (n_dofs,).

        Returns
        -------
        float
            Displacement correction error. Dimensionless, in [0, ∞).
        """
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


class CompositeCriterion(ConvergenceCriterion):
    """
    Composite criterion that combines multiple criteria with AND/OR logic.

    Allows flexible composition of independent criteria based on physics
    and requirements. Users specify:
    1. Which criteria to use
    2. How to combine them (AND = all must converge, OR = any can converge)

    **Use Cases**:

    1. **Robust implicit (AND)**:
       - Both residual AND displacement must converge
       - For nonlinear, large-deformation problems
       ```python
       CompositeCriterion(
           [ForceResidualCriterion(1e-3), DisplacementIncrementCriterion(1e-2)],
           combine_logic='AND'
       )
       ```

    2. **Elastic with equilibrium check (AND)**:
       - Both strain AND residual converge
       - For elastic with complex BC
       ```python
       CompositeCriterion(
           [StrainBasedCriterion(1e-7), ForceResidualCriterion(1e-3)],
           combine_logic='AND'
       )
       ```

    3. **Progressive convergence (OR)**:
       - Any criterion can trigger exit (flexible for research)
       ```python
       CompositeCriterion(
           [StrainBasedCriterion(1e-5), ForceResidualCriterion(1e-4)],
           combine_logic='OR'
       )
       ```

    Parameters
    ----------
    criteria : list of ConvergenceCriterion
        List of criteria to combine.
    combine_logic : {'AND', 'OR'}, optional
        - 'AND': All criteria must be converged (default)
        - 'OR': Any criterion can be converged
    name : str, optional
        Descriptive name. Default: auto-generated from criteria.
    """

    def __init__(
        self,
        criteria: List[ConvergenceCriterion],
        combine_logic: str = "AND",
        name: Optional[str] = None,
    ):
        assert (
            combine_logic.upper() in ["AND", "OR"]
        ), "combine_logic must be 'AND' or 'OR'"
        assert len(criteria) > 0, "Must provide at least one criterion"
        
        # Prevent nesting of CompositeCriterion for clarity
        for i, crit in enumerate(criteria):
            assert not isinstance(crit, CompositeCriterion), (
                f"Nested CompositeCriterion at index {i} not supported. "
                "Use flat composition with all criteria at the same level. "
                "Example: CompositeCriterion([Residual(), Displacement()], 'AND')"
            )

        self.criteria = criteria
        self.combine_logic = combine_logic.upper()

        # Auto-generate name if not provided
        if name is None:
            crit_names = "+".join(c.name for c in criteria)
            op = "·" if self.combine_logic == "AND" else "∪"
            name = f"({crit_names})[{op}]"

        # Use max tolerance for display
        max_tol = max(c.tolerance for c in criteria)
        super().__init__(tolerance=max_tol, name=name)

    def initialize(self, momentum_eq: LinearMomentumBase) -> None:
        """
        Initialize all child criteria.

        Parameters
        ----------
        momentum_eq : LinearMomentumBase
            Momentum equation solver.

        Returns
        -------
        None
        """
        for criterion in self.criteria:
            criterion.initialize(momentum_eq)
        self.reset_history()

    def compute_error(
        self,
        momentum_eq: LinearMomentumBase,
        stress: to.Tensor,
        displacement_new: to.Tensor,
        displacement_old: to.Tensor,
    ) -> Tuple[float, ...]:
        """
        Compute error metrics from all child criteria.

        Evaluates each criterion in sequence and returns tuple of error values.
        Handles cases where strain data may not be available (e.g., residual-only).

        Parameters
        ----------
        momentum_eq : LinearMomentumBase
            Momentum equation solver.
        stress : torch.Tensor
            Current stress field, shape (n_elems, 3, 3).
        displacement_new : torch.Tensor
            Displacement at iteration k+1, shape (n_dofs,).
        displacement_old : torch.Tensor
            Displacement at iteration k, shape (n_dofs,).

        Returns
        -------
        tuple of float
            Error values from each child criterion (in order).
        """
        error_values = []

        for criterion in self.criteria:
            if isinstance(criterion, StrainBasedCriterion):
                # Strain criterion needs strain fields from momentum_eq
                strain_prev = getattr(momentum_eq, "_eps_tot_k", None)
                strain_curr = getattr(momentum_eq, "_eps_tot", None)
                if strain_prev is not None and strain_curr is not None:
                    error_val = criterion.compute_error(strain_prev, strain_curr)
                else:
                    # Graceful degradation: skip if strains unavailable
                    error_val = 0.0

            elif isinstance(criterion, ForceResidualCriterion):
                error_val = criterion.compute_error(momentum_eq, stress)

            elif isinstance(criterion, DisplacementIncrementCriterion):
                error_val = criterion.compute_error(
                    momentum_eq, displacement_new, displacement_old
                )

            else:
                # Unknown criterion type (shouldn't reach here)
                error_val = 0.0

            error_values.append(error_val)

        # Store composite error snapshot with metadata
        self.history.append(
            {
                "criterion_names": [c.name for c in self.criteria],
                "errors": error_values,
            }
        )

        return tuple(error_values)

    def is_converged(self, *errors: float) -> bool:
        """
        Check composite convergence based on logical combination rule.

        Parameters
        ----------
        *errors : float
            Error values from compute_error() (one per child criterion).

        Returns
        -------
        bool
            - AND logic: All errors ≤ respective tolerances
            - OR logic: Any error ≤ respective tolerance

        Raises
        ------
        ValueError
            If number of errors doesn't match number of criteria.
        """
        if len(errors) != len(self.criteria):
            raise ValueError(
                f"Expected {len(self.criteria)} error values, "
                f"got {len(errors)}"
            )

        convergence_flags = [
            errors[i] <= self.criteria[i].tolerance
            for i in range(len(self.criteria))
        ]

        if self.combine_logic == "AND":
            return all(convergence_flags)  # All criteria converged
        else:  # OR
            return any(convergence_flags)  # At least one converged

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
            "combine_logic": self.combine_logic,
            "criteria": [
                {
                    "name": c.name,
                    "tolerance": c.tolerance,
                    "history": c.get_history(),
                }
                for c in self.criteria
            ],
            "composite_history": self.history,
        }

