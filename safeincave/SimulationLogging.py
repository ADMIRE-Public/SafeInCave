# Copyright (c) 2026, The SafeInCave Developers
#
# SPDX-License-Identifier: BSD-3-Clause

"""Simulation report logger for element-level monitoring.

Provides centralized logging infrastructure to track element-level variables during
nonlinear iterations. Uses a flexible variable registry system that allows easy
extension without modifying the core logger.
"""

import torch as to
import numpy as np
from typing import Callable, Dict, List, Optional


# ============================================================================
# Variable Registry System
# ============================================================================

# Global registry for variable definitions
VARIABLE_REGISTRY: Dict[str, tuple] = {}


def register_variable(name: str, header: str, extractor: Callable) -> None:
    """
    Register a new variable in the global registry.
    
    This allows users to add custom variables without modifying the logger code.
    
    Parameters
    ----------
    name : str
        Variable name (used in variables_to_track list).
    header : str
        Column header for CSV output (e.g., 'Stress XX (Pa)').
    extractor : Callable
        Function that extracts the variable value. Signature:
        def extractor(elem_id: int, **kwargs) -> float:
            ...
    
    Examples
    --------
    >>> def extract_temperature(elem_id, **kwargs):
    ...     T = kwargs.get('temperature')
    ...     return T[elem_id].item() if T is not None else 0.0
    >>> register_variable('T', 'Temperature (K)', extract_temperature)
    """
    VARIABLE_REGISTRY[name] = (header, extractor)


def get_variable(name: str) -> Optional[tuple]:
    """
    Get variable definition from registry.
    
    Parameters
    ----------
    name : str
        Variable name.
    
    Returns
    -------
    tuple or None
        (header, extractor) tuple if found, None otherwise.
    """
    return VARIABLE_REGISTRY.get(name)


def list_registered_variables() -> List[str]:
    """
    List all registered variable names.
    
    Returns
    -------
    list of str
        Names of all registered variables.
    """
    return list(VARIABLE_REGISTRY.keys())


# ============================================================================
# Built-in Variable Extractors
# ============================================================================

def _extract_stress_component(component: str):
    """Create an extractor for stress tensor components."""
    component_map = {
        'xx': (0, 0), 'yy': (1, 1), 'zz': (2, 2),
        'xy': (0, 1), 'yx': (0, 1), 'yz': (1, 2), 'zy': (1, 2),
        'xz': (0, 2), 'zx': (0, 2)
    }
    
    if component not in component_map:
        raise ValueError(f"Unknown stress component: {component}")
    
    i, j = component_map[component]
    
    def extractor(elem_id: int, stress: to.Tensor, **kwargs) -> float:
        """Extract stress component for given element."""
        return stress[elem_id, i, j].item()
    
    return extractor


def _extract_strain_component(component: str):
    """Create an extractor for strain tensor components."""
    component_map = {
        'xx': (0, 0), 'yy': (1, 1), 'zz': (2, 2),
        'xy': (0, 1), 'yx': (0, 1), 'yz': (1, 2), 'zy': (1, 2),
        'xz': (0, 2), 'zx': (0, 2)
    }
    
    if component not in component_map:
        raise ValueError(f"Unknown strain component: {component}")
    
    i, j = component_map[component]
    
    def extractor(elem_id: int, strain: Optional[to.Tensor], **kwargs) -> float:
        """Extract strain component for given element."""
        if strain is None:
            return 0.0
        return strain[elem_id, i, j].item()
    
    return extractor


def _extract_mean_stress():
    """Create extractor for mean stress (hydrostatic)."""
    def extractor(elem_id: int, stress: to.Tensor, **kwargs) -> float:
        """Extract mean stress for given element."""
        return (stress[elem_id, 0, 0] + stress[elem_id, 1, 1] + stress[elem_id, 2, 2]).item() / 3.0
    return extractor


def _extract_von_mises_stress():
    """Create extractor for von Mises stress."""
    def extractor(elem_id: int, stress: to.Tensor, **kwargs) -> float:
        """Extract von Mises stress for given element."""
        s = stress[elem_id]
        s_xx, s_yy, s_zz = s[0, 0], s[1, 1], s[2, 2]
        s_xy, s_yz, s_xz = s[0, 1], s[1, 2], s[0, 2]
        
        # Compute mean stress (hydrostatic part)
        s_m = (s_xx + s_yy + s_zz) / 3.0
        
        # Compute deviatoric stress components
        s_xx_dev = s_xx - s_m
        s_yy_dev = s_yy - s_m
        s_zz_dev = s_zz - s_m
        
        # Compute second deviatoric invariant J2 using deviatoric components
        j2 = ((s_xx_dev - s_yy_dev)**2 + (s_yy_dev - s_zz_dev)**2 + 
              (s_zz_dev - s_xx_dev)**2 + 6.0 * (s_xy**2 + s_yz**2 + s_xz**2)) / 6.0
        
        # Von Mises stress: σ_vm = √(3·J2)
        # Use to.sqrt to ensure tensor operations, then convert to float
        von_mises = to.sqrt(3.0 * j2)
        return von_mises.item()
    return extractor


def _extract_yield_function():
    """Create extractor for yield function value."""
    def extractor(elem_id: int, **kwargs) -> float:
        """Extract yield function value for given element."""
        F = kwargs.get('yield_function')
        if F is None:
            return 0.0
        return F[elem_id].item()
    return extractor


def _extract_plastic_multiplier():
    """Create extractor for plastic multiplier."""
    def extractor(elem_id: int, **kwargs) -> float:
        """Extract plastic multiplier for given element."""
        delta_lambda = kwargs.get('plastic_multiplier')
        if delta_lambda is None:
            return 0.0
        return delta_lambda[elem_id].item()
    return extractor


def _extract_yield_indicator():
    """Create extractor for yield indicator (boolean)."""
    def extractor(elem_id: int, **kwargs) -> int:
        """Extract yield indicator for given element."""
        yield_indicator = kwargs.get('yield_indicator')
        if yield_indicator is None:
            return 0
        return int(yield_indicator[elem_id].item())
    return extractor


# ============================================================================
# Register Default Variables
# ============================================================================

def _register_default_variables():
    """Register the default set of variables."""
    register_variable('sxx', 'sxx (Pa)', _extract_stress_component('xx'))
    register_variable('syy', 'syy (Pa)', _extract_stress_component('yy'))
    register_variable('szz', 'szz (Pa)', _extract_stress_component('zz'))
    register_variable('sm', 'sm (Pa)', _extract_mean_stress())
    register_variable('svM', 'svM (Pa)', _extract_von_mises_stress())
    register_variable('exx', 'exx (-)', _extract_strain_component('xx'))
    register_variable('eyy', 'eyy (-)', _extract_strain_component('yy'))
    register_variable('ezz', 'ezz (-)', _extract_strain_component('zz'))
    register_variable('F', 'F (Pa)', _extract_yield_function())
    register_variable('dl', 'dl (-)', _extract_plastic_multiplier())
    register_variable('YieldDP', 'YieldDP', _extract_yield_indicator())


# Register default variables when module is loaded
_register_default_variables()


# ============================================================================
# Simulation Logging Class
# ============================================================================

class SimulationLogging:
    """
    Central simulation report logger for element-level monitoring.
    
    Tracks user-specified variables at a monitored element, writing results to CSV
    with full solver context (iteration count, nonlinear error, time step).
    Uses a flexible variable registry system for easy extension.
    
    Parameters
    ----------
    target_point : list or np.ndarray
        Target point [x, y, z] for element selection.
    variables_to_track : list of str, optional
        Names of variables to track (e.g., ['sxx', 'syy', 'szz', 'exx', 'eyy', 'ezz']).
        Variables must be registered in the VARIABLE_REGISTRY.
        Default tracks all registered variables.
    time_conversion : float, optional
        Conversion factor for time units (e.g., 3600 for hours to seconds).
        Default 1.0 (no conversion).
    
    Examples
    --------
    >>> # Simple usage with default variables
    >>> logger = sf.SimulationLogging(
    ...     target_point=[0.0, 0.0, 1.0],
    ...     variables_to_track=['sxx', 'syy', 'szz', 'exx', 'eyy', 'ezz']
    ... )
    """
    
    def __init__(
        self,
        target_point: list | np.ndarray,
        variables_to_track: Optional[List[str]] = None,
        time_conversion: float = 1.0,
    ):
        self.target_point = np.array(target_point)
        self.time_conversion = time_conversion

        # Internal state - set by configure()
        self._monitored_elem_id = None
        self._log_file = None
        self.bulk_modulus = None
        self.shear_modulus = None
        
        # Set variables to track
        if variables_to_track is None:
            # Track all registered variables
            self.variables_to_track = list_registered_variables()
        else:
            self.variables_to_track = variables_to_track
        
        # Build variables dictionary from registry
        self.variables = {}
        for var_name in self.variables_to_track:
            var_def = get_variable(var_name)
            if var_def is not None:
                self.variables[var_name] = var_def
            else:
                print(f"WARNING: Variable '{var_name}' not found in registry. Skipping.")
    
    def configure(
        self,
        centroids: np.ndarray,
        log_file: str | None = None,
        time_conversion: float | None = None,
    ) -> None:
        """
        Configure logger with grid and file information.
        
        Parameters
        ----------
        centroids : np.ndarray
            Element centroids, shape (N, 3).
        log_file : str, optional
            Path to output CSV file. If None, logging is disabled.
        time_conversion : float, optional
            Time conversion factor. If provided, overrides the value from __init__.
        """
        if log_file is None:
            self._monitored_elem_id = None
            self._log_file = None
            return
        
        # Update time_conversion if provided
        if time_conversion is not None:
            self.time_conversion = time_conversion
        
        # Select element closest to target point
        distances = np.linalg.norm(centroids - self.target_point, axis=1)
        self._monitored_elem_id = np.argmin(distances)
        self._log_file = log_file

        # Print initialization message
        print(f"{log_file}")

        # Initialize CSV file with header
        self._setup_file()
    
    def _setup_file(self) -> None:
        """Write CSV header to log file (MPI rank 0 only)."""
        try:
            from mpi4py import MPI
            rank = MPI.COMM_WORLD.rank
        except Exception:
            rank = 0
        
        if rank != 0:
            return
        
        try:
            # Create parent directory if it doesn't exist
            import os
            os.makedirs(os.path.dirname(self._log_file), exist_ok=True)
            
            # Use CSV writer for proper formatting
            import csv
            with open(self._log_file, 'w', newline='') as f:
                # Build header as list of column names
                header = ['Step', 'dt (h)', 't/t_final', 'Iters', 'NL_Error']
                
                # Add variable columns in order
                for var_name, (var_header, _) in self.variables.items():
                    header.append(var_header)

                writer = csv.writer(f)
                writer.writerow(header)
        except Exception as e:
            import sys
            print(f"ERROR initializing log file {self._log_file}: {e}", file=sys.stderr)
    
    def log_initial_state(
        self,
        time: float = 0.0,
        time_final: float = 1.0,
        stress: Optional[to.Tensor] = None,
        strain: Optional[to.Tensor] = None,
        **kwargs
    ) -> None:
        """
        Log the initial state (step 0) with zero dt and zero iterations.

        Parameters
        ----------
        time : float, optional
            Initial time (raw, SI seconds). Default 0.0.
        time_final : float, optional
            Final time endpoint. Default 1.0.
        stress : torch.Tensor, optional
            Full stress tensor, shape (N, 3, 3).
        strain : torch.Tensor, optional
            Total strain tensor, shape (N, 3, 3).
        **kwargs : dict
            Additional variables for extractors.
        """
        self.log_step(
            step=0,
            time=time,
            time_step=0.0,
            time_final=time_final,
            iteration=-1,          # displays as Iters=0
            nonlinear_error=0.0,
            stress=stress,
            strain=strain,
            **kwargs
        )

    def log_step(
        self,
        step: int,
        time: float = 0.0,
        time_step: float = 0.0,
        time_final: float = 1.0,
        iteration: int = 0,
        nonlinear_error: float = 0.0,
        stress: Optional[to.Tensor] = None,
        strain: Optional[to.Tensor] = None,
        **kwargs
    ) -> None:
        """
        Log simulation state at monitored element.
        
        Parameters
        ----------
        step : int
            Time step number (raw).
        time : float, optional
            Current absolute time (raw, SI seconds). Default 0.0.
        time_step : float, optional
            Time step size (raw, SI seconds). Default 0.0.
        time_final : float, optional
            Final time endpoint. Default 1.0.
        iteration : int, optional
            Newton iteration number (raw 0-indexed). Default 0.
        nonlinear_error : float, optional
            Nonlinear error estimate. Default 0.0.
        stress : torch.Tensor, optional
            Full stress tensor, shape (N, 3, 3).
        strain : torch.Tensor, optional
            Total strain tensor, shape (N, 3, 3).
        **kwargs : dict
            Additional variables for extractors (e.g., yield_function, 
            plastic_multiplier, yield_indicator, temperature, etc.).
        """
        if self._monitored_elem_id is None or self._log_file is None:
            return
        
        try:
            elem_id = self._monitored_elem_id
            
            # Convert dt from physical units (SI seconds) to display units
            time_step_converted = time_step / self.time_conversion if self.time_conversion != 1.0 else time_step
            
            # Increment iteration by 1 for display (show as 1, 2, 3... instead of 0, 1, 2...)
            iteration_display = iteration + 1
            
            # Compute time ratio for display
            time_ratio = 0.0 if time_final == 0.0 else time / time_final
            
            # Build base row with standard columns
            row = [step, time_step_converted, time_ratio, iteration_display, nonlinear_error]
            
            # Extract and add each variable
            for var_name, (var_header, extractor) in self.variables.items():
                try:
                    value = extractor(elem_id, stress=stress, strain=strain, **kwargs)
                    row.append(value)
                except Exception:
                    # If extraction fails, log 0.0 and continue
                    row.append(0.0)

            # Write to log file (rank 0 only)
            try:
                from mpi4py import MPI
                rank = MPI.COMM_WORLD.rank
            except Exception:
                rank = 0
            
            if rank == 0:
                import csv
                with open(self._log_file, 'a', newline='') as f:
                    writer = csv.writer(f)
                    writer.writerow(row)
        except Exception as e:
            import sys
            print(f"ERROR logging state: {e}", file=sys.stderr)
    
    @staticmethod
    def get_centroids(grid):
        """
        Compute element centroids from mesh geometry.

        Parameters
        ----------
        grid : GridHandlerGMSH
            Grid object with mesh and topology information.

        Returns
        -------
        numpy.ndarray
            Element centroids, shape (n_elems, 3), where each row is the
            arithmetic mean of the four vertices of the corresponding tetrahedral cell.
        """
        conn = grid.mesh.topology.connectivity(3, 0).array.reshape(-1, 4)
        coord = grid.mesh.geometry.x
        centroids = np.zeros((grid.n_elems, 3))
        for i in range(grid.n_elems):
            nodes = conn[i]
            centroids[i, :] = (
                coord[nodes[0], :] + coord[nodes[1], :] + coord[nodes[2], :] + coord[nodes[3], :]
            ) / 4.0
        return centroids
    
    @staticmethod
    def auto_configure_from_simulator(logger, simulator, outputs) -> None:
        """
        Auto-configure logger from simulator and outputs.
        
        This is a convenience method that extracts configuration from a simulator
        and configures the logger. Use this in simulator __init__ methods.
        
        Parameters
        ----------
        logger : SimulationLogging
            Logger to configure.
        simulator : Simulator
            Simulator instance (Simulator_M, Simulator_TM, etc.).
        outputs : list of SaveFields
            Output handlers to extract output folder from.
        """
        if logger is None:
            return
        
        # Get centroids from grid
        try:
            grid = getattr(simulator, 'eq_mom', None)
            if grid is None:
                grid = getattr(simulator, 'eq_heat', None)
            if grid is None:
                return
            
            centroids = SimulationLogging.get_centroids(grid.grid)
        except Exception as e:
            print(f"WARNING: Could not extract centroids for logger: {e}")
            return
        
        # Auto-detect log file path from outputs
        log_file = None
        for output in outputs:
            if hasattr(output, 'output_folder') and output.output_folder:
                import os
                log_file = os.path.join(output.output_folder, "log.csv")
                break
        
        # Get time_conversion from time controller
        time_conversion = getattr(simulator, 't_control', None)
        time_conversion = getattr(time_conversion, 'time_conversion', 1.0) if time_conversion else 1.0
        
        # Configure logger
        if log_file is not None:
            logger.configure(
                centroids=centroids,
                log_file=log_file,
                time_conversion=time_conversion,
            )
    
    @staticmethod
    def extract_yield_variables(material) -> dict:
        """
        Extract yield-related variables from material for logging.
        
        Parameters
        ----------
        material : Material
            Material object with non-elastic elements.
        
        Returns
        -------
        dict
            Dictionary with yield_function, plastic_multiplier, and yield_indicator
            if available, otherwise empty dict.
        """
        log_kwargs = {}
        
        if not hasattr(material, 'elems_ne'):
            return log_kwargs
        
        # Try to extract yield indicators from non-elastic elements
        for elem_ne in material.elems_ne:
            if hasattr(elem_ne, 'F') and hasattr(elem_ne, 'delta_lambda'):
                log_kwargs['yield_function'] = elem_ne.F
                log_kwargs['plastic_multiplier'] = elem_ne.delta_lambda
                if hasattr(elem_ne, 'YieldDP'):
                    log_kwargs['yield_indicator'] = elem_ne.YieldDP
                elif hasattr(elem_ne, 'is_plastic'):
                    log_kwargs['yield_indicator'] = elem_ne.is_plastic
                break
        
        return log_kwargs


# ============================================================================
# Backward Compatibility Wrapper
# ============================================================================

class SimulationLoggingLegacy(SimulationLogging):
    """
    Legacy interface for SimulationLogging with backward compatibility.
    
    This class maintains the old API including auto_configure() and log_state().
    Deprecated: Use SimulationLogging instead.
    """
    
    def __init__(
        self,
        target_point: list | np.ndarray,
        variables_to_track: Optional[List[str]] = None,
    ):
        """Initialize with legacy interface."""
        super().__init__(target_point=target_point, variables_to_track=variables_to_track)
    
    def configure(
        self,
        centroids: np.ndarray,
        log_file: str | None = None,
        mom_eq=None,
        bulk_modulus: Optional[to.Tensor] = None,
        shear_modulus: Optional[to.Tensor] = None,
        time_conversion: Optional[float] = None,
    ) -> None:
        """
        Configure with legacy interface.
        
        .. deprecated::
            Use the new interface without mom_eq parameter.
        """
        import warnings
        warnings.warn(
            "SimulationLoggingLegacy is deprecated. Use SimulationLogging instead.",
            DeprecationWarning,
            stacklevel=2
        )
        
        # Update time_conversion if provided
        if time_conversion is not None:
            self.time_conversion = time_conversion
        
        # Call parent configure
        super().configure(centroids, log_file)
        
        # Store bulk/shear modulus if provided
        if bulk_modulus is not None:
            self.bulk_modulus = bulk_modulus
        if shear_modulus is not None:
            self.shear_modulus = shear_modulus
    
    def auto_configure(
        self,
        grid,
        material,
        outputs: Optional[List] = None,
        log_file: Optional[str] = None,
        time_conversion: Optional[float] = None,
    ) -> None:
        """
        Legacy auto-configure method.
        
        .. deprecated::
            Use configure() directly with centroids and log_file instead.
        """
        import warnings
        warnings.warn(
            "auto_configure() is deprecated. Use configure() with centroids and log_file instead.",
            DeprecationWarning,
            stacklevel=2
        )
        
        # Extract centroids from grid
        centroids = self.get_centroids(grid)
        
        # Auto-detect log file from outputs if not provided
        actual_log_file = log_file
        if actual_log_file is None and outputs:
            for output in outputs:
                if hasattr(output, 'output_folder') and output.output_folder:
                    import os
                    actual_log_file = os.path.join(output.output_folder, "log.csv")
                    break
        
        # Extract bulk/shear modulus from material if available
        bulk_mod = None
        shear_mod = None
        if material is not None:
            if hasattr(material, 'K_elems'):
                bulk_mod = material.K_elems
                shear_mod = material.G_elems
            elif hasattr(material, 'C') and material.C is not None and material.C.shape[0] > 0:
                C = material.C
                bulk_mod = (C[:, 0, 0] + C[:, 1, 1] + C[:, 2, 2] + 2*(C[:, 0, 1] + C[:, 0, 2] + C[:, 1, 2])) / 9.0
                shear_mod = C[:, 3, 3]
        
        # Call configure
        self.configure(
            centroids=centroids,
            log_file=actual_log_file,
            bulk_modulus=bulk_mod,
            shear_modulus=shear_mod,
            time_conversion=time_conversion,
        )
    
    def log_state(
        self,
        step: int,
        stress: to.Tensor,
        F: to.Tensor,
        delta_lambda: to.Tensor,
        YieldDP: to.Tensor,
        iteration: int = 0,
        nonlinear_error: float = 0.0,
        time: float = 0.0,
        time_step: float = 0.0,
        time_final: float = 1.0,
        strain: Optional[to.Tensor] = None,
        sm: Optional[to.Tensor] = None,
        svM: Optional[to.Tensor] = None,
    ) -> None:
        """
        Legacy log_state method.
        
        .. deprecated::
            Use log_step() with new interface instead.
        """
        import warnings
        warnings.warn(
            "log_state() is deprecated. Use log_step() with new interface.",
            DeprecationWarning,
            stacklevel=2
        )
        
        # Convert legacy call to new interface
        kwargs = {
            'yield_function': F,
            'plastic_multiplier': delta_lambda,
            'yield_indicator': YieldDP,
        }
        
        # If sm and svM are provided, override the extractors temporarily
        original_vars = self.variables.copy()
        if sm is not None and 'sm' in self.variables:
            def make_sm_extractor(sm_tensor):
                def extractor(elem_id, **kwargs):
                    return sm_tensor[elem_id].item()
                return extractor
            self.variables['sm'] = ('sm (Pa)', make_sm_extractor(sm))
        
        if svM is not None and 'svM' in self.variables:
            def make_svm_extractor(svm_tensor):
                def extractor(elem_id, **kwargs):
                    return svm_tensor[elem_id].item()
                return extractor
            self.variables['svM'] = ('svM (Pa)', make_svm_extractor(svM))
        
        # Call new log_step
        self.log_step(
            step=step,
            time=time,
            time_step=time_step,
            time_final=time_final,
            iteration=iteration,
            nonlinear_error=nonlinear_error,
            stress=stress,
            strain=strain,
            **kwargs
        )
        
        # Restore original variables
        self.variables = original_vars