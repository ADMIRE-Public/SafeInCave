# Copyright (c) 2026, The SafeInCave Developers
#
# SPDX-License-Identifier: BSD-3-Clause

"""Base simulator class with shared step scaffolding and bisection/cutback logic."""

from abc import ABC, abstractmethod
import copy

class Simulator(ABC):
    """
    Abstract simulation driver interface.

    Subclasses implement a concrete `run()` method that advances one or more
    coupled PDE solvers in time, handles I/O, and updates material/internal
    variables as needed.
    """

    @abstractmethod
    def run(self):
        """
        Execute the simulation.

        Returns
        -------
        None
        """
        pass

    def _final_consuming_solve(self, stress_to, t: float, dt: float):
        """
        Advance eps_ne_k with the converged inelastic rate before committing.

        The nonlinear loop's order is solve -> stress -> compute_eps_ne_rate,
        so the rate computed in the final iteration is never pushed into eps_ne_k
        before update_eps_ne_old commits it. This call advances eps_ne_k without
        resolving the momentum equation (displacement is already converged).

        Note: does NOT call increment_internal_variables again here — the main
        loop's last iteration already incremented hardening variables (e.g.
        ViscoplasticDesai.alpha) using self.r/self.h/self.P linearized at that
        iteration's stress. Calling it again here reapplies that same stale
        linearization against a second, different stress delta, double-counting
        the hardening update every step. Harmless for models without persistent
        hardening state via this hook, but on models that have it
        (ViscoplasticDesai, MunsonDawsonCreep, ModifiedCamClayViscoplastic) it
        drove alpha to drift until compute_CT's consistent tangent went singular.
        """
        # Commit the plain (un-relaxed) rate, not the over-relaxed iterate.
        for elem_ne in getattr(self.eq_mom.mat, "elems_ne", []):
            fn = getattr(elem_ne, "finalize_step_rate", None)
            if callable(fn):
                fn()
        # Advance eps_ne_k with the final converged rate, then update stress.
        stress_k_to = stress_to.clone()
        self.eq_mom.compute_eps_ne_k(dt)
        eps_tot_to = self.eq_mom.compute_total_strain()
        stress_to = self.eq_mom.compute_stress(eps_tot_to)
        return stress_to, stress_k_to, eps_tot_to

    def _compute_error(self) -> float:
        """
        Compute current raw convergence error.

        Criterion is responsible for fetching required state from momentum_eq.
        Enables any error metric (strain, residual, displacement, composite)
        without loop coupling to specific implementations.

        Returns
        -------
        float
            Raw criterion error value for convergence checks and reporting.
        """
        return self.convergence_handler.compute_error(self.eq_mom)

    def _initialize_convergence_criterion(self) -> None:
        """
        Initialize active convergence strategy at time-step start.

        Parameters
        ----------
        None
        """
        self.convergence_handler.initialize_step(self.eq_mom)

    def _finalize_outputs(self) -> None:
        """
        Close and finalize all output handlers.

        Closes all XDMF file handles and copies mesh files to the output
        directory for provenance. This should be called at the end of a
        simulation's `run()` method.

        Returns
        -------
        None
        """
        for output in self.outputs:
            output.close()
            output.save_mesh()

    @staticmethod
    def _is_array_like(value) -> bool:
        """Return True for array-like values that expose copy/shape/dtype."""
        return (
            hasattr(value, "copy")
            and hasattr(value, "shape")
            and hasattr(value, "dtype")
        )

    @staticmethod
    def _clone_value(value):
        """Best-effort cloning utility for tensors/arrays/containers/scalars."""
        if hasattr(value, "clone"):
            try:
                return value.clone()
            except Exception:
                pass
        if hasattr(value, "copy"):
            try:
                return value.copy()
            except Exception:
                pass
        try:
            return copy.deepcopy(value)
        except Exception:
            return value

    def _snapshot_function_arrays(self, obj) -> dict:
        """Capture all dolfinx-like Function arrays found in object attributes."""
        snapshot = {}
        for attr_name, attr_value in getattr(obj, "__dict__", {}).items():
            if hasattr(attr_value, "x") and hasattr(attr_value.x, "array"):
                try:
                    snapshot[attr_name] = attr_value.x.array.copy()
                except Exception:
                    continue
        return snapshot

    def _restore_function_arrays(self, obj, snapshot: dict) -> None:
        """Restore dolfinx-like Function arrays captured by _snapshot_function_arrays."""
        for attr_name, saved_array in snapshot.items():
            attr_value = getattr(obj, attr_name, None)
            if hasattr(attr_value, "x") and hasattr(attr_value.x, "array"):
                attr_value.x.array[:] = saved_array

    def _snapshot_object_state(self, obj) -> dict:
        """Capture mutable object attributes using best-effort cloning."""
        snapshot = {}
        for attr_name, attr_value in getattr(obj, "__dict__", {}).items():
            if callable(attr_value):
                continue
            if hasattr(attr_value, "x") and hasattr(attr_value.x, "array"):
                continue
            if isinstance(attr_value, (int, float, bool, str, type(None), tuple, list, dict)):
                snapshot[attr_name] = self._clone_value(attr_value)
            elif hasattr(attr_value, "clone") or self._is_array_like(attr_value):
                snapshot[attr_name] = self._clone_value(attr_value)
        return snapshot

    def _restore_object_state(self, obj, snapshot: dict) -> None:
        """Restore object attributes captured by _snapshot_object_state."""
        for attr_name, saved_value in snapshot.items():
            setattr(obj, attr_name, self._clone_value(saved_value))

    def _snapshot_material_internal_state(self, eq_mom) -> list:
        """Capture non-elastic internal variable state from material elements."""
        material = getattr(eq_mom, "mat", None)
        if material is None or not hasattr(material, "elems_ne"):
            return []
        return [self._snapshot_object_state(elem) for elem in material.elems_ne]

    def _restore_material_internal_state(self, eq_mom, snapshot: list) -> None:
        """Restore non-elastic internal variable state for material elements."""
        material = getattr(eq_mom, "mat", None)
        if material is None or not hasattr(material, "elems_ne"):
            return
        for elem, elem_snapshot in zip(material.elems_ne, snapshot):
            self._restore_object_state(elem, elem_snapshot)

    def _snapshot_caverns_state(self) -> dict:
        """Capture cavern model mutable states."""
        caverns = getattr(self, "caverns", None)
        if caverns is None:
            return {}
        snapshot = {}
        for cavern_list_name in ("caverns_T", "caverns_PT", "caverns_MFlux"):
            cavern_list = getattr(caverns, cavern_list_name, [])
            snapshot[cavern_list_name] = [
                self._snapshot_object_state(cavern) for cavern in cavern_list
            ]
        return snapshot

    def _restore_caverns_state(self, snapshot: dict) -> None:
        """Restore cavern model mutable states."""
        caverns = getattr(self, "caverns", None)
        if caverns is None:
            return
        for cavern_list_name in ("caverns_T", "caverns_PT", "caverns_MFlux"):
            cavern_list = getattr(caverns, cavern_list_name, [])
            saved_list = snapshot.get(cavern_list_name, [])
            for cavern, cavern_snapshot in zip(cavern_list, saved_list):
                self._restore_object_state(cavern, cavern_snapshot)

    def _capture_step_state(self, include_heat: bool = False, include_caverns: bool = False) -> dict:
        """Capture state needed to rollback a failed nonlinear step attempt."""
        state = {
            "time": {
                "t": self.t_control.t,
                "dt": self.t_control.dt,
                "step_counter": self.t_control.step_counter,
            },
            "eq_mom_functions": self._snapshot_function_arrays(self.eq_mom),
            "material_state": self._snapshot_material_internal_state(self.eq_mom),
        }
        if include_heat and hasattr(self, "eq_heat"):
            state["eq_heat_functions"] = self._snapshot_function_arrays(self.eq_heat)
        if include_caverns:
            state["caverns_state"] = self._snapshot_caverns_state()
        return state

    def _restore_step_state(self, state: dict, include_heat: bool = False, include_caverns: bool = False) -> None:
        """Restore state captured by _capture_step_state."""
        time_state = state["time"]
        self.t_control.t = time_state["t"]
        self.t_control.dt = time_state["dt"]
        self.t_control.step_counter = time_state["step_counter"]

        self._restore_function_arrays(self.eq_mom, state["eq_mom_functions"])
        self._restore_material_internal_state(self.eq_mom, state["material_state"])

        if include_heat and hasattr(self, "eq_heat") and "eq_heat_functions" in state:
            self._restore_function_arrays(self.eq_heat, state["eq_heat_functions"])
        if include_caverns and "caverns_state" in state:
            self._restore_caverns_state(state["caverns_state"])
