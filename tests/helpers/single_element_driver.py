# Copyright (c) 2026, The SafeInCave Developers
#
# SPDX-License-Identifier: BSD-3-Clause

"""
Single-element constitutive driver.

Drives one (or a small batch of) plastic constitutive model(s) through a
prescribed total-strain path with no mesh, no FE solve and no global
iteration — the fastest possible plastic-model debugger. Because the strain
is prescribed, each step is strain-driven and the return mapping gives the
exact incremental response in a single evaluation.

The driver uses the models' exact-trial path: it passes ``eps_tot`` into
``compute_eps_ne_rate`` so each model reconstructs the elastic trial stress
``C : (eps_tot − Σ eps_ne_old)`` itself (including a linked partner's
committed strain, which is how the DP–Rankine pair shares one trial), then
commits the resulting increment into ``eps_ne_old``.
"""

from __future__ import annotations
import csv
import torch as to
from safeincave.Utils import dotdot_torch

_VOIGT = [(0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2)]
_VOIGT_NAMES = ["xx", "yy", "zz", "xy", "xz", "yz"]


class SingleElementDriver:
    """
    Parameters
    ----------
    elems : list
        Constitutive models (NonElasticElement subclasses) sharing the same
        batch size. For a coupled DP–Rankine pair pass both (in any order);
        their ``cede_to`` link must already be set up.
    C_voigt : torch.Tensor, optional
        Elastic stiffness (N, 6, 6) in tensorial Voigt form. Defaults to the
        first model's ``_C_voigt`` (present on the plastic models and
        asserted identical across linked models).
    """

    def __init__(self, elems, C_voigt: to.Tensor | None = None):
        self.elems = list(elems)
        if C_voigt is None:
            C_voigt = self.elems[0]._C_voigt
        self.C = C_voigt
        for elem in self.elems:
            # Over-relaxation is a device of the global staggered iteration;
            # the driver commits every increment directly, so it must be off
            # for the committed strain to be the plain return-mapping result.
            if getattr(elem, "acceleration", None) not in (None, 1.0):
                elem.acceleration = None

    def total_inelastic_strain(self) -> to.Tensor:
        eps_ne = to.zeros_like(self.elems[0].eps_ne_old)
        for elem in self.elems:
            eps_ne = eps_ne + elem.eps_ne_old
        return eps_ne

    def stress(self, eps_tot: to.Tensor) -> to.Tensor:
        """Equilibrium stress at the committed state, C:(ε_tot − Σ ε_ne)."""
        return dotdot_torch(self.C, eps_tot - self.total_inelastic_strain())

    def step(self, eps_tot: to.Tensor, Temp: to.Tensor | None = None) -> dict:
        """
        Advance one strain-driven step to total strain ``eps_tot`` (N, 3, 3):
        return-map every model on the shared exact trial, commit the
        increments, and report the committed stress and diagnostics.
        """
        n = eps_tot.shape[0]
        T = Temp if Temp is not None else to.zeros(n, dtype=to.float64)
        stress_prev = self.stress(eps_tot)
        for elem in self.elems:
            # phi1 = 1 with no cached phi2 makes the stored rate equal the
            # plastic strain increment (rate-independent models pre-divide
            # by phi2).
            elem._phi2 = None
            elem.compute_eps_ne_rate(stress_prev, 1.0, T, eps_tot=eps_tot)
        for elem in self.elems:
            elem.eps_ne_old = elem.eps_ne_old + elem.eps_ne_rate
            elem.eps_ne_rate = to.zeros_like(elem.eps_ne_rate)
            elem.eps_ne_rate_old = to.zeros_like(elem.eps_ne_rate_old)
        sigma = self.stress(eps_tot)
        record = {"eps_tot": eps_tot.clone(), "sigma": sigma}
        for elem in self.elems:
            name = getattr(elem, "name", type(elem).__name__)
            record[f"{name}_delta_lambda"] = elem.delta_lambda.clone()
            mode = getattr(elem, "yield_mode", None)
            if mode is not None:
                record[f"{name}_yield_mode"] = mode.clone()
        return record

    def run(
        self,
        eps_path: to.Tensor,
        Temp: to.Tensor | None = None,
        csv_path: str | None = None,
    ) -> list[dict]:
        """
        Run a full strain path, shape (n_steps, N, 3, 3). Optionally write a
        CSV (element 0 only) with strain/stress Voigt components and each
        model's Δλ and yield mode per step.
        """
        records = [self.step(eps_t, Temp) for eps_t in eps_path]
        if csv_path is not None:
            self._write_csv(records, csv_path)
        return records

    @staticmethod
    def _write_csv(records: list[dict], csv_path: str) -> None:
        header = ["step"]
        header += [f"eps_{c}" for c in _VOIGT_NAMES]
        header += [f"sig_{c}" for c in _VOIGT_NAMES]
        extra_keys = [
            k for k in records[0] if k not in ("eps_tot", "sigma")
        ]
        header += extra_keys
        with open(csv_path, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(header)
            for step, rec in enumerate(records):
                row = [step]
                row += [float(rec["eps_tot"][0, i, j]) for i, j in _VOIGT]
                row += [float(rec["sigma"][0, i, j]) for i, j in _VOIGT]
                row += [float(rec[k][0]) for k in extra_keys]
                writer.writerow(row)
