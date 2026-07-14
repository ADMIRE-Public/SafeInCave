# Copyright (c) 2026, The SafeInCave Developers
#
# SPDX-License-Identifier: BSD-3-Clause

"""
Reusable finite-difference probe for algorithmic (consistent) tangents.

The probe compares a constitutive model's analytic tangent
``G_alg = d(Δε_p)/dσ_trial`` against a central finite difference of its
return mapping. A silent Voigt-convention error in the tangent degrades the
global iteration's convergence rate without producing wrong-looking results,
so this probe is the merge gate for any tangent change.

Voigt convention (must match ``Utils.dotdot_torch`` operators): order
``[xx, yy, zz, xy, xz, yz]``; rows carry plain tensor components; columns
carry doubled shear entries. The doubling arises because both symmetric
stress entries vary together, which the probe reproduces by perturbing
``σ[i, j]`` and ``σ[j, i]`` simultaneously (this also keeps the perturbed
tensor symmetric — required by eigendecomposition-based return mappings).
"""

from __future__ import annotations
import torch as to

# Voigt order [xx, yy, zz, xy, xz, yz]
_ROW_I = [0, 1, 2, 0, 0, 1]
_ROW_J = [0, 1, 2, 1, 2, 2]
_VOIGT_PAIRS = [(0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2)]


def to_voigt_rows(tensor: to.Tensor) -> to.Tensor:
    """Plain-component Voigt row vector of an (N, 3, 3) tensor, shape (N, 6)."""
    return tensor[:, _ROW_I, _ROW_J]


def fd_tangent(
    return_map, stress_trial: to.Tensor, rel_step: float = 1e-6
) -> to.Tensor:
    """
    Central-difference tangent d(Δε_p)/dσ_trial of a return mapping.

    Parameters
    ----------
    return_map : Callable[[torch.Tensor], torch.Tensor]
        Maps trial stress (N, 3, 3) to plastic strain increment (N, 3, 3).
        Must be stateless with respect to repeated calls.
    stress_trial : torch.Tensor
        Trial stresses to probe at, shape (N, 3, 3). Each state must sit in
        the interior of its return-mapping region — finite differences
        across a region boundary (yield-surface kink) are meaningless.
    rel_step : float
        Perturbation relative to the largest stress magnitude.

    Returns
    -------
    torch.Tensor
        Tangent (N, 6, 6): rows plain components, columns doubled shear.
    """
    n = stress_trial.shape[0]
    scale = float(to.amax(stress_trial.abs()))
    h = rel_step * max(scale, 1e-30)
    G = to.zeros((n, 6, 6), dtype=to.float64)
    sp = stress_trial.clone()
    for k, (i, j) in enumerate(_VOIGT_PAIRS):
        sp[:, i, j] += h
        if i != j:
            sp[:, j, i] += h
        eps_a = to_voigt_rows(return_map(sp))
        sp[:, i, j] -= 2.0 * h
        if i != j:
            sp[:, j, i] -= 2.0 * h
        eps_b = to_voigt_rows(return_map(sp))
        sp[:, i, j] += h
        if i != j:
            sp[:, j, i] += h
        G[:, :, k] = (eps_a - eps_b) / (2.0 * h)
    return G


def element_return_map(elem, Temp: to.Tensor | None = None):
    """
    Wrap a NonElasticElement's return mapping as Δε_p(σ_trial).

    Uses ``compute_eps_ne_rate(..., return_eps_ne=True)`` on the stale-rate
    fallback path with a zeroed stored rate, so the internal trial stress
    equals the input, and rescales by the cached φ₂ so the result is the
    plastic strain *increment* regardless of prior ``compute_G_B`` calls.
    """

    def _return_map(sig: to.Tensor) -> to.Tensor:
        n = sig.shape[0]
        T = Temp if Temp is not None else to.zeros(n, dtype=to.float64)
        # The stale-rate fallback adds C:(phi2*eps_ne_rate) to the input to
        # form its trial; zero the stored rate so the trial IS the input
        # (a prior stateful call may have left a nonzero rate behind).
        elem.eps_ne_rate = to.zeros_like(elem.eps_ne_rate)
        phi2 = elem._phi2 if elem._phi2 is not None else 1.0
        return phi2 * elem.compute_eps_ne_rate(sig, 1.0, T, return_eps_ne=True)

    return _return_map


def element_analytic_tangent(
    elem, stress_trial: to.Tensor, Temp: to.Tensor | None = None
) -> to.Tensor:
    """
    The model's own algorithmic tangent G_alg at ``stress_trial``.

    Stateful: runs the return mapping (to populate the model's trial cache)
    followed by ``compute_G_B`` with dt=2, θ=0.5 so φ₂ = 1 and the stored
    ``G`` equals ``tangent_penalty · G_alg``; the penalty is divided out.
    The stored rate is zeroed first so the fallback trial equals the input.
    """
    n = stress_trial.shape[0]
    T = Temp if Temp is not None else to.zeros(n, dtype=to.float64)
    elem.eps_ne_rate = to.zeros_like(elem.eps_ne_rate)
    elem._phi2 = None
    elem.compute_eps_ne_rate(stress_trial, 1.0, T)
    elem.compute_G_B(stress_trial, 2.0, 0.5, T)
    return elem.G / elem.tangent_penalty


def element_stress_map(elem, Temp: to.Tensor | None = None, dt: float = 1.0):
    """
    Wrap the Newton interface as σ(eps_el_trial) for FD probing of D_ep.
    compute_stress_and_tangent is stateless across calls by contract.
    """

    def _stress_map(eps_el_trial: to.Tensor) -> to.Tensor:
        n = eps_el_trial.shape[0]
        T = Temp if Temp is not None else to.zeros(n, dtype=to.float64)
        sigma, _ = elem.compute_stress_and_tangent(eps_el_trial, T, dt)
        return sigma

    return _stress_map


def compare_element_dep(
    elem,
    eps_el_trial: to.Tensor,
    Temp: to.Tensor | None = None,
    dt: float = 1.0,
    rel_step: float = 1e-6,
) -> to.Tensor:
    """
    Probe the Newton interface: per-element relative error between the FD
    derivative dσ/dε of compute_stress_and_tangent and its returned D_ep.
    This is the merge gate for any change to the consistent tangent.
    """
    n = eps_el_trial.shape[0]
    T = Temp if Temp is not None else to.zeros(n, dtype=to.float64)
    _, D_ep = elem.compute_stress_and_tangent(eps_el_trial, T, dt)
    D_fd = fd_tangent(element_stress_map(elem, Temp, dt), eps_el_trial, rel_step)
    return tangent_error(D_fd, D_ep)


def tangent_error(G_fd: to.Tensor, G_analytic: to.Tensor) -> to.Tensor:
    """
    Per-element relative Frobenius error between two (N, 6, 6) tangents,
    normalized by the batch-wide largest tangent norm so that zero-tangent
    (elastic) elements report the error on an absolute scale.
    """
    diff = to.linalg.matrix_norm(G_fd - G_analytic)
    ref = to.maximum(to.linalg.matrix_norm(G_fd), to.linalg.matrix_norm(G_analytic))
    ref_floor = max(float(ref.max()), 1e-300)
    return diff / to.clamp(ref, min=ref_floor)


def compare_element_tangent(
    elem,
    stress_trial: to.Tensor,
    Temp: to.Tensor | None = None,
    rel_step: float = 1e-6,
) -> to.Tensor:
    """
    Probe a model at the given trial states: per-element relative error
    between the FD tangent of its return mapping and its analytic tangent.
    """
    G_analytic = element_analytic_tangent(elem, stress_trial, Temp)
    G_fd = fd_tangent(element_return_map(elem, Temp), stress_trial, rel_step)
    return tangent_error(G_fd, G_analytic)
