"""Bayesian evidence utilities for synodic-only Nordtvedt models (Steps 016, 073)."""

from __future__ import annotations

import numpy as np
from scipy.stats import gaussian_kde

from scripts.utils.llr_constants import ETA_SCALE_FACTOR

B_INTERCEPT_LO = -0.1
B_INTERCEPT_HI = 0.1
B_PRIOR_WIDTH = B_INTERCEPT_HI - B_INTERCEPT_LO

# Documented η priors for Savage–Dickey sensitivity (Step 016).
ETA_PRIOR_SPECS = [
    {
        "prior_id": "uniform_reference",
        "prior_label": "Uniform reference (pipeline default)",
        "eta_lo": -0.01,
        "eta_hi": 0.01,
        "reference": True,
    },
    {
        "prior_id": "uniform_tight",
        "prior_label": "Uniform tight (±10⁻³)",
        "eta_lo": -0.001,
        "eta_hi": 0.001,
        "reference": False,
    },
    {
        "prior_id": "uniform_wide",
        "prior_label": "Uniform wide (±10⁻¹)",
        "eta_lo": -0.1,
        "eta_hi": 0.1,
        "reference": False,
    },
    {
        "prior_id": "uniform_half_cent",
        "prior_label": "Uniform half-cent (±5×10⁻⁴)",
        "eta_lo": -5e-4,
        "eta_hi": 5e-4,
        "reference": False,
    },
]

KDE_BANDWIDTH_METHODS = (
    ("scott_default", "scott"),
    ("silverman", "silverman"),
    ("narrow_0.5x", 0.5),
    ("wide_2.0x", 2.0),
)


def uniform_prior_density_at_zero(eta_lo: float, eta_hi: float) -> float:
    """Uniform prior π(η) evaluated at η = 0 (interior point)."""
    if not (eta_lo < 0 < eta_hi):
        raise ValueError(f"η = 0 must lie strictly inside [{eta_lo}, {eta_hi}]")
    return 1.0 / (eta_hi - eta_lo)


def savage_dickey_bf(
    eta_samples: np.ndarray,
    eta_lo: float,
    eta_hi: float,
    bw_method: str | float = "scott",
) -> float:
    """BF₁₀ = π(η=0) / p(η=0|data) under the stated uniform η prior."""
    prior_at_zero = uniform_prior_density_at_zero(eta_lo, eta_hi)
    kde = gaussian_kde(eta_samples, bw_method=bw_method)
    post_at_zero = float(kde.evaluate(0.0)[0])
    if post_at_zero <= 0:
        return np.inf
    return prior_at_zero / post_at_zero


def savage_dickey_sensitivity_row(
    eta_samples: np.ndarray,
    eta_lo: float,
    eta_hi: float,
    prior_id: str,
    prior_label: str,
) -> dict:
    """One row of the documented Savage–Dickey prior × bandwidth table."""
    prior_at_zero = uniform_prior_density_at_zero(eta_lo, eta_hi)
    bandwidth_bfs = {}
    for label, bw in KDE_BANDWIDTH_METHODS:
        bandwidth_bfs[label] = float(
            savage_dickey_bf(eta_samples, eta_lo, eta_hi, bw_method=bw)
        )
    bf_vals = list(bandwidth_bfs.values())
    bf_positive = [v for v in bf_vals if v > 0 and np.isfinite(v)]
    bf_gmean = float(np.exp(np.mean(np.log(bf_positive)))) if bf_positive else np.inf
    return {
        "prior_id": prior_id,
        "prior_label": prior_label,
        "eta_prior_lo": float(eta_lo),
        "eta_prior_hi": float(eta_hi),
        "intercept_prior": f"b ∈ [{B_INTERCEPT_LO}, {B_INTERCEPT_HI}] m (uniform)",
        "prior_density_at_eta_zero": float(prior_at_zero),
        "bandwidth_bf": bandwidth_bfs,
        "bf_geometric_mean": bf_gmean,
        "bf_min": float(min(bf_vals)),
        "bf_max": float(max(bf_vals)),
        "bf_range_ratio": float(max(bf_vals) / min(bf_vals))
        if min(bf_vals) > 0
        else np.inf,
    }


def build_prior_sensitivity_table(
    posterior_by_prior: dict[str, np.ndarray],
    prior_specs: list[dict],
) -> list[dict]:
    """Combine per-prior MCMC chains into the documented sensitivity table."""
    rows = []
    for spec in prior_specs:
        key = spec["prior_id"]
        if key not in posterior_by_prior:
            raise KeyError(f"Missing MCMC samples for prior_id={key!r}")
        rows.append(
            savage_dickey_sensitivity_row(
                posterior_by_prior[key],
                spec["eta_lo"],
                spec["eta_hi"],
                spec["prior_id"],
                spec["prior_label"],
            )
        )
    return rows


def log_likelihood_synodic(
    eta: np.ndarray,
    b: np.ndarray,
    x: np.ndarray,
    y: np.ndarray,
    y_err: np.ndarray,
) -> np.ndarray:
    """
    Gaussian log-likelihood for r = 13 η cos(D) + b.

    Supports (n_samples,) η and b, or 2-D grids with matching shape.
    """
    eta = np.asarray(eta, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    log2pi = np.log(2 * np.pi)

    if eta.ndim == 1 and b.ndim == 1:
        model = ETA_SCALE_FACTOR * eta[:, np.newaxis] * x + b[:, np.newaxis]
        resid = y - model
        return -0.5 * np.sum((resid / y_err) ** 2 + log2pi + 2 * np.log(y_err), axis=1)

    if eta.shape != b.shape:
        raise ValueError("η and b must share the same shape for grid evaluation")
    model = (
        ETA_SCALE_FACTOR * eta[..., np.newaxis] * x
        + b[..., np.newaxis]
    )
    resid = y - model
    return -0.5 * np.sum(
        (resid / y_err) ** 2 + log2pi + 2 * np.log(y_err), axis=-1
    )


def log_prior_box(
    eta: np.ndarray,
    b: np.ndarray,
    eta_lo: float,
    eta_hi: float,
) -> np.ndarray:
    """Log uniform prior on (η, b) with independent boxes."""
    shape = np.broadcast(eta, b).shape
    out = np.full(shape, -np.inf, dtype=np.float64)
    ok = (
        (eta >= eta_lo)
        & (eta <= eta_hi)
        & (b >= B_INTERCEPT_LO)
        & (b <= B_INTERCEPT_HI)
    )
    log_p = -np.log((eta_hi - eta_lo) * B_PRIOR_WIDTH)
    out[ok] = log_p
    return out


def log_evidence_grid(
    x: np.ndarray,
    y: np.ndarray,
    y_err: np.ndarray,
    eta_lo: float,
    eta_hi: float,
    n_eta: int = 121,
    n_b: int = 81,
) -> dict:
    """
    Log marginal likelihood ∫ p(y|η,b) π(η,b) dη db via log-sum-exp quadrature.
    Bandwidth-free cross-check for the synodic-only model (chunked for memory).
    """
    eta_grid = np.linspace(eta_lo, eta_hi, n_eta)
    b_grid = np.linspace(B_INTERCEPT_LO, B_INTERCEPT_HI, n_b)
    deta = eta_grid[1] - eta_grid[0]
    db = b_grid[1] - b_grid[0]
    log_prior_const = -np.log((eta_hi - eta_lo) * B_PRIOR_WIDTH)

    log_posts = []
    for eta in eta_grid:
        eta_row = np.full(len(b_grid), eta)
        log_like = log_likelihood_synodic(eta_row, b_grid, x, y, y_err)
        log_posts.append(log_like + log_prior_const)
    log_z1 = float(np.logaddexp.reduce(np.concatenate(log_posts)) + np.log(deta * db))

    # M₀: η = 0, integrate b only
    eta0 = np.zeros_like(b_grid)
    log_like0 = log_likelihood_synodic(eta0, b_grid, x, y, y_err)
    log_pr0 = -np.log(B_PRIOR_WIDTH)
    log_z0 = float(np.logaddexp.reduce(log_like0 + log_pr0) + np.log(db))

    return {
        "log_evidence_m1": log_z1,
        "log_evidence_m0": log_z0,
        "log_bf10": log_z1 - log_z0,
        "bf10": float(np.exp(log_z1 - log_z0)),
        "n_eta": n_eta,
        "n_b": n_b,
    }


def _log_posterior_chunks(
    eta: np.ndarray,
    b: np.ndarray,
    x: np.ndarray,
    y: np.ndarray,
    y_err: np.ndarray,
    eta_lo: float,
    eta_hi: float,
    chunk: int = 512,
) -> np.ndarray:
    """Chunked M₁ log posterior to avoid O(N_samples × N_obs) memory spikes."""
    out = np.empty(len(eta), dtype=np.float64)
    for s in range(0, len(eta), chunk):
        e = min(s + chunk, len(eta))
        out[s:e] = log_likelihood_synodic(eta[s:e], b[s:e], x, y, y_err) + log_prior_box(
            eta[s:e], b[s:e], eta_lo, eta_hi
        )
    return out


def bridge_sampling_synodic(
    x: np.ndarray,
    y: np.ndarray,
    y_err: np.ndarray,
    eta_lo: float,
    eta_hi: float,
    mcmc_samples_m1: np.ndarray,
    n_iter: int = 5,
    max_samples: int = 4000,
) -> dict:
    """
    Iterative bridge sampling (Meng & Wong 1996) between M₁ (η,b) and M₀ (η=0,b).
    Subsamples MCMC draws for stable memory use on large N.
    """
    if len(mcmc_samples_m1) > max_samples:
        rng = np.random.default_rng(42)
        idx = rng.choice(len(mcmc_samples_m1), size=max_samples, replace=False)
        samples = mcmc_samples_m1[idx]
    else:
        samples = mcmc_samples_m1

    eta_s = samples[:, 0]
    b_s = samples[:, 1]
    n = len(eta_s)

    w = 1.0 / y_err**2
    b_mle = np.sum(w * y) / np.sum(w)
    b_spread = np.std(y - b_mle) * 0.05
    rng = np.random.default_rng(43)
    b_m0 = np.clip(
        rng.normal(b_mle, max(b_spread, 1e-5), size=n),
        B_INTERCEPT_LO,
        B_INTERCEPT_HI,
    )

    log_p1 = _log_posterior_chunks(eta_s, b_s, x, y, y_err, eta_lo, eta_hi)
    log_p0 = log_likelihood_synodic(np.zeros(n), b_m0, x, y, y_err) - np.log(B_PRIOR_WIDTH)

    c = 0.5
    for _ in range(n_iter):
        log_m1 = np.logaddexp(np.log(1 - c) + log_p1, np.log(c) + log_p0)
        log_m0 = np.logaddexp(np.log(c) + log_p1, np.log(1 - c) + log_p0)
        log_num = np.logaddexp.reduce(log_p1 - log_m1) - np.log(n)
        log_den = np.logaddexp.reduce(log_p0 - log_m0) - np.log(n)
        log_bf = log_num - log_den
        c = 1.0 / (1.0 + np.exp(log_bf))

    return {
        "log_bf10": float(log_bf),
        "bf10": float(np.exp(log_bf)),
        "bridge_mixture_c": float(c),
        "n_samples": n,
        "n_iterations": n_iter,
        "max_samples_used": max_samples,
    }
