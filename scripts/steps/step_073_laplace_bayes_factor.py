#!/usr/bin/env python3
"""Step 073: Laplace BF, grid/nested/bridge cross-checks, posterior sign probability."""
import argparse
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

import emcee
import numpy as np
import pandas as pd

from scripts.utils.bayesian_evidence import bridge_sampling_synodic, log_evidence_grid
from scripts.utils.config import get_config
from scripts.utils.llr_constants import ETA_SCALE_FACTOR
from scripts.utils.logger import TEPLogger, print_status, set_step_logger, set_verbose_mode
from scripts.utils.statistical_utils import detect_outliers_sigma

TEP_CONFIG = get_config()
ETA_LO, ETA_HI = -0.01, 0.01


def log_prob(theta, x, y, y_err):
    eta, b = theta
    if not (ETA_LO <= eta <= ETA_HI and -0.1 <= b <= 0.1):
        return -np.inf
    model = ETA_SCALE_FACTOR * eta * x + b
    return -0.5 * np.sum((y - model) ** 2 / y_err**2 + np.log(2 * np.pi) + 2 * np.log(y_err))


def log_prob_batch(theta, x, y, y_err):
    """Vector log-density for emcee (rows of theta are walkers)."""
    eta = theta[:, 0]
    b = theta[:, 1]
    ll = np.full(eta.shape[0], -np.inf, dtype=np.float64)
    ok = (ETA_LO <= eta) & (eta <= ETA_HI) & (-0.1 <= b) & (b <= 0.1)
    if not np.any(ok):
        return ll
    model = ETA_SCALE_FACTOR * np.outer(eta[ok], x) + b[ok, np.newaxis]
    ll[ok] = -0.5 * np.sum((y - model) ** 2 / y_err**2 + np.log(2 * np.pi) + 2 * np.log(y_err), axis=1)
    return ll


def run_laplace_bf(n_walkers=32, n_steps=5000, burn_in=1000, verbose=False):
    print_status("BAYESIAN EVIDENCE CROSS-CHECKS (Step 073)", "INFO")
    df = pd.read_csv(
        PROJECT_ROOT / "data" / "processed" / "INPOP19a_all_stations_residuals.csv"
    )
    mask = ~detect_outliers_sigma(df["residual_m"].values, 6.0)
    df = df[mask].copy()
    x = df["cos_elong_rad"].values
    y = df["residual_m"].values
    y_err = np.clip(df["sigma_m"].values.astype(float), 1e-6, None)

    np.random.seed(TEP_CONFIG.get("RANDOM_SEED", 42))
    X = np.column_stack([ETA_SCALE_FACTOR * x, np.ones_like(x)])
    ols, _, _, _ = np.linalg.lstsq(
        X / y_err[:, np.newaxis], y / y_err, rcond=None
    )
    pos = np.array([ols[0], ols[1]]) + 1e-4 * np.random.randn(n_walkers, 2)

    sampler = emcee.EnsembleSampler(n_walkers, 2, log_prob, args=(x, y, y_err))
    sampler.run_mcmc(pos, n_steps, progress=verbose)
    samples = sampler.get_chain(discard=burn_in, flat=True)
    eta_s = samples[:, 0]

    p_negative = float(np.mean(eta_s < 0))
    p_positive = float(np.mean(eta_s > 0))

    row_chunk = 256
    logp = np.concatenate(
        [
            log_prob_batch(samples[s : s + row_chunk], x, y, y_err)
            for s in range(0, len(samples), row_chunk)
        ]
    )

    # Deterministic MAP via weighted least squares (exact for uniform prior + Gaussian likelihood)
    # Design matrix X = [13*cos(D), 1], weights w = 1/y_err^2
    # MAP = (X^T W X)^{-1} X^T W y  (computed without forming full diagonal W)
    X = np.column_stack([ETA_SCALE_FACTOR * x, np.ones_like(x)])
    w_vec = 1.0 / y_err**2
    with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
        XtWX = X.T @ (w_vec[:, np.newaxis] * X)
        XtWy = X.T @ (w_vec * y)
    map_eta, map_b = np.linalg.solve(XtWX, XtWy)

    # Inverse Hessian = (X^T W X)^{-1} (exact curvature of the log-posterior at MAP)
    inv_hessian = np.linalg.inv(XtWX)
    det_inv_hessian = np.linalg.det(inv_hessian)

    # Log-likelihood at MAP (with full Gaussian normalisation, matching log_likelihood_synodic)
    model_map = ETA_SCALE_FACTOR * map_eta * x + map_b
    loglik_map = -0.5 * np.sum((y - model_map) ** 2 / y_err**2 + np.log(2 * np.pi) + 2 * np.log(y_err))

    d = 2
    log_ml_laplace = loglik_map + 0.5 * d * np.log(2 * np.pi) + 0.5 * np.log(det_inv_hessian)

    w = 1.0 / y_err**2
    b_mle_m0 = np.sum(w * y) / np.sum(w)
    logp_m0 = -0.5 * np.sum((y - b_mle_m0) ** 2 / y_err**2 + np.log(2 * np.pi) + 2 * np.log(y_err))
    log_prior_m1 = -np.log(0.02 * 0.2)
    log_prior_m0 = -np.log(0.2)
    # Laplace volume correction for M0 (1 parameter: intercept b).
    # Posterior variance sigma_b^2 = 1 / sum(w_i) for weighted mean with known errors.
    sigma_b_sq = 1.0 / np.sum(w)
    log_ml_m0 = logp_m0 + log_prior_m0 + 0.5 * np.log(2 * np.pi) + 0.5 * np.log(sigma_b_sq)
    log_bf = log_ml_laplace + log_prior_m1 - log_ml_m0
    bf_laplace = float(np.exp(log_bf))

    # BIC from heteroskedastic log-likelihood (consistent with Laplace/grid/bridge models)
    n = len(y)
    bic_m1 = -2 * loglik_map + 2 * np.log(n)
    bic_m0 = -2 * logp_m0 + 1 * np.log(n)
    bf_bic = float(np.exp(0.5 * (bic_m0 - bic_m1)))

    grid_evidence = log_evidence_grid(x, y, y_err, ETA_LO, ETA_HI)
    bridge_evidence = bridge_sampling_synodic(
        x, y, y_err, ETA_LO, ETA_HI, samples, n_iter=8
    )

    results = {
        "step_id": "step_073",
        "status": "PASS",
        "model": "synodic_only: r = 13 η cos(D) + b",
        "prior_eta": f"uniform on [{ETA_LO}, {ETA_HI}]",
        "prior_intercept": "uniform on [-0.1, 0.1] m",
        "posterior_prob_negative_eta": p_negative,
        "posterior_prob_positive_eta": p_positive,
        "map_eta": float(map_eta),
        "map_b": float(map_b),
        "map_method": "weighted_least_squares",
        "inv_hessian_det": float(det_inv_hessian),
        "log_ml_laplace": float(log_ml_laplace),
        "log_bf_laplace": float(log_bf),
        "bf_laplace": bf_laplace,
        "bf_bic": bf_bic,
        "evidence_cross_checks": {
            "grid_quadrature": grid_evidence,
            "bridge_sampling": bridge_evidence,
            "note": (
                "Grid quadrature is the primary bandwidth-free cross-check. Laplace "
                "uses the exact WLS MAP and inverse Hessian (X^T W X)^{-1}. BIC uses the "
                "heteroskedastic log-likelihood for consistency. Bridge sampling is "
                "unstable for this nested model (M0 ⊂ M1) with BF ≫ 1 because the M1 "
                "posterior (η ≈ -3e-4) and M0 posterior (η = 0) barely overlap; the "
                "Savage-Dickey ratio in Step 016 is the recommended bridge estimator."
            ),
        },
        "n_samples": int(len(samples)),
        "n_obs": int(n),
    }

    print_status(f"P(η<0|data) = {p_negative:.6f}", "CALC")
    print_status(f"Laplace BF₁₀ = {bf_laplace:.3e}", "CALC")
    print_status(f"BIC BF₁₀ = {bf_bic:.3e}", "CALC")
    print_status(f"Grid BF₁₀ = {grid_evidence['bf10']:.3e}", "CALC")
    print_status(f"Bridge BF₁₀ = {bridge_evidence['bf10']:.3e}", "CALC")

    log_dir = PROJECT_ROOT / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    logger = TEPLogger("step_073", str(log_dir / "step_073_laplace_bayes_factor.log"))
    set_step_logger(logger)
    logger.save_step_results(results, PROJECT_ROOT, "step_073_laplace_bayes_factor")
    print_status("✓ Complete", "SUCCESS")
    return results


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()
    set_verbose_mode(args.verbose)
    run_laplace_bf(
        n_walkers=TEP_CONFIG["MCMC_STANDARD_WALKERS"],
        n_steps=TEP_CONFIG["MCMC_STANDARD_STEPS"],
        burn_in=TEP_CONFIG["MCMC_BURN_IN"],
        verbose=args.verbose,
    )
