#!/usr/bin/env python3
"""
Step 062: False-Positive Rate Simulation (Parametric Bootstrap under GR Null)
==============================================================================

Computes the exact false-positive rate for the primary full-systematic η
estimand by parametric bootstrap under the GR null hypothesis (η = 0).

Method:
1. Fit the full-systematic model (cosD + annual + monthly + thermal cos2D)
   to the real data to obtain the null residuals (with cosD partialled out).
2. Fit an AR(1) model to the null residuals to capture temporal correlation.
3. Generate 10,000 synthetic datasets by simulating AR(1) noise with the
   fitted ρ and residual variance, adding back nuisance structure only.
4. For each synthetic dataset, fit the full-systematic model and record η.
5. The false-positive rate is the fraction of simulations where |η| ≥ |η_obs|.

This provides an exact, data-driven p-value that accounts for temporal
autocorrelation (ρ ≈ 0.413) and the full nuisance design matrix.
"""

import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

import argparse
import numpy as np
import pandas as pd
from scipy import linalg as spla

from scripts.utils.logger import TEPLogger, set_step_logger, print_status
from scripts.utils.llr_constants import ETA_SCALE_FACTOR
from scripts.utils.statistical_utils import detect_outliers_sigma, robust_regression
from scripts.utils.config import get_config
from scripts.utils.numerics import suppress_scipy_array_api_matmul_runtime_warning

TEP_CONFIG = get_config()


def fit_full_systematic_model(df):
    """Build design matrix and fit the full-systematic model."""
    t = df['date_julian'].values
    cosD = df['cos_elong_rad'].values
    y = df['residual_m'].values

    cos2D = np.cos(2.0 * np.arccos(np.clip(cosD, -1.0, 1.0)))
    monthly_sin = np.sin(2.0 * np.pi * t / 27.32)
    monthly_cos = np.cos(2.0 * np.pi * t / 27.32)
    annual_sin = np.sin(2.0 * np.pi * t / 365.25)
    annual_cos = np.cos(2.0 * np.pi * t / 365.25)

    X = np.column_stack([cosD, cos2D, monthly_sin, monthly_cos,
                         annual_sin, annual_cos, np.ones(len(y))])
    names = ['cosD', 'cos2D', 'monthly_sin', 'monthly_cos',
             'annual_sin', 'annual_cos', 'const']

    reg = robust_regression(y, X, weights=None, scale_errors_by_birge=False)
    coeffs = reg['coefficients']
    with suppress_scipy_array_api_matmul_runtime_warning(), \
         np.errstate(over='ignore', divide='ignore', invalid='ignore'):
        resid = y - X @ coeffs

    return {
        'X': X, 'y': y, 'names': names, 'coeffs': coeffs,
        'resid': resid, 'cosD': cosD, 't': t,
        'eta_obs': coeffs[0] / ETA_SCALE_FACTOR,
        'eta_err_obs': reg['errors'][0] / ETA_SCALE_FACTOR,
    }


def estimate_ar1(residuals):
    """Estimate AR(1) parameter from residuals."""
    rho = np.sum(residuals[1:] * residuals[:-1]) / np.sum(residuals[:-1]**2)
    sigma_eps = np.std(residuals[1:] - rho * residuals[:-1], ddof=1)
    return rho, sigma_eps


def simulate_ar1_noise(n, rho, sigma_eps, seed=None):
    """Generate AR(1) noise sequence."""
    rng = np.random.default_rng(seed)
    eps = rng.normal(0, sigma_eps, size=n)
    noise = np.zeros(n)
    noise[0] = eps[0] / np.sqrt(1 - rho**2) if abs(rho) < 1 else eps[0]
    for i in range(1, n):
        noise[i] = rho * noise[i-1] + eps[i]
    return noise


def _simulate_ar1_noise_matrix(n: int, rho: float, sigma_eps: float, n_sim: int,
                                seed: int = None) -> np.ndarray:
    """
    n_sim independent AR(1) columns with the same (rho, sigma_eps).
    Matches simulate_ar1_noise column-wise: column j uses default_rng(seed + j).
    """
    base_seed = seed if seed is not None else 42
    noise = np.zeros((n, n_sim), dtype=np.float64)
    eps = np.empty((n, n_sim), dtype=np.float64)
    for j in range(n_sim):
        rng = np.random.default_rng(base_seed + j)
        eps[:, j] = rng.normal(0, sigma_eps, size=n)
    if abs(rho) < 1.0:
        scale0 = 1.0 / np.sqrt(1.0 - rho**2)
        noise[0] = eps[0] * scale0
    else:
        noise[0] = eps[0]
    for t in range(1, n):
        noise[t] = rho * noise[t - 1] + eps[t]
    return noise


def _ols_first_coef_from_qr(
    y_matrix: np.ndarray,
    q: np.ndarray,
    r: np.ndarray,
) -> np.ndarray:
    """Return first regression coefficient for each column of y (same fixed design)."""
    qty = q.T @ y_matrix
    beta = spla.solve_triangular(r, qty, check_finite=False)
    return beta[0, :]


def run_false_positive_simulation(n_simulations=10000, verbose=False):
    """Run the false-positive rate simulation."""
    print_status("=" * 60, "INFO")
    print_status("FALSE-POSITIVE RATE SIMULATION (Step 062)", "TITLE")
    print_status("=" * 60, "INFO")

    input_path = PROJECT_ROOT / "data" / "processed" / \
        "INPOP19a_all_stations_residuals.csv"
    df = pd.read_csv(input_path)

    outlier_mask = detect_outliers_sigma(df['residual_m'].values, sigma_threshold=6.0)
    n_outliers = int(np.sum(outlier_mask))
    df_clean = df[~outlier_mask].sort_values(
        ["date_julian", "station"], kind="mergesort"
    ).reset_index(drop=True)

    print_status(f"Data: {len(df_clean):,} observations after "
                 f"6σ MAD cleaning ({n_outliers} outliers removed)", "INFO")

    fit = fit_full_systematic_model(df_clean)
    eta_obs = fit['eta_obs']
    eta_err_obs = fit['eta_err_obs']
    snr_obs = abs(eta_obs) / eta_err_obs

    print_status(f"Observed η = {eta_obs:.6e} ± {eta_err_obs:.6e} "
                 f"(SNR = {snr_obs:.2f}σ)", "CALC")

    X_null = fit['X'][:, 1:]
    coeffs_null = fit['coeffs'][1:]
    with suppress_scipy_array_api_matmul_runtime_warning(), \
         np.errstate(over='ignore', divide='ignore', invalid='ignore'):
        null_residuals = fit['y'] - X_null @ coeffs_null

    rho, sigma_eps = estimate_ar1(null_residuals)
    print_status(f"Null-residual AR(1): ρ = {rho:.4f}, σε = {sigma_eps:.4f} m", "CALC")

    n = len(df_clean)
    X_full = np.asarray(fit['X'], dtype=np.float64)

    sqrt_w = np.ones(n, dtype=np.float64)
    xw = X_full * sqrt_w[:, np.newaxis]
    q, r_mat = spla.qr(xw, mode="economic", check_finite=False)
    s = spla.svdvals(r_mat)
    cond = float(s[0] / s[-1]) if s.size and s[-1] > 0 else np.inf
    if cond > 1e12:
        raise RuntimeError(
            f"Step 062: full-systematic design ill-conditioned (κ={cond:.2e}); "
            "refuse batch shortcut — inspect data / design."
        )

    print_status(f"Running {n_simulations:,} simulations (vectorized AR + QR batch)...", "PROCESS")
    eta_sim = np.empty(n_simulations, dtype=float)
    y0 = (X_null @ coeffs_null).astype(np.float64, copy=False)
    chunk = min(2048, max(256, n_simulations // 8))
    done = 0
    while done < n_simulations:
        mloc = min(chunk, n_simulations - done)
        noise_blk = _simulate_ar1_noise_matrix(
            n, rho, sigma_eps, mloc,
            seed=TEP_CONFIG.get("RANDOM_SEED", 42))
        with suppress_scipy_array_api_matmul_runtime_warning(), \
             np.errstate(over='ignore', divide='ignore', invalid='ignore'):
            y_blk = y0[:, np.newaxis] + noise_blk
        eta_sim[done : done + mloc] = _ols_first_coef_from_qr(y_blk, q, r_mat) / ETA_SCALE_FACTOR
        done += mloc
        if verbose and done % max(1000, chunk) == 0:
            print_status(f"  Completed {done:,}/{n_simulations:,}", "INFO")

    n_extreme = int(np.sum(np.abs(eta_sim) >= np.abs(eta_obs)))
    p_value = n_extreme / n_simulations
    # Clopper-Pearson exact 95% upper bound for binomial proportion
    # For 0 successes in n trials: p_upper = 1 - alpha^(1/n)
    p_value_upper = 1.0 - (0.05 ** (1.0 / n_simulations))

    mean_null = float(np.mean(eta_sim))
    std_null = float(np.std(eta_sim, ddof=1))
    ci_95_null = [float(np.percentile(eta_sim, 2.5)),
                  float(np.percentile(eta_sim, 97.5))]
    ci_99_null = [float(np.percentile(eta_sim, 0.5)),
                  float(np.percentile(eta_sim, 99.5))]

    sigma_equiv = abs(eta_obs) / std_null if std_null > 0 else float('inf')

    print_status("", "INFO")
    print_status("RESULTS", "TITLE")
    print_status(f"  Simulations: {n_simulations:,}", "CALC")
    print_status(f"  |η| ≥ |η_obs| in {n_extreme} / {n_simulations:,} trials", "CALC")
    print_status(f"  False-positive rate (exact p): {p_value:.6e}", "CALC")
    print_status(f"  Conservative upper bound:     {p_value_upper:.6e}", "CALC")
    print_status(f"  Null η mean:  {mean_null:.6e}", "CALC")
    print_status(f"  Null η std:   {std_null:.6e}", "CALC")
    print_status(f"  Null 95% CI:  [{ci_95_null[0]:.6e}, {ci_95_null[1]:.6e}]", "CALC")
    print_status(f"  Null 99% CI:  [{ci_99_null[0]:.6e}, {ci_99_null[1]:.6e}]", "CALC")
    print_status(f"  Equivalent Gaussian σ: {sigma_equiv:.2f}", "CALC")

    if p_value == 0.0:
        print_status("  ZERO false positives in {:,} trials".format(n_simulations),
                     "SUCCESS")
        print_status("  p < 1/{:,} = {:.2e}".format(n_simulations, 1.0/n_simulations),
                     "SUCCESS")

    results = {
        "step_id": "step_062",
        "n_simulations": n_simulations,
        "observed_eta": float(eta_obs),
        "observed_eta_error": float(eta_err_obs),
        "observed_snr": float(snr_obs),
        "null_ar1_rho": float(rho),
        "null_ar1_sigma_eps": float(sigma_eps),
        "n_extreme": int(n_extreme),
        "false_positive_rate": float(p_value),
        "false_positive_rate_upper_bound": float(p_value_upper),
        "null_eta_mean": mean_null,
        "null_eta_std": std_null,
        "null_eta_ci95": ci_95_null,
        "null_eta_ci99": ci_99_null,
        "equivalent_gaussian_sigma": float(sigma_equiv),
        "n_observations": int(n),
        "n_outliers_removed": int(n_outliers),
        "method": "Parametric bootstrap under GR null with AR(1) noise model",
        "status": "PASS" if p_value < 0.05 else "WARNING"
    }

    return results


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Step 062: False-Positive Rate Simulation")
    parser.add_argument("--n-simulations", type=int, default=10000,
                        help="Number of bootstrap iterations (default: 10000)")
    parser.add_argument("--verbose", action="store_true",
                        help="Print progress every 1000 simulations")
    args = parser.parse_args()

    log_dir = PROJECT_ROOT / "logs"
    logger = TEPLogger("step_062", str(
        log_dir / "step_062_false_positive_simulation.log"))
    set_step_logger(logger)

    summary = run_false_positive_simulation(
        n_simulations=args.n_simulations, verbose=args.verbose)

    logger.save_step_results(summary, PROJECT_ROOT,
                             "step_062_false_positive_simulation")
    print_status("False-Positive Simulation Complete.", "SUCCESS")
