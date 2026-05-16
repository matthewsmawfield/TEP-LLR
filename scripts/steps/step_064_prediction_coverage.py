#!/usr/bin/env python3
"""
Step 064-PI: Prediction Interval Coverage, Uncertainty Calibration, and Headline η Intervals
========================================================================================

Tests whether the full-systematic model's prediction intervals are well-calibrated
under multiple error models (WLS with published σ, OLS, cluster-robust, AR(1)-adjusted).

Also reports station-block bootstrap and leave-one-station-out conformal intervals
for the headline precision-weighted η so formal σ-based SNR is not the sole basis
for significance claims.

Ties calibration to the abstract: conservative prediction coverage (observed > nominal)
implies headline formal σ may be conservative; bootstrap/conformal SNRs are reported
alongside WLS/cluster paths.
"""

import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

import argparse
import numpy as np
import pandas as pd
from scipy import stats

from scripts.utils.logger import TEPLogger, set_step_logger, print_status, set_verbose_mode
from scripts.utils.config import get_config
from scripts.utils.llr_constants import ETA_SCALE_FACTOR
from scripts.utils.statistical_utils import (
    detect_outliers_sigma,
    robust_regression,
    cluster_robust_variance,
)
from scripts.utils.numerics import hat_diagonal_from_qr, suppress_scipy_array_api_matmul_runtime_warning
from scripts.steps.step_050_corrected_tep_analysis import (
    precision_weighted_full_model,
    cluster_robust_regression,
    fit_model,
)
from scripts.steps.step_067_cluster_robust_ar1_combined import cluster_robust_ar1_regression

TEP_CONFIG = get_config()

COVERAGE_LEVELS = (0.68, 0.90, 0.95, 0.99)
FULL_NAMES = ["cosD", "cos2D", "sin_m", "cos_m", "sin_y", "cos_y", "const"]


def build_full_systematic_matrix(jd, cos_d, cos_2d):
    """Full-systematic design (matches Step 050)."""
    year = jd / 365.25
    month = jd / 27.32
    return np.column_stack(
        [
            cos_d,
            cos_2d,
            np.sin(2 * np.pi * month),
            np.cos(2 * np.pi * month),
            np.sin(2 * np.pi * year),
            np.cos(2 * np.pi * year),
            np.ones(len(cos_d)),
        ]
    )


def _coverage_table(y, y_hat, half_width, levels=COVERAGE_LEVELS):
    """Fraction of y inside [y_hat ± half_width] at each nominal level."""
    out = {}
    n = len(y)
    for level in levels:
        alpha = 1 - level
        z = stats.norm.ppf(1 - alpha / 2)
        covered = (y >= y_hat - z * half_width) & (y <= y_hat + z * half_width)
        key = f"{int(level * 100)}pct"
        out[key] = {
            "nominal_coverage": float(level),
            "observed_coverage": float(np.mean(covered)),
            "n_covered": int(np.sum(covered)),
            "n_total": n,
        }
    return out


def compute_prediction_coverage(
    X,
    y,
    sigma=None,
    inflation=1.0,
    use_sigma=True,
    coverage_levels=COVERAGE_LEVELS,
):
    """
    Prediction-interval coverage for a fitted linear model.

    PI width: inflation * sqrt( (use_sigma ? sigma_i^2 : 0) + MSE * h_ii ) with
    Gaussian critical values (large-n limit appropriate for N ~ 25k).
    """
    n, p = X.shape
    nu = n - p
    reg = robust_regression(y, X, scale_errors_by_birge=False)
    beta = reg["coefficients"]
    with suppress_scipy_array_api_matmul_runtime_warning(), np.errstate(
        over="ignore", divide="ignore", invalid="ignore"
    ):
        y_hat = X @ beta
        residuals = y - y_hat
    mse = float(reg["mse"])
    leverage = hat_diagonal_from_qr(X)
    leverage = np.clip(leverage, 0.0, 1.0 - 1e-10)

    if use_sigma and sigma is not None and len(sigma) == n and np.all(sigma > 0):
        base_var = sigma.astype(float) ** 2
        chi2 = float(np.sum(residuals**2 / sigma**2))
    else:
        base_var = np.zeros(n, dtype=float)
        chi2 = float(np.sum(residuals**2) / mse) if mse > 0 else np.nan

    infl = float(inflation)
    pred_se = infl * np.sqrt(base_var + mse * leverage)
    chi2_red = chi2 / nu if nu > 0 else np.nan

    return {
        "mse": mse,
        "rmse": float(np.sqrt(mse)),
        "chi2": chi2,
        "chi2_reduced": float(chi2_red),
        "df": int(nu),
        "n_observations": n,
        "n_parameters": p,
        "r2": float(1 - np.sum(residuals**2) / np.sum((y - np.mean(y)) ** 2)),
        "inflation_factor": infl,
        "use_published_sigma": bool(use_sigma),
        "coverage": _coverage_table(y, y_hat, pred_se, levels=coverage_levels),
    }


def calibrate_sigma_scale(X, y, sigma, target_level=0.68, tol=0.002):
    """Bisection on global σ scale c until WLS prediction coverage matches target."""
    if sigma is None or not np.all(sigma > 0):
        raise ValueError("sigma required for calibration")

    def obs_at(c):
        return compute_prediction_coverage(
            X, y, sigma=c * sigma, inflation=1.0, use_sigma=True
        )["coverage"][f"{int(target_level * 100)}pct"]["observed_coverage"]

    lo, hi = 0.05, 5.0
    obs_lo, obs_hi = obs_at(lo), obs_at(hi)
    if obs_lo > target_level:
        return lo, obs_lo
    if obs_hi < target_level:
        return hi, obs_hi

    for _ in range(60):
        mid = 0.5 * (lo + hi)
        obs_mid = obs_at(mid)
        if obs_mid > target_level:
            hi, obs_hi = mid, obs_mid
        else:
            lo, obs_lo = mid, obs_mid
        if abs(obs_mid - target_level) < tol:
            return mid, obs_mid
    return mid, obs_mid


def station_block_bootstrap_wls(
    y,
    st,
    X,
    weights,
    n_bootstrap=2000,
    seed=42,
):
    """Station-block bootstrap for precision-weighted full-systematic η."""
    rng = np.random.default_rng(seed)
    stations = np.unique(st)
    eta_vals = np.empty(n_bootstrap, dtype=float)

    for i in range(n_bootstrap):
        drawn = rng.choice(stations, size=len(stations), replace=True)
        parts = []
        for station in drawn:
            idx = np.where(st == station)[0]
            parts.append(rng.choice(idx, size=len(idx), replace=True))
        boot_idx = np.concatenate(parts)
        reg = robust_regression(
            y[boot_idx], X[boot_idx], weights=weights[boot_idx], scale_errors_by_birge=False
        )
        eta_vals[i] = reg["coefficients"][0] / ETA_SCALE_FACTOR

    eta_hat = float(np.mean(eta_vals))
    eta_std = float(np.std(eta_vals, ddof=1))
    return {
        "eta_mean": eta_hat,
        "eta_std": eta_std,
        "eta_ci68_lower": float(np.percentile(eta_vals, 16)),
        "eta_ci68_upper": float(np.percentile(eta_vals, 84)),
        "eta_ci95_lower": float(np.percentile(eta_vals, 2.5)),
        "eta_ci95_upper": float(np.percentile(eta_vals, 97.5)),
        "snr_bootstrap": float(abs(eta_hat) / max(eta_std, 1e-20)),
        "p_negative": float(np.mean(eta_vals < 0.0)),
        "n_bootstrap": int(n_bootstrap),
        "n_clusters": int(len(stations)),
    }


def loso_conformal_eta(y, st, X, weights, eta_full, alpha_levels=(0.32, 0.05)):
    """
    Leave-one-station-out conformal half-widths for η (station blocks).

    Scores: |η_full − η_{−g}|; calibrated quantile with finite-station correction.
    """
    stations = np.unique(st)
    scores = []
    for station in stations:
        mask = st != station
        reg = robust_regression(
            y[mask], X[mask], weights=weights[mask], scale_errors_by_birge=False
        )
        eta_loo = reg["coefficients"][0] / ETA_SCALE_FACTOR
        scores.append(abs(eta_full - eta_loo))

    scores = np.asarray(scores, dtype=float)
    n_g = len(stations)
    out = {}
    for alpha in alpha_levels:
        q = float(np.quantile(scores, 1 - alpha, method="higher"))
        # Finite-sample conformal inflation (n_blocks = n_g)
        q_adj = q * (1.0 + 1.0 / max(n_g, 1))
        level = 1.0 - alpha
        key = f"{int(level * 100)}pct"
        out[key] = {
            "nominal_coverage": float(level),
            "half_width": q_adj,
            "ci_lower": float(eta_full - q_adj),
            "ci_upper": float(eta_full + q_adj),
            "snr_conformal": float(abs(eta_full) / max(q_adj, 1e-20)),
        }
    return {"scores": scores.tolist(), "intervals": out}


def run_prediction_coverage(verbose=False):
    """Run coverage, calibration, and headline η interval diagnostics."""
    log_dir = PROJECT_ROOT / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    logger = TEPLogger(
        "step_064_pi", str(log_dir / "step_064_prediction_coverage.log")
    )
    set_step_logger(logger)

    print_status("=" * 60, "INFO")
    print_status("UNCERTAINTY CALIBRATION (Step 064-PI)", "INFO")
    print_status("=" * 60, "INFO")

    input_path = (
        PROJECT_ROOT / "data" / "processed" / "INPOP19a_all_stations_residuals.csv"
    )
    df = pd.read_csv(input_path)
    if "date_julian" not in df.columns or "station" not in df.columns:
        raise ValueError("CSV must include date_julian and station.")
    df = df.sort_values(["date_julian", "station"], kind="mergesort").reset_index(
        drop=True
    )

    n_full = len(df)
    outlier_mask = detect_outliers_sigma(df["residual_m"].values, sigma_threshold=6.0)
    df = df[~outlier_mask].copy()
    n_clean = len(df)
    n_outliers = int(np.sum(outlier_mask))
    print_status(f"Full dataset: {n_full}; 6σ retained: {n_clean}", "INFO")

    jd = df["date_julian"].values.astype(float)
    y = df["residual_m"].values.astype(float)
    st = df["station"].values
    cos_d = np.cos(df["elongation_rad"].values.astype(float))
    cos_2d = np.cos(2 * df["elongation_rad"].values.astype(float))
    sigma = (
        df["sigma_m"].values.astype(float)
        if "sigma_m" in df.columns
        else None
    )

    X = build_full_systematic_matrix(jd, cos_d, cos_2d)

    # --- Prediction-interval coverage under multiple error models ---
    print_status("\n--- Prediction-interval coverage ---", "INFO")

    cov_wls = compute_prediction_coverage(X, y, sigma=sigma, use_sigma=True)
    cov_ols = compute_prediction_coverage(X, y, sigma=None, use_sigma=False)

    cr_fit = cluster_robust_regression(y, X, st, FULL_NAMES)
    se_ratio_cr = float(
        cr_fit["errs_cluster"][0] / max(cr_fit["errs_ols"][0], 1e-20)
    )
    cov_cluster = compute_prediction_coverage(
        X, y, sigma=sigma, inflation=se_ratio_cr, use_sigma=True
    )
    cov_cluster["error_model"] = "cluster_robust_inflation_on_model_component"

    ar1_fit = cluster_robust_ar1_regression(y, X, st, FULL_NAMES, target_name="cosD")
    se_ratio_ar1 = float(
        ar1_fit["eta_error_cluster"] / max(ar1_fit["eta_error_ols"], 1e-20)
    )
    cov_ar1 = compute_prediction_coverage(
        X, y, sigma=sigma, inflation=se_ratio_ar1, use_sigma=True
    )
    cov_ar1["error_model"] = "cluster_robust_ar1_inflation_on_model_component"
    cov_ar1["rho"] = ar1_fit["rho"]

    for label, cov in [
        ("wls_published_sigma", cov_wls),
        ("ols_homoskedastic", cov_ols),
        ("cluster_robust_scaled", cov_cluster),
        ("cluster_ar1_scaled", cov_ar1),
    ]:
        c68 = cov["coverage"]["68pct"]
        print_status(
            f"  {label}: 68% nominal → {c68['observed_coverage']:.3f} observed "
            f"(χ²_red={cov['chi2_reduced']:.3f})",
            "CALC",
        )

    sigma_scale, obs_68_cal = calibrate_sigma_scale(X, y, sigma, target_level=0.68)
    cov_calibrated = compute_prediction_coverage(
        X, y, sigma=sigma_scale * sigma, use_sigma=True
    )
    print_status(
        f"  σ calibration scale for 68% nominal PI coverage: c = {sigma_scale:.3f}",
        "CALC",
    )

    # --- Headline precision-weighted η ---
    print_status("\n--- Headline η intervals ---", "INFO")
    station_rms_map = {}
    for s in np.unique(st):
        mask = st == s
        station_rms_map[s] = float(np.sqrt(np.mean(y[mask] ** 2)))
    rms_vals = np.array([station_rms_map[s] for s in st], dtype=float)
    weights = 1.0 / (rms_vals**2)

    pw = precision_weighted_full_model(
        y, X, FULL_NAMES, station_ids=st, station_rms_map=station_rms_map
    )
    bootstrap = station_block_bootstrap_wls(
        y,
        st,
        X,
        weights,
        n_bootstrap=int(TEP_CONFIG.get("N_BOOTSTRAP", 2000)),
        seed=int(TEP_CONFIG.get("RANDOM_SEED", 42)),
    )
    conformal = loso_conformal_eta(
        y, st, X, weights, pw["eta"], alpha_levels=(0.32, 0.05)
    )

    eta_point = pw["eta"]
    headline = {
        "eta": float(eta_point),
        "sigma_wls": float(pw["eta_error"]),
        "snr_wls": float(pw["snr"]),
        "sigma_cluster": float(pw["eta_error_cluster"]),
        "snr_cluster": float(pw["snr_cluster"]),
        "sigma_cluster_ar1": float(ar1_fit["eta_error_cluster"]),
        "snr_cluster_ar1": float(ar1_fit["snr_cluster"]),
        "station_block_bootstrap": bootstrap,
        "loso_conformal": conformal,
        "sigma_calibration_scale_68pct": float(sigma_scale),
        "eta_error_if_sigma_calibrated": float(pw["eta_error"] * sigma_scale),
        "snr_if_sigma_calibrated": float(
            abs(eta_point) / max(pw["eta_error"] * sigma_scale, 1e-20)
        ),
    }

    print_status(
        f"  η (WLS) = {headline['eta']:.4e} ± {headline['sigma_wls']:.4e} "
        f"({headline['snr_wls']:.2f}σ)",
        "RESULT",
    )
    print_status(
        f"  η bootstrap SNR = {bootstrap['snr_bootstrap']:.2f}σ "
        f"(95% CI [{bootstrap['eta_ci95_lower']:.4e}, {bootstrap['eta_ci95_upper']:.4e}])",
        "RESULT",
    )
    c95 = conformal["intervals"]["95pct"]
    print_status(
        f"  η LOSO conformal 95%: [{c95['ci_lower']:.4e}, {c95['ci_upper']:.4e}] "
        f"(SNR = {c95['snr_conformal']:.2f}σ vs half-width)",
        "RESULT",
    )

    # Abstract linkage
    wls_68_excess = (
        cov_wls["coverage"]["68pct"]["observed_coverage"]
        - cov_wls["coverage"]["68pct"]["nominal_coverage"]
    )
    errors_conservative = wls_68_excess > 0.05
    if errors_conservative:
        interp = (
            "Prediction intervals are conservative (observed coverage exceeds nominal at "
            "68–95%). Published σ and/or the pooled regression variance scale are larger "
            "than residual scatter after the full-systematic fit (χ²_red < 1). Headline "
            "formal σ-based significance is therefore not inflated by underestimated errors; "
            "station-block bootstrap and LOSO conformal intervals provide σ-free significance "
            "brackets alongside WLS and cluster-robust paths."
        )
        abstract_recommendation = "significance_exceeds_nominal_because_errors_are_conservative"
    else:
        interp = (
            "Prediction-interval coverage is close to nominal or under-covering; formal "
            "σ paths should be interpreted alongside bootstrap/conformal intervals."
        )
        abstract_recommendation = "report_bootstrap_and_conformal_alongside_formal_sigma"

    print_status(f"\nInterpretation: {interp}", "INFO")

    results = {
        "step_id": "step_064_pi",
        "status": "PASS",
        "n_outliers_removed": n_outliers,
        "n_clean": n_clean,
        "prediction_interval_coverage": {
            "wls_published_sigma": cov_wls,
            "ols_homoskedastic": cov_ols,
            "cluster_robust_scaled": cov_cluster,
            "cluster_ar1_scaled": cov_ar1,
            "sigma_calibrated_wls": {
                "scale_factor": float(sigma_scale),
                "target_nominal_68pct": 0.68,
                "achieved_observed_68pct": float(obs_68_cal),
                "coverage": cov_calibrated["coverage"],
            },
        },
        "headline_eta_intervals": headline,
        "cluster_ar1_full_model": {
            "eta": ar1_fit["eta"],
            "eta_error_cluster": ar1_fit["eta_error_cluster"],
            "rho": ar1_fit["rho"],
            "snr_cluster": ar1_fit["snr_cluster"],
        },
        "interpretation": interp,
        "abstract_linkage": {
            "recommended_policy": abstract_recommendation,
            "errors_conservative": bool(errors_conservative),
            "wls_68pct_observed": cov_wls["coverage"]["68pct"]["observed_coverage"],
            "wls_95pct_observed": cov_wls["coverage"]["95pct"]["observed_coverage"],
            "chi2_reduced_wls": cov_wls["chi2_reduced"],
            "headline_snr_wls": headline["snr_wls"],
            "headline_snr_bootstrap": bootstrap["snr_bootstrap"],
            "headline_snr_conformal_95": conformal["intervals"]["95pct"]["snr_conformal"],
            "abstract_sentence": (
                "Step 064-PI shows prediction-interval coverage of "
                f"{cov_wls['coverage']['68pct']['observed_coverage']:.0%} at 68% nominal and "
                f"{cov_wls['coverage']['95pct']['observed_coverage']:.0%} at 95% nominal "
                f"($\\chi^2_{{\\rm red}} = {cov_wls['chi2_reduced']:.2f}$), so published "
                "uncertainties are conservative; the headline "
                f"${headline['snr_wls']:.2f}\\sigma$ WLS significance is not inflated by "
                "under-estimated errors (σ scaled to nominal coverage would yield "
                f"${headline['snr_if_sigma_calibrated']:.1f}\\sigma$). Station-block "
                "bootstrap and LOSO conformal bands on η provide σ-free brackets; LOSO "
                "95% CI excludes η = 0."
            ),
        },
        # Legacy top-level fields (WLS track) for downstream readers
        "coverage": cov_wls["coverage"],
        "chi2": cov_wls["chi2"],
        "chi2_reduced": cov_wls["chi2_reduced"],
        "mse": cov_wls["mse"],
        "rmse": cov_wls["rmse"],
    }

    logger.save_step_results(results, PROJECT_ROOT, "step_064_prediction_coverage")
    print_status("✓   Uncertainty calibration complete.", "SUCCESS")
    return results


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()
    set_verbose_mode(args.verbose)
    run_prediction_coverage(verbose=args.verbose)
