#!/usr/bin/env python3
"""
Step 074: Wild Cluster Bootstrap and Station-Era Block Bootstrap
=================================================================

Addresses the small-G vulnerability of cluster-robust standard errors.
With only five stations, analytical cluster-robust SEs (even with
G/(G-1) finite-cluster corrections) are known to have poor coverage
and can be anti-conservative.  This step treats the wild cluster
bootstrap with Webb weights as a co-primary uncertainty companion
to the headline η estimate, not a downstream sensitivity check.

Methods:
1. Wild cluster bootstrap (Webb 6-point weights) on the full-systematic
   OLS model.  Webb weights are specifically recommended for very small
   G (Cameron & Miller 2015) because they match the first four moments
   of the standard normal, giving more reliable inference than
   Rademacher weights when G < 10.

2. Station-era block bootstrap: define blocks by station × era
   (pre-1990 / 1990s / post-2000) and resample blocks with replacement.
   This raises the effective number of clusters from 5 to ~15,
   mitigating the small-station concern through design rather than
   correction.

3. Contrast with analytical cluster-robust SE to quantify the
   degree of anti-conservatism.

Output is written to the evidence ledger as main_results and
referenced directly in the manuscript abstract and results table.
"""

import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

import argparse
import numpy as np
import pandas as pd

from scripts.utils.logger import TEPLogger, set_step_logger, print_status, set_verbose_mode
from scripts.utils.config import get_config
from scripts.utils.llr_constants import ETA_SCALE_FACTOR
from scripts.utils.statistical_utils import (
    detect_outliers_sigma,
    robust_regression,
    wild_cluster_bootstrap,
    block_bootstrap_station_era,
    cluster_robust_variance,
)

TEP_CONFIG = get_config()
log_dir = PROJECT_ROOT / "logs"
log_dir.mkdir(parents=True, exist_ok=True)
logger = TEPLogger("step_074", str(log_dir / "step_074_wild_cluster_bootstrap.log"))
set_step_logger(logger)

FULL_NAMES = ["cosD", "cos2D", "sin_m", "cos_m", "sin_y", "cos_y", "const"]


def build_full_systematic_matrix(jd, el):
    """Build full-systematic design matrix (matches Step 050)."""
    year = jd / 365.25
    month = jd / 27.32
    cos_c = np.cos(el)
    cos2d = np.cos(2 * el)
    return np.column_stack(
        [
            cos_c,
            cos2d,
            np.sin(2 * np.pi * month),
            np.cos(2 * np.pi * month),
            np.sin(2 * np.pi * year),
            np.cos(2 * np.pi * year),
            np.ones(len(cos_c)),
        ]
    )


def assign_era(jd):
    """Assign each Julian date to a coarse era for station-era blocking."""
    eras = np.empty(len(jd), dtype=object)
    # Use calendar year boundaries for physical interpretability
    years = jd / 365.25 + 1858.0  # rough conversion; actual thresholds in Julian days
    # More robust: use fixed Julian day thresholds calibrated to 1990-01-01 and 2000-01-01
    jd_1990 = 2447892.5
    jd_2000 = 2451545.0
    eras[jd < jd_1990] = "pre1990"
    eras[(jd >= jd_1990) & (jd < jd_2000)] = "1990s"
    eras[jd >= jd_2000] = "post2000"
    return eras


def run_wild_cluster_bootstrap(verbose=False):
    print_status("═══ Step 073: Wild Cluster Bootstrap (small-G robustness) ═══", "TITLE")

    input_path = PROJECT_ROOT / "data" / "processed" / "INPOP19a_all_stations_residuals.csv"
    df = pd.read_csv(input_path)
    df = df.sort_values(["date_julian", "station"], kind="mergesort").reset_index(drop=True)

    # 6σ cleaning (canonical)
    outlier_mask = detect_outliers_sigma(df["residual_m"].values, sigma_threshold=6.0)
    n_outliers = int(np.sum(outlier_mask))
    df = df[~outlier_mask].copy()
    n_clean = len(df)
    print_status(f"Dataset: N={n_clean:,} (removed {n_outliers} outliers)", "DATA")

    jd = df["date_julian"].values.astype(float)
    y = df["residual_m"].values.astype(float)
    st = df["station"].values
    el = df["elongation_rad"].values.astype(float)

    X = build_full_systematic_matrix(jd, el)

    # --- Reference OLS fit ---
    reg_ols = robust_regression(y, X, weights=None, scale_errors_by_birge=False)
    beta_ols = reg_ols["coefficients"]
    with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
        resid_ols = y - X @ beta_ols

    eta_ols = float(beta_ols[0] / ETA_SCALE_FACTOR)
    se_ols = float(reg_ols["errors"][0] / ETA_SCALE_FACTOR)
    snr_ols = float(abs(eta_ols) / max(se_ols, 1e-20))

    # --- Analytical cluster-robust (for contrast) ---
    cr = cluster_robust_variance(X, resid_ols, st, small_sample_correction=True)
    se_cr = float(cr["se_cluster"][0] / ETA_SCALE_FACTOR)
    snr_cr = float(abs(eta_ols) / max(se_cr, 1e-20))
    n_clusters = int(cr["n_clusters"])

    print_status(f"  η (OLS)            = {eta_ols:.4e} ± {se_ols:.4e} ({snr_ols:.2f}σ)", "RESULT")
    print_status(f"  η (cluster-robust) = {eta_ols:.4e} ± {se_cr:.4e} ({snr_cr:.2f}σ)  G={n_clusters}", "RESULT")

    # --- 1. Wild cluster bootstrap (Webb weights) ---
    print_status("--- Wild cluster bootstrap (Webb 6-point) ---", "INFO")
    n_boot = int(TEP_CONFIG.get("N_BOOTSTRAP", 10000))
    seed = int(TEP_CONFIG.get("RANDOM_SEED", 42))

    wild_webb = wild_cluster_bootstrap(
        y, X, st,
        n_bootstrap=n_boot,
        weight_scheme="webb",
        seed=seed,
        target_idx=0,
        scale_factor=ETA_SCALE_FACTOR,
    )
    wild_scaled = wild_webb["scaled"]

    print_status(
        f"  η = {wild_scaled['beta']:.4e} ± {wild_scaled['se']:.4e} "
        f"({wild_scaled['snr']:.2f}σ wild bootstrap)",
        "RESULT",
    )
    print_status(
        f"  95% CI: [{wild_scaled['ci_95_lower']:.4e}, {wild_scaled['ci_95_upper']:.4e}] "
        f"(p = {wild_webb['p_value_two_tailed']:.4f})",
        "RESULT",
    )

    # --- 2. Wild cluster bootstrap (Rademacher, for comparison) ---
    wild_rad = wild_cluster_bootstrap(
        y, X, st,
        n_bootstrap=n_boot,
        weight_scheme="rademacher",
        seed=seed + 1,
        target_idx=0,
        scale_factor=ETA_SCALE_FACTOR,
    )
    rad_scaled = wild_rad["scaled"]
    print_status(
        f"  η = {rad_scaled['beta']:.4e} ± {rad_scaled['se']:.4e} "
        f"({rad_scaled['snr']:.2f}σ Rademacher)",
        "RESULT",
    )

    # --- 3. Station-era block bootstrap ---
    print_status("--- Station-era block bootstrap ---", "INFO")
    era_ids = assign_era(jd)
    era_boot = block_bootstrap_station_era(
        y, X, st, era_ids,
        n_bootstrap=5000,
        seed=seed + 2,
        target_idx=0,
        scale_factor=ETA_SCALE_FACTOR,
    )
    era_scaled = era_boot["scaled"]
    print_status(
        f"  η = {era_scaled['beta']:.4e} ± {era_scaled['se']:.4e} "
        f"({era_scaled['snr']:.2f}σ, {era_boot['n_blocks']} blocks)",
        "RESULT",
    )
    print_status(
        f"  95% CI: [{era_scaled['ci_95_lower']:.4e}, {era_scaled['ci_95_upper']:.4e}]",
        "RESULT",
    )

    # --- Comparative diagnostic ---
    se_ratio_wild_to_cr = wild_scaled["se"] / max(se_cr, 1e-20)
    se_ratio_era_to_cr = era_scaled["se"] / max(se_cr, 1e-20)

    print_status("--- Small-G diagnostic ---", "INFO")
    print_status(f"  wild_SE / cluster_SE  = {se_ratio_wild_to_cr:.3f}", "RESULT")
    print_status(f"  era_SE / cluster_SE   = {se_ratio_era_to_cr:.3f}", "RESULT")

    if se_ratio_wild_to_cr > 1.1:
        diag = (
            "Wild bootstrap standard error exceeds analytical cluster-robust; "
            "the headline cluster-robust significance may be anti-conservative."
        )
    elif se_ratio_wild_to_cr < 0.9:
        diag = (
            "Wild bootstrap standard error is smaller than analytical cluster-robust; "
            "analytical SEs appear conservative."
        )
    else:
        diag = (
            "Wild bootstrap and analytical cluster-robust standard errors agree "
            "to within 10%; small-G bias in analytical SEs is modest."
        )
    print_status(f"  {diag}", "INFO")

    # --- Result packaging ---
    results = {
        "step_id": "step_074",
        "status": "PASS",
        "main_results": True,
        "dataset": {
            "n_obs": n_clean,
            "n_outliers_removed": n_outliers,
            "n_clusters_station": n_clusters,
        },
        "reference_ols": {
            "eta": eta_ols,
            "eta_error": se_ols,
            "snr": snr_ols,
        },
        "reference_cluster_robust": {
            "eta": eta_ols,
            "eta_error": se_cr,
            "snr": snr_cr,
            "n_clusters": n_clusters,
        },
        "wild_cluster_bootstrap": {
            "weight_scheme": wild_webb["weight_scheme"],
            "n_bootstrap": wild_webb["n_bootstrap"],
            "n_clusters": wild_webb["n_clusters"],
            "eta": wild_scaled["beta"],
            "eta_error": wild_scaled["se"],
            "snr": wild_scaled["snr"],
            "ci_95_lower": wild_scaled["ci_95_lower"],
            "ci_95_upper": wild_scaled["ci_95_upper"],
            "ci_99_lower": wild_scaled["ci_99_lower"],
            "ci_99_upper": wild_scaled["ci_99_upper"],
            "p_value_two_tailed": wild_webb["p_value_two_tailed"],
            "se_ratio_to_cluster_robust": se_ratio_wild_to_cr,
        },
        "wild_cluster_rademacher": {
            "weight_scheme": wild_rad["weight_scheme"],
            "eta": rad_scaled["beta"],
            "eta_error": rad_scaled["se"],
            "snr": rad_scaled["snr"],
            "ci_95_lower": rad_scaled["ci_95_lower"],
            "ci_95_upper": rad_scaled["ci_95_upper"],
            "p_value_two_tailed": wild_rad["p_value_two_tailed"],
        },
        "station_era_block_bootstrap": {
            "n_bootstrap": era_boot["n_bootstrap"],
            "n_blocks": era_boot["n_blocks"],
            "eta": era_scaled["beta"],
            "eta_error": era_scaled["se"],
            "snr": era_scaled["snr"],
            "ci_95_lower": era_scaled["ci_95_lower"],
            "ci_95_upper": era_scaled["ci_95_upper"],
            "p_value_two_tailed": era_boot["p_value_two_tailed"],
            "se_ratio_to_cluster_robust": se_ratio_era_to_cr,
        },
        "interpretation": {
            "small_g_diagnostic": diag,
            "headline_recommendation": (
                "Report wild cluster bootstrap (Webb) as co-primary uncertainty "
                "companion to analytical cluster-robust SE.  The Webb-weight wild "
                "bootstrap is specifically designed for G ≈ 5 and provides "
                "asymptotically valid inference under heteroskedasticity and "
                "intra-cluster correlation without relying on large-G normality."
            ),
        },
    }

    logger.save_step_results(results, PROJECT_ROOT, "step_074_wild_cluster_bootstrap")
    print_status("Step 074 complete.", "SUCCESS")
    return results


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()
    set_verbose_mode(args.verbose)
    run_wild_cluster_bootstrap(verbose=args.verbose)
