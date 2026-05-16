#!/usr/bin/env python3
"""
Step 067: Cluster-Robust + AR(1) Combined Standard Errors
===========================================================

The primary estimand in Step 050 treats AR(1) GLS and cluster-robust SEs
as separate checks. This is methodologically incomplete: the data exhibit
both strong autocorrelation (ρ ≈ 0.425, DW ≈ 1.15) AND station-level
clustering. A defensible primary estimator must correct for BOTH
simultaneously.

This step implements a two-way error model:
  1. Cochrane-Orcutt pre-whitening to remove AR(1) structure
  2. Cluster-robust sandwich estimator on the pre-whitened data

The cluster-robust Newey-West (CR-NW) variance is:
    V_cr_nw = (X'X)^-1  [ Σ_g X_g' u_g u_g' X_g ]  (X'X)^-1

where u_g are the pre-whitened residuals within cluster g. This is
asymptotically valid because pre-whitening makes the within-cluster
residuals approximately uncorrelated across time, restoring the
Liang-Zeger sandwich conditions.

Scientific rationale:
- Cluster-robust alone ignores AR(1), inflating significance
- AR(1) GLS alone ignores clustering, understating significance when
  between-station variance is small (as it is: cluster-robust SEs are
  *smaller* than vanilla in Step 050, confirming low between-station
  heterogeneity)
- The combined estimator is conservative and correct
"""

import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

import json
import numpy as np
import pandas as pd
from scipy import stats

from scripts.utils.logger import TEPLogger, set_step_logger, print_status
from scripts.utils.config import get_config
from scripts.utils.llr_constants import ETA_SCALE_FACTOR
from scripts.steps.step_050_corrected_tep_analysis import ar1_gls_regression
from scripts.utils.statistical_utils import (
    robust_regression, detect_outliers_sigma, cluster_robust_variance
)
from scripts.utils.numerics import suppress_scipy_array_api_matmul_runtime_warning

TEP_CONFIG = get_config()
log_dir = PROJECT_ROOT / "logs"
log_dir.mkdir(parents=True, exist_ok=True)
logger = TEPLogger("step_067", str(log_dir / "step_067_cluster_robust_ar1_combined.log"))
set_step_logger(logger)


def cluster_robust_ar1_regression(y, X, cluster_ids, names, target_name='cosD'):
    """
    Combined cluster-robust + AR(1) estimator.

    Steps:
    1. Fit OLS on raw data, extract residuals
    2. Estimate ρ from lag-1 autocorrelation of residuals
    3. Apply Cochrane-Orcutt quasi-differencing to y and X
    4. Re-fit OLS on transformed data
    5. Compute cluster-robust SEs on transformed residuals

    Returns dict with eta, combined SE, rho, and significance.
    """
    n = len(y)

    # Step 1: OLS on full model
    reg_ols = robust_regression(y, X, weights=None, scale_errors_by_birge=False)
    c = reg_ols['coefficients']
    with np.errstate(over='ignore', divide='ignore', invalid='ignore'):
        resid = y - X @ c

    # Step 2: Estimate ρ (must sort by time first; caller guarantees this)
    rho = np.sum(resid[1:] * resid[:-1]) / np.sum(resid[:-1]**2)
    rho = float(np.clip(rho, -0.99, 0.99))
    dw_stat = float(np.sum(np.diff(resid)**2) / np.sum(resid**2))

    # Step 3: Cochrane-Orcutt transform (lose first observation)
    y_star = y[1:] - rho * y[:-1]
    X_star = X[1:] - rho * X[:-1]
    cluster_ids_star = cluster_ids[1:]

    # Step 4: OLS on transformed data
    reg_star = robust_regression(y_star, X_star, weights=None, scale_errors_by_birge=False)
    coeffs = reg_star['coefficients']
    errs_ols = reg_star['errors']

    # Step 5: Cluster-robust on transformed residuals
    with np.errstate(over='ignore', divide='ignore', invalid='ignore'):
        u_star = y_star - X_star @ coeffs
    cr = cluster_robust_variance(X_star, u_star, cluster_ids_star, small_sample_correction=True)
    errs_cluster = cr['se_cluster']

    # Extract target
    target_idx = names.index(target_name)
    eta = coeffs[target_idx] / ETA_SCALE_FACTOR
    eta_err_ols = errs_ols[target_idx] / ETA_SCALE_FACTOR
    eta_err_cluster = errs_cluster[target_idx] / ETA_SCALE_FACTOR

    # Effective sample size under AR(1)
    n_eff = n * (1 - rho) / (1 + rho) if abs(rho) < 1 else n

    return {
        'coeffs': coeffs,
        'errs_ols': errs_ols,
        'errs_cluster': errs_cluster,
        'eta': float(eta),
        'eta_error_ols': float(eta_err_ols),
        'eta_error_cluster': float(eta_err_cluster),
        'snr_ols': float(abs(eta) / max(eta_err_ols, 1e-20)),
        'snr_cluster': float(abs(eta) / max(eta_err_cluster, 1e-20)),
        'rho': float(rho),
        'durbin_watson': float(dw_stat),
        'n_obs': n,
        'n_eff': float(n_eff),
        'n_clusters': cr['n_clusters'],
        'names': names,
    }


def main():
    print_status("═══ Step 067: Cluster-Robust + AR(1) Combined SEs ═══", "TITLE")

    input_path = PROJECT_ROOT / "data" / "processed" / "INPOP19a_all_stations_residuals.csv"
    df = pd.read_csv(input_path)

    # Time sort for AR(1)
    df = df.sort_values(['date_julian', 'station'], kind='mergesort').reset_index(drop=True)
    res = df['residual_m'].values
    st = df['station'].values
    el = df['elongation_rad'].values
    jd = df['date_julian'].values

    # Outlier cleaning (same as step_050)
    outlier_mask = detect_outliers_sigma(res, 6.0)
    n_outliers = int(np.sum(outlier_mask))
    res_c = res[~outlier_mask]
    st_c = st[~outlier_mask]
    el_c = el[~outlier_mask]
    jd_c = jd[~outlier_mask]

    print_status(f"Dataset: N={len(res_c):,} (removed {n_outliers} outliers)", "DATA")

    # Build physics terms
    year = jd_c / 365.25
    sin_y = np.sin(2 * np.pi * year)
    cos_y = np.cos(2 * np.pi * year)
    month = jd_c / 27.32
    sin_m = np.sin(2 * np.pi * month)
    cos_m = np.cos(2 * np.pi * month)
    cos_c = np.cos(el_c)
    cos2d = np.cos(2 * el_c)

    # Full systematic design matrix (Model 5 from step_050)
    X_full = np.column_stack([cos_c, cos2d, sin_m, cos_m, sin_y, cos_y, np.ones(len(cos_c))])
    names = ['cosD', 'cos2D', 'sin_m', 'cos_m', 'sin_y', 'cos_y', 'const']

    # -----------------------------------------------------------------
    # Primary estimand: cluster-robust + AR(1) combined
    # -----------------------------------------------------------------
    print_status("--- Full systematic model with cluster-robust + AR(1) ---", "INFO")
    result = cluster_robust_ar1_regression(res_c, X_full, st_c, names, target_name='cosD')

    print_status(f"  ρ = {result['rho']:.4f}, DW = {result['durbin_watson']:.3f}", "RESULT")
    print_status(f"  η (AR(1)-GLS OLS SE)  = {result['eta']:.4e} ± {result['eta_error_ols']:.4e} "
                 f"({result['snr_ols']:.2f}σ)", "RESULT")
    print_status(f"  η (cluster-robust)    = {result['eta']:.4e} ± {result['eta_error_cluster']:.4e} "
                 f"({result['snr_cluster']:.2f}σ)", "RESULT")
    print_status(f"  Effective N (AR(1))   = {result['n_eff']:.0f}", "RESULT")

    # -----------------------------------------------------------------
    # For comparison: pure cluster-robust (no AR(1)) on same data
    # -----------------------------------------------------------------
    print_status("--- Pure cluster-robust (no AR(1)) for comparison ---", "INFO")
    reg_ols = robust_regression(res_c, X_full, weights=None, scale_errors_by_birge=False)
    c_ols = reg_ols['coefficients']
    with np.errstate(over='ignore', divide='ignore', invalid='ignore'):
        resid_ols = res_c - X_full @ c_ols
    cr_ols = cluster_robust_variance(X_full, resid_ols, st_c, small_sample_correction=True)
    eta_cr_only = c_ols[0] / ETA_SCALE_FACTOR
    se_cr_only = cr_ols['se_cluster'][0] / ETA_SCALE_FACTOR
    snr_cr_only = abs(eta_cr_only) / max(se_cr_only, 1e-20)
    print_status(f"  η (cluster-robust only) = {eta_cr_only:.4e} ± {se_cr_only:.4e} "
                 f"({snr_cr_only:.2f}σ)", "RESULT")

    # -----------------------------------------------------------------
    # For comparison: pure AR(1) GLS (no clustering) on same data
    # -----------------------------------------------------------------
    print_status("--- Pure AR(1) GLS (no clustering) for comparison ---", "INFO")
    gls_only = ar1_gls_regression(res_c, X_full, names, target_name='cosD', cluster_ids=None)
    print_status(f"  η (AR(1) GLS only) = {gls_only['eta']:.4e} ± {gls_only['eta_error']:.4e} "
                 f"({abs(gls_only['eta'])/max(gls_only['eta_error'], 1e-20):.2f}σ)", "RESULT")

    # -----------------------------------------------------------------
    # cosD-only model for direct comparability with step_003
    # -----------------------------------------------------------------
    print_status("--- cosD-only model with cluster-robust + AR(1) ---", "INFO")
    X_cosd = np.column_stack([cos_c, np.ones(len(cos_c))])
    names_cosd = ['cosD', 'const']
    result_cosd = cluster_robust_ar1_regression(res_c, X_cosd, st_c, names_cosd, target_name='cosD')
    print_status(f"  η = {result_cosd['eta']:.4e} ± {result_cosd['eta_error_cluster']:.4e} "
                 f"({result_cosd['snr_cluster']:.2f}σ)", "RESULT")

    # -----------------------------------------------------------------
    # Save results
    # -----------------------------------------------------------------
    output = {
        "step_id": "step_067",
        "status": "PASS",
        "dataset": {"n_obs": int(len(res_c)), "n_outliers_removed": int(n_outliers)},
        "full_systematic_cluster_ar1": {
            "eta": result['eta'],
            "eta_error_ols": result['eta_error_ols'],
            "eta_error_cluster": result['eta_error_cluster'],
            "snr_ols": result['snr_ols'],
            "snr_cluster": result['snr_cluster'],
            "rho": result['rho'],
            "durbin_watson": result['durbin_watson'],
            "n_eff": result['n_eff'],
            "n_clusters": result['n_clusters'],
        },
        "cosd_only_cluster_ar1": {
            "eta": result_cosd['eta'],
            "eta_error_cluster": result_cosd['eta_error_cluster'],
            "snr_cluster": result_cosd['snr_cluster'],
            "rho": result_cosd['rho'],
        },
        "comparison": {
            "cluster_robust_only": {
                "eta": float(eta_cr_only),
                "eta_error": float(se_cr_only),
                "snr": float(snr_cr_only),
            },
            "ar1_gls_only": {
                "eta": float(gls_only['eta']),
                "eta_error": float(gls_only['eta_error']),
                "snr": float(abs(gls_only['eta']) / max(gls_only['eta_error'], 1e-20)),
            },
        },
        "interpretation": (
            "The combined cluster-robust + AR(1) row is the most conservative "
            "frequentist sensitivity bound on the full-systematic design; the headline "
            "remains precision-weighted WLS (Step 050) because it retains all observations "
            "while down-weighting heterogeneous station noise."
        ),
    }

    logger.save_step_results(output, PROJECT_ROOT, "step_067_cluster_robust_ar1_combined")
    print_status("Step 067 complete.", "SUCCESS")
    return output


if __name__ == "__main__":
    main()
