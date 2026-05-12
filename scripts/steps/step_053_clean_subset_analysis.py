#!/usr/bin/env python3
"""
Step 053: Clean Subset High-SNR Analysis
=========================================

Option A from critical assessment: restrict to the two highest-quality stations
(Grasse C-SPAD era 2010+ and APO) to test whether the TEP signal strengthens
when phase coverage, precision, and temporal span are optimised.

Rationale:
- Grasse 2010+: C-SPAD detector, ~2 cm RMS, continuous coverage, good phase sampling
- APO: 3.2 cm RMS, excellent phase coverage (mean|cosD|=0.517), independent station
- Excludes: Haleakala (opposite sign, underpowered), Matera (extreme phase truncation,
  N=346), McDonald2 (moderate truncation, larger noise)

This tests whether the signal is a genuine physical modulation detectable in the
best data, or an artifact of pooling heterogeneous stations.
"""

import sys
from pathlib import Path
import json

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from scripts.utils.logger import TEPLogger, set_step_logger, print_status
from scripts.utils.statistical_utils import detect_outliers_sigma, robust_regression, cluster_robust_variance
from scripts.utils.llr_constants import ETA_SCALE_FACTOR
import numpy as np
import pandas as pd
from scipy import stats

logger = TEPLogger("step_053")
set_step_logger(logger)


def fit_model(y, X, names):
    """Fit OLS via QR, return dict with coeffs, errs, rss, aic, resid."""
    result = robust_regression(y, X, scale_errors_by_birge=False)
    c = result['coefficients']
    errs = result['errors']
    n, k = X.shape
    resid = y - X @ c
    mse = result['mse']
    rss = np.sum(resid**2)
    aic = 2 * k + n * np.log(rss / max(n, 1))
    return {
        'coeffs': c, 'errs': errs, 'snrs': np.abs(c) / np.maximum(errs, 1e-20),
        'rss': rss, 'aic': aic, 'mse': mse, 'resid': resid,
        'n': n, 'k': k, 'names': names
    }


def cluster_robust_regression(y, X, cluster_ids, names):
    """Fit OLS with cluster-robust standard errors."""
    result = robust_regression(y, X, scale_errors_by_birge=False)
    c = result['coefficients']
    errs = result['errors']
    n, k = X.shape
    resid = y - X @ c
    mse = result['mse']
    rss = np.sum(resid**2)
    aic = 2 * k + n * np.log(rss / max(n, 1))

    # Cluster-robust SEs
    cr = cluster_robust_variance(X, resid, cluster_ids, small_sample_correction=True)
    se_cluster = cr['se_cluster']

    return {
        'coeffs': c, 'errs': errs, 'se_cluster': se_cluster,
        'snrs': np.abs(c) / np.maximum(se_cluster, 1e-20),
        'rss': rss, 'aic': aic, 'mse': mse, 'resid': resid,
        'n': n, 'k': k, 'names': names, 'n_clusters': cr['n_clusters']
    }


def temporal_cv(df, models, split_year=2013.0):
    """Temporal hold-out: train pre-split, test post-split."""
    train = df[df['date_julian_year'] < split_year]
    test = df[df['date_julian_year'] >= split_year]
    results = {}
    for name, (X_train, y_train, X_test, y_test) in models.items():
        reg = robust_regression(y_train, X_train, scale_errors_by_birge=False)
        c = reg['coefficients']
        y_pred = X_test @ c
        ss_res = np.sum((y_test - y_pred)**2)
        ss_tot = np.sum((y_test - np.mean(y_test))**2)
        r2 = 1 - ss_res / ss_tot if ss_tot > 0 else float('nan')
        results[name] = {
            'n_train': len(y_train), 'n_test': len(y_test),
            'r2_pred': float(r2), 'rmse': float(np.sqrt(ss_res / len(y_test))),
            'rss_test': float(ss_res)
        }
    return results


def main():
    print_status("═══ Step 053: Clean Subset High-SNR Analysis ═══", "TITLE")

    # Load data
    input_path = PROJECT_ROOT / "data" / "processed" / "INPOP19a_all_stations_residuals.csv"
    df = pd.read_csv(input_path)

    # Subset: Grasse 2010+ and APO all years
    grasse_modern = df[(df['station'] == 'Grasse') & (df['date_julian_year'] >= 2010.0)]
    apo_all = df[df['station'] == 'APO']
    df_clean = pd.concat([grasse_modern, apo_all], ignore_index=True)

    print_status(f"Clean subset: Grasse 2010+ (N={len(grasse_modern)}) + APO (N={len(apo_all)}) = {len(df_clean)}", "INFO")

    # Outlier cleaning (6 sigma)
    res = df_clean['residual_m'].values
    outlier_mask = detect_outliers_sigma(res, 6.0)
    df_sub = df_clean[~outlier_mask].copy()
    n_outliers = int(np.sum(outlier_mask))
    print_status(f"Removed {n_outliers} outliers, N={len(df_sub)}", "INFO")

    # Extract variables
    res_c = df_sub['residual_m'].values
    st_c = df_sub['station'].values
    cos_c = np.cos(df_sub['elongation_rad'].values)
    jd_c = df_sub['date_julian'].values
    el_c = df_sub['elongation_rad'].values
    year_c = jd_c / 365.25

    # Physics terms
    sin_y = np.sin(2 * np.pi * year_c)
    cos_y = np.cos(2 * np.pi * year_c)
    month = jd_c / 27.32
    sin_m = np.sin(2 * np.pi * month)
    cos_m = np.cos(2 * np.pi * month)
    cos2d = np.cos(2 * el_c)

    # Cluster IDs for station-based robust SEs
    cluster_map = {s: i for i, s in enumerate(np.unique(st_c))}
    cluster_ids = np.array([cluster_map[s] for s in st_c])

    results = {
        'step': '053_clean_subset_analysis',
        'subset_description': 'Grasse 2010+ and APO',
        'n_total': len(df_sub),
        'n_grasse': int(np.sum(st_c == 'Grasse')),
        'n_apo': int(np.sum(st_c == 'APO')),
        'date_range': [float(df_sub['date_julian_year'].min()), float(df_sub['date_julian_year'].max())],
        'mean_abs_cosD': float(np.mean(np.abs(cos_c))),
        'residual_rms_m': float(np.sqrt(np.mean(res_c**2)))
    }

    # --- MODEL 1: cosD only ---
    print_status("--- Model 1: cosD only ---", "INFO")
    m1 = fit_model(res_c, np.column_stack([cos_c, np.ones(len(cos_c))]), ['cosD', 'const'])
    eta1 = m1['coeffs'][0] / ETA_SCALE_FACTOR
    err1 = m1['errs'][0] / ETA_SCALE_FACTOR
    snr1 = abs(eta1) / err1
    print_status(f"  η = {eta1:.4e} ± {err1:.4e} ({snr1:.2f}σ)", "RESULT")
    results['m1_cosD_only'] = {'eta': float(eta1), 'eta_err': float(err1), 'snr': float(snr1)}

    # --- MODEL 4: cosD + cos(2D) + annual (full systematic) ---
    print_status("--- Model 4: cosD + cos(2D) + annual + monthly (full) ---", "INFO")
    m4 = fit_model(res_c, np.column_stack([cos_c, cos2d, sin_y, cos_y, sin_m, cos_m, np.ones(len(cos_c))]),
                   ['cosD', 'cos2D', 'sin_y', 'cos_y', 'sin_m', 'cos_m', 'const'])
    eta4 = m4['coeffs'][0] / ETA_SCALE_FACTOR
    err4 = m4['errs'][0] / ETA_SCALE_FACTOR
    snr4 = abs(eta4) / err4
    print_status(f"  η = {eta4:.4e} ± {err4:.4e} ({snr4:.2f}σ)", "RESULT")
    results['m4_full'] = {'eta': float(eta4), 'eta_err': float(err4), 'snr': float(snr4)}

    # --- MODEL 4 with cluster-robust SEs ---
    print_status("--- Model 4: cluster-robust SEs ---", "INFO")
    m4_cr = cluster_robust_regression(res_c, np.column_stack([cos_c, cos2d, sin_y, cos_y, sin_m, cos_m, np.ones(len(cos_c))]),
                                       cluster_ids, ['cosD', 'cos2D', 'sin_y', 'cos_y', 'sin_m', 'cos_m', 'const'])
    eta4_cr = m4_cr['coeffs'][0] / ETA_SCALE_FACTOR
    err4_cr = m4_cr['se_cluster'][0] / ETA_SCALE_FACTOR
    snr4_cr = abs(eta4_cr) / err4_cr
    print_status(f"  η = {eta4_cr:.4e} ± {err4_cr:.4e} ({snr4_cr:.2f}σ, cluster-robust, {m4_cr['n_clusters']} clusters)", "RESULT")
    results['m4_full_cluster_robust'] = {
        'eta': float(eta4_cr), 'eta_err': float(err4_cr), 'snr': float(snr4_cr),
        'n_clusters': int(m4_cr['n_clusters'])
    }

    # --- Temporal CV: pre-2013 / post-2013 ---
    print_status("--- Temporal CV (pre-2013 / post-2013) ---", "INFO")
    split_year = 2013.0
    train = df_sub[df_sub['date_julian_year'] < split_year]
    test = df_sub[df_sub['date_julian_year'] >= split_year]
    print_status(f"  Train: {len(train)}, Test: {len(test)}", "INFO")

    def build_X(df):
        yr = df['date_julian'].values / 365.25
        mo = df['date_julian'].values / 27.32
        return np.column_stack([
            np.cos(df['elongation_rad'].values),
            np.cos(2 * df['elongation_rad'].values),
            np.sin(2 * np.pi * yr), np.cos(2 * np.pi * yr),
            np.sin(2 * np.pi * mo), np.cos(2 * np.pi * mo),
            np.ones(len(df))
        ])

    X_tr = build_X(train)
    X_te = build_X(test)
    y_tr = train['residual_m'].values
    y_te = test['residual_m'].values

    reg = robust_regression(y_tr, X_tr, scale_errors_by_birge=False)
    c = reg['coefficients']
    y_pred = X_te @ c
    ss_res = np.sum((y_te - y_pred)**2)
    ss_tot = np.sum((y_te - np.mean(y_te))**2)
    r2_pred = 1 - ss_res / ss_tot if ss_tot > 0 else float('nan')
    eta_train = c[0] / ETA_SCALE_FACTOR
    eta_test = None  # Not defined for test set
    print_status(f"  R²_pred = {r2_pred:.4f}, RMSE = {np.sqrt(ss_res/len(y_te)):.4f}", "RESULT")
    results['temporal_cv'] = {
        'split_year': split_year,
        'n_train': len(y_tr), 'n_test': len(y_te),
        'r2_pred': float(r2_pred), 'rmse': float(np.sqrt(ss_res / len(y_te))),
        'eta_train': float(eta_train)
    }

    # --- Per-station within subset ---
    print_status("--- Per-station cosD-only within clean subset ---", "INFO")
    per_station = {}
    for stn in ['Grasse', 'APO']:
        sub = df_sub[df_sub['station'] == stn]
        y = sub['residual_m'].values
        X = np.column_stack([np.cos(sub['elongation_rad'].values), np.ones(len(sub))])
        reg = robust_regression(y, X, scale_errors_by_birge=False)
        eta = reg['coefficients'][0] / ETA_SCALE_FACTOR
        err = reg['errors'][0] / ETA_SCALE_FACTOR
        snr = abs(eta) / max(err, 1e-20)
        per_station[stn] = {
            'eta': float(eta), 'eta_err': float(err), 'snr': float(snr), 'n': len(sub)
        }
        print_status(f"  {stn}: η = {eta:.4e} ± {err:.4e} ({snr:.2f}σ), N={len(sub)}", "RESULT")
    results['per_station'] = per_station

    # Save
    output_path = PROJECT_ROOT / "results/outputs/step_053_clean_subset_analysis.json"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, 'w') as f:
        json.dump(results, f, indent=2)
    print_status(f"Saved to {output_path}", "INFO")


if __name__ == '__main__':
    main()
