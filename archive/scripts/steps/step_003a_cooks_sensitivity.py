#!/usr/bin/env python3
"""
Step 003a: Cook's D Threshold Sensitivity Analysis
Test different Cook's D excision thresholds to evaluate impact on OLS:Theil-Sen ratio.
"""

import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

import pandas as pd
import numpy as np
from scipy import stats

from scripts.utils.logger import TEPLogger, set_step_logger, set_verbose_mode, print_status
from scripts.utils.llr_constants import ETA_SCALE_FACTOR
from scripts.utils.statistical_utils import linear_regression, detect_outliers_sigma


def cooks_d_threshold_sensitivity():
    """
    Test different Cook's D thresholds and evaluate impact on:
    - OLS eta and SNR
    - Theil-Sen eta and SNR
    - OLS:Theil-Sen ratio
    - Number of excised points
    """
    print_status("═══ Cook's D Threshold Sensitivity Analysis", "TITLE")
    print_status("═══ STEP PURPOSE: Evaluate impact of different Cook's D excision thresholds", "INFO")
    
    # Load data
    input_path = PROJECT_ROOT / "data" / "processed" / "INPOP19a_all_stations_residuals.csv"
    df = pd.read_csv(input_path)
    
    # Apply 6σ MAD outlier cleaning (standard)
    outlier_mask = detect_outliers_sigma(df['residual_m'].values, sigma_threshold=6.0)
    df_clean = df[~outlier_mask].copy()
    
    # Prepare data
    x = np.cos(df_clean['elongation_rad'].values)
    y = df_clean['residual_m'].values
    n = len(y)
    
    # Design matrix
    X = np.column_stack([x, np.ones_like(x)])
    
    # Full sample OLS
    reg_ols = linear_regression(y, x, weights=None)
    eta_ols_full = reg_ols['eta']
    snr_ols_full = abs(eta_ols_full) / reg_ols['eta_error']
    
    # Full sample Theil-Sen
    from scipy.stats import theilslopes
    slope_ts, intercept_ts, _, _ = theilslopes(y, x)
    eta_ts_full = slope_ts / ETA_SCALE_FACTOR
    # Theil-Sen error estimation (bootstrap)
    n_boot = 1000
    eta_ts_boot = []
    for _ in range(n_boot):
        idx = np.random.choice(n, n, replace=True)
        slope_boot, _, _, _ = theilslopes(y[idx], x[idx])
        eta_ts_boot.append(slope_boot / ETA_SCALE_FACTOR)
    eta_ts_err_full = np.std(eta_ts_boot)
    snr_ts_full = abs(eta_ts_full) / eta_ts_err_full
    
    ratio_full = abs(eta_ols_full / eta_ts_full)
    
    print_status("═══ Full Sample (no Cook's D excision)", "INFO")
    print_status(f"  OLS η = {eta_ols_full:.8e} ± {reg_ols['eta_error']:.8e} ({snr_ols_full:.2f}σ)", "CALC")
    print_status(f"  Theil-Sen η = {eta_ts_full:.8e} ± {eta_ts_err_full:.8e} ({snr_ts_full:.2f}σ)", "CALC")
    print_status(f"  OLS:Theil-Sen ratio = {ratio_full:.2f}", "CALC")
    print_status("", "INFO")
    
    # Compute Cook's D once from the FULL sample (standard for threshold sweep).
    # The influence of each point is measured relative to the full-sample fit.
    # PERFORMANCE FIX: O(n) leverage computation avoids materialising the n×n hat
    # matrix (26000² ≈ 676M entries ≈ 5.4 GB).
    XtX_inv = np.linalg.inv(X.T @ X)
    leverage = np.sum((X @ XtX_inv) * X, axis=1)
    residuals_full = y - (reg_ols['eta'] * 13.0 * x + reg_ols['intercept'])
    mse_full = np.sum(residuals_full**2) / (n - 2)
    D_full = (residuals_full**2 / (2 * mse_full)) * (leverage / (1 - leverage)**2)

    # Test different thresholds
    thresholds = [2/n, 4/n, 8/n, 16/n, 32/n, 64/n, 1e-4, 1e-3, 1e-2]
    results = []

    for thresh in thresholds:
        # Excise points above threshold (using full-sample Cook's D)
        mask = D_full <= thresh
        n_excised = np.sum(~mask)
        
        # Excised OLS
        reg_exc = linear_regression(y[mask], x[mask], weights=None)
        eta_exc = reg_exc['eta']
        snr_exc = abs(eta_exc) / reg_exc['eta_error']
        
        # Excised Theil-Sen
        slope_exc, _, _, _ = theilslopes(y[mask], x[mask])
        eta_ts_exc = slope_exc / ETA_SCALE_FACTOR
        eta_ts_boot_exc = []
        for _ in range(n_boot):
            idx = np.random.choice(np.sum(mask), np.sum(mask), replace=True)
            slope_boot, _, _, _ = theilslopes(y[mask][idx], x[mask][idx])
            eta_ts_boot_exc.append(slope_boot / ETA_SCALE_FACTOR)
        eta_ts_err_exc = np.std(eta_ts_boot_exc)
        snr_ts_exc = abs(eta_ts_exc) / eta_ts_err_exc
        
        ratio_exc = abs(eta_exc / eta_ts_exc) if eta_ts_exc != 0 else np.inf
        
        results.append({
            'threshold': thresh,
            'n_excised': n_excised,
            'pct_excised': 100 * n_excised / n,
            'eta_ols': eta_exc,
            'eta_err_ols': reg_exc['eta_error'],
            'snr_ols': snr_exc,
            'eta_ts': eta_ts_exc,
            'eta_err_ts': eta_ts_err_exc,
            'snr_ts': snr_ts_exc,
            'ratio': ratio_exc
        })
        
        print_status(f"  Threshold {thresh:.2e}: {n_excised} excised ({100*n_excised/n:.1f}%)", "CALC")
        print_status(f"    OLS η = {eta_exc:.8e} ({snr_exc:.2f}σ)", "CALC")
        print_status(f"    Theil-Sen η = {eta_ts_exc:.8e} ({snr_ts_exc:.2f}σ)", "CALC")
        print_status(f"    Ratio = {ratio_exc:.2f}", "CALC")
    
    # Find optimal threshold (ratio closest to 1.0)
    ratios = [r['ratio'] for r in results]
    optimal_idx = np.argmin([abs(r - 1.0) for r in ratios])
    optimal = results[optimal_idx]
    
    print_status("", "INFO")
    print_status("═══ OPTIMAL THRESHOLD (OLS:Theil-Sen ratio ≈ 1.0)", "INFO")
    print_status(f"  Threshold = {optimal['threshold']:.2e}", "CALC")
    print_status(f"  Excised = {optimal['n_excised']} ({optimal['pct_excised']:.1f}%)", "CALC")
    print_status(f"  OLS η = {optimal['eta_ols']:.8e} ({optimal['snr_ols']:.2f}σ)", "CALC")
    print_status(f"  Theil-Sen η = {optimal['eta_ts']:.8e} ({optimal['snr_ts']:.2f}σ)", "CALC")
    print_status(f"  Ratio = {optimal['ratio']:.2f}", "CALC")
    
    return results, optimal


if __name__ == "__main__":
    set_verbose_mode(True)
    results, optimal = cooks_d_threshold_sensitivity()
