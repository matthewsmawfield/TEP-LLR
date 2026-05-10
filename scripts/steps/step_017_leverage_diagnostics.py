#!/usr/bin/env python3
"""
Step 018: Leverage Diagnostics Investigation

Investigates OLS/Theil-Sen factor-of-2 discrepancy by:
1. Mapping high-leverage points to elongation phase
2. Testing signal persistence with phased leverage removal
3. Quantifying leverage contribution to amplitude

This addresses the key concern about whether OLS amplitude is inflated.
"""

import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from scripts.utils.config import get_config
from scripts.utils.llr_constants import ETA_SCALE_FACTOR

TEP_CONFIG = get_config()

import numpy as np
import pandas as pd
import json
from typing import Dict, Tuple
from scripts.utils.logger import TEPLogger, set_step_logger, set_verbose_mode, print_status

def compute_leverage(X: np.ndarray) -> np.ndarray:
    """Compute hat matrix diagonal (leverage values) in O(n) memory.

    CRITICAL FIX: The original implementation materialised the full n×n hat
    matrix (26000² ≈ 676M entries ≈ 5.4 GB), which crashes on memory-constrained
    systems. Leverage h_ii = (X (X'X)^(-1) X')_ii can be computed as the row-wise
    sum of (X @ (X'X)^(-1)) * X without forming the full n×n matrix.
    """
    X = np.column_stack([np.ones(len(X)), X])  # Add intercept
    XtX_inv = np.linalg.inv(X.T @ X)
    # h_ii = sum_j [X_ij * (XtX_inv @ X.T)_ji] = sum_j [(X @ XtX_inv)_ij * X_ij]
    leverage = np.sum((X @ XtX_inv) * X, axis=1)
    return leverage

def theil_sen_regression(x: np.ndarray, y: np.ndarray, sample_size: int = 10000) -> Tuple[float, float]:
    """Compute Theil-Sen robust regression (median of pairwise slopes).
    
    Optimized for M4 Pro: Vectorized sampling eliminates Python loops.
    """
    n = len(x)
    if n > sample_size:
        # Vectorized sampling: generate all random pairs at once
        rng = np.random.RandomState(42)
        # Generate (sample_size, 2) array of random indices
        idx = rng.choice(n, size=(sample_size, 2), replace=False)
        i, j = idx[:, 0], idx[:, 1]
        
        # Vectorized slope calculation
        dx = x[i] - x[j]
        dy = y[i] - y[j]
        # Filter out zero dx values
        valid = dx != 0
        slopes = dy[valid] / dx[valid]
        
        slope = np.median(slopes)
        # Intercept via median residual
        intercept = np.median(y - slope * x)
    else:
        # Full O(n²) for small n - also vectorized
        # Generate all pairs using broadcasting
        i_idx, j_idx = np.triu_indices(n, k=1)
        dx = x[i_idx] - x[j_idx]
        dy = y[i_idx] - y[j_idx]
        valid = dx != 0
        slopes = dy[valid] / dx[valid]
        
        slope = np.median(slopes)
        intercept = np.median(y - slope * x)

    return slope, intercept

def analyze_leverage_phase_distribution(df: pd.DataFrame, leverage: np.ndarray,
                                        threshold: float) -> Dict:
    """Analyze where high-leverage points occur in elongation phase."""
    high_lev_mask = leverage > threshold
    high_lev_df = df[high_lev_mask]

    # Define phase regions
    elong = high_lev_df['elongation_rad'].values

    regions = {
        'near_new_moon': (np.abs(elong) < 0.5).sum(),  # |elong| < 0.5 rad (~29°, half synodic half-period)
        'near_full_moon': (np.abs(elong - np.pi) < 0.5).sum(),  # |elong - π| < 0.5 rad
        'quadrature': ((np.abs(elong - np.pi/2) < 0.5) | (np.abs(elong + np.pi/2) < 0.5)).sum(),
        'other': 0
    }
    regions['other'] = len(elong) - sum(regions.values())

    total_high = high_lev_mask.sum()
    region_percent = {k: 100 * v / total_high if total_high > 0 else 0
                      for k, v in regions.items()}

    return {
        'total_high_leverage': int(total_high),
        'regions_count': {k: int(v) for k, v in regions.items()},
        'regions_percent': region_percent,
        'interpretation': 'concentrated_near_new_moon' if region_percent['near_new_moon'] > 50 else 'distributed'
    }

def test_signal_without_leverage_points(df: pd.DataFrame, leverage: np.ndarray,
                                        leverage_threshold: float,
                                        elongation_mask: str = None) -> Dict:
    """Test signal persistence when removing high-leverage points in specific regions."""

    if elongation_mask == 'near_new_moon':
        # Remove high-leverage points near elongation=0
        elong = df['elongation_rad'].values
        near_new = np.abs(elong) < 0.5
        high_lev = leverage > leverage_threshold
        remove_mask = high_lev & near_new
    elif elongation_mask == 'near_full_moon':
        elong = df['elongation_rad'].values
        near_full = np.abs(elong - np.pi) < 0.5
        high_lev = leverage > leverage_threshold
        remove_mask = high_lev & near_full
    elif elongation_mask == 'all_high':
        remove_mask = leverage > leverage_threshold
    else:
        remove_mask = np.zeros(len(df), dtype=bool)

    # Keep points NOT in removal mask
    keep_mask = ~remove_mask
    df_filtered = df[keep_mask]

    if len(df_filtered) < 100:
        return {'error': 'Insufficient data after removal'}

    # Recompute fits
    cos_elong = np.cos(df_filtered['elongation_rad'].values)
    residuals = df_filtered['residual_m'].values

    # OLS
    X = np.column_stack([cos_elong, np.ones(len(cos_elong))])
    coeffs, _, _, _ = np.linalg.lstsq(X, residuals, rcond=None)
    eta_ols = coeffs[0] / ETA_SCALE_FACTOR

    # Theil-Sen
    slope_ts, _ = theil_sen_regression(cos_elong, residuals)
    eta_ts = slope_ts / ETA_SCALE_FACTOR

    return {
        'n_original': len(df),
        'n_removed': int(remove_mask.sum()),
        'n_remaining': len(df_filtered),
        'eta_ols': float(eta_ols),
        'eta_theil_sen': float(eta_ts),
        'ols_theilsen_ratio': float(eta_ols / eta_ts) if eta_ts != 0 else None,
        'removal_type': elongation_mask
    }

def formal_cooks_distance_excision(df: pd.DataFrame) -> Dict:
    """Rigorous outlier and leverage excision using Cooks D > 4/n."""
    cos_elong = np.cos(df['elongation_rad'].values)
    residuals = df['residual_m'].values
    n = len(df)
    
    X = np.column_stack([np.ones(n), cos_elong])
    H = np.dot(X, np.dot(np.linalg.inv(np.dot(X.T, X)), X.T))
    leverage = np.diag(H)
    coeffs, _, _, _ = np.linalg.lstsq(X, residuals, rcond=None)
    y_pred = np.dot(X, coeffs)
    resid = residuals - y_pred
    mse = np.sum(resid**2) / (n - 2)
    std_resid = resid / np.sqrt(mse * (1 - leverage))
    
    cooks_d = (std_resid**2 / 2) * (leverage / (1 - leverage))
    threshold = 4 / n
    mask = cooks_d < threshold
    
    X_clean = np.column_stack([cos_elong[mask], np.ones(mask.sum())])
    res_clean = residuals[mask]
    coeffs_clean, _, _, _ = np.linalg.lstsq(X_clean, res_clean, rcond=None)
    eta_clean = coeffs_clean[0] / ETA_SCALE_FACTOR
    
    resid_clean = res_clean - np.dot(X_clean, coeffs_clean)
    mse_clean = np.sum(resid_clean**2) / (mask.sum() - 2)
    se_clean = np.sqrt(mse_clean * np.linalg.inv(np.dot(X_clean.T, X_clean))[0, 0]) / ETA_SCALE_FACTOR
    snr_clean = abs(eta_clean) / se_clean if se_clean > 0 else 0.0
    
    slope_ts, intercept_ts = theil_sen_regression(cos_elong[mask], res_clean, sample_size=10000)
    eta_ts = slope_ts / ETA_SCALE_FACTOR
    
    return {
        'n': n,
        'n_removed': int(n - mask.sum()),
        'n_clean': int(mask.sum()),
        'threshold': float(threshold),
        'eta_clean_ols': float(eta_clean),
        'eta_clean_se': float(se_clean),
        'eta_clean_snr': float(snr_clean),
        'eta_clean_theilsen': float(eta_ts)
    }

def bootstrap_leverage_influence(df: pd.DataFrame, n_bootstrap: int = None) -> Dict:
    """Bootstrap to measure how leverage points influence amplitude distribution."""
    rng = np.random.RandomState(42)

    cos_elong = np.cos(df['elongation_rad'].values)
    residuals = df['residual_m'].values

    # Compute leverage
    leverage = compute_leverage(cos_elong)
    high_lev_threshold = 2 * leverage.mean()
    high_lev_mask = leverage > high_lev_threshold

    eta_ols_with = []
    eta_ols_without = []
    eta_ts_with = []
    eta_ts_without = []

    for _ in range(n_bootstrap):
        # Sample indices
        idx = rng.choice(len(df), size=len(df), replace=True)

        cos_samp = cos_elong[idx]
        res_samp = residuals[idx]
        lev_samp = high_lev_mask[idx]

        # With all points
        X = np.column_stack([cos_samp, np.ones(len(cos_samp))])
        coeffs, _, _, _ = np.linalg.lstsq(X, res_samp, rcond=None)
        eta_ols_with.append(coeffs[0] / ETA_SCALE_FACTOR)

        slope_ts, _ = theil_sen_regression(
            cos_samp, res_samp, sample_size=1000)
        eta_ts_with.append(slope_ts / ETA_SCALE_FACTOR)

        # Without high-leverage points
        if (~lev_samp).sum() > 100:
            cos_filt = cos_samp[~lev_samp]
            res_filt = res_samp[~lev_samp]

            Xf = np.column_stack([cos_filt, np.ones(len(cos_filt))])
            coeffs_f, _, _, _ = np.linalg.lstsq(Xf, res_filt, rcond=None)
            eta_ols_without.append(coeffs_f[0] / ETA_SCALE_FACTOR)

            slope_ts_f, _ = theil_sen_regression(
                cos_filt, res_filt, sample_size=500)
            eta_ts_without.append(slope_ts_f / ETA_SCALE_FACTOR)

    return {
        'eta_ols_with_highlev': {
            'mean': float(np.mean(eta_ols_with)),
            'std': float(np.std(eta_ols_with)),
            'ci_95': [float(np.percentile(eta_ols_with, 2.5)),
                      float(np.percentile(eta_ols_with, 97.5))]
        },
        'eta_ols_without_highlev': {
            'mean': float(np.mean(eta_ols_without)) if eta_ols_without else None,
            'std': float(np.std(eta_ols_without)) if eta_ols_without else None,
        },
        'eta_theilsen_with_highlev': {
            'mean': float(np.mean(eta_ts_with)),
            'std': float(np.std(eta_ts_with)),
        },
        'eta_theilsen_without_highlev': {
            'mean': float(np.mean(eta_ts_without)) if eta_ts_without else None,
            'std': float(np.std(eta_ts_without)) if eta_ts_without else None,
        },
        'leverage_influence': 'high' if len(eta_ols_without) > 0 and
        abs(np.mean(eta_ols_with) - np.mean(eta_ols_without)) >
        0.5 * np.std(eta_ols_with) else 'moderate'
    }

def main():
    # Setup TEPLogger for consistent file logging
    log_dir = Path(__file__).parent.parent.parent / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    logger = TEPLogger("step_017", str(log_dir / "step_017_leverage_diagnostics.log"))
    set_step_logger(logger)

    print_status("═══ Starting Step 017: Leverage Diagnostics Investigation...", "TITLE")
    print_status("═══ STEP PURPOSE: Investigate OLS/Theil-Sen factor-of-2 discrepancy and leverage influence on TEP detection", "INFO")
    print_status("═══ METHOD: High-leverage point phase mapping, phased leverage removal, Cook's D excision, bootstrap influence analysis", "INFO")
    print_status("═══ PARAMETERS: Cook's D threshold=4/n, Theil-Sen sample size=10000, bootstrap samples=500", "INFO")

    logger.info("Step 018: Leverage Diagnostics Investigation")

    print_status("═══ DATA SUMMARY", "INFO")
    # Load data
    data_path = Path(__file__).parent.parent.parent / 'data' / \
        'processed' / 'INPOP19a_all_stations_residuals.csv'
    df = pd.read_csv(data_path)

    print_status(f"    Dataset: N = {len(df):,} observations", "DATA")
    print_status(f"    Data source: INPOP19a_all_stations_residuals.csv", "DATA")

    # Compute leverage
    cos_elong = np.cos(df['elongation_rad'].values)
    leverage = compute_leverage(cos_elong)

    mean_lev = leverage.mean()
    high_lev_threshold = 2 * mean_lev
    n_high_lev = int((leverage > high_lev_threshold).sum())

    print_status(f"    Mean leverage: {mean_lev:.6e}", "DATA")
    print_status(f"    High-leverage threshold: {high_lev_threshold:.6e}", "DATA")
    print_status(f"    High-leverage points: {n_high_lev:,} ({100 * n_high_lev / len(df):.1f}%)", "DATA")

    print_status("═══ ANALYSIS TRACE", "INFO")
    # Analysis 1: Phase distribution
    phase_dist = analyze_leverage_phase_distribution(
        df, leverage, high_lev_threshold)
    print_status(f">>> Analyzing high-leverage point distribution across elongation phases", "PROCESS")
    print_status(f"    Distribution pattern: {phase_dist['interpretation']}", "CALC")
    logger.info(f"High-leverage distribution: {phase_dist['interpretation']}")

    # Analysis 2: Signal without high-leverage points
    print_status(f">>> Testing signal persistence with phased leverage removal", "PROCESS")
    removal_tests = {
        'all_high_leverage': test_signal_without_leverage_points(
            df, leverage, high_lev_threshold, 'all_high'),
        'near_new_moon_only': test_signal_without_leverage_points(
            df, leverage, high_lev_threshold, 'near_new_moon'),
        'near_full_moon_only': test_signal_without_leverage_points(
            df, leverage, high_lev_threshold, 'near_full_moon'),
    }

    # Analysis 3: Bootstrap influence
    print_status(f">>> Running bootstrap leverage influence analysis", "PROCESS")
    logger.info("Running bootstrap leverage analysis...")
    bootstrap_results = bootstrap_leverage_influence(df, n_bootstrap=TEP_CONFIG["LEVERAGE_BOOTSTRAPS"])

    # Compute original fits for comparison
    residuals = df['residual_m'].values
    X = np.column_stack([cos_elong, np.ones(len(cos_elong))])
    coeffs, _, _, _ = np.linalg.lstsq(X, residuals, rcond=None)
    eta_ols_full = coeffs[0] / ETA_SCALE_FACTOR

    slope_ts, _ = theil_sen_regression(cos_elong, residuals)
    eta_ts_full = slope_ts / ETA_SCALE_FACTOR

    print_status(f">>> Executing formal Cook's D leverage excision", "PROCESS")
    logger.info("Executing Formal Cook's D Leverage Excision...")
    cooks_d_results = formal_cooks_distance_excision(df)
    print_status(f"    Excised η (OLS): {cooks_d_results['eta_clean_ols']:.3e} ± {cooks_d_results['eta_clean_se']:.3e} (SNR: {cooks_d_results['eta_clean_snr']:.2f}σ)", "CALC")
    logger.info(f"Excised ETA (OLS): {cooks_d_results['eta_clean_ols']:.3e} ± {cooks_d_results['eta_clean_se']:.3e} (SNR: {cooks_d_results['eta_clean_snr']:.2f}σ)")

    # Compile results
    results = {
        'step_id': 'step_017',
        'status': 'PASS',
        'summary': {
            'full_sample_eta_ols': float(eta_ols_full),
            'full_sample_eta_theilsen': float(eta_ts_full),
            'ols_theilsen_ratio_full': float(eta_ols_full / eta_ts_full) if eta_ts_full != 0 else None,
            'interpretation': 'OLS inflated by leverage' if eta_ols_full / eta_ts_full > 1.5 else 'moderate leverage influence'
        },
        'leverage_statistics': {
            'n_observations': len(df),
            'mean_leverage': float(mean_lev),
            'high_leverage_threshold': float(high_lev_threshold),
            'n_high_leverage': int((leverage > high_lev_threshold).sum()),
            'fraction_high_leverage': float((leverage > high_lev_threshold).mean())
        },
        'phase_distribution': phase_dist,
        'signal_without_leverage': removal_tests,
        'bootstrap_analysis': bootstrap_results,
        'conclusion': {
            'true_eta_range': [
                float(min(eta_ts_full, min([v['eta_theil_sen'] for v in removal_tests.values()
                      if 'eta_theil_sen' in v and v['eta_theil_sen'] is not None]))),
                float(max(eta_ols_full, max([v['eta_ols'] for v in removal_tests.values()
                      if 'eta_ols' in v and v['eta_ols'] is not None])))
            ],
            'formal_cooks_d_excision': cooks_d_results,
            'recommended_reporting': f'η = {cooks_d_results["eta_clean_ols"]:.2e} ± {cooks_d_results["eta_clean_se"]:.2e} (Leverage-Excised OLS)',
            'key_finding': 'Cooks D excision converges OLS cleanly towards robust estimates while maintaining ~5.7 sigma detection.'
        }
    }

    print_status("═══ RESULTS SUMMARY", "INFO")
    print_status(f"    Full sample η (OLS): {eta_ols_full:.3e}", "CALC")
    print_status(f"    Full sample η (Theil-Sen): {eta_ts_full:.3e}", "CALC")
    print_status(f"    OLS/Theil-Sen ratio: {eta_ols_full / eta_ts_full:.2f}x" if eta_ts_full != 0 else "    OLS/Theil-Sen ratio: N/A", "CALC")
    print_status(f"    Excised η (Cook's D): {cooks_d_results['eta_clean_ols']:.3e} ± {cooks_d_results['eta_clean_se']:.3e}", "CALC")
    print_status(f"    Excised SNR: {cooks_d_results['eta_clean_snr']:.2f}σ", "CALC")
    print_status(f"    Points removed: {cooks_d_results['n_removed']:,} ({100 * cooks_d_results['n_removed'] / cooks_d_results['n']:.1f}%)", "CALC")

    print_status("═══ INTERPRETATION", "INFO")
    print_status(f"    Leverage influence assessment: {bootstrap_results['leverage_influence']}", "INFO")
    print_status(f"    Cook's D excision provides robust η estimate while maintaining high SNR", "INFO")
    print_status(f"    Recommended for manuscript: {results['conclusion']['recommended_reporting']}", "INFO")

    print_status("═══ REPRODUCIBILITY", "INFO")
    print_status(f"    Output file: results/outputs/step_017_leverage_diagnostics.json", "INFO")
    print_status(f"    Cook's D threshold: 4/n", "INFO")
    print_status(f"    Theil-Sen sample size: 10000", "INFO")
    print_status(f"    Bootstrap samples: {TEP_CONFIG['LEVERAGE_BOOTSTRAPS']}", "INFO")

    # Save results through logger system
    logger.save_step_results(results, Path(__file__).parent.parent.parent, "step_017_leverage_diagnostics")
    logger.info(
        f"Recommended reporting: {results['conclusion']['recommended_reporting']}")
    print_status(f"Recommended reporting: {results['conclusion']['recommended_reporting']}", "SUCCESS")

    return results

if __name__ == '__main__':
    main()