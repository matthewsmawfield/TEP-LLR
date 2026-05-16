#!/usr/bin/env python3
"""
Step 063: INPOP19a Outlier Threshold Sensitivity Sweep
=======================================================

Mirrors step_006b (DE430 outlier robustness) for the primary INPOP19a dataset.
Computes:

1. Threshold sweep (3σ–10σ MAD) to verify the η signal is robust to
   outlier-removal threshold choice, not an artefact of a single cutoff.
2. Phase-bin chi-square test on 6σ outliers to test whether they cluster
   non-uniformly in elongation.
3. Bootstrap resampling (1,000 draws) on 6σ-cleaned data for 95% CI.
4. Permutation null test (10,000 shuffles) for data-driven p-value.

All numbers extracted from live pipeline outputs; no hard-coded strings.
"""

import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

import argparse
import json
import numpy as np
import pandas as pd
from scipy import stats

from scripts.utils.config import get_config
from scripts.utils.numerics import suppress_scipy_array_api_matmul_runtime_warning
from scripts.utils.llr_constants import ETA_SCALE_FACTOR
from scripts.utils.logger import TEPLogger, set_step_logger, print_status
from scripts.utils.statistical_utils import linear_regression, detect_outliers_sigma

TEP_CONFIG = get_config()


def _fit_regression(residuals, cos_elong):
    """Fit OLS regression and return correlation + regression metrics."""
    with suppress_scipy_array_api_matmul_runtime_warning():
        r, p = stats.pearsonr(residuals, cos_elong)
    reg = linear_regression(residuals, cos_elong)
    snr = abs(reg['eta']) / reg['eta_error'] if reg['eta_error'] > 0 else 0.0
    return {
        'r': float(r),
        'p_value': float(p),
        'eta': float(reg['eta']),
        'eta_error': float(reg['eta_error']),
        'snr': float(snr),
        'n': int(len(residuals)),
    }


def threshold_sweep(residuals, cos_elong, elongation, thresholds):
    """Test signal stability across MAD outlier thresholds."""
    results = {}
    for sigma in thresholds:
        outlier_mask = detect_outliers_sigma(residuals, sigma_threshold=sigma)
        n_outliers = int(np.sum(outlier_mask))
        n_kept = int(np.sum(~outlier_mask))

        if n_kept < 100:
            results[f"{sigma:.0f}sigma"] = {
                'sigma_threshold': sigma,
                'n_outliers': n_outliers,
                'n_kept': n_kept,
                'error': 'Insufficient data after outlier removal',
            }
            continue

        fit = _fit_regression(residuals[~outlier_mask], cos_elong[~outlier_mask])
        results[f"{sigma:.0f}sigma"] = {
            'sigma_threshold': sigma,
            'n_outliers': n_outliers,
            'n_kept': n_kept,
            **fit,
        }
    return results


def phase_bin_chi2(residuals, elongation, outlier_mask, n_bins=8):
    """Chi-square test for non-uniform outlier distribution across elongation."""
    elongation_deg = np.degrees(elongation) % 360
    bin_edges = np.linspace(0, 360, n_bins + 1)
    all_counts, _ = np.histogram(elongation_deg, bins=bin_edges)
    outlier_counts, _ = np.histogram(elongation_deg[outlier_mask], bins=bin_edges)

    expected = all_counts * (np.sum(outlier_mask) / len(residuals))

    # Mask bins with zero expected counts to avoid NaN
    valid = expected > 0
    chi2_terms = np.zeros_like(expected, dtype=float)
    with np.errstate(divide='ignore', invalid='ignore'):
        chi2_terms[valid] = (outlier_counts[valid] - expected[valid])**2 / expected[valid]
    chi2 = float(np.sum(chi2_terms))
    n_valid_bins = int(np.sum(valid))
    p_chi2 = 1 - stats.chi2.cdf(chi2, max(1, n_valid_bins - 1)) if n_valid_bins > 1 else float('nan')

    return {
        'chi2': chi2,
        'p_value': float(p_chi2),
        'n_bins': n_bins,
        'n_valid_bins': n_valid_bins,
        'all_counts': all_counts.tolist(),
        'outlier_counts': outlier_counts.tolist(),
        'expected_counts': expected.tolist(),
    }


def bootstrap_ci(residuals, cos_elong, n_bootstrap=1000, seed=42):
    """Bootstrap 95% CI for correlation and η."""
    rng = np.random.default_rng(seed)
    n = len(residuals)
    r_boot = np.empty(n_bootstrap)
    eta_boot = np.empty(n_bootstrap)

    for i in range(n_bootstrap):
        idx = rng.integers(0, n, size=n)
        fit = _fit_regression(residuals[idx], cos_elong[idx])
        r_boot[i] = fit['r']
        eta_boot[i] = fit['eta']

    return {
        'r_ci95': [float(np.percentile(r_boot, 2.5)),
                    float(np.percentile(r_boot, 97.5))],
        'eta_ci95': [float(np.percentile(eta_boot, 2.5)),
                     float(np.percentile(eta_boot, 97.5))],
        'r_mean': float(np.mean(r_boot)),
        'eta_mean': float(np.mean(eta_boot)),
        'n_bootstrap': n_bootstrap,
    }


def permutation_test(residuals, cos_elong, n_permutations=10000, seed=42):
    """Permutation test for correlation significance."""
    rng = np.random.default_rng(seed)
    r_obs, _ = stats.pearsonr(residuals, cos_elong)
    r_perm = np.empty(n_permutations)

    for i in range(n_permutations):
        shuffled = rng.permutation(residuals)
        r_perm[i], _ = stats.pearsonr(shuffled, cos_elong)

    n_extreme = int(np.sum(np.abs(r_perm) >= np.abs(r_obs)))
    p_value = n_extreme / n_permutations

    return {
        'observed_r': float(r_obs),
        'n_extreme': n_extreme,
        'n_permutations': n_permutations,
        'p_value': float(p_value),
        'perm_r_mean': float(np.mean(r_perm)),
        'perm_r_std': float(np.std(r_perm, ddof=1)),
    }


def run_outlier_sensitivity():
    """Run the full outlier sensitivity analysis."""
    print_status("=" * 60, "INFO")
    print_status("INPOP19a OUTLIER THRESHOLD SENSITIVITY (Step 063)", "TITLE")
    print_status("=" * 60, "INFO")

    input_path = PROJECT_ROOT / "data" / "processed" / \
        "INPOP19a_all_stations_residuals.csv"
    df = pd.read_csv(input_path)

    residuals = df['residual_m'].values
    cos_elong = df['cos_elong_rad'].values
    elongation = df['elongation_rad'].values

    print_status(f"Full dataset: {len(df):,} observations", "INFO")

    # 1. Threshold sweep
    thresholds = [3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]
    print_status("Running threshold sweep (3σ–10σ MAD)...", "PROCESS")
    sweep_results = threshold_sweep(residuals, cos_elong, elongation, thresholds)

    print_status("Threshold sweep results:", "CALC")
    for key, val in sweep_results.items():
        if 'error' in val:
            print_status(f"  {key}: {val['error']}", "WARNING")
        else:
            print_status(f"  {key}: η={val['eta']:.6e} ± {val['eta_error']:.6e} "
                         f"(SNR={val['snr']:.2f}σ, N={val['n_kept']:,})", "CALC")

    # 2. Phase-bin chi-square on 6σ outliers
    outlier_mask_6 = detect_outliers_sigma(residuals, sigma_threshold=6.0)
    n_outliers_6 = int(np.sum(outlier_mask_6))
    print_status(f"\n6σ outliers: {n_outliers_6} ({100*n_outliers_6/len(df):.1f}%)", "INFO")
    chi2_result = phase_bin_chi2(residuals, elongation, outlier_mask_6)
    print_status(f"Phase-bin χ² = {chi2_result['chi2']:.2f}, "
                 f"p = {chi2_result['p_value']:.4f}", "CALC")

    # 3. Bootstrap on 6σ-cleaned data
    clean_mask = ~outlier_mask_6
    print_status(f"\nRunning bootstrap ({TEP_CONFIG.get('BOOTSTRAP_ITERATIONS', 1000)} draws)...",
                 "PROCESS")
    boot_result = bootstrap_ci(residuals[clean_mask], cos_elong[clean_mask],
                               n_bootstrap=TEP_CONFIG.get('BOOTSTRAP_ITERATIONS', 1000))
    print_status(f"Bootstrap 95% CI for η: [{boot_result['eta_ci95'][0]:.6e}, "
                 f"{boot_result['eta_ci95'][1]:.6e}]", "CALC")

    # 4. Permutation test on 6σ-cleaned data
    print_status(f"\nRunning permutation test "
                 f"({TEP_CONFIG.get('PERMUTATION_ITERATIONS', 10000)} shuffles)...",
                 "PROCESS")
    perm_result = permutation_test(residuals[clean_mask], cos_elong[clean_mask],
                                   n_permutations=TEP_CONFIG.get('PERMUTATION_ITERATIONS', 10000))
    print_status(f"Permutation p = {perm_result['p_value']:.6e} "
                 f"({perm_result['n_extreme']}/{perm_result['n_permutations']:,})", "CALC")

    # 5. Compute η stability metric
    eta_values = [sweep_results[k]['eta'] for k in sweep_results if 'eta' in sweep_results[k]]
    eta_range = max(eta_values) - min(eta_values) if eta_values else 0.0
    eta_mean = float(np.mean(eta_values)) if eta_values else 0.0
    print_status(f"\nη range across thresholds: {eta_range:.6e}", "CALC")
    print_status(f"η mean across thresholds:  {eta_mean:.6e}", "CALC")

    results = {
        "step_id": "step_063",
        "threshold_sweep": sweep_results,
        "phase_bin_chi2": chi2_result,
        "bootstrap": boot_result,
        "permutation": perm_result,
        "eta_range_across_thresholds": float(eta_range),
        "eta_mean_across_thresholds": eta_mean,
        "n_total": int(len(df)),
        "status": "PASS"
    }

    return results


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Step 063: INPOP19a Outlier Threshold Sensitivity")
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()

    log_dir = PROJECT_ROOT / "logs"
    logger = TEPLogger("step_063", str(
        log_dir / "step_063_outlier_sensitivity.log"))
    set_step_logger(logger)

    summary = run_outlier_sensitivity()
    logger.save_step_results(summary, PROJECT_ROOT,
                             "step_063_outlier_sensitivity")
    print_status("Outlier Sensitivity Analysis Complete.", "SUCCESS")
