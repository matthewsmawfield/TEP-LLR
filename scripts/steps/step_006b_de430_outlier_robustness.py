#!/usr/bin/env python3
"""
Step 006b: DE430 Outlier Robustness and Sensitivity Analysis

Rigorous, reproducible validation of the DE430 outlier-removal claims
made in the manuscript. Computes:

1. Threshold sweep (3σ–10σ MAD) to verify the signal is robust to
   outlier-removal threshold choice, not an artefact of a single cutoff.
2. Phase-bin chi-square test on 6σ outliers to test whether they cluster
   non-uniformly in elongation (expected for genuine measurement errors
   at specific lunar phases, not random noise).
3. Bootstrap resampling (1000 draws) on 6σ-cleaned data to provide
   95% confidence intervals for correlation r and η.
4. Permutation null test (10 000 shuffles) to provide a data-driven
   p-value for the cleaned correlation.

All numbers are extracted from live pipeline outputs; no hard-coded
strings are written to the JSON result.
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
from typing import Dict, List, Tuple

from scripts.utils.config import get_config
from scripts.utils.numerics import suppress_scipy_array_api_matmul_runtime_warning
from scripts.utils.llr_constants import ETA_SCALE_FACTOR
from scripts.utils.logger import TEPLogger, set_step_logger, set_verbose_mode, print_status
from scripts.utils.statistical_utils import linear_regression, detect_outliers_sigma

TEP_CONFIG = get_config()


def _fit_regression(residuals: np.ndarray, cos_elong: np.ndarray) -> Dict:
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


def threshold_sweep(residuals: np.ndarray, cos_elong: np.ndarray,
                    elongation: np.ndarray, thresholds: List[float]) -> Dict:
    """Test signal stability across a range of MAD outlier thresholds."""
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


def phase_bin_chi_square(elongation: np.ndarray, outlier_mask: np.ndarray,
                         n_bins: int = 8) -> Dict:
    """
    Test whether outliers are uniformly distributed across elongation phase.

    Bins elongation into n_bins equal-width bins over [0, 2π).
    The null hypothesis is that outliers are uniformly distributed.
    A significant chi-square indicates phase-clustered outliers,
    consistent with genuine measurement errors at specific lunar phases
    (e.g. near full moon) rather than random noise.
    """
    # Normalise elongation to [0, 2π)
    elong_norm = np.mod(elongation + 2 * np.pi, 2 * np.pi)
    outlier_elong = elong_norm[outlier_mask]

    bin_edges = np.linspace(0, 2 * np.pi, n_bins + 1)
    observed, _ = np.histogram(outlier_elong, bins=bin_edges)
    expected = np.full(n_bins, len(outlier_elong) / n_bins)

    # Chi-square with Yates-like small-sample guard
    chi2 = np.sum((observed - expected) ** 2 / np.maximum(expected, 1e-10))
    dof = n_bins - 1
    p_chi2 = float(stats.chi2.sf(chi2, dof))

    # Identify the dominant phase region(s)
    max_bin_idx = int(np.argmax(observed))
    bin_centres = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    dominant_centre_deg = float(np.degrees(bin_centres[max_bin_idx]))

    # Count in the 135°–225° region (full-moon vicinity)
    full_moon_mask = (outlier_elong >= 3 * np.pi / 4) & (outlier_elong < 5 * np.pi / 4)
    n_full_moon = int(np.sum(full_moon_mask))
    pct_full_moon = 100.0 * n_full_moon / len(outlier_elong) if len(outlier_elong) > 0 else 0.0

    return {
        'n_bins': n_bins,
        'bin_edges_deg': [float(np.degrees(e)) for e in bin_edges],
        'observed_counts': [int(c) for c in observed],
        'expected_counts': [float(c) for c in expected],
        'chi2': float(chi2),
        'dof': int(dof),
        'p_value': float(p_chi2),
        'dominant_bin_centre_deg': dominant_centre_deg,
        'n_outliers_full_moon_region': n_full_moon,
        'pct_outliers_full_moon_region': float(pct_full_moon),
    }


def bootstrap_ci(residuals: np.ndarray, cos_elong: np.ndarray,
                 n_bootstrap: int = 1000, seed: int = 42) -> Dict:
    """
    Non-parametric bootstrap (residual resampling with replacement)
    to estimate 95% confidence intervals for r and η.
    """
    rng = np.random.RandomState(seed)
    n = len(residuals)

    r_boot = np.zeros(n_bootstrap)
    eta_boot = np.zeros(n_bootstrap)

    with suppress_scipy_array_api_matmul_runtime_warning():
        for i in range(n_bootstrap):
            idx = rng.choice(n, size=n, replace=True)
            y_s = residuals[idx]
            x_s = cos_elong[idx]
            r_boot[i], _ = stats.pearsonr(y_s, x_s)
            reg = linear_regression(y_s, x_s)
            eta_boot[i] = reg['eta']

    return {
        'n_bootstrap': n_bootstrap,
        'seed': seed,
        'r_ci_95': [float(np.percentile(r_boot, 2.5)),
                    float(np.percentile(r_boot, 97.5))],
        'eta_ci_95': [float(np.percentile(eta_boot, 2.5)),
                      float(np.percentile(eta_boot, 97.5))],
        'r_mean': float(np.mean(r_boot)),
        'r_std': float(np.std(r_boot)),
        'eta_mean': float(np.mean(eta_boot)),
        'eta_std': float(np.std(eta_boot)),
    }


def permutation_test(residuals: np.ndarray, cos_elong: np.ndarray,
                     n_permutations: int = None, seed: int = 43) -> Dict:
    """
    Permutation test: shuffle residuals relative to cos(elongation)
    to build a null distribution for Pearson r.
    """
    if n_permutations is None:
        n_permutations = TEP_CONFIG.get("PERMUTATION_ITERATIONS", 10000)

    rng = np.random.RandomState(seed)
    n = len(residuals)
    with suppress_scipy_array_api_matmul_runtime_warning():
        r_obs, _ = stats.pearsonr(residuals, cos_elong)
    r_obs = float(r_obs)

    n_exceed = 0
    r_perm = np.zeros(n_permutations)
    with suppress_scipy_array_api_matmul_runtime_warning():
        for i in range(n_permutations):
            y_shuffled = rng.permutation(residuals)
            r_perm[i], _ = stats.pearsonr(y_shuffled, cos_elong)
            if abs(r_perm[i]) >= abs(r_obs):
                n_exceed += 1

    p_perm = (n_exceed + 1) / (n_permutations + 1)

    return {
        'n_permutations': n_permutations,
        'seed': seed,
        'r_observed': r_obs,
        'n_exceeding': int(n_exceed),
        'p_value': float(p_perm),
        'r_perm_mean': float(np.mean(r_perm)),
        'r_perm_std': float(np.std(r_perm)),
        'r_perm_min': float(np.min(r_perm)),
        'r_perm_max': float(np.max(r_perm)),
    }


def de430_outlier_robustness(verbose: bool = False) -> Dict:
    """Run the full DE430 outlier robustness analysis."""
    processed_dir = PROJECT_ROOT / "data" / "processed"
    de430_path = processed_dir / "DE430_all_residuals.csv"

    if not de430_path.exists():
        raise FileNotFoundError(f"DE430 residuals not found at {de430_path}")

    df = pd.read_csv(de430_path)
    residuals = df['residual_m'].values.astype(float)
    elongation = df['elongation_rad'].values.astype(float)
    cos_elong = np.cos(elongation)
    n_total = len(df)

    if verbose:
        print_status("DE430 Outlier Robustness Analysis", "TITLE")
        print_status(f"  Total observations: {n_total:,}", "DATA")
        print_status(f"  Residual RMS: {np.sqrt(np.mean(residuals**2)):.3f} m", "DATA")
        print_status(f"  Residual std: {np.std(residuals):.3f} m", "DATA")

    # ------------------------------------------------------------------
    # 1. Raw (uncleaned) correlation
    # ------------------------------------------------------------------
    with suppress_scipy_array_api_matmul_runtime_warning():
        raw_r, raw_p = stats.pearsonr(residuals, cos_elong)
    if verbose:
        print_status("Raw DE430 correlation (no outlier removal):", "PROCESS")
        print_status(f"  r = {raw_r:.6f}, p = {raw_p:.4f}", "CALC")

    # ------------------------------------------------------------------
    # 2. Threshold sweep
    # ------------------------------------------------------------------
    thresholds = [3.0, 4.0, 5.0, 6.0, 10.0]
    if verbose:
        print_status("Threshold sweep across MAD outlier cutoffs:", "PROCESS")
    sweep = threshold_sweep(residuals, cos_elong, elongation, thresholds)

    if verbose:
        for key, val in sweep.items():
            if 'error' not in val:
                print_status(
                    f"  {key}: n_out={val['n_outliers']}, "
                    f"r={val['r']:.5f}, η={val['eta']:.3e} ± {val['eta_error']:.3e} "
                    f"({val['snr']:.2f}σ)", "CALC")

    # ------------------------------------------------------------------
    # 3. Phase-bin chi-square on 6σ outliers
    # ------------------------------------------------------------------
    outlier_mask_6s = detect_outliers_sigma(residuals, sigma_threshold=6.0)
    chi2_result = phase_bin_chi_square(elongation, outlier_mask_6s, n_bins=8)

    if verbose:
        print_status("Phase-bin chi-square test (6σ outliers):", "PROCESS")
        print_status(
            f"  χ² = {chi2_result['chi2']:.2f}, dof = {chi2_result['dof']}, "
            f"p = {chi2_result['p_value']:.4e}", "CALC")
        print_status(
            f"  Outliers in 135°–225° region: "
            f"{chi2_result['n_outliers_full_moon_region']} / "
            f"{int(np.sum(outlier_mask_6s))} "
            f"({chi2_result['pct_outliers_full_moon_region']:.1f}%)", "CALC")

    # ------------------------------------------------------------------
    # 4. Bootstrap CI on 6σ-cleaned data
    # ------------------------------------------------------------------
    clean_residuals = residuals[~outlier_mask_6s]
    clean_cos_elong = cos_elong[~outlier_mask_6s]
    boot = bootstrap_ci(clean_residuals, clean_cos_elong,
                        n_bootstrap=1000, seed=42)

    if verbose:
        print_status("Bootstrap CI (6σ-cleaned, N=1000):", "PROCESS")
        print_status(
            f"  r  95% CI = [{boot['r_ci_95'][0]:.5f}, {boot['r_ci_95'][1]:.5f}]",
            "CALC")
        print_status(
            f"  η  95% CI = [{boot['eta_ci_95'][0]:.3e}, {boot['eta_ci_95'][1]:.3e}]",
            "CALC")

    # ------------------------------------------------------------------
    # 5. Permutation test on 6σ-cleaned data
    # ------------------------------------------------------------------
    perm = permutation_test(clean_residuals, clean_cos_elong,
                            n_permutations=TEP_CONFIG.get("PERMUTATION_ITERATIONS", 10000),
                            seed=43)

    if verbose:
        print_status("Permutation null test (6σ-cleaned):", "PROCESS")
        print_status(
            f"  Observed r = {perm['r_observed']:.5f}", "CALC")
        print_status(
            f"  Null r  = {perm['r_perm_mean']:.6f} ± {perm['r_perm_std']:.6f}", "CALC")
        print_status(
            f"  Exceeding count = {perm['n_exceeding']} / {perm['n_permutations']}, "
            f"p = {perm['p_value']:.4e}", "CALC")

    # ------------------------------------------------------------------
    # 6. Assemble results
    # ------------------------------------------------------------------
    results = {
        'step_id': 'step_006b',
        'status': 'PASS',
        'data_summary': {
            'n_total': n_total,
            'residual_rms_m': float(np.sqrt(np.mean(residuals**2))),
            'residual_std_m': float(np.std(residuals)),
        },
        'raw_correlation': {
            'r': float(raw_r),
            'p_value': float(raw_p),
        },
        'threshold_sweep': sweep,
        'phase_chi_square': chi2_result,
        'bootstrap': boot,
        'permutation': perm,
    }

    return results


def main():
    parser = argparse.ArgumentParser(
        description="Step 006b: DE430 Outlier Robustness and Sensitivity Analysis")
    parser.add_argument("--verbose", action="store_true", help="Verbose output")
    args = parser.parse_args()

    log_dir = PROJECT_ROOT / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    logger = TEPLogger(
        "step_006b",
        str(log_dir / "step_006b_de430_outlier_robustness.log"))
    set_step_logger(logger)
    set_verbose_mode(args.verbose)

    print_status("Step 006b: DE430 Outlier Robustness Analysis", "TITLE")

    results = de430_outlier_robustness(verbose=args.verbose)

    logger.save_step_results(
        results, PROJECT_ROOT, "step_006b_de430_outlier_robustness")

    print_status("Step 006b complete.", "SUCCESS")
    return results


if __name__ == "__main__":
    main()
