#!/usr/bin/env python3
"""
Step 007: Meta-Analysis of Ephemeris Results
Combines INPOP19a and DE430 results using inverse-variance weighting,
baseline-aware weights, and an ephemeris-difference systematic term.

Author: TEP-LLR Analysis Pipeline
Date: 2026-05-10
"""

import sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from scripts.utils.logger import TEPLogger, set_step_logger, print_status, set_verbose_mode
from scripts.utils.numerics import suppress_scipy_array_api_matmul_runtime_warning
from scripts.utils.statistical_utils import detect_outliers_sigma

import argparse


def meta_analysis_ephemerides(verbose=False):
    """
    Combine INPOP19a and DE430 η estimates with baseline-aware weighting.

    Loads η and formal errors from step_003 and step_006 outputs. Pearson r and
    p-values for each ephemeris are computed from the corresponding processed
    residual CSVs using the same 6σ MAD outlier mask as Step 003 for INPOP19a.
    """
    processed_dir = PROJECT_ROOT / "data" / "processed"
    outputs_dir = PROJECT_ROOT / "results" / "outputs"

    step_003_file = outputs_dir / "step_003_statistical_analysis.json"
    if not step_003_file.exists():
        print_status("step_003_statistical_analysis.json not found", "ERROR")
        return None

    import json
    with open(step_003_file, 'r') as f:
        step_003_results = json.load(f)

    # Load DE430 results from step_006
    step_006_file = outputs_dir / "step_006_multi_ephemeris_comparison.json"
    if not step_006_file.exists():
        print_status("Step 006 results not found", "ERROR")
        return None

    with open(step_006_file, 'r') as f:
        step_006_results = json.load(f)

    inpop_csv = processed_dir / "INPOP19a_all_stations_residuals.csv"
    if not inpop_csv.exists():
        raise FileNotFoundError(f"Missing processed INPOP file: {inpop_csv}")

    df_inpop = pd.read_csv(inpop_csv)
    out_inpop = detect_outliers_sigma(df_inpop["residual_m"].values, sigma_threshold=6.0)
    df_in_clean = df_inpop[~out_inpop]
    y_in = df_in_clean["residual_m"].values
    if "cos_elong_rad" in df_in_clean.columns:
        x_in = df_in_clean["cos_elong_rad"].values
    else:
        x_in = np.cos(df_in_clean["elongation_rad"].values)
    with suppress_scipy_array_api_matmul_runtime_warning():
        inpop_r, inpop_p = stats.pearsonr(y_in, x_in)
    baseline_inpop = float(
        df_in_clean["date_julian_year"].max() - df_in_clean["date_julian_year"].min()
    )

    # Extract INPOP19a statistics from step_003 (η and errors; r from data above)
    stats_inpop = {
        'name': 'INPOP19a',
        'n_obs': step_003_results['regression_metrics']['n_obs'],
        'eta': step_003_results['eta_ols'],
        'eta_error': step_003_results['eta_ols_error'],
        'snr': abs(step_003_results['eta_ols']) / step_003_results['eta_ols_error'],
        'baseline_years': baseline_inpop,
        'r': float(inpop_r),
        'p_value': float(inpop_p)
    }

    # Extract DE430 statistics from step_006
    de430_data = step_006_results['comparisons']['DE430']

    # Compute correlation for DE430 from actual data
    de430_file = processed_dir / "DE430_all_residuals.csv"
    if not de430_file.exists():
        raise FileNotFoundError(
            f"DE430 processed residuals missing: {de430_file}. "
            "Run DE430 preprocessing (step 002) before meta-analysis."
        )

    df_de430 = pd.read_csv(de430_file)
    de430_residuals = df_de430['residual_m'].values
    de430_cos_elong = np.cos(df_de430['elongation_rad'].values)
    with suppress_scipy_array_api_matmul_runtime_warning():
        de430_r, de430_p = stats.pearsonr(de430_residuals, de430_cos_elong)
    baseline_de430 = float(
        df_de430["date_julian_year"].max() - df_de430["date_julian_year"].min()
    )

    stats_de430 = {
        'name': 'DE430',
        'n_obs': de430_data['n_used'],
        'eta': de430_data['eta'],
        'eta_error': de430_data['eta_error'],
        'snr': de430_data['snr'],
        'baseline_years': baseline_de430,
        'r': float(de430_r) if de430_r is not None else None,
        'p_value': float(de430_p) if de430_p is not None else None
    }

    if verbose:
        print_status(
            f"INPOP19a (from step_003): η = {stats_inpop['eta']:.3e} ± {stats_inpop['eta_error']:.3e} "
            f"({stats_inpop['snr']:.2f}σ), r = {stats_inpop['r']:.5f}",
            "CALC",
        )
        print_status(
            f"DE430 (from step_006): η = {stats_de430['eta']:.3e} ± {stats_de430['eta_error']:.3e} "
            f"({stats_de430['snr']:.2f}σ), r = {stats_de430['r']:.5f}",
            "CALC",
        )

    # Meta-analysis: Bayesian combination
    # Weight by inverse variance, but also account for baseline and sign consistency
    
    # Check sign consistency
    sign_consistent = (stats_inpop['eta'] < 0 and stats_de430['eta'] < 0) or \
                     (stats_inpop['eta'] > 0 and stats_de430['eta'] > 0)
    
    if verbose:
        print_status(f"Sign consistency: {sign_consistent}", "CALC")

    # Variance weights
    var_inpop = stats_inpop['eta_error']**2
    var_de430 = stats_de430['eta_error']**2
    
    # Baseline weighting (longer baseline gets more weight)
    baseline_inpop = stats_inpop['baseline_years']
    baseline_de430 = stats_de430['baseline_years']
    total_baseline = baseline_inpop + baseline_de430
    
    weight_inpop = (1/var_inpop) * (baseline_inpop / total_baseline)
    weight_de430 = (1/var_de430) * (baseline_de430 / total_baseline)
    
    # Normalize weights
    total_weight = weight_inpop + weight_de430
    weight_inpop_norm = weight_inpop / total_weight
    weight_de430_norm = weight_de430 / total_weight
    
    if verbose:
        print_status(f"INPOP19a weight: {weight_inpop_norm:.3f} (baseline: {baseline_inpop:.1f} years)", "CALC")
        print_status(f"DE430 weight:   {weight_de430_norm:.3f} (baseline: {baseline_de430:.1f} years)", "CALC")

    # Combined estimate
    eta_combined = weight_inpop_norm * stats_inpop['eta'] + weight_de430_norm * stats_de430['eta']
    
    # Combined error (propagation)
    var_combined = (weight_inpop_norm**2 * var_inpop) + (weight_de430_norm**2 * var_de430)
    eta_error_combined = np.sqrt(var_combined)
    
    # Combined SNR
    if eta_error_combined > 0:
        snr_combined = abs(eta_combined) / eta_error_combined
    else:
        snr_combined = 0
    
    # Systematic uncertainty from ephemeris differences
    # The difference between INPOP19a and DE430 eta estimates represents ephemeris-level systematic uncertainty
    # Standard metrological practice for combining measurements with unknown systematic biases:
    # Use the difference between measurements as the systematic uncertainty
    # (BIPM GUM: Guide to the Expression of Uncertainty in Measurement)
    systematic_diff = abs(stats_inpop['eta'] - stats_de430['eta'])
    systematic_uncertainty = systematic_diff  # Ephemeris difference as systematic uncertainty
    
    # Total uncertainty (statistical + systematic)
    eta_error_total = np.sqrt(eta_error_combined**2 + systematic_uncertainty**2)
    snr_total = abs(eta_combined) / eta_error_total if eta_error_total > 0 else 0
    
    if verbose:
        print_status("Meta-Analysis Results:", "CALC")
        print_status(f"  Combined η (statistical only): {eta_combined:.3e} ± {eta_error_combined:.3e} ({snr_combined:.2f}σ)", "CALC")
        print_status(f"  Systematic uncertainty: {systematic_uncertainty:.3e}", "CALC")
        print_status(f"  Combined η (total uncertainty): {eta_combined:.3e} ± {eta_error_total:.3e} ({snr_total:.2f}σ)", "CALC")

    # Sign consistency check (qualitative only — not used for quantitative weighting)
    # Both ephemerides showing the same sign is evidence against random fluctuation,
    # but applying a multiplicative factor to both weights and renormalising is a no-op.
    # The proper way to incorporate sign consistency would be Bayesian model averaging
    # with a sign-consistency prior, which is beyond the scope of this meta-analysis.
    if sign_consistent and verbose:
        print_status("  Sign consistency: both ephemerides show same sign (qualitative check)", "CALC")

    # Results summary
    results = {
        "step_id": "step_007",
        "individual_results": {
            "INPOP19a": {
                "eta": float(stats_inpop['eta']),
                "eta_error": float(stats_inpop['eta_error']),
                "snr": float(stats_inpop['snr']),
                "n_obs": int(stats_inpop['n_obs']),
                "baseline_years": float(stats_inpop['baseline_years']),
                "r": float(stats_inpop['r']),
                "p_value": float(stats_inpop['p_value'])
            },
            "DE430": {
                "eta": float(stats_de430['eta']),
                "eta_error": float(stats_de430['eta_error']),
                "snr": float(stats_de430['snr']),
                "n_obs": int(stats_de430['n_obs']),
                "baseline_years": float(stats_de430['baseline_years']),
                "r": float(stats_de430['r']),
                "p_value": float(stats_de430['p_value'])
            }
        },
        "meta_analysis": {
            "sign_consistent": bool(sign_consistent),
            "sign_consistent_note": "Qualitative check only; not used for quantitative weighting (identical weights after renormalisation)",
            "eta_combined_statistical": float(eta_combined),
            "eta_error_statistical": float(eta_error_combined),
            "snr_statistical": float(snr_combined),
            "systematic_uncertainty": float(systematic_uncertainty),
            "eta_combined_total": float(eta_combined),
            "eta_error_total": float(eta_error_total),
            "snr_total": float(snr_total),
            "weights": {
                "INPOP19a": float(weight_inpop_norm),
                "DE430": float(weight_de430_norm)
            },
            "interpretation_notes": {
                "headline_estimands": (
                    "Primary inference is anchored in Step 040 (unified table) and "
                    "Step 050 (corrected full-systematic / cluster-robust models), "
                    "not in snr_total from this ephemeris-combination bookkeeping step."
                ),
                "snr_total": (
                    "snr_total inflates the denominator by a conservative ephemeris-difference "
                    "systematic term; it is not an independent discovery statistic and can be "
                    "much smaller than snr_statistical even when both ephemerides are sign-consistent."
                ),
            },
        },
        "interpretation": {
            "primary_result": f"η = {eta_combined:.3e} ± {eta_error_total:.3e} ({snr_total:.2f}σ)",
            "conclusion": "Sign-consistent ephemerides strengthen evidence" if sign_consistent else "Ephemerides show inconsistent signs"
        },
        "status": "PASS"
    }

    print_status("Meta-Analysis Complete", "SUCCESS")
    print_status(f"Primary result: {results['interpretation']['primary_result']}", "INFO")

    return results


def main():
    parser = argparse.ArgumentParser(
        description="Step 007: Meta-Analysis of Ephemeris Results")
    args = parser.parse_args()

    log_dir = PROJECT_ROOT / "logs"
    logger = TEPLogger("step_007", str(log_dir / "step_007_meta_analysis.log"))
    set_step_logger(logger)
    set_verbose_mode(True)

    print_status("Meta-Analysis of INPOP19a and DE430 Results", "TITLE")
    print_status("="*80, "TITLE")

    results = meta_analysis_ephemerides(verbose=True)

    if results:
        outputs_dir = PROJECT_ROOT / "results" / "outputs"
        outputs_dir.mkdir(parents=True, exist_ok=True)
        
        logger.save_step_results(results, PROJECT_ROOT, "step_007_meta_analysis")
        print_status("Step 007 Complete", "SUCCESS")
    else:
        print_status("Step 007 Failed", "ERROR")


if __name__ == "__main__":
    main()
