#!/usr/bin/env python3
"""
Step 005: Multi-Ephemeris Comparison for TEP Nordtvedt Signal Detection
"""

import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

import argparse
import pandas as pd
import numpy as np

from scripts.utils.statistical_utils import linear_regression, detect_outliers_sigma
from scripts.utils.logger import TEPLogger, set_step_logger, set_verbose_mode, print_status

def compare_ephemerides(verbose=False):
    processed_dir = PROJECT_ROOT / "data" / "processed"
    ephemeris_files = {
        "INPOP19a": processed_dir / "INPOP19a_all_stations_residuals.csv",
        "DE430": processed_dir / "DE430_all_residuals.csv",
        "DE440": processed_dir / "DE440_all_residuals.csv",
        "INPOP21a": processed_dir / "INPOP21a_all_stations_residuals.csv"  # Added INPOP21a for comprehensive comparison
    }

    comparisons = {}
    for name, path in ephemeris_files.items():
        if path.exists():
            df = pd.read_csv(path)
            residuals = df['residual_m'].values
            cos_elong = np.cos(df['elongation_rad'].values)

            # Robust pre-filter for severe outliers before regression
            outlier_mask = detect_outliers_sigma(
                residuals, sigma_threshold=6.0)
            kept_mask = ~outlier_mask
            n_outliers = int(np.sum(outlier_mask))

            if np.sum(kept_mask) >= 100:
                residuals_fit = residuals[kept_mask]
                cos_elong_fit = cos_elong[kept_mask]
            else:
                # Fail loudly: outlier filtering should not remove >99% of data
                raise ValueError(f"Outlier filtering too aggressive: only {np.sum(kept_mask)}/{len(residuals)} points remain. Check data quality.")

            reg = linear_regression(residuals_fit, cos_elong_fit)
            snr = abs(reg['eta']) / \
                reg['eta_error'] if reg['eta_error'] > 0 else np.nan
            comparisons[name] = {
                "n_obs": len(df),
                "n_used": int(len(residuals_fit)),
                "n_outliers_removed": n_outliers,
                "eta": reg['eta'],
                "eta_error": reg['eta_error'],
                "snr": snr
            }
            if verbose:
                print_status(f"  Ephemeris {name} summary:", "CALC")
                print_status(
                    f"    η   = {reg['eta']:.4e} ± {reg['eta_error']:.4e}", "CALC")
                print_status(
                    f"    SNR = {comparisons[name]['snr']:.2f}σ", "CALC")
                print_status(
                    f"    N   = {len(df)} (used={len(residuals_fit)}, outliers_removed={n_outliers})", "CALC")

    return comparisons

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Step 005: Multi-Ephemeris Comparison")
    args = parser.parse_args()

    log_dir = PROJECT_ROOT / "logs"
    logger = TEPLogger("step_005", str(
        log_dir / "step_005_multi_ephemeris_comparison.log"))
    set_step_logger(logger)
    set_verbose_mode(True)

    print_status("Starting Multi-Ephemeris Comparison...", "TITLE")

    comp_results = compare_ephemerides(verbose=True)

    results = {
        "step_id": "step_005",
        "comparisons": comp_results,
        "status": "PASS" if comp_results else "FAIL"
    }

    logger.save_step_results(results, PROJECT_ROOT,
                             "step_005_multi_ephemeris_comparison")
    print_status("Multi-Ephemeris Comparison Complete.", "SUCCESS")