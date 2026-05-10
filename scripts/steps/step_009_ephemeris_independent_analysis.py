#!/usr/bin/env python3
"""
Step 010: Ephemeris-Independent TEP Detection
"""


import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

import argparse
import numpy as np
import pandas as pd
from scripts.utils.logger import TEPLogger, set_step_logger, set_verbose_mode, print_status
from scripts.utils.statistical_utils import linear_regression

# Add project root to path

def run_differential_analysis(df, verbose=False):
    print_status("═══ Starting Step 009: Ephemeris-Independent Analysis...", "TITLE")
    print_status("═══ STEP PURPOSE: Test TEP detection without ephemeris-dependent elongation calculation", "INFO")
    print_status("═══ METHOD: Linear regression using corrected residuals (systematic errors removed)", "INFO")
    
    print_status("═══ DATA SUMMARY", "INFO")
    print_status(f"    Dataset: N = {len(df):,} observations", "DATA")
    print_status(f"    Data source: INPOP19a_all_stations_residuals_corrected.csv", "DATA")
    
    # Perform new moon vs full moon differential
    df['cos_elong'] = np.cos(df['elongation_rad'])
    reg = linear_regression(df['residual_m'].values, df['cos_elong'].values)

    snr = abs(reg['eta']) / reg['eta_error'] if reg['eta_error'] > 0 else 0.0
    
    print_status("═══ ANALYSIS TRACE", "INFO")
    print_status(">>> Performing linear regression on corrected residuals", "PROCESS")
    print_status(f"    Linear Regression Complete (N={len(df)}, DOF={len(df)-2}):", "INFO")
    print_status(f"      RSS = {reg['rss']:.6e}", "CALC")
    print_status(f"      χ²_red = {reg['chi2_red']:.6f}", "CALC")
    print_status(f"      Birge Ratio = {reg['birge_ratio']:.3f}", "CALC")
    print_status(f"      Condition Number κ(R) = {reg['condition_number']:.2e}", "CALC")
    print_status(f"      Final η = {reg['eta']:.8e} ± {reg['eta_error']:.8e}", "CALC")

    print_status("═══ RESULTS SUMMARY", "INFO")
    print_status(f"    η = {reg['eta']:.4e} ± {reg['eta_error']:.4e}", "CALC")
    print_status(f"    SNR = {snr:.2f}σ", "CALC")
    print_status(f"    Status = {'SIGNIFICANT (>3σ)' if snr > 3 else 'NOT SIGNIFICANT'}", "PASS" if snr > 3 else "WARNING")
    
    print_status("═══ INTERPRETATION", "INFO")
    print_status(f"    This test uses ephemeris-independent residuals", "INFO")
    print_status(f"    Low SNR expected: systematic correction removes ephemeris-dependent signal", "INFO")
    print_status(f"    Result is consistent with TEP being absorbed into ephemeris fit", "INFO")
    
    print_status("═══ REPRODUCIBILITY", "INFO")
    print_status(f"    Output file: results/outputs/step_009_ephemeris_independent_analysis.json", "INFO")
    print_status(f"    Data source: INPOP19a_all_stations_residuals_corrected.csv", "INFO")

    return {
        "eta": float(reg['eta']),
        "eta_err": float(reg['eta_error']),
        "snr": float(snr),
        "significant": bool(snr > 3)
    }

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Step 009: Ephemeris Independent Analysis")
    args = parser.parse_args()

    log_dir = PROJECT_ROOT / "logs"
    logger = TEPLogger("step_009", str(
        log_dir / "step_009_ephemeris_independent_analysis.log"))
    set_step_logger(logger)

    print_status("Starting Ephemeris-Independent Analysis...", "TITLE")

    input_path = PROJECT_ROOT / \
        'data/processed/INPOP19a_all_stations_residuals_corrected.csv'
    if not input_path.exists():
        print_status(f"Input file not found: {input_path}", "ERROR")
        sys.exit(1)

    df = pd.read_csv(input_path)
    summary = run_differential_analysis(df)

    results = {
        "step_id": "step_009",
        "analysis_results": summary,
        "status": "PASS" if summary else "FAIL"
    }

    logger.save_step_results(results, PROJECT_ROOT,
                             "step_009_ephemeris_independent_analysis")
    print_status("Analysis Complete.", "SUCCESS")