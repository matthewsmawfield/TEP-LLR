#!/usr/bin/env python3
"""
Step 009: Corrected Data Analysis for TEP-LLR
"""


import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

import argparse
import pandas as pd
import numpy as np
from scripts.utils.logger import TEPLogger, set_step_logger, set_verbose_mode, print_status
from scripts.utils.statistical_utils import linear_regression

def analyze_corrected_data(verbose=False):
    processed_dir = PROJECT_ROOT / "data" / "processed"
    input_file = processed_dir / "INPOP19a_all_stations_residuals_corrected.csv"

    if not input_file.exists():
        print_status(f"Corrected data not found: {input_file}", "WARNING")
        return None

    df = pd.read_csv(input_file)
    reg = linear_regression(df['residual_m'].values,
                            np.cos(df['elongation_rad'].values))

    if verbose:
        print_status(f"Analyzed {len(df)} corrected observations", "INFO")
        print_status("Corrected data summary:", "CALC")
        print_status(
            f"  η   = {reg['eta']:.4e} ± {reg['eta_error']:.4e}", "CALC")
        print_status(
            f"  SNR = {abs(reg['eta']) / reg['eta_error']:.2f}σ", "CALC")
        print_status(f"  N   = {len(df)}", "CALC")

    return {
        "n_obs": len(df),
        "eta": float(reg['eta']),
        "eta_err": float(reg['eta_error']),
        "snr": float(abs(reg['eta']) / reg['eta_error']) if reg['eta_error'] > 0 else 0.0
    }

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Step 009: Corrected Data Analysis")
    args = parser.parse_args()

    log_dir = PROJECT_ROOT / "logs"
    logger = TEPLogger("step_009", str(
        log_dir / "step_009_corrected_data_analysis.log"))
    set_step_logger(logger)
    set_verbose_mode(True)

    print_status("Analyzing Corrected LLR Data...", "TITLE")
    summary = analyze_corrected_data(verbose=True)

    results = {
        "step_id": "step_009",
        "analysis_results": summary,
        "status": "PASS" if summary else "FAIL"
    }

    logger.save_step_results(results, PROJECT_ROOT,
                             "step_009_corrected_data_analysis")