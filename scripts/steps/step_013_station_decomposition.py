#!/usr/bin/env python3
"""
Step 014: Station Decomposition for TEP-LLR
"""


import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

import argparse
import pandas as pd
import numpy as np
from scripts.utils.statistical_utils import linear_regression
from scripts.utils.logger import TEPLogger, set_step_logger, set_verbose_mode, print_status

# Add the project root to the Python path

def run_station_decomposition(df, verbose=False):
    stations = df['station'].unique()
    station_results = {}
    for s in stations:
        station_df = df[df['station'] == s]
        reg = linear_regression(station_df['residual_m'].values, np.cos(
            station_df['elongation_rad'].values))
        station_results[s] = {"eta": float(reg['eta']), "snr": float(
            abs(reg['eta'])/reg['eta_error'])}
        if verbose:
            print_status(f"  Station {s} decomposition summary:", "CALC")
            print_status(
                f"    η   = {reg['eta']:.4e} ± {reg['eta_error']:.4e}", "CALC")
            print_status(f"    SNR = {station_results[s]['snr']:.2f}σ", "CALC")
            print_status(f"    N   = {len(station_df)}", "CALC")

    return station_results

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Step 013: Station Decomposition")
    args = parser.parse_args()

    log_dir = PROJECT_ROOT / "logs"
    logger = TEPLogger("step_013", str(
        log_dir / "step_013_station_decomposition.log"))
    set_step_logger(logger)
    set_verbose_mode(True)

    print_status("Starting Station Decomposition Analysis...", "TITLE")

    input_path = PROJECT_ROOT / 'data/processed/INPOP19a_all_stations_residuals.csv'
    if not input_path.exists():
        print_status(f"Input file not found: {input_path}", "ERROR")
        sys.exit(1)

    df = pd.read_csv(input_path)
    summary = run_station_decomposition(df)

    results = {
        "step_id": "step_013",
        "station_results": summary,
        "status": "PASS"
    }

    logger.save_step_results(results, PROJECT_ROOT,
                             "step_013_station_decomposition")
    print_status("Decomposition Complete.", "SUCCESS")