#!/usr/bin/env python3
"""
Step 001b: DE430 Ephemeris Independent Cross-Check
Processes JPL DE430 residuals for comparison with INPOP19a.
"""

import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from scripts.utils.astronomical_utils import compute_elongation
from scripts.utils.logger import TEPLogger, set_step_logger, print_status, set_verbose_mode
from scripts.utils.parse_de430 import parse_de430_file

import argparse
import numpy as np

# Add project root to path

def process_de430(verbose=False):
    raw_dir = PROJECT_ROOT / "data" / "raw"
    processed_dir = PROJECT_ROOT / "data" / "processed"

    # Use original DE430 data (2014-2018) only
    # Investigation revealed:
    # 1. Geoazur data has opposite sign convention (corrected by negation)
    # 2. Original DE430 shows NO correlation with cos(elongation) (~0)
    # 3. Geoazur shows negative correlation after sign correction
    # 4. Systematic differences between datasets (offset, std)
    # 5. Combining them doesn't strengthen evidence due to dominant original dataset
    original_file = raw_dir / "DE430_2014-2018_residuals.dat"

    if original_file.exists():
        print_status(f"Processing DE430 data from {original_file.name}", "PROCESS")
        df = parse_de430_file(original_file)
    else:
        print_status(f"DE430 residuals not found", "WARNING")
        return None

    if verbose:
        print_status(f"Parsed {len(df)} observations", "INFO")

    print_status("Computing Moon-Sun elongation angles...", "INFO")
    df = compute_elongation(df)

    # Cleaning data
    print_status("Cleaning data...", "INFO")
    initial_count = len(df)
    df = df.dropna(subset=['residual_m', 'elongation_rad'])
    cleaned_count = len(df)

    if verbose:
        print_status(
            f"Cleaned dataset: {cleaned_count} observations (removed {initial_count - cleaned_count})", "INFO")
        print_status("Residual stats for DE430:", "CALC")
        print_status(f"  Mean: {df['residual_m'].mean():.6e} m", "CALC")
        print_status(f"  Std:  {df['residual_m'].std():.6e} m", "CALC")
        print_status(
            f"  RMS:  {np.sqrt(np.mean(df['residual_m']**2)):.6e} m", "CALC")

    output_path = processed_dir / "DE430_all_residuals.csv"
    df.to_csv(output_path, index=False)

    stats = {
        "n_obs": cleaned_count,
        "residual_rms_m": float(df['residual_m'].std()),
        "date_range": [float(df['date_julian_year'].min()), float(df['date_julian_year'].max())],
        "output_file": str(output_path.relative_to(PROJECT_ROOT))
    }

    print_status(
        f"DE430 Preprocessing Complete: {cleaned_count} observations", "SUCCESS")
    return stats

def main():
    parser = argparse.ArgumentParser(
        description="Step 002: DE430 Preprocessing")
    log_dir = PROJECT_ROOT / "logs"

    logger = TEPLogger("step_002", str(
        log_dir / "step_002_de430_preprocessing.log"))
    set_step_logger(logger)
    set_verbose_mode(True)

    summary = process_de430(verbose=True)

    results = {
        "step_id": "step_002",
        "summary": summary,
        "status": "PASS" if summary else "FAIL"
    }

    logger.save_step_results(results, PROJECT_ROOT,
                             "step_002_de430_preprocessing")

if __name__ == "__main__":
    main()