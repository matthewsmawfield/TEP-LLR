#!/usr/bin/env python3
"""
Step 001e: Extended DE430 Data Ingestion
Processes extended DE430 residuals (2019-2024) for improved statistical power.

This step is designed to handle extended DE430 baseline data when available.
The current DE430 dataset covers 2014-2018 (4,597 observations, 0.33σ significance).
Extending to 2019-2024 could increase N to ~10,000+ and improve significance to ~0.5-0.7σ.

Author: TEP-LLR Analysis Pipeline
Date: 2026-05-09
"""

import sys
from pathlib import Path
import pandas as pd
import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from scripts.utils.astronomical_utils import compute_elongation
from scripts.utils.parse_de430 import parse_de430_file
from scripts.utils.logger import TEPLogger, set_step_logger, print_status, set_verbose_mode

import argparse


def process_extended_de430(verbose=False):
    """
    Process extended DE430 data from 2019-2024.

    This function looks for extended DE430 residual files and processes them
    to combine with the existing 2014-2018 dataset.

    Expected file naming convention:
    - DE430_2014-2018_residuals.dat (existing)
    - DE430_2019-2024_residuals.dat (to be added)
    - Or combined: DE430_extended_residuals.dat

    Args:
        verbose: Enable verbose output

    Returns:
        Dictionary with processing statistics or None if no extended data found
    """
    raw_dir = PROJECT_ROOT / "data" / "raw"
    processed_dir = PROJECT_ROOT / "data" / "processed"

    # Look for extended DE430 files
    potential_files = [
        "DE430_2019-2024_residuals.dat",
        "DE430_extended_residuals.dat",
        "DE430_2019_residuals.dat",
        "DE430_2020_residuals.dat",
        "DE430_2021_residuals.dat",
        "DE430_2022_residuals.dat",
        "DE430_2023_residuals.dat",
        "DE430_2024_residuals.dat",
    ]

    extended_files = []
    for filename in potential_files:
        filepath = raw_dir / filename
        if filepath.exists():
            extended_files.append(filepath)
            print_status(f"Found extended DE430 file: {filename}", "INFO")

    if not extended_files:
        print_status("No extended DE430 files found", "WARNING")
        print_status("Expected files:", "INFO")
        for filename in potential_files:
            print_status(f"  - {filename}", "INFO")
        print_status("Place extended DE430 data in data/raw/ directory", "INFO")
        return None

    print_status(f"Processing {len(extended_files)} extended DE430 file(s)", "PROCESS")

    # Process existing DE430 data
    existing_file = raw_dir / "DE430_2014-2018_residuals.dat"
    if existing_file.exists():
        print_status("Loading existing DE430 data (2014-2018)...", "INFO")
        df_existing = parse_de430_file(existing_file)
        df_existing = compute_elongation(df_existing)
        print_status(f"Loaded {len(df_existing)} existing observations", "INFO")
    else:
        print_status("No existing DE430 data found", "WARNING")
        df_existing = pd.DataFrame()

    # Process extended files
    df_extended_list = []
    for filepath in extended_files:
        print_status(f"Processing {filepath.name}...", "PROCESS")
        try:
            df = parse_de430_file(filepath)
            if len(df) > 0:
                df = compute_elongation(df)
                df_extended_list.append(df)
                print_status(f"Loaded {len(df)} observations from {filepath.name}", "SUCCESS")
            else:
                print_status(f"No valid data in {filepath.name}", "WARNING")
        except Exception as e:
            print_status(f"Error processing {filepath.name}: {e}", "ERROR")

    if not df_extended_list:
        print_status("No valid extended data processed", "WARNING")
        return None

    # Combine extended data
    df_extended = pd.concat(df_extended_list, ignore_index=True)
    print_status(f"Total extended observations: {len(df_extended)}", "INFO")

    # Combine with existing data
    if len(df_existing) > 0:
        df_combined = pd.concat([df_existing, df_extended], ignore_index=True)
        print_status(f"Combined dataset: {len(df_combined)} total observations", "SUCCESS")
    else:
        df_combined = df_extended
        print_status(f"Using extended data only: {len(df_combined)} observations", "INFO")

    # Clean combined data
    initial_count = len(df_combined)
    df_combined = df_combined.dropna(subset=['residual_m', 'elongation_rad'])
    cleaned_count = len(df_combined)

    if verbose:
        print_status(f"Cleaned dataset: {cleaned_count} observations (removed {initial_count - cleaned_count})", "INFO")
        print_status("Extended DE430 residual stats:", "CALC")
        print_status(f"  Mean: {df_combined['residual_m'].mean():.6e} m", "CALC")
        print_status(f"  Std:  {df_combined['residual_m'].std():.6e} m", "CALC")
        print_status(f"  RMS:  {np.sqrt(np.mean(df_combined['residual_m']**2)):.6e} m", "CALC")
        print_status(f"  Date range: {df_combined['date_julian_year'].min():.2f} to {df_combined['date_julian_year'].max():.2f}", "CALC")

    # Save combined dataset
    output_path = processed_dir / "DE430_extended_residuals.csv"
    df_combined.to_csv(output_path, index=False)
    print_status(f"Saved extended dataset to {output_path}", "SUCCESS")

    # Also update the main DE430 file for consistency
    main_output_path = processed_dir / "DE430_all_residuals.csv"
    df_combined.to_csv(main_output_path, index=False)
    print_status(f"Updated main DE430 file: {main_output_path}", "INFO")

    stats = {
        "n_obs": cleaned_count,
        "n_existing": len(df_existing) if len(df_existing) > 0 else 0,
        "n_extended": len(df_extended),
        "residual_rms_m": float(df_combined['residual_m'].std()),
        "date_range": [float(df_combined['date_julian_year'].min()), float(df_combined['date_julian_year'].max())],
        "baseline_years": float(df_combined['date_julian_year'].max() - df_combined['date_julian_year'].min()),
        "output_file": str(output_path.relative_to(PROJECT_ROOT))
    }

    print_status(f"Extended DE430 Preprocessing Complete", "SUCCESS")
    print_status(f"Baseline extended to {stats['baseline_years']:.1f} years", "INFO")
    print_status(f"Expected significance improvement: ~0.33σ → ~0.5-0.7σ (based on √N scaling)", "INFO")

    return stats


def main():
    parser = argparse.ArgumentParser(
        description="Step 001e: Extended DE430 Preprocessing")
    args = parser.parse_args()

    log_dir = PROJECT_ROOT / "logs"
    logger = TEPLogger("step_001e", str(log_dir / "step_001e_extended_de430.log"))
    set_step_logger(logger)
    set_verbose_mode(True)

    print_status("Extended DE430 Data Ingestion", "TITLE")
    print_status("="*80, "TITLE")

    summary = process_extended_de430(verbose=True)

    results = {
        "step_id": "step_001e",
        "summary": summary,
        "status": "PASS" if summary else "FAIL"
    }

    logger.save_step_results(results, PROJECT_ROOT, "step_001e_extended_de430")

    if summary:
        print_status("Step 001e Complete", "SUCCESS")
    else:
        print_status("Step 001e Skipped - No extended data available", "INFO")


if __name__ == "__main__":
    main()
