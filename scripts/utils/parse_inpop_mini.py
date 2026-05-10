#!/usr/bin/env python3
"""
Parse INPOP MINI format LLR residuals and convert to TEP analysis format.
"""

from scripts.utils.astronomical_utils import compute_elongation
from scripts.utils.llr_constants import SIGMA_UNCERTAINTY_FLOOR_MM
from scripts.utils.logger import TEPLogger, set_step_logger, print_status, get_verbose_mode
import sys
from pathlib import Path
from typing import Dict, Any

import pandas as pd
import numpy as np
from astropy.time import Time

# Add the project root to the Python path
PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))


def parse_mini_file(filepath: Path, station_name: str) -> pd.DataFrame:
    """
    Parse INPOP MINI format file.
    """
    data = []
    line_num = -1  # CRITICAL FIX: initialise before loop to prevent NameError on empty files

    with open(filepath, 'r') as f:
        for line_num, line in enumerate(f):
            parts = line.strip().split()
            if len(parts) < 4:
                continue

            try:
                date_julian_year = float(parts[0])
                residual_m = float(parts[1])
                # Handle newer format (2018+) where reflector is concatenated
                # with the following ID field, e.g. '5<timestamp...>'
                reflector = int(parts[2][0])
                timestamp = parts[3]

                # Convert decimal year to Julian date using astropy
                time = Time(date_julian_year, format='decimalyear')

                # Do not extract uncertainty from parts[6].
                # Empirical investigation across all five INPOP19a station files
                # shows parts[6] values are 10^5–10^6 times larger than the
                # residual RMS and do not represent measurement uncertainties:
                #   APO:     mean ~ 20,700  (ratio_to_rms ~ 656,000x)
                #   Grasse:  mean ~ 342,500 (ratio_to_rms ~ 3,194,000x)
                #   Matera:  mean ~ 2,581   (ratio_to_rms ~ 40,800x)
                #   McDonald2: mean ~ 22,859 (ratio_to_rms ~ 239,000x)
                #   Haleakala: mean ~ 7,616 (ratio_to_rms ~ 55,000x)
                # These values appear to be return-rate counts or other metadata.
                # Using them as sigma_m would produce unphysical uncertainties.
                # Fall back to the data-driven station-RMS estimate below.
                sigma_m = None

                data.append({
                    'date_julian': time.jd,
                    'date_julian_year': date_julian_year,
                    'residual_m': residual_m,
                    'sigma_m': sigma_m,
                    'reflector': reflector,
                    'timestamp': timestamp,
                    'station': station_name
                })
            except (ValueError, IndexError) as e:
                # Use debug for individual line errors
                print_status(
                    f"Line {line_num} parse error in {filepath.name}: {e}", "VERBOSE")
                continue

    df = pd.DataFrame(data)

    # Handle missing uncertainties with data-driven estimation
    # Group by station and use station-specific RMS for observations without sigma_m
    if 'sigma_m' in df.columns and df['sigma_m'].isna().any():
        stations_with_missing = df[df['sigma_m'].isna()]['station'].unique()
        for station in stations_with_missing:
            # Use RMS of residuals from this station as uncertainty estimate
            station_data = df[df['station'] == station]
            station_rms = station_data['residual_m'].std()
            # Apply floor
            station_rms = max(station_rms, SIGMA_UNCERTAINTY_FLOOR_MM / 1000.0)
            df.loc[(df['station'] == station) & (df['sigma_m'].isna()), 'sigma_m'] = station_rms
            
            if get_verbose_mode():
                print_status(f"  Estimated sigma_m = {station_rms:.4f} m for {len(station_data[station_data['sigma_m'].isna()])} observations from station {station}", "CALC")
    
    # If still missing (all observations missing), use global RMS
    if df['sigma_m'].isna().any():
        global_rms = df['residual_m'].std()
        global_rms = max(global_rms, SIGMA_UNCERTAINTY_FLOOR_MM / 1000.0)
        df['sigma_m'] = df['sigma_m'].fillna(global_rms)
        if get_verbose_mode():
            print_status(f"  Using global RMS = {global_rms:.4f} m for remaining observations", "CALC")
    
    n_parsed = len(df)
    n_failed = line_num + 1 - n_parsed
    
    if n_parsed > 0 and (n_failed / (n_parsed + n_failed)) > 0.001:
        raise RuntimeError(
            f"Parsing failure threshold exceeded in {filepath.name}: "
            f"{n_failed} failed of {n_parsed + n_failed} total lines (>0.1%). "
            "Verify file format integrity."
        )

    if get_verbose_mode() and len(df) > 0:
        print_status(f"Parsed file {filepath.name}:", "CALC")
        print_status(f"  Rows: {len(df)}", "CALC")
        print_status(
            f"  Residual Mean: {df['residual_m'].mean():.6e} m", "CALC")
        print_status(
            f"  Residual Std:  {df['residual_m'].std():.6e} m", "CALC")

    return df


def process_station_file(input_path: Path, station_name: str, output_path: Path) -> Dict[str, Any]:
    """
    Process a single station's MINI file and output to TEP format.
    """
    # Use relative path for cleaner logs
    input_rel = input_path.relative_to(PROJECT_ROOT) if input_path.is_relative_to(PROJECT_ROOT) else input_path
    print_status(
        f"Processing {station_name} data from {input_rel}", "PROCESS")

    # Parse MINI format
    df = parse_mini_file(input_path, station_name)

    if len(df) == 0:
        print_status(f"No data found in {input_rel}", "WARNING")
        return {"n_obs": 0}

    print_status(f"Parsed {len(df)} observations", "INFO")

    # Compute elongation
    print_status(
        f"Computing Moon elongation angles for {station_name}...", "INFO")
    df = compute_elongation(df)

    # PERFORMANCE FIX: Pre-compute cos(elongation_rad) to avoid redundant trigonometric computations
    # across all analysis steps
    df['cos_elong_rad'] = np.cos(df['elongation_rad'])

    # Select output columns
    output_columns = ['date_julian', 'date_julian_year',
                      'residual_m', 'sigma_m', 'elongation_rad', 'cos_elong_rad', 'station', 'reflector']
    df_output = df[output_columns].copy()

    if get_verbose_mode():
        print_status(
            f"Preparation for station {station_name} complete:", "CALC")
        print_status(f"  Columns: {list(df_output.columns)}", "CALC")
        print_status(
            f"  Memory:  {df_output.memory_usage(deep=True).sum() / 1024:.2f} KB", "CALC")

    # Save to CSV
    output_path.parent.mkdir(parents=True, exist_ok=True)
    df_output.to_csv(output_path, index=False)

    stats = {
        "n_obs": len(df_output),
        "residual_min_m": float(df_output['residual_m'].min()),
        "residual_max_m": float(df_output['residual_m'].max()),
        "residual_rms_m": float(df_output['residual_m'].std()),
        "date_min": float(df_output['date_julian_year'].min()),
        "date_max": float(df_output['date_julian_year'].max())
    }

    print_status(
        f"Saved {stats['n_obs']} observations to {output_path.name}", "SUCCESS")
    print_status(f"  Residual RMS: {stats['residual_rms_m']:.3f} m", "INFO")

    return stats


def process_all_stations() -> Dict[str, Any]:
    """Process all downloaded INPOP station files."""
    print_status("Processing INPOP LLR residuals for TEP analysis", "TITLE")

    data_dir = PROJECT_ROOT / "data" / "raw"
    output_dir = PROJECT_ROOT / "data" / "processed"
    output_dir.mkdir(parents=True, exist_ok=True)

    # Station files downloaded
    # All stations included in preprocessing to avoid signal encoding
    # Haleakala will be evaluated post-detection for systematic effects
    station_files = {
        'APO': data_dir / "INPOP19a_APO_residuals.txt",
        'Grasse': data_dir / "INPOP19a_Grasse_residuals.txt",
        'Matera': data_dir / "INPOP19a_Matera_residuals.txt",
        'McDonald2': data_dir / "INPOP19a_McDonald2_residuals.txt",
        'Haleakala': data_dir / "INPOP19a_Haleakala_residuals.txt"
    }

    all_data = []
    station_stats = {}

    for station_name, input_path in station_files.items():
        if not input_path.exists():
            print_status(f"File not found: {input_path.name}", "WARNING")
            continue

        output_path = output_dir / f"{station_name}_residuals.csv"
        stats = process_station_file(input_path, station_name, output_path)
        station_stats[station_name] = stats

        # Read back for combining
        if output_path.exists() and output_path.stat().st_size > 0:
            df = pd.read_csv(output_path)
            if len(df) > 0:
                all_data.append(df)

    summary = {
        "individual_stations": station_stats,
        "combined": {}
    }

    # Combine all stations into single file
    if all_data:
        combined_df = pd.concat(all_data, ignore_index=True)
        combined_path = output_dir / "INPOP19a_all_stations_residuals.csv"
        combined_df.to_csv(combined_path, index=False)
        print_status(
            f"Combined {len(combined_df)} observations from all stations", "SUCCESS")

        summary["combined"] = {
            "total_obs": len(combined_df),
            "stations": list(combined_df['station'].unique()),
            "date_range": [float(combined_df['date_julian_year'].min()), float(combined_df['date_julian_year'].max())],
            "residual_rms_m": float(combined_df['residual_m'].std()),
            "output_file": str(combined_path.relative_to(PROJECT_ROOT))
        }

        # Print combined statistics
        print_status("\nCombined statistics:", "TITLE")
        print_status(
            f"Total observations: {summary['combined']['total_obs']}", "INFO")
        print_status(
            f"Residual RMS: {summary['combined']['residual_rms_m']:.3f} m", "INFO")
    else:
        print_status("No station data processed", "ERROR")

    return summary


if __name__ == "__main__":
    log_dir = PROJECT_ROOT / "logs"
    logger = TEPLogger("parse_inpop_mini", str(
        log_dir / "parse_inpop_mini.log"))
    set_step_logger(logger)

    process_all_stations()
