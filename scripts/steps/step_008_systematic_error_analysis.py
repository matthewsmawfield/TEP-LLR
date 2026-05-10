#!/usr/bin/env python3
"""
Step 006: Systematic Error Analysis for TEP Nordtvedt Signal Detection
Enhanced with comprehensive systematic error budget table
"""

import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

import argparse
import pandas as pd
import numpy as np

from scripts.utils.logger import TEPLogger, set_step_logger, set_verbose_mode, print_status

def analyze_systematics(df, verbose=False, logger=None):
    stations = df['station'].unique()
    systematics = {}

    print_status(f"Analyzing {len(stations)} stations for systematic errors...", "INFO")
    print_status(f"Total observations: {len(df):,}", "INFO")
    print_status(f"Residual threshold for 'clean' status: 0.15 m", "INFO")
    print_status("", "INFO")

    clean_count = 0
    noisy_count = 0

    for s in stations:
        station_data = df[df['station'] == s]
        rms = station_data['residual_m'].std()
        is_clean = bool(rms < 0.15)
        systematics[s] = {
            "n_obs": len(station_data),
            "residual_rms": float(rms),
            "clean_profile": is_clean
        }

        # Always log station analysis (not just verbose)
        status_icon = "✓" if is_clean else "⚠"
        status_label = "CLEAN" if is_clean else "NOISY"
        print_status(f"{status_icon} Station {s}:", "PASS" if is_clean else "WARNING")
        print_status(f"    N obs     = {len(station_data):,}", "INFO")
        print_status(f"    RMS       = {rms:.4f} m", "INFO")
        print_status(f"    Status    = {status_label}", "PASS" if is_clean else "WARNING")

        if is_clean:
            clean_count += 1
        else:
            noisy_count += 1

    print_status("", "INFO")
    print_status(f"Systematic Analysis Summary:", "TITLE")
    print_status(f"  Clean stations:  {clean_count}/{len(stations)}", "PASS" if clean_count == len(stations) else "WARNING")
    print_status(f"  Noisy stations:  {noisy_count}/{len(stations)}", "INFO" if noisy_count == 0 else "WARNING")
    print_status(f"  Overall status:  {'PASS' if all(v.get('clean_profile', False) for v in systematics.values()) else 'WARNING'}", "PASS" if all(v.get('clean_profile', False) for v in systematics.values()) else "WARNING")

    return systematics

def generate_systematic_error_budget(df, systematics, verbose=False, logger=None):
    """
    Generate comprehensive systematic error budget table
    Quantifies contributions from various error sources
    """
    print_status("", "INFO")
    print_status("Generating Systematic Error Budget...", "TITLE")
    print_status("", "INFO")

    # Calculate global statistics
    global_rms = df['residual_m'].std()
    global_mean = df['residual_m'].mean()
    n_total = len(df)

    # Error budget components (quantitative estimates based on LLR literature)
    error_budget = {
        "ephemeris_modeling": {
            "source": "Ephemeris modeling (INPOP19a)",
            "magnitude_cm": 1.5,  # ~1.5 cm from INPOP19a orbit determination
            "description": "Uncertainty in lunar and planetary ephemeris fitting",
            "reference": "Fienga et al. 2019"
        },
        "atmospheric_delay": {
            "source": "Atmospheric delay modeling",
            "magnitude_cm": 0.5,  # ~0.5 cm typical for modern corrections
            "description": "Tropospheric delay correction uncertainties",
            "reference": "Degnan 1993"
        },
        "instrumental": {
            "source": "Instrumental systematic",
            "magnitude_cm": 0.3,  # ~0.3 cm for modern SPAD detectors
            "description": "Detector calibration and timing electronics",
            "reference": "Murphy et al. 2014"
        },
        "center_of_mass": {
            "source": "Retroreflector center-of-mass",
            "magnitude_cm": 0.2,  # ~0.2 cm uncertainty in reflector position
            "description": "Uncertainty in retroreflector array center-of-mass position",
            "reference": "Murphy et al. 2010"
        },
        "tidal_modeling": {
            "source": "Tidal modeling",
            "magnitude_cm": 0.4,  # ~0.4 cm from solid Earth and ocean tides
            "description": "Solid Earth and ocean tide model uncertainties",
            "reference": "Williams et al. 2013"
        },
        "thermal_expansion": {
            "source": "Thermal expansion",
            "magnitude_cm": 0.1,  # ~0.1 cm from telescope thermal effects
            "description": "Telescope and retroreflector thermal expansion",
            "reference": "Murphy 2012"
        }
    }

    # Calculate contributions as percentages of total RMS
    total_budget_cm = sum(v["magnitude_cm"] for v in error_budget.values())
    for key in error_budget:
        error_budget[key]["percentage"] = (error_budget[key]["magnitude_cm"] / total_budget_cm) * 100
        error_budget[key]["variance_contribution"] = (error_budget[key]["magnitude_cm"] / 100.0) ** 2  # Convert to meters

    # Add observed residual RMS
    observed_variance = (global_rms / 100.0) ** 2  # Convert cm to m
    budget_variance = sum(v["variance_contribution"] for v in error_budget.values())
    unexplained_variance = max(0, observed_variance - budget_variance)
    unexplained_rms_cm = np.sqrt(unexplained_variance) * 100 if unexplained_variance > 0 else 0

    error_budget["unexplained"] = {
        "source": "Unexplained residual variance",
        "magnitude_cm": unexplained_rms_cm,
        "description": "Residual variance not accounted for by known systematics",
        "percentage": (unexplained_rms_cm / total_budget_cm) * 100 if total_budget_cm > 0 else 0,
        "variance_contribution": unexplained_variance
    }

    # Print error budget table
    print_status("Systematic Error Budget Table:", "TITLE")
    print_status(f"{'Source':<35} {'Magnitude (cm)':<18} {'% Contribution':<15}", "INFO")
    print_status("-" * 70, "INFO")
    
    for key, value in error_budget.items():
        print_status(f"{value['source']:<35} {value['magnitude_cm']:<18.2f} {value['percentage']:<15.1f}", "INFO")
    
    print_status("-" * 70, "INFO")
    print_status(f"{'Total Budget':<35} {total_budget_cm:<18.2f} {'100.0':<15}", "INFO")
    print_status(f"{'Observed RMS':<35} {global_rms * 100:<18.2f} {'-':<15}", "INFO")
    print_status("", "INFO")

    return error_budget

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Step 008: Systematic Error Analysis with Error Budget")
    args = parser.parse_args()

    log_dir = PROJECT_ROOT / "logs"
    logger = TEPLogger("step_008", str(log_dir / "step_008_systematic_error_analysis.log"))
    set_step_logger(logger)

    print_status("Starting Systematic Error Analysis...", "TITLE")

    input_path = PROJECT_ROOT / 'data/processed/INPOP19a_all_stations_residuals.csv'
    if not input_path.exists():
        print_status(f"Input file not found: {input_path}", "ERROR")
        sys.exit(1)

    df = pd.read_csv(input_path)
    print_status(f"Loaded {len(df):,} observations from {input_path.relative_to(PROJECT_ROOT)}", "INFO")
    print_status("", "INFO")

    sys_results = analyze_systematics(df)
    error_budget = generate_systematic_error_budget(df, sys_results)
    
    all_clean = all(v.get("clean_profile", False)
                    for v in sys_results.values())

    results = {
        "step_id": "step_008",
        "systematics": sys_results,
        "error_budget": error_budget,
        "status": "PASS" if all_clean else "WARNING"
    }

    logger.save_step_results(results, PROJECT_ROOT,
                             "step_008_systematic_error_analysis")
    print_status("Systematic Error Analysis Complete.", "SUCCESS")