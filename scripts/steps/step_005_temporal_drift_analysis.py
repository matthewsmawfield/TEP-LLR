#!/usr/bin/env python3
"""
Step 005: Temporal Drift Analysis for TEP Nordtvedt Signal Detection
Enhanced with temporal autocorrelation analysis
"""

import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from scripts.utils.logger import TEPLogger, set_step_logger, print_status, set_verbose_mode
from scripts.utils.statistical_utils import linear_regression
from scripts.utils.llr_constants import TEMPORAL_DRIFT_MAX_LAG, TEMPORAL_DRIFT_ERA_SPLIT_YEAR

import argparse
import pandas as pd
import numpy as np
from scipy import stats

# Add the project root to the Python path

def analyze_temporal_autocorrelation(df, station, max_lag=TEMPORAL_DRIFT_MAX_LAG, verbose=False):
    """
    Analyze temporal autocorrelation in residuals to check for temporal dependencies
    """
    station_data = df[df['station'] == station].sort_values('date_julian_year')
    if len(station_data) < 50:
        return {"error": "Insufficient data for autocorrelation analysis"}

    residuals = station_data['residual_m'].values

    # Calculate autocorrelation function
    autocorr = []
    lags = []
    for lag in range(1, min(max_lag + 1, len(residuals) // 2)):
        if lag >= len(residuals):
            break
        corr = np.corrcoef(residuals[:-lag], residuals[lag:])[0, 1]
        autocorr.append(corr)
        lags.append(lag)

    # Calculate Durbin-Watson statistic
    dw_stat = 0
    if len(residuals) > 1:
        diff_residuals = np.diff(residuals)
        dw_stat = np.sum(diff_residuals ** 2) / np.sum(residuals ** 2)

    # Test for significant autocorrelation
    significant_lags = []
    for i, (lag, corr) in enumerate(zip(lags, autocorr)):
        # Standard error for autocorrelation under null hypothesis of no autocorrelation
        se = 1 / np.sqrt(len(residuals))
        z_score = corr / se
        p_value = 2 * (1 - stats.norm.cdf(abs(z_score)))
        if p_value < 0.05:
            significant_lags.append({"lag": lag, "autocorr": corr, "p_value": p_value})

    if verbose:
        print_status(f"  {station} Temporal Autocorrelation Analysis:", "CALC")
        print_status(f"    Durbin-Watson statistic: {dw_stat:.3f}", "CALC")
        print_status(f"    Significant autocorrelations at 5% level: {len(significant_lags)}", "CALC")
        if significant_lags:
            for sig in significant_lags[:3]:  # Show first 3 significant lags
                print_status(f"      Lag {sig['lag']}: r={sig['autocorr']:.3f}, p={sig['p_value']:.3f}", "CALC")

    return {
        "durbin_watson": float(dw_stat),
        "autocorrelations": autocorr,
        "lags": lags,
        "significant_autocorrelations": significant_lags,
        "interpretation": "No significant temporal autocorrelation" if len(significant_lags) == 0 else f"Significant autocorrelation at {len(significant_lags)} lags"
    }

def era_analysis(df, station, split_year=TEMPORAL_DRIFT_ERA_SPLIT_YEAR, verbose=False):
    station_data = df[df['station'] == station]  # PERFORMANCE FIX: Removed unnecessary .copy()
    if len(station_data) < 100:
        return {"error": "Insufficient data"}

    early = station_data[station_data['date_julian_year'] < split_year]
    late = station_data[station_data['date_julian_year'] >= split_year]

    res = {"early": {}, "late": {}}
    for label, data in [("early", early), ("late", late)]:
        if len(data) >= 50:
            reg = linear_regression(data['residual_m'].values, np.cos(
                data['elongation_rad'].values))
            res[label] = {"n_obs": len(data), "eta": reg['eta'], "snr": abs(
                reg['eta'])/reg['eta_error']}
            if verbose:
                print_status(
                    f"  {station} {label} era (Split={split_year}):", "CALC")
                print_status(
                    f"    η   = {reg['eta']:.4e} ± {reg['eta_error']:.4e}", "CALC")
                print_status(f"    SNR = {res[label]['snr']:.2f}σ", "CALC")
                print_status(f"    N   = {len(data)}", "CALC")

    return res

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Step 005: Temporal Drift Analysis with Autocorrelation")
    args = parser.parse_args()

    log_dir = PROJECT_ROOT / "logs"
    logger = TEPLogger("step_005", str(
        log_dir / "step_005_temporal_drift_analysis.log"))
    set_step_logger(logger)
    set_verbose_mode(True)

    print_status("Starting Temporal Drift Analysis...", "TITLE")

    input_path = PROJECT_ROOT / 'data/processed/INPOP19a_all_stations_residuals.csv'
    if not input_path.exists():
        print_status(f"Input file not found: {input_path}", "ERROR")
        sys.exit(1)

    df = pd.read_csv(input_path)
    stations = df['station'].unique()

    print_status("", "INFO")
    print_status("Temporal Autocorrelation Analysis:", "TITLE")
    print_status("", "INFO")

    autocorr_res = {}
    era_res = {}
    for s in stations:
        autocorr_res[s] = analyze_temporal_autocorrelation(df, s)
        era_res[s] = era_analysis(df, s)

    print_status("", "INFO")
    print_status("Autocorrelation Summary:", "TITLE")
    stations_with_autocorr = sum(1 for v in autocorr_res.values() if "significant_autocorrelations" in v and v["significant_autocorrelations"])
    print_status(f"  Stations with significant autocorrelation: {stations_with_autocorr}/{len(stations)}", "PASS" if stations_with_autocorr == 0 else "WARNING")
    print_status("", "INFO")

    results = {
        "step_id": "step_005",
        "autocorrelation_analysis": autocorr_res,
        "era_analysis": era_res,
        "status": "PASS"
    }

    logger.save_step_results(results, PROJECT_ROOT,
                             "step_005_temporal_drift_analysis")
    print_status("Temporal Drift Analysis Complete.", "SUCCESS")
