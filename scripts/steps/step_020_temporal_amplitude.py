#!/usr/bin/env python3
"""
Step 020: Temporal Amplitude Evolution Analysis

Investigates whether the TEP signal amplitude varies over time:
1. Sliding window correlation analysis
2. Tests for secular trends in amplitude
3. Correlates amplitude with known instrumental changes

Addresses: Is early Grasse η = -20.6×10⁻⁴ vs late η = -3.56×10⁻⁴ physical or systematic?
"""


import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

import numpy as np
from scripts.utils.numerics import stable_lstsq
import pandas as pd
import json
from typing import Dict, List
from scipy import stats
from scripts.utils.logger import TEPLogger, set_step_logger, set_verbose_mode, print_status
from scripts.utils.llr_constants import ETA_SCALE_FACTOR

def sliding_window_eta(df: pd.DataFrame, window_years: int = 5,
                       step_years: int = 2) -> List[Dict]:
    """Compute eta in sliding time windows."""
    df = df.copy()  # Keep this copy since we add a column
    df['year'] = df['date_julian_year'].astype(int)

    min_year = int(df['year'].min())
    max_year = int(df['year'].max())

    results = []
    for start_year in range(min_year, max_year - window_years + 1, step_years):
        end_year = start_year + window_years

        mask = (df['year'] >= start_year) & (df['year'] < end_year)
        window_df = df[mask]

        if len(window_df) < 500:  # Require minimum statistics
            continue

        cos_elong = np.cos(window_df['elongation_rad'].values)
        residuals = window_df['residual_m'].values

        # OLS fit
        X = np.column_stack([cos_elong, np.ones(len(cos_elong))])
        coeffs, _, _, _ = stable_lstsq(X, residuals)
        eta = coeffs[0] / ETA_SCALE_FACTOR

        # Error estimate (unbiased MSE: divide by n-2, not n)
        with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
            resid_fit = residuals - X @ coeffs
        n_win = len(resid_fit)
        mse = np.sum(resid_fit**2) / (n_win - 2) if n_win > 2 else np.nan
        cos_centered = cos_elong - np.mean(cos_elong)
        var_A = mse / np.sum(cos_centered ** 2) if np.sum(cos_centered**2) > 0 else np.nan
        eta_error = np.sqrt(var_A) / ETA_SCALE_FACTOR if np.isfinite(var_A) and var_A > 0 else np.nan

        # Correlation
        r = np.corrcoef(residuals, cos_elong)[0, 1]

        is_significant = False
        if np.isfinite(eta_error) and eta_error > 0:
            is_significant = bool(abs(eta / eta_error) > 3) if eta_error > 0 else False

        results.append({
            'window_start': start_year,
            'window_end': end_year,
            'window_center': (start_year + end_year) / 2,
            'n_observations': len(window_df),
            'eta': float(eta),
            'eta_error': float(eta_error),
            'snr': float(eta / eta_error) if np.isfinite(eta_error) and eta_error > 0 else 0,
            'correlation_r': float(r),
            'significant': is_significant
        })

    return results

def test_amplitude_trend(window_results: List[Dict]) -> Dict:
    """Test for secular trend in amplitude over time."""
    if len(window_results) < 3:
        return {'error': 'Insufficient windows'}

    years = np.array([w['window_center'] for w in window_results])
    etas = np.array([w['eta'] for w in window_results])
    errors = np.array([w['eta_error'] for w in window_results])

    # Linear regression of eta vs time
    slope, intercept, r_value, p_value, std_err = stats.linregress(years, etas)

    # Weighted regression (accounting for error bars)
    weights = 1 / errors**2
    X = np.column_stack([years, np.ones(len(years))])
    W = np.diag(weights)
    XtWX = X.T @ W @ X
    XtWy = X.T @ W @ etas
    try:
        coeffs_w = np.linalg.solve(XtWX, XtWy)
        slope_weighted = coeffs_w[0]
        cov_w = np.linalg.pinv(XtWX, rcond=1e-10, hermitian=True)
        slope_err_weighted = np.sqrt(cov_w[0, 0])
    except (np.linalg.LinAlgError, ValueError):
        slope_weighted = slope
        slope_err_weighted = std_err

    # Chi-squared test for consistency with constant amplitude
    if len(etas) > 1:
        mean_eta = np.average(etas, weights=weights)
        chi2 = np.sum(weights * (etas - mean_eta)**2)
        dof = len(etas) - 1
        chi2_per_dof = chi2 / dof if dof > 0 else 0
    else:
        chi2_per_dof = 0

    return {
        'linear_trend_slope_per_year': float(slope),
        'linear_trend_pvalue': float(p_value),
        'linear_trend_significant': bool(p_value < 0.05),
        'weighted_slope_per_year': float(slope_weighted),
        'weighted_slope_error': float(slope_err_weighted),
        'chi2_consistency': float(chi2_per_dof),
        'interpretation': 'amplitude_varies_significantly' if chi2_per_dof > 3 else 'consistent_with_constant'
    }


def correlate_with_instrumental_changes(window_results: List[Dict],
                                        instrument_timeline: Dict) -> Dict:
    """Test if amplitude correlates with known instrumental changes."""
    # Timeline of major LLR events
    # These would need to be verified against actual observatory records

    # Example timeline (to be replaced with actual data if available)
    events = {
        1984: 'Haleakala_start',
        1987: 'McDonald2_upgrade',
        1994: 'Grasse_laser_upgrade',
        2000: 'APO_start',
        2005: 'Grasse_detector_replacement',
        2012: 'APO_timing_upgrade',
    }

    # Find windows near events
    results = {'events_analyzed': []}

    for year, event in events.items():
        nearby_windows = [w for w in window_results
                          if abs(w['window_center'] - year) <= 3]
        if nearby_windows:
            avg_eta = np.mean([w['eta'] for w in nearby_windows])
            results['events_analyzed'].append({
                'year': year,
                'event': event,
                'nearby_eta': float(avg_eta),
                'n_windows': len(nearby_windows)
            })

    # Simple test: does amplitude jump near upgrade years?
    if len(results['events_analyzed']) >= 2:
        etas_by_event = [e['nearby_eta'] for e in results['events_analyzed']]
        etas_other = [w['eta'] for w in window_results
                      if not any(abs(w['window_center'] - e['year']) <= 3
                                 for e in results['events_analyzed'])]

        if etas_other:
            # Two-sample test
            t_stat, p_val = stats.ttest_ind(etas_by_event, etas_other)
            results['instrumental_correlation'] = {
                't_statistic': float(t_stat),
                'p_value': float(p_val),
                'correlated': bool(p_val < 0.05)
            }

    return results

def station_specific_temporal_analysis(df: pd.DataFrame, station: str,
                                       window_years: int = 7) -> List[Dict]:
    """Temporal evolution for a single station."""
    station_df = df[df['station'] == station].copy()

    if len(station_df) < 1000:
        return []

    station_df['year'] = station_df['date_julian_year'].astype(int)

    min_year = int(station_df['year'].min())
    max_year = int(station_df['year'].max())

    results = []
    for start_year in range(min_year, max_year - window_years + 1, 3):
        end_year = start_year + window_years
        mask = (station_df['year'] >= start_year) & (
            station_df['year'] < end_year)
        window_df = station_df[mask]

        if len(window_df) < 200:
            continue

        cos_elong = np.cos(window_df['elongation_rad'].values)
        residuals = window_df['residual_m'].values

        X = np.column_stack([cos_elong, np.ones(len(cos_elong))])
        coeffs, _, _, _ = stable_lstsq(X, residuals)
        eta = coeffs[0] / ETA_SCALE_FACTOR

        results.append({
            'station': station,
            'window_start': start_year,
            'window_end': end_year,
            'n_obs': len(window_df),
            'eta': float(eta),
            'rms_cm': float(np.sqrt(np.mean(residuals**2)) * 100)
        })

    return results

def main():
    # Setup TEPLogger for consistent file logging
    log_dir = Path(__file__).parent.parent.parent / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    logger = TEPLogger("step_020", str(log_dir / "step_020_temporal_amplitude.log"))
    set_step_logger(logger)

    print_status("═══ Starting Step 020: Temporal Amplitude Evolution Analysis...", "TITLE")
    print_status("═══ STEP PURPOSE: Investigate whether TEP signal amplitude varies over time (secular trends, instrumental changes)", "INFO")
    print_status("═══ METHOD: Sliding window correlation analysis, secular trend testing, instrumental timeline correlation", "INFO")
    print_status("═══ PARAMETERS: Window size=5 years, step size=2 years, station window=7 years", "INFO")

    logger.info("Step 020: Temporal Amplitude Evolution Analysis")

    print_status("═══ DATA SUMMARY", "INFO")
    # Load data
    data_path = Path(__file__).parent.parent.parent / 'data' / \
        'processed' / 'INPOP19a_all_stations_residuals.csv'
    df = pd.read_csv(data_path)

    print_status(f"    Dataset: N = {len(df):,} observations", "DATA")
    print_status(f"    Data source: INPOP19a_all_stations_residuals.csv", "DATA")

    print_status("═══ ANALYSIS TRACE", "INFO")
    # Global sliding window analysis
    print_status(f">>> Running sliding window analysis", "PROCESS")
    logger.info("Running sliding window analysis...")
    window_results = sliding_window_eta(df, window_years=5, step_years=2)

    # Test for trend
    print_status(f">>> Testing for secular trends", "PROCESS")
    logger.info("Testing for secular trends...")
    trend_analysis = test_amplitude_trend(window_results)

    # Instrumental correlation
    print_status(f">>> Correlating with instrumental timeline", "PROCESS")
    logger.info("Correlating with instrumental timeline...")
    instrumental = correlate_with_instrumental_changes(window_results, {})

    # Station-specific evolution (focus on Grasse and APO)
    print_status(f">>> Station-specific temporal analysis (Grasse, APO)", "PROCESS")
    logger.info("Station-specific temporal analysis...")
    grasse_evolution = station_specific_temporal_analysis(
        df, 'Grasse', window_years=7)
    apo_evolution = station_specific_temporal_analysis(
        df, 'APO', window_years=7)

    # Compile results
    results = {
        'step_id': 'step_020',
        'status': 'PASS',
        'summary': {
            'n_windows': len(window_results),
            'significant_windows': len([w for w in window_results if w['significant']]),
            'amplitude_trend_significant': trend_analysis.get('linear_trend_significant', False),
            'consistency_assessment': trend_analysis.get('interpretation', 'unknown')
        },
        'sliding_window_results': window_results,
        'trend_analysis': trend_analysis,
        'instrumental_correlation': instrumental,
        'station_evolution': {
            'Grasse': grasse_evolution,
            'APO': apo_evolution
        },
        'conclusions': {
            'physical_variation_likely': bool(trend_analysis.get('chi2_consistency', 0) < 3 and
                                              not trend_analysis.get('linear_trend_significant', True)),
            'systematic_drift_likely': bool(trend_analysis.get('chi2_consistency', 0) > 3 or
                                            trend_analysis.get('linear_trend_significant', False)),
            'key_finding': 'amplitude_stable' if trend_analysis.get('chi2_consistency', 0) < 3
            else 'amplitude_varies_instrumental'
        }
    }

    print_status("═══ RESULTS SUMMARY", "INFO")
    print_status(f"    Windows analyzed: {len(window_results)}", "CALC")
    print_status(f"    Significant windows: {results['summary']['significant_windows']}", "CALC")
    print_status(f"    Amplitude trend significant: {results['summary']['amplitude_trend_significant']}", "CALC")
    print_status(f"    Consistency assessment: {results['summary']['consistency_assessment']}", "CALC")

    print_status("═══ INTERPRETATION", "INFO")
    print_status(f"    Trend interpretation: {results['conclusions']['key_finding']}", "INFO")
    print_status(f"    Physical variation likely: {results['conclusions']['physical_variation_likely']}", "INFO")
    print_status(f"    Systematic drift likely: {results['conclusions']['systematic_drift_likely']}", "INFO")

    print_status("═══ REPRODUCIBILITY", "INFO")
    print_status(f"    Output file: results/outputs/step_020_temporal_amplitude.json", "INFO")
    print_status(f"    Window size: 5 years", "INFO")
    print_status(f"    Step size: 2 years", "INFO")

    # Save results
    project_root = Path(__file__).parent.parent.parent
    output_path = project_root / 'results' / 'outputs' / 'step_020_temporal_amplitude.json'
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, 'w') as f:
        json.dump(results, f, indent=2)

    output_rel = output_path.relative_to(project_root) if output_path.is_relative_to(project_root) else output_path
    logger.info(f"Results saved to {output_rel}")
    logger.info(f"Windows analyzed: {len(window_results)}")
    logger.info(
        f"Significant windows: {results['summary']['significant_windows']}")
    logger.info(
        f"Trend interpretation: {results['conclusions']['key_finding']}")

    return results

if __name__ == '__main__':
    main()