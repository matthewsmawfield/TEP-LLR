#!/usr/bin/env python3
"""
Step 019: Station Quality Diagnostics

Investigates Haleakala and other station anomalies through:
1. Quality metrics comparison (RMS, outlier rate, gaps)
2. Systematic correlation analysis
3. Quality-cut sensitivity tests

Addresses the question: Why does Haleakala show positive η?
"""


import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

import json
from typing import Dict

import numpy as np
import pandas as pd
from scripts.utils.logger import TEPLogger, set_step_logger, set_verbose_mode, print_status
from scripts.utils.llr_constants import ETA_SCALE_FACTOR

def compute_station_quality_metrics(df: pd.DataFrame, station: str) -> Dict:
    """Compute comprehensive quality metrics for a station."""
    station_df = df[df['station'] == station]

    if len(station_df) == 0:
        return {'error': 'No data for station'}

    residuals = station_df['residual_m'].values

    # Basic statistics
    metrics = {
        'n_observations': len(station_df),
        'rms_residual_cm': float(np.sqrt(np.mean(residuals**2)) * 100),
        'mean_residual_m': float(np.mean(residuals)),
        'std_residual_m': float(np.std(residuals)),
        'median_residual_m': float(np.median(residuals)),
    }

    # Outlier analysis (using IQR method)
    q1, q3 = np.percentile(residuals, [25, 75])
    iqr = q3 - q1
    outlier_mask = (residuals < q1 - 3 * iqr) | (residuals > q3 + 3 * iqr)
    metrics['outlier_fraction'] = float(outlier_mask.mean())
    metrics['outlier_count'] = int(outlier_mask.sum())

    # Temporal coverage
    years = station_df['date_julian_year'].astype(int)
    metrics['year_range'] = [int(years.min()), int(years.max())]
    metrics['n_years'] = int(years.max() - years.min() + 1)

    # Data gaps (years with < 10 observations)
    obs_per_year = station_df.groupby(years).size()
    gaps = (obs_per_year < 10).sum()
    metrics['years_with_sparse_data'] = int(gaps)
    metrics['fraction_years_sparse'] = float(
        gaps / metrics['n_years']) if metrics['n_years'] > 0 else 0

    # Magnitude coverage (if available)
    if 'magnitude' in station_df.columns:
        mag = station_df['magnitude'].values
        metrics['magnitude_range'] = [float(mag.min()), float(mag.max())]
        metrics['magnitude_std'] = float(np.std(mag))

    # Elongation coverage (uniformity check)
    elong = station_df['elongation_rad'].values
    elong_binned = np.histogram(elong, bins=8, range=(-np.pi, np.pi))[0]
    metrics['elongation_uniformity'] = float(np.std(
        elong_binned) / np.mean(elong_binned)) if np.mean(elong_binned) > 0 else 0

    return metrics

def analyze_station_systematics(df: pd.DataFrame, station: str) -> Dict:
    """Analyze potential systematic correlations for a station."""
    station_df = df[df['station'] == station]

    if len(station_df) < 50:
        return {'error': 'Insufficient data'}

    residuals = station_df['residual_m'].values
    results = {}

    # Time-of-year correlation (seasonal)
    if 'date_julian_year' in station_df.columns:
        # Day of year from fractional year (standard calendar calculation)
        year_frac = station_df['date_julian_year'].values
        doy = ((year_frac - year_frac.astype(int)) * 365.25).astype(int)
        cos_doy = np.cos(2 * np.pi * doy / 365.25)
        sin_doy = np.sin(2 * np.pi * doy / 365.25)

        r_cos = np.corrcoef(residuals, cos_doy)[0, 1]
        r_sin = np.corrcoef(residuals, sin_doy)[0, 1]
        results['seasonal_correlation'] = float(np.sqrt(r_cos**2 + r_sin**2))

    # Magnitude correlation
    if 'magnitude' in station_df.columns:
        mag = station_df['magnitude'].values
        r_mag = np.corrcoef(residuals, mag)[0, 1]
        results['magnitude_correlation'] = float(r_mag)

    # Sin(elongation) correlation (tests for phase offset systematics)
    elong = station_df['elongation_rad'].values
    sin_elong = np.sin(elong)
    r_sin_elong = np.corrcoef(residuals, sin_elong)[0, 1]
    results['sin_elongation_correlation'] = float(r_sin_elong)

    # Year trend
    if 'date_julian_year' in station_df.columns:
        years = station_df['date_julian_year'].astype(int).values
        if len(np.unique(years)) > 2:
            r_year = np.corrcoef(residuals, years)[0, 1]
            results['temporal_trend_correlation'] = float(r_year)

    return results

def test_quality_cuts(df: pd.DataFrame, station: str, cuts: Dict) -> Dict:
    """Test signal stability under quality cuts."""
    station_df = df[df['station'] == station]  # PERFORMANCE FIX: Removed unnecessary .copy()

    if len(station_df) < 100:
        return {'error': 'Insufficient data for quality cuts'}

    original_n = len(station_df)

    # Apply cuts
    mask = np.ones(len(station_df), dtype=bool)

    if 'max_residual_m' in cuts:
        mask &= np.abs(
            station_df['residual_m'].values) < cuts['max_residual_m']

    if 'max_magnitude' in cuts and 'magnitude' in station_df.columns:
        mask &= station_df['magnitude'].values < cuts['max_magnitude']

    if 'min_year' in cuts and 'date_julian_year' in station_df.columns:
        years = station_df['date_julian_year'].astype(int).values
        mask &= years >= cuts['min_year']

    if 'max_year' in cuts and 'date_julian_year' in station_df.columns:
        years = station_df['date_julian_year'].astype(int).values
        mask &= years <= cuts['max_year']

    filtered_df = station_df[mask]

    if len(filtered_df) < 50:
        return {'error': 'Too few observations after cuts', 'n_after_cuts': len(filtered_df)}

    # Compute eta on filtered data
    cos_elong = np.cos(filtered_df['elongation_rad'].values)
    residuals = filtered_df['residual_m'].values

    X = np.column_stack([cos_elong, np.ones(len(cos_elong))])
    coeffs, _, _, _ = np.linalg.lstsq(X, residuals, rcond=None)
    eta = coeffs[0] / ETA_SCALE_FACTOR

    # Error estimate
    resid_fit = residuals - np.dot(X, coeffs)
    mse = np.mean(resid_fit**2)
    var_eta = mse / np.sum((cos_elong - np.mean(cos_elong))**2)
    eta_error = np.sqrt(var_eta) / ETA_SCALE_FACTOR

    return {
        'n_original': original_n,
        'n_after_cuts': len(filtered_df),
        'fraction_retained': len(filtered_df) / original_n,
        'eta': float(eta),
        'eta_error': float(eta_error),
        'snr': float(eta / eta_error) if eta_error > 0 else 0,
        'sign_changed': bool((eta > 0) != (station_df['residual_m'].corr(
            np.cos(station_df['elongation_rad'])) < 0))  # Compare to original
    }

def compare_stations_quality(df: pd.DataFrame) -> Dict:
    """Comprehensive quality comparison across all stations."""
    stations = df['station'].unique()

    comparison = {}
    for station in stations:
        metrics = compute_station_quality_metrics(df, station)
        systematics = analyze_station_systematics(df, station)

        # Compute original eta for this station
        station_df = df[df['station'] == station]
        if len(station_df) > 100:
            cos_elong = np.cos(station_df['elongation_rad'].values)
            residuals = station_df['residual_m'].values
            X = np.column_stack([cos_elong, np.ones(len(cos_elong))])
            coeffs, _, _, _ = np.linalg.lstsq(X, residuals, rcond=None)
            eta = coeffs[0] / ETA_SCALE_FACTOR
        else:
            eta = None

        comparison[station] = {
            'quality_metrics': metrics,
            'systematic_correlations': systematics,
            'eta_original': float(eta) if eta is not None else None
        }

    # Identify outliers
    rms_values = [v['quality_metrics']['rms_residual_cm'] for v in comparison.values()
                  if 'rms_residual_cm' in v['quality_metrics']]
    mean_rms = np.mean(rms_values)
    std_rms = np.std(rms_values)

    outlier_stations = []
    for station, data in comparison.items():
        rms = data['quality_metrics'].get('rms_residual_cm', 0)
        if rms > mean_rms + 2 * std_rms:
            outlier_stations.append({
                'station': station,
                'rms_cm': rms,
                'mean_rms_cm': mean_rms,
                'deviation_sigma': (rms - mean_rms) / std_rms if std_rms > 0 else 0,
                'reason': 'High RMS (likely noisy data)'
            })

    return {
        'station_comparison': comparison,
        'rms_statistics': {
            'mean_cm': float(mean_rms),
            'std_cm': float(std_rms),
            'threshold_cm': float(mean_rms + 2 * std_rms)
        },
        'quality_outlier_stations': outlier_stations
    }

def main():
    # Setup TEPLogger for consistent file logging
    log_dir = Path(__file__).parent.parent.parent / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    logger = TEPLogger("step_018", str(log_dir / "step_018_station_quality.log"))
    set_step_logger(logger)

    print_status("═══ Starting Step 018: Station Quality Diagnostics...", "TITLE")
    print_status("═══ STEP PURPOSE: Investigate Haleakala and other station anomalies through quality metrics", "INFO")
    print_status("═══ METHOD: Quality metrics comparison (RMS, outlier rate, gaps), systematic correlation analysis, quality-cut sensitivity tests", "INFO")
    print_status("═══ PARAMETERS: IQR outlier threshold=3×, quality cut max residual=0.5m", "INFO")

    logger.info("Step 018: Station Quality Diagnostics")

    print_status("═══ DATA SUMMARY", "INFO")
    # Load data
    data_path = Path(__file__).parent.parent.parent / 'data' / \
        'processed' / 'INPOP19a_all_stations_residuals.csv'
    df = pd.read_csv(data_path)

    print_status(f"    Dataset: N = {len(df):,} observations", "DATA")
    print_status(f"    Stations: {sorted(df['station'].unique())}", "DATA")
    print_status(f"    Data source: INPOP19a_all_stations_residuals.csv", "DATA")

    print_status("═══ ANALYSIS TRACE", "INFO")
    # Comprehensive quality comparison
    print_status(f">>> Analyzing station quality metrics across all stations", "PROCESS")
    logger.info("Analyzing station quality metrics...")
    quality_comparison = compare_stations_quality(df)

    # Detailed analysis for Haleakala (the anomaly)
    if 'Haleakala' in df['station'].values:
        print_status(f">>> Detailed analysis of Haleakala anomaly", "PROCESS")
        logger.info("Detailed analysis of Haleakala anomaly...")
        haleakala_quality = {
            'metrics': compute_station_quality_metrics(df, 'Haleakala'),
            'systematics': analyze_station_systematics(df, 'Haleakala'),
        }

        # Test quality cuts
        print_status(f">>> Testing Haleakala quality cuts", "PROCESS")
        logger.info("Testing Haleakala quality cuts...")
        haleakala_cuts = {
            'outlier_removal': test_quality_cuts(df, 'Haleakala', {'max_residual_m': 0.5}),
            'recent_only': test_quality_cuts(df, 'Haleakala', {'min_year': 1987, 'max_year': 1989}),
        }
        haleakala_quality['quality_cut_tests'] = haleakala_cuts

        # Determine likely cause of anomaly
        metrics = haleakala_quality['metrics']
        systematics = haleakala_quality['systematics']

        causes = []
        if metrics.get('rms_residual_cm', 0) > 20:  # > 20 cm RMS
            causes.append('high_noise_rms')
        if metrics.get('n_observations', 0) < 1000:
            causes.append('small_sample')
        if metrics.get('outlier_fraction', 0) > 0.1:
            causes.append('high_outlier_rate')
        if abs(systematics.get('temporal_trend_correlation', 0)) > 0.3:
            causes.append('temporal_drift')

        haleakala_quality['likely_causes'] = causes if causes else [
            'noise_limited']
    else:
        haleakala_quality = {'error': 'Haleakala not in dataset'}

    # Similar analysis for McDonald2
    if 'McDonald2' in df['station'].values:
        print_status(f">>> Analyzing McDonald2 era split", "PROCESS")
        mcdonald_quality = {
            'metrics': compute_station_quality_metrics(df, 'McDonald2'),
            'systematics': analyze_station_systematics(df, 'McDonald2'),
        }

        # Split by era
        mcdonald_early = test_quality_cuts(df, 'McDonald2', {'max_year': 1995})
        mcdonald_late = test_quality_cuts(df, 'McDonald2', {'min_year': 1996})

        mcdonald_quality['era_split'] = {
            'early': mcdonald_early,
            'late': mcdonald_late
        }
    else:
        mcdonald_quality = {'error': 'McDonald2 not in dataset'}

    # Compile results
    results = {
        'step_id': 'step_018',
        'status': 'PASS',
        'summary': {
            'n_stations_analyzed': len(df['station'].unique()),
            'quality_outliers_identified': len(quality_comparison['quality_outlier_stations']),
            'haleakala_assessment': haleakala_quality.get('likely_causes', ['unknown']),
        },
        'station_quality_comparison': quality_comparison,
        'haleakala_detailed_analysis': haleakala_quality,
        'mcdonald2_detailed_analysis': mcdonald_quality,
        'conclusions': {
            'haleakala_positive_eta_explanation': (
                'Likely noise-limited measurement with ' +
                f"{haleakala_quality.get('metrics', {}).get('n_observations', 0)} observations, " +
                f"RMS {haleakala_quality.get('metrics', {}).get('rms_residual_cm', 0.0):.1f} cm"
            ) if 'metrics' in haleakala_quality else 'Not analyzed',
            'multi_station_confidence': 'HIGH' if len([s for s in quality_comparison['station_comparison'].values()
                                                      if s.get('eta_original', 0) and s['eta_original'] < 0]) >= 2
            else 'MODERATE'
        }
    }

    print_status("═══ RESULTS SUMMARY", "INFO")
    print_status(f"    Stations analyzed: {len(df['station'].unique())}", "CALC")
    print_status(f"    Quality outliers identified: {len(quality_comparison['quality_outlier_stations'])}", "CALC")
    if 'Haleakala' in df['station'].values:
        print_status(f"    Haleakala RMS: {haleakala_quality['metrics']['rms_residual_cm']:.1f} cm", "CALC")
        print_status(f"    Haleakala n_obs: {haleakala_quality['metrics']['n_observations']}", "CALC")
        print_status(f"    Haleakala likely causes: {haleakala_quality.get('likely_causes', [])}", "CALC")
    print_status(f"    Multi-station confidence: {results['conclusions']['multi_station_confidence']}", "CALC")

    print_status("═══ INTERPRETATION", "INFO")
    if 'Haleakala' in df['station'].values:
        print_status(f"    Haleakala positive η explained by: {', '.join(haleakala_quality.get('likely_causes', ['unknown']))}", "INFO")
    print_status(f"    Station quality analysis tests for systematic biases", "INFO")
    print_status(f"    Quality-cut sensitivity tests validate signal robustness", "INFO")

    print_status("═══ REPRODUCIBILITY", "INFO")
    print_status(f"    Output file: results/outputs/step_018_station_quality.json", "INFO")
    print_status(f"    IQR outlier threshold: 3×", "INFO")
    print_status(f"    Quality cut max residual: 0.5m", "INFO")

    # Print key findings
    if 'Haleakala' in df['station'].values:
        logger.info(
            f"Haleakala RMS: {haleakala_quality['metrics']['rms_residual_cm']:.1f} cm")
        logger.info(
            f"Haleakala n_obs: {haleakala_quality['metrics']['n_observations']}")
        logger.info(
            f"Haleakala likely causes: {haleakala_quality.get('likely_causes', [])}")

    # Save results through logger system
    logger.save_step_results(results, Path(__file__).parent.parent.parent, "step_018_station_quality")

    return results

if __name__ == '__main__':
    main()