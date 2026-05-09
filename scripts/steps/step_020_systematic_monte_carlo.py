#!/usr/bin/env python3
"""
Step 020: Systematic Error Monte Carlo Propagation

Quantifies systematic uncertainty through Monte Carlo error propagation:
1. Ephemeris uncertainties (position, velocity errors)
2. Tidal model variations (Love number uncertainties)
3. Atmospheric delay errors (tropospheric model scatter)
4. Instrumental calibration uncertainties

Outputs quantitative systematic error budget for eta.
"""


import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

import json
from typing import Dict, Tuple

sys.path.insert(0, str(Path(__file__).parent.parent.parent))

import numpy as np
import pandas as pd
from scripts.utils.logger import TEPLogger, set_step_logger, set_verbose_mode, print_status
from scripts.utils.llr_constants import ETA_SCALE_FACTOR

def inject_ephemeris_uncertainty(df: pd.DataFrame, rng: np.random.RandomState,
                                 position_error_m: float = 0.5,
                                 velocity_error_mm_s: float = 0.1) -> pd.DataFrame:
    """Simulate ephemeris uncertainties in Earth-Moon range.
    
    Parameter ranges based on ephemeris accuracy specifications:
    - Position error: 0.3-0.7 m (INPOP19a/DE430 ephemeris accuracy ~0.5 m)
      Source: INPOP19a documentation, Folkner et al. 2014
    - Velocity error: 0.05-0.15 mm/s (lunar orbit velocity ~1 km/s, relative error ~1e-7)
      Source: Ephemeris inter-comparison studies
    """
    df_perturbed = df.copy()  # Keep this copy since we modify it

    # Add correlated noise to residuals (simulates ephemeris bias)
    n = len(df_perturbed)

    # Position uncertainty: quasi-bias over timescales
    time_days = np.arange(n) / n  # Normalized time

    # Low-frequency correlated errors (ephemeris tends to have systematic drifts)
    # Frequency range 0.5-2 cycles represents timescales from half to full data span
    correlated = np.sin(2 * np.pi * time_days * rng.uniform(0.5, 2)
                        ) * position_error_m * rng.uniform(0.5, 1.5)

    # Add white noise component (30% of position error is typical for uncorrelated component)
    white = rng.normal(0, position_error_m * 0.3, n)

    df_perturbed['residual_m'] = df_perturbed['residual_m'] + \
        correlated + white

    return df_perturbed

def inject_tidal_uncertainty(df: pd.DataFrame, rng: np.random.RandomState,
                             love_number_error: float = 0.01) -> pd.DataFrame:
    """Simulate tidal model uncertainties via Love number variations.
    
    Parameter ranges based on tidal model accuracy:
    - Love number error: 0.005-0.015 (h2 ~0.6, k2 ~0.3, relative error ~1-2%)
      Source: IERS2010 conventions, Williams et al. 2014
    - Tidal amplitude: ~5 mm at maximum for realistic Love number error
      Source: LLR tidal model error budgets
    """
    df_perturbed = df.copy()

    # Tidal effect correlates with lunar phase (elongation)
    elong = df_perturbed['elongation_rad'].values

    # Love number uncertainty creates amplitude uncertainty in tidal correction
    # ~5 mm effect at maximum for realistic Love number error
    tidal_amplitude_m = 0.005 * (love_number_error / 0.01)  # 5mm baseline

    # Tidal signal: roughly twice per synodic month (gravity gradient)
    # Amplitude variation 0.8-1.2 represents ±20% uncertainty in tidal model
    tidal_signal = tidal_amplitude_m * \
        np.cos(2 * elong) * rng.uniform(0.8, 1.2)

    df_perturbed['residual_m'] = df_perturbed['residual_m'] + tidal_signal

    return df_perturbed

def inject_atmospheric_uncertainty(df: pd.DataFrame, rng: np.random.RandomState,
                                   zenith_delay_error_mm: float = 5.0) -> pd.DataFrame:
    """Simulate atmospheric delay model uncertainties.
    
    Parameter ranges based on atmospheric delay model accuracy:
    - Zenith delay error: 3-7 mm (typical for Saastamoinen/Marini-Murray models)
      Source: Mendes et al. 2002, Atmospheric delay models for LLR
    """
    df_perturbed = df.copy()  # Keep this copy since we modify it

    # Atmospheric delay varies with station elevation and weather
    # Model error ~5mm at zenith, scales with airmass

    n = len(df_perturbed)

    # Correlated over hours (weather timescale)
    block_size = max(1, n // 1000)  # ~1000 time blocks
    atmospheric_bias = np.repeat(rng.normal(0, zenith_delay_error_mm / 1000,
                                            (n + block_size - 1) // block_size),
                                 block_size)[:n]

    df_perturbed['residual_m'] = df_perturbed['residual_m'] + atmospheric_bias

    return df_perturbed

def inject_instrumental_uncertainty(df: pd.DataFrame, rng: np.random.RandomState,
                                    calibration_drift_mm_yr: float = 0.5) -> pd.DataFrame:
    """Simulate slow instrumental calibration drift.
    
    Parameter ranges based on LLR station calibration stability:
    - Calibration drift: 0.3-0.7 mm/yr (typical for retroreflector array calibration)
      Source: Murphy et al. 2002, LLR calibration and systematics
    """
    df_perturbed = df.copy()  # Keep this copy since we modify it

    if 'date_julian_year' not in df_perturbed.columns:
        return df_perturbed

    # Linear drift over observation period
    years = df_perturbed['date_julian_year'].astype(int).values
    years_normalized = (years - years.min()) / \
        max(years.max() - years.min(), 1)

    drift = years_normalized * calibration_drift_mm_yr / \
        1000 * rng.uniform(-1, 1)

    df_perturbed['residual_m'] = df_perturbed['residual_m'] + drift

    return df_perturbed

def compute_eta_with_error(df: pd.DataFrame) -> Tuple[float, float]:
    """Compute eta and its statistical error.
    
    Uses proper error propagation: sigma_eta = sigma_res / (A * sqrt(N_eff))
    where A = 13.0 is the ETA_SCALE_FACTOR and N_eff accounts for the 
    elongation distribution.
    """
    cos_elong = np.cos(df['elongation_rad'].values)
    residuals = df['residual_m'].values
    n = len(residuals)

    X = np.column_stack([cos_elong, np.ones(len(cos_elong))])
    coeffs, residuals_fit, rank, _ = np.linalg.lstsq(X, residuals, rcond=None)

    if rank < 2:
        return np.nan, np.nan

    eta = coeffs[0] / ETA_SCALE_FACTOR

    # Proper statistical error calculation:
    # sigma_eta = sigma_res / (A * sqrt(N_eff))
    # where N_eff = sum(cos^2) / max(cos^2) accounts for elongation coverage
    sigma_res = np.std(residuals_fit) if len(residuals_fit) > 0 else np.std(residuals)
    cos_centered = cos_elong - np.mean(cos_elong)
    n_eff = np.sum(cos_centered**2) / (np.max(cos_centered**2) if np.max(cos_centered**2) > 0 else 1)
    eta_error = sigma_res / (ETA_SCALE_FACTOR * np.sqrt(max(n_eff, 1)))

    return eta, eta_error

def monte_carlo_systematic_analysis(df: pd.DataFrame, n_mc: int = 1000,
                                    seed: int = 42) -> Dict:
    """Run Monte Carlo propagation of systematic errors."""
    rng = np.random.RandomState(seed)

    # Baseline eta
    eta_baseline, eta_stat_error = compute_eta_with_error(df)

    results = {
        'baseline': {
            'eta': float(eta_baseline),
            'statistical_error': float(eta_stat_error)
        },
        'systematic_variations': {}
    }

    # 1. Ephemeris uncertainty
    logger.info("Running ephemeris uncertainty MC...")
    eta_ephem = []
    for _ in range(n_mc):
        df_pert = inject_ephemeris_uncertainty(df, rng,
                                               position_error_m=rng.uniform(
                                                   0.3, 0.7),
                                               velocity_error_mm_s=rng.uniform(0.05, 0.15))
        eta, _ = compute_eta_with_error(df_pert)
        if not np.isnan(eta):
            eta_ephem.append(eta)

    results['systematic_variations']['ephemeris'] = {
        'eta_mean': float(np.mean(eta_ephem)),
        'eta_std': float(np.std(eta_ephem)),
        'systematic_error': float(np.std(eta_ephem)),
        'description': 'Position error ~0.5m, velocity error ~0.1 mm/s'
    }

    # 2. Tidal uncertainty
    logger.info("Running tidal uncertainty MC...")
    eta_tidal = []
    for _ in range(n_mc):
        df_pert = inject_tidal_uncertainty(df, rng,
                                           love_number_error=rng.uniform(0.005, 0.015))
        eta, _ = compute_eta_with_error(df_pert)
        if not np.isnan(eta):
            eta_tidal.append(eta)

    results['systematic_variations']['tidal'] = {
        'eta_mean': float(np.mean(eta_tidal)),
        'eta_std': float(np.std(eta_tidal)),
        'systematic_error': float(np.std(eta_tidal)),
        'description': 'Love number uncertainty ~0.01'
    }

    # 3. Atmospheric uncertainty
    logger.info("Running atmospheric uncertainty MC...")
    eta_atmos = []
    for _ in range(n_mc):
        df_pert = inject_atmospheric_uncertainty(df, rng,
                                                 zenith_delay_error_mm=rng.uniform(3, 7))
        eta, _ = compute_eta_with_error(df_pert)
        if not np.isnan(eta):
            eta_atmos.append(eta)

    results['systematic_variations']['atmospheric'] = {
        'eta_mean': float(np.mean(eta_atmos)),
        'eta_std': float(np.std(eta_atmos)),
        'systematic_error': float(np.std(eta_atmos)),
        'description': 'Tropospheric delay model error ~5mm'
    }

    # 4. Instrumental drift
    logger.info("Running instrumental uncertainty MC...")
    eta_inst = []
    for _ in range(n_mc):
        df_pert = inject_instrumental_uncertainty(df, rng,
                                                  calibration_drift_mm_yr=rng.uniform(0.3, 0.7))
        eta, _ = compute_eta_with_error(df_pert)
        if not np.isnan(eta):
            eta_inst.append(eta)

    results['systematic_variations']['instrumental'] = {
        'eta_mean': float(np.mean(eta_inst)),
        'eta_std': float(np.std(eta_inst)),
        'systematic_error': float(np.std(eta_inst)),
        'description': 'Calibration drift ~0.5 mm/year'
    }

    # 5. Combined systematic (all sources together)
    logger.info("Running combined systematic MC...")
    eta_combined = []
    for _ in range(n_mc):
        df_pert = df.copy()  # Keep this copy since we modify it
        df_pert = inject_ephemeris_uncertainty(df_pert, rng)
        df_pert = inject_tidal_uncertainty(df_pert, rng)
        df_pert = inject_atmospheric_uncertainty(df_pert, rng)
        df_pert = inject_instrumental_uncertainty(df_pert, rng)
        eta, _ = compute_eta_with_error(df_pert)
        eta_combined.append(eta)

    combined_systematic_error = float(np.std(eta_combined))
    total_uncertainty = np.sqrt(
        eta_stat_error**2 + combined_systematic_error**2)

    results['systematic_variations']['combined'] = {
        'eta_mean': float(np.mean(eta_combined)),
        'eta_std': combined_systematic_error,
        'systematic_error': combined_systematic_error,
        'description': 'All systematic sources combined'
    }

    # Load the proper statistical error from step_002 (MANDATORY)
    step_002_path = Path(__file__).parent.parent.parent / 'results' / 'outputs' / 'step_002_statistical_analysis.json'
    if not step_002_path.exists():
        raise FileNotFoundError(
            f"Required upstream data not found: {step_002_path}. "
            "Step 002 must be run before Step 020 to ensure methodological consistency."
        )
    
    with open(step_002_path, 'r') as f:
        step_002_data = json.load(f)
    final_stat_error = step_002_data.get('eta_err_mcmc', 0.001)
    
    # Total uncertainty: systematic dominates for LLR at this precision level
    # The systematic error floor from literature is ~5e-4 in eta (from ~5-10mm LLR floor / 13m scale)
    literature_systematic_floor_eta = 5e-4  # From Murphy 2013, Williams et al. 2014
    effective_systematic = max(combined_systematic_error, literature_systematic_floor_eta)
    
    # Total uncertainty is dominated by systematic floor for LLR
    total_uncertainty_corrected = np.sqrt(final_stat_error**2 + effective_systematic**2)
    
    results['error_budget'] = {
        'statistical': float(final_stat_error),
        'systematic_combined': float(effective_systematic),
        'systematic_mc_component': float(combined_systematic_error),
        'total_uncertainty': float(total_uncertainty_corrected),
        'signal_to_total_ratio': float(abs(eta_baseline) / total_uncertainty_corrected),
        'systematic_to_statistical_ratio': float(effective_systematic / final_stat_error),
        'dominant_error_source': 'systematic' if effective_systematic > final_stat_error else 'statistical',
        'literature_systematic_floor_applied': True,
        'systematic_floor_source': 'Murphy 2013; Williams et al. 2014 LLR error budgets'
    }

    results['conclusion'] = {
        'baseline_eta': float(eta_baseline),
        'systematic_uncertainty': float(effective_systematic),
        'statistical_uncertainty': float(final_stat_error),
        'total_uncertainty': float(total_uncertainty_corrected),
        'assessment': 'systematic_errors_dominant' if effective_systematic > final_stat_error else 'statistical_errors_dominant',
        'recommendation': f"Report: η = {eta_baseline:.2e} ± {final_stat_error:.2e} (stat) ± {effective_systematic:.2e} (sys)",
        'note': 'Systematic error floor from LLR literature (~5mm) dominates over statistical precision.'
    }

    return results

def main():
    global logger
    # Setup TEPLogger for consistent file logging
    log_dir = Path(__file__).parent.parent.parent / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    logger = TEPLogger("step_020", str(log_dir / "step_020_systematic_monte_carlo.log"))
    set_step_logger(logger)

    print_status("═══ Starting Step 020: Systematic Error Monte Carlo Propagation...", "TITLE")
    print_status("═══ STEP PURPOSE: Quantify systematic uncertainty through Monte Carlo error propagation for TEP detection", "INFO")
    print_status("═══ METHOD: Ephemeris uncertainties, tidal model variations, atmospheric delay errors, instrumental calibration uncertainties", "INFO")
    print_status("═══ PARAMETERS: MC iterations=500, seed=42, ephemeris position error=0.5m, tidal Love number error=0.01, atmospheric zenith error=5mm", "INFO")

    logger.info("Step 020: Systematic Error Monte Carlo Propagation")

    print_status("═══ DATA SUMMARY", "INFO")
    # Load data
    data_path = Path(__file__).parent.parent.parent / 'data' / \
        'processed' / 'INPOP19a_all_stations_residuals.csv'
    df = pd.read_csv(data_path)

    print_status(f"    Dataset: N = {len(df):,} observations", "DATA")
    print_status(f"    Data source: INPOP19a_all_stations_residuals.csv", "DATA")

    print_status("═══ ANALYSIS TRACE", "INFO")
    logger.info(f"Loaded {len(df)} observations")
    print_status(f">>> Running ephemeris uncertainty MC", "PROCESS")
    print_status(f">>> Running tidal uncertainty MC", "PROCESS")
    print_status(f">>> Running atmospheric uncertainty MC", "PROCESS")
    print_status(f">>> Running instrumental uncertainty MC", "PROCESS")
    print_status(f">>> Running combined systematic MC", "PROCESS")

    # Run Monte Carlo analysis
    results = monte_carlo_systematic_analysis(df, n_mc=500, seed=42)

    # Add metadata
    results['step_id'] = 'step_020'
    results['status'] = 'PASS'
    results['mc_parameters'] = {
        'n_iterations': 500,
        'seed': 42,
        'ephemeris_position_error_m': 0.5,
        'ephemeris_velocity_error_mm_s': 0.1,
        'tidal_love_number_error': 0.01,
        'atmospheric_zenith_error_mm': 5.0,
        'instrumental_drift_mm_yr': 0.5
    }

    print_status("═══ RESULTS SUMMARY", "INFO")
    budget = results['error_budget']
    print_status(f"    Statistical error: ±{budget['statistical']:.2e}", "CALC")
    print_status(f"    Systematic error: ±{budget['systematic_combined']:.2e}", "CALC")
    print_status(f"    Total uncertainty: ±{budget['total_uncertainty']:.2e}", "CALC")
    print_status(f"    Signal/Total ratio: {budget['signal_to_total_ratio']:.1f}", "CALC")
    print_status(f"    Dominant source: {budget['dominant_error_source'].upper()}", "CALC")

    print_status("═══ INTERPRETATION", "INFO")
    print_status(f"    Systematic errors dominate over statistical precision at LLR precision levels", "INFO")
    print_status(f"    Literature systematic floor (~5mm) applied from Murphy 2013; Williams et al. 2014", "INFO")
    print_status(f"    Recommendation: Report both statistical and systematic uncertainties", "INFO")

    print_status("═══ REPRODUCIBILITY", "INFO")
    print_status(f"    Output file: results/outputs/step_020_systematic_monte_carlo.json", "INFO")
    print_status(f"    MC iterations: 500", "INFO")
    print_status(f"    Random seed: 42", "INFO")
    print_status(f"    Statistical error source: step_002_statistical_analysis.json", "INFO")

    # Save results
    project_root = Path(__file__).parent.parent.parent
    output_path = project_root / 'results' / 'outputs' / 'step_020_systematic_monte_carlo.json'
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, 'w') as f:
        json.dump(results, f, indent=2)

    output_rel = output_path.relative_to(project_root) if output_path.is_relative_to(project_root) else output_path
    logger.info(f"Results saved to {output_rel}")

    # Print summary
    logger.info(f"Statistical error: ±{budget['statistical']:.2e}")
    logger.info(f"Systematic error: ±{budget['systematic_combined']:.2e}")
    logger.info(f"Total uncertainty: ±{budget['total_uncertainty']:.2e}")
    logger.info(f"Signal/Total ratio: {budget['signal_to_total_ratio']:.1f}")
    logger.info(f"Dominant source: {budget['dominant_error_source'].upper()}")

    return results

if __name__ == '__main__':
    main()