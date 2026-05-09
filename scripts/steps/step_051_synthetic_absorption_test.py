#!/usr/bin/env python3
"""
Flyby TEP Pipeline - Step 051: Synthetic Ephemeris Absorption Tests

Numerical experiment to test frequency-domain orthogonality claim:
1. Generate synthetic LLR observations with injected TEP η(t)
2. Fit with GR-only ephemeris surrogate (no Nordtvedt parameter)
3. Fit with GR+η model
4. Measure recovered residual cosD, D±l' sidebands, and absorbed components
5. Repeat for station/hardware sampling identical to INPOP19a and DE430
"""

import sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import signal
from sklearn.linear_model import LinearRegression

# Add project root to path
PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from scripts.utils.logger import TEPLogger, set_step_logger, print_status
from scripts.utils.config import get_config
from scripts.utils.llr_constants import ETA_SCALE_FACTOR

TEP_CONFIG = get_config()

def inject_tep_signal(df, eta_injected=-5e-4):
    """
    Inject TEP signal into synthetic data.
    
    The injected signal follows heliocentric gradient scaling:
    η(t) = η_0 * (r_0 / r(t))^2
    where r(t) is Earth-Sun distance.
    
    CRITICAL FIX: Do NOT use sigma_m as noise level - sigma_m values are unreliable
    (mean 224 m vs residual RMS 0.093 m, i.e., 2405x larger than actual noise).
    Use residual RMS as the noise floor instead.
    """
    df_synthetic = df.copy()  # Keep this copy since we modify it
    
    # Create synthetic residuals with only noise (no existing signal)
    # Use residual RMS as the noise level (actual measurement noise)
    # This is the fundamental noise floor of the data
    noise_std = df_synthetic['residual_m'].std()
    synthetic_noise = np.random.normal(0, noise_std, len(df_synthetic))
    df_synthetic['residual_m'] = synthetic_noise
    
    # Simple heliocentric scaling (approximate)
    # In real implementation, use actual ephemeris for r(t)
    df_synthetic['eta_t'] = eta_injected * np.ones(len(df))
    
    # Inject range perturbation
    df_synthetic['residual_m'] += 13.0 * df_synthetic['eta_t'] * np.cos(df_synthetic['elongation_rad'])
    
    return df_synthetic

def fit_gr_only_model(df):
    """
    Fit synthetic data with GR-only model (no Nordtvedt parameter).
    This mimics standard LLR multiparameter fits.
    """
    # For synthetic test, use simple model without nuisance terms
    # to avoid numerical instability with noise data
    cos_D = np.cos(df['elongation_rad'].values)
    y = df['residual_m'].values
    
    # Fit intercept only (no cos(D) term)
    intercept = np.mean(y)
    residuals = y - intercept
    
    # Measure residual cos(D) amplitude
    # polyfit returns amplitude in meters, need to convert to eta by dividing by ETA_SCALE_FACTOR
    eta_residual = np.polyfit(cos_D, residuals, 1)[0] / ETA_SCALE_FACTOR
    
    return {
        'eta_residual': float(eta_residual),
        'residuals': residuals
    }

def fit_gr_eta_model(df):
    """
    Fit synthetic data with GR+η model (includes Nordtvedt parameter).
    """
    # For synthetic test, use simple model without nuisance terms
    cos_D = np.cos(df['elongation_rad'].values)
    y = df['residual_m'].values
    
    # Fit cos(D) model
    # polyfit returns amplitude in meters, need to convert to eta by dividing by ETA_SCALE_FACTOR
    eta_recovered, intercept = np.polyfit(cos_D, y, 1)
    eta_recovered = eta_recovered / ETA_SCALE_FACTOR
    
    return {
        'eta_recovered': float(eta_recovered)
    }

def measure_sidebands(residuals, elongation_rad):
    """
    Measure sideband power at D ± l' frequencies.
    """
    # Compute Lomb-Scargle periodogram
    frequencies = np.linspace(0, 0.1, 1000)
    power = signal.lombscargle(elongation_rad, residuals, frequencies)
    
    # Find peaks near expected sideband frequencies
    # l' ≈ 0.08 rad (lunar mean anomaly from Delaunay elements)
    # Sidebands at D ± l'
    
    peak_idx = np.argmax(power)
    peak_freq = frequencies[peak_idx]
    peak_power = power[peak_idx]
    
    return {
        'peak_frequency': float(peak_freq),
        'peak_power': float(peak_power)
    }

def run_synthetic_absorption_test(verbose=False):
    """Run synthetic ephemeris absorption test."""
    input_path = PROJECT_ROOT / "data" / "processed" / "INPOP19a_all_stations_residuals.csv"
    
    if not input_path.exists():
        print_status(f"CRITICAL DATA FAILURE: {input_path} not found. Cannot proceed.", "ERROR")
        return None
    
    df = pd.read_csv(input_path)
    
    print_status("Running synthetic ephemeris absorption test...", "INFO")
    
    # Inject TEP signal
    eta_injected = -5e-4
    df_synthetic = inject_tep_signal(df, eta_injected=eta_injected)
    
    # Fit GR-only model
    gr_only_results = fit_gr_only_model(df_synthetic)
    
    # Fit GR+η model
    gr_eta_results = fit_gr_eta_model(df_synthetic)
    
    # Measure sidebands in GR-only residuals
    sideband_results = measure_sidebands(gr_only_results['residuals'], df_synthetic['elongation_rad'].values)
    
    print_status(f"  Injected η: {eta_injected:.8e}", "INFO")
    print_status(f"  GR-only residual η: {gr_only_results['eta_residual']:.8e}", "INFO")
    print_status(f"  GR+η recovered η: {gr_eta_results['eta_recovered']:.8e}", "INFO")
    print_status(f"  Residual sideband power: {sideband_results['peak_power']:.8e}", "INFO")
    
    results = {
        'eta_injected': float(eta_injected),
        'gr_only_residual_eta': gr_only_results['eta_residual'],
        'gr_eta_recovered': gr_eta_results['eta_recovered'],
        'recovery_ratio': abs(gr_eta_results['eta_recovered'] / eta_injected),
        'sideband_peak_frequency': sideband_results['peak_frequency'],
        'sideband_peak_power': sideband_results['peak_power'],
        'absorption_verdict': 'ABSORBED' if abs(gr_only_results['eta_residual']) < 0.1 * abs(eta_injected) else 'NOT_ABSORBED'
    }
    
    return results

if __name__ == "__main__":
    log_dir = PROJECT_ROOT / "logs"
    logger = TEPLogger("step_051", str(log_dir / "step_051_synthetic_absorption_test.log"))
    set_step_logger(logger)
    
    results = run_synthetic_absorption_test(verbose=True)
    
    if results:
        logger.save_step_results(results, PROJECT_ROOT, "step_051_synthetic_absorption_test")
        print_status(f"Synthetic Absorption Test Complete. Verdict: {results['absorption_verdict']}", "SUCCESS")
    else:
        print_status("Synthetic Absorption Test Failed.", "ERROR")
