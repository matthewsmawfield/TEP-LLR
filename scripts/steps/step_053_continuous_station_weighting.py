#!/usr/bin/env python3
"""
Flyby TEP Pipeline - Step 053: Continuous Station Weighting

Replaces binary "powered/underpowered" classification with continuous station weights.

Weight formula:
w_s = (σ_r,s^2 N_s / σ_cosD,s^2) × Q_phase,s × Q_hardware,s

Global η as weighted average:
η_global = Σ_s w_s η_s / Σ_s w_s
"""

import sys
from pathlib import Path
import numpy as np
import pandas as pd

# Add project root to path
PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from scripts.utils.logger import TEPLogger, set_step_logger, print_status
from scripts.utils.config import get_config
from scripts.utils.llr_constants import ETA_SCALE_FACTOR

TEP_CONFIG = get_config()

def calculate_station_weight(df_station):
    """
    Calculate continuous station weight.
    
    Args:
        df_station: DataFrame for a single station
    
    Returns:
        dict with weight components and total weight
    """
    N = len(df_station)
    sigma_r = df_station['residual_m'].std()
    cos_D = np.cos(df_station['elongation_rad'].values)
    sigma_cosD = cos_D.std()
    
    # Phase coverage quality factor (based on uniformity)
    # Simple metric: ratio of actual std to theoretical max for uniform distribution
    # Standard deviation of cos(U) for U ~ Uniform(-π, π) is exactly 1/√2 ≈ 0.707
    # Derivation: std(cos(U)) = sqrt(Var(cos(U))) = sqrt(0.5) = 1/√2
    phase_uniformity = sigma_cosD / 0.707  # 0.707 = 1/√2 is exact std for uniform phase
    Q_phase = min(phase_uniformity, 1.0)
    
    # Hardware quality factor (based on year range and era)
    years = df_station['year'].values if 'year' in df_station.columns else pd.to_datetime(df_station['jd'] - 2440587.5, unit='D').dt.year.values
    modern_fraction = np.mean(years >= 2009)  # Fraction of C-SPAD era
    Q_hardware = 0.5 + 0.5 * modern_fraction  # 0.5 for legacy, 1.0 for modern
    
    # Calculate weight
    # CORRECTED: Higher weight for stations with lower noise (sigma_r), more observations (N),
    # better phase coverage (sigma_cosD), and better hardware (Q_hardware)
    weight = (N * sigma_cosD**2 / sigma_r**2) * Q_phase * Q_hardware
    
    # Estimate η for this station
    # polyfit returns amplitude in meters, need to convert to eta by dividing by ETA_SCALE_FACTOR
    eta_station = np.polyfit(cos_D, df_station['residual_m'].values, 1)[0] / ETA_SCALE_FACTOR
    
    return {
        'N': int(N),
        'sigma_r': float(sigma_r),
        'sigma_cosD': float(sigma_cosD),
        'Q_phase': float(Q_phase),
        'Q_hardware': float(Q_hardware),
        'weight': float(weight),
        'eta_station': float(eta_station)
    }

def run_continuous_weighting_analysis(verbose=False):
    """Run continuous station weighting analysis."""
    input_path = PROJECT_ROOT / "data" / "processed" / "INPOP19a_all_stations_residuals.csv"
    
    if not input_path.exists():
        print_status(f"CRITICAL DATA FAILURE: {input_path} not found. Cannot proceed.", "ERROR")
        return None
    
    df = pd.read_csv(input_path)
    
    # Add year column if not present
    if 'year' not in df.columns:
        df['year'] = df['date_julian_year'].values
    
    print_status("Running continuous station weighting analysis...", "INFO")
    
    # Calculate weights for each station
    stations = df['station'].unique()
    station_results = {}
    
    total_weight = 0.0
    weighted_eta_sum = 0.0
    
    for station in stations:
        df_station = df[df['station'] == station]
        result = calculate_station_weight(df_station)
        station_results[station] = result
        
        total_weight += result['weight']
        weighted_eta_sum += result['weight'] * result['eta_station']
        
        print_status(f"  {station}:", "INFO")
        print_status(f"    N = {result['N']:,}, weight = {result['weight']:.6e}", "INFO")
        print_status(f"    η = {result['eta_station']:.8e}", "INFO")
    
    # Calculate global weighted η
    eta_global = weighted_eta_sum / total_weight if total_weight > 0 else 0.0
    
    # Calculate weighted standard error
    weighted_variance_sum = 0.0
    for station in stations:
        result = station_results[station]
        eta_station = result['eta_station']
        weighted_variance_sum += result['weight'] * (eta_station - eta_global)**2
    
    eta_err_global = np.sqrt(weighted_variance_sum / total_weight) if total_weight > 0 else 0.0
    snr_global = abs(eta_global) / eta_err_global if eta_err_global > 0 else 0.0
    
    print_status(f"  Global weighted η = {eta_global:.8e} ± {eta_err_global:.8e}", "INFO")
    print_status(f"  SNR = {snr_global:.2f}σ", "INFO")
    
    results = {
        'eta_global': float(eta_global),
        'eta_err_global': float(eta_err_global),
        'snr_global': float(snr_global),
        'station_weights': station_results,
        'total_weight': float(total_weight)
    }
    
    return results

if __name__ == "__main__":
    log_dir = PROJECT_ROOT / "logs"
    logger = TEPLogger("step_053", str(log_dir / "step_053_continuous_station_weighting.log"))
    set_step_logger(logger)
    
    results = run_continuous_weighting_analysis(verbose=True)
    
    if results:
        logger.save_step_results(results, PROJECT_ROOT, "step_053_continuous_station_weighting")
        print_status(f"Continuous Station Weighting Complete. η_global = {results['eta_global']:.8e}", "SUCCESS")
    else:
        print_status("Continuous Station Weighting Failed.", "ERROR")
