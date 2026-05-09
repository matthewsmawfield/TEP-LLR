#!/usr/bin/env python3
"""
Flyby TEP Pipeline - Step 054: Reflector-Level Controls and Lunar Libration

WARNING: This script contains SIMULATED/FABRICATED DATA and should NOT be used
for production analysis or publication results.

The actual INPOP19a dataset does not include reflector identifiers. This script
simulates reflector assignments with random probabilities as a placeholder for
future implementation when actual reflector-level data becomes available.

DO NOT USE RESULTS FROM THIS SCRIPT FOR SCIENTIFIC CONCLUSIONS.

Analyzes LLR residuals by reflector array (Apollo 11, 14, 15, Lunokhod 1, 2) to account for
reflector geometry and lunar libration effects.

Includes frequency controls: sin(D), 2D, annual, anomalistic, draconic terms.
"""

import sys
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.linear_model import LinearRegression

# Add project root to path
PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from scripts.utils.logger import TEPLogger, set_step_logger, print_status
from scripts.utils.config import get_config
from scripts.utils.llr_constants import ETA_SCALE_FACTOR

TEP_CONFIG = get_config()

def analyze_reflector(df_reflector):
    """
    Analyze residuals for a single reflector.
    """
    N = len(df_reflector)
    sigma_r = df_reflector['residual_m'].std()
    
    # Phase coverage
    elongation_rad = df_reflector['elongation_rad'].values
    phase_min = elongation_rad.min()
    phase_max = elongation_rad.max()
    phase_coverage = f"({phase_min/np.pi:.1f}π-{phase_max/np.pi:.1f}π)"
    
    # Estimate η
    # polyfit returns amplitude in meters, need to convert to eta by dividing by ETA_SCALE_FACTOR
    cos_D = np.cos(elongation_rad)
    eta_reflector = np.polyfit(cos_D, df_reflector['residual_m'].values, 1)[0] / ETA_SCALE_FACTOR
    
    # Estimate libration coupling (correlation with elongation)
    libration_coupling = np.abs(np.corrcoef(elongation_rad, df_reflector['residual_m'].values)[0, 1])
    
    # Categorize libration coupling
    if libration_coupling > 0.3:
        coupling_category = "High"
    elif libration_coupling > 0.15:
        coupling_category = "Moderate"
    else:
        coupling_category = "Low"
    
    return {
        'N': int(N),
        'phase_coverage': phase_coverage,
        'eta': float(eta_reflector),
        'RMS_cm': float(sigma_r * 100),
        'libration_coupling': float(libration_coupling),
        'coupling_category': coupling_category
    }

def fit_frequency_controls(df):
    """
    Fit model with frequency controls: sin(D), 2D, annual, anomalistic, draconic.
    """
    elongation_rad = df['elongation_rad'].values
    jd = df['date_julian'].values
    
    # Frequency controls
    sin_D = np.sin(elongation_rad)
    cos_2D = np.cos(2 * elongation_rad)
    sin_2D = np.sin(2 * elongation_rad)
    
    # Annual terms (1 year = 365.25 days)
    year_phase = 2 * np.pi * (jd - jd.min()) / 365.25
    annual_sin = np.sin(year_phase)
    annual_cos = np.cos(year_phase)
    
    # Anomalistic terms (anomalistic year ≈ 365.26 days)
    anomalistic_phase = 2 * np.pi * (jd - jd.min()) / 365.26
    anomalistic_sin = np.sin(anomalistic_phase)
    anomalistic_cos = np.cos(anomalistic_phase)
    
    # Draconic terms (draconic month ≈ 27.21 days)
    draconic_phase = 2 * np.pi * (jd - jd.min()) / 27.21
    draconic_sin = np.sin(draconic_phase)
    draconic_cos = np.cos(draconic_phase)
    
    # Assemble design matrix
    X_dict = {
        'cos_D': np.cos(elongation_rad),
        'sin_D': sin_D,
        'cos_2D': cos_2D,
        'sin_2D': sin_2D,
        'annual_sin': annual_sin,
        'annual_cos': annual_cos,
        'anomalistic_sin': anomalistic_sin,
        'anomalistic_cos': anomalistic_cos,
        'draconic_sin': draconic_sin,
        'draconic_cos': draconic_cos,
    }
    
    X = pd.DataFrame(X_dict)
    y = df['residual_m'].values
    
    model = LinearRegression(fit_intercept=True)
    model.fit(X, y)
    
    # Extract η coefficient (first feature is cos_D)
    # The coefficient from LinearRegression is the amplitude in meters
    # Need to convert to dimensionless eta by dividing by ETA_SCALE_FACTOR
    eta_with_controls = model.coef_[0] / ETA_SCALE_FACTOR
    
    # Get coefficient names
    coefficients = dict(zip(X.columns.tolist(), model.coef_.tolist()))
    
    return {
        'eta_with_controls': float(eta_with_controls),
        'coefficients': {k: float(v) for k, v in coefficients.items()},
        'intercept': float(model.intercept_)
    }

def run_reflector_analysis(verbose=False):
    """Run reflector-level analysis."""
    print_status("CRITICAL WARNING: This step uses simulated reflector data.", "ERROR")
    print_status("Actual INPOP19a dataset does not include reflector identifiers.", "ERROR")
    print_status("This step is DISABLED to prevent fabricated data from entering results.", "ERROR")
    return None
    
    input_path = PROJECT_ROOT / "data" / "processed" / "INPOP19a_all_stations_residuals.csv"
    
    if not input_path.exists():
        print_status(f"CRITICAL DATA FAILURE: {input_path} not found. Cannot proceed.", "ERROR")
        return None
    
    df = pd.read_csv(input_path)
    
    print_status("Running reflector-level analysis...", "INFO")
    
    # Note: The actual dataset may not have reflector information
    # For this implementation, we'll simulate reflector assignment based on station
    # In a real implementation, this would use actual reflector data
    
    # Simulate reflector assignment (this is a placeholder)
    # In reality, each observation would have a reflector identifier
    np.random.seed(42)
    reflectors = ['Apollo 11', 'Apollo 14', 'Apollo 15', 'Lunokhod 1', 'Lunokhod 2']
    df['reflector'] = np.random.choice(reflectors, size=len(df), p=[0.2, 0.3, 0.3, 0.1, 0.1])
    
    # Analyze each reflector
    reflector_results = {}
    
    for reflector in reflectors:
        df_reflector = df[df['reflector'] == reflector]
        if len(df_reflector) > 0:
            result = analyze_reflector(df_reflector)
            reflector_results[reflector] = result
            
            print_status(f"  {reflector}:", "INFO")
            print_status(f"    N = {result['N']:,}, η = {result['eta']:.8e}", "INFO")
            print_status(f"    RMS = {result['RMS_cm']:.1f} cm, coupling = {result['coupling_category']}", "INFO")
    
    # Fit frequency controls on full dataset
    frequency_results = fit_frequency_controls(df)
    
    print_status(f"  η with frequency controls: {frequency_results['eta_with_controls']:.8e}", "INFO")
    
    results = {
        'reflector_analysis': reflector_results,
        'frequency_controls': frequency_results,
        'note': 'Reflector assignment simulated - actual reflector data required for production'
    }
    
    return results

if __name__ == "__main__":
    log_dir = PROJECT_ROOT / "logs"
    logger = TEPLogger("step_054", str(log_dir / "step_054_reflector_level_analysis.log"))
    set_step_logger(logger)
    
    results = run_reflector_analysis(verbose=True)
    
    if results:
        logger.save_step_results(results, PROJECT_ROOT, "step_054_reflector_level_analysis")
        print_status("Reflector-Level Analysis Complete.", "SUCCESS")
    else:
        print_status("Reflector-Level Analysis Failed.", "ERROR")
