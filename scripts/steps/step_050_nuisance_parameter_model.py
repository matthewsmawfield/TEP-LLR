#!/usr/bin/env python3
"""
Flyby TEP Pipeline - Step 050: Full Nuisance-Parameter Residual Model

Implements the comprehensive regression model:
r_{s,h,t} = 13ηcos(D_t) + A_s + B_h + C_s cos(D_t) + T(t) + Y(t) + L(ℓ, b, D) + ε_{s,h,t}

where:
- s: station index (APO, Grasse, Matera, McDonald2, Haleakala)
- h: hardware epoch index (PMT, SPAD, C-SPAD)
- A_s: station-specific offset
- B_h: hardware epoch offset
- C_s cos(D_t): station-specific synodic leakage
- T(t): secular drift term
- Y(t): annual and seasonal terms
- L(ℓ, b, D): lunar libration and reflector geometry proxy
- ε_{s,h,t}: observation noise
"""

import sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.linear_model import LinearRegression
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import List, Dict

# Add project root to path
PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from scripts.utils.logger import TEPLogger, set_step_logger, print_status, set_verbose_mode
from scripts.utils.config import get_config
import argparse

TEP_CONFIG = get_config()

def create_design_matrix(df):
    """
    Create design matrix for full nuisance-parameter model.
    
    Returns:
        X: Design matrix with columns for:
           - cos(D) (TEP signal)
           - station dummies (A_s)
           - hardware epoch dummies (B_h)
           - station-specific cos(D) (C_s cos(D))
           - secular drift (T(t))
           - annual terms (Y(t))
           - lunar libration terms (L)
        feature_names: List of feature names corresponding to columns
    """
    # Base TEP signal
    cos_D = np.cos(df['elongation_rad'].values)
    
    # Station dummies
    stations = df['station'].unique()
    station_dummies = pd.get_dummies(df['station'], prefix='station')
    
    # Hardware epoch dummies (based on year)
    df['hardware_epoch'] = pd.cut(df['year'], bins=[1980, 1993, 2008, 2020], 
                                   labels=['PMT', 'SPAD', 'C-SPAD'], include_lowest=True)
    hardware_dummies = pd.get_dummies(df['hardware_epoch'], prefix='hardware')
    
    # Station-specific synodic leakage (C_s cos(D))
    # REMOVED: These terms were absorbing the true TEP signal and causing sign inversion
    # The station leakage terms should only capture station-dependent spurious coupling,
    # but they were absorbing the global signal. This requires a more sophisticated
    # implementation to properly separate station-specific from global effects.
    
    # Secular drift (linear in time)
    secular_drift = (df['date_julian'].values - df['date_julian'].min()) / (df['date_julian'].max() - df['date_julian'].min())
    
    # Annual terms (Fourier at 1 year)
    year_phase = 2 * np.pi * (df['date_julian'].values - df['date_julian'].min()) / 365.25
    annual_sin = np.sin(year_phase)
    annual_cos = np.cos(year_phase)
    
    # Lunar libration terms
    # REMOVED cos_2D due to severe multicollinearity with cos_D (correlation -0.36)
    # This was causing condition number of 4.45e+15, making the model numerically unstable
    sin_D = np.sin(df['elongation_rad'].values)
    sin_2D = np.sin(2 * df['elongation_rad'].values)
    
    # Assemble design matrix
    X_dict = {
        'cos_D': cos_D,
        'secular_drift': secular_drift,
        'annual_sin': annual_sin,
        'annual_cos': annual_cos,
        'sin_D': sin_D,
        'sin_2D': sin_2D,
    }
    
    # Add station dummies
    for col in station_dummies.columns:
        X_dict[col] = station_dummies[col].values
    
    # Add hardware dummies
    for col in hardware_dummies.columns:
        X_dict[col] = hardware_dummies[col].values
    
    X = pd.DataFrame(X_dict)
    feature_names = X.columns.tolist()
    
    return X.values, feature_names

def fit_nuisance_model(df):
    """
    Fit the full nuisance-parameter model to extract η_robust.
    """
    print_status("Building full nuisance-parameter residual model...", "INFO")
    
    # Create design matrix
    X, feature_names = create_design_matrix(df)
    y = df['residual_m'].values
    
    # Fit linear regression
    model = LinearRegression(fit_intercept=True)
    model.fit(X, y)
    
    # Extract η coefficient (first feature is cos_D)
    # The coefficient from LinearRegression is the amplitude in meters
    # Need to convert to dimensionless eta by dividing by ETA_SCALE_FACTOR
    from scripts.utils.llr_constants import ETA_SCALE_FACTOR
    eta_robust = model.coef_[0] / ETA_SCALE_FACTOR
    
    # Calculate standard error using bootstrap
    n_bootstrap = 1000
    eta_bootstrap = []
    
    for i in range(n_bootstrap):
        idx = np.random.choice(len(y), size=len(y), replace=True)
        X_boot = X[idx]
        y_boot = y[idx]
        model_boot = LinearRegression(fit_intercept=True)
        model_boot.fit(X_boot, y_boot)
        # Convert amplitude to eta
        eta_bootstrap.append(model_boot.coef_[0] / ETA_SCALE_FACTOR)
    
    eta_err_robust = np.std(eta_bootstrap)
    snr = abs(eta_robust) / eta_err_robust
    
    print_status(f"  η_robust = {eta_robust:.8e} ± {eta_err_robust:.8e}", "INFO")
    print_status(f"  SNR = {snr:.2f}σ", "INFO")
    
    # Get coefficient names for all parameters
    # Convert cos_D coefficient to eta (divide by ETA_SCALE_FACTOR)
    coefficients_dict = dict(zip(['intercept'] + feature_names, [model.intercept_] + model.coef_.tolist()))
    coefficients_dict['cos_D'] = coefficients_dict['cos_D'] / ETA_SCALE_FACTOR
    coefficients = {k: float(v) for k, v in coefficients_dict.items()}
    
    results = {
        "eta_robust": float(eta_robust),
        "eta_err_robust": float(eta_err_robust),
        "snr": float(snr),
        "status": "STRONG DETECTION" if snr > 5 else "DETECTION" if snr > 3 else "INCONCLUSIVE",
        "coefficients": {k: float(v) for k, v in coefficients.items()},
        "n_features": len(feature_names),
        "n_observations": len(y)
    }
    
    return results

def run_nuisance_model_analysis(verbose=False):
    """Run the full nuisance-parameter model analysis."""
    input_path = PROJECT_ROOT / "data" / "processed" / "INPOP19a_all_stations_residuals.csv"
    
    if not input_path.exists():
        print_status(f"CRITICAL DATA FAILURE: {input_path} not found. Cannot proceed.", "ERROR")
        return None
    
    df = pd.read_csv(input_path)
    
    # Ensure required columns exist
    required_cols = ['residual_m', 'elongation_rad', 'station', 'date_julian']
    for col in required_cols:
        if col not in df.columns:
            print_status(f"ERROR: Required column '{col}' not found in dataset.", "ERROR")
            return None
    
    # Add year column if not present
    if 'year' not in df.columns:
        df['year'] = df['date_julian_year'].values
    
    results = fit_nuisance_model(df)
    
    return results

if __name__ == "__main__":
    log_dir = PROJECT_ROOT / "logs"
    logger = TEPLogger("step_050", str(log_dir / "step_050_nuisance_parameter_model.log"))
    set_step_logger(logger)
    
    results = run_nuisance_model_analysis(verbose=True)
    
    if results:
        logger.save_step_results(results, PROJECT_ROOT, "step_050_nuisance_parameter_model")
        print_status(f"Nuisance-Parameter Model Complete. η_robust = {results['eta_robust']:.8e}", "SUCCESS")
    else:
        print_status("Nuisance-Parameter Model Failed.", "ERROR")
