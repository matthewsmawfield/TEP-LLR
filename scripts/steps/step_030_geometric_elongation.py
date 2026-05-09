#!/usr/bin/env python3
"""
Step 030: True Geometric Elongation vs Mean Synodic Phase Fault Nullifier
Tests whether the TEP signal correlates strictly with mathematical periodic mean phase, or the TRUE geometric elongation defined by J2000 planetary vectors.
"""


import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

import pandas as pd
import numpy as np
from astropy.time import Time
from astropy.coordinates import get_sun, get_body, EarthLocation
import astropy.units as u
import statsmodels.api as sm
import json
from scripts.utils.logger import TEPLogger, set_step_logger, set_verbose_mode, print_status

# Add project root to path

def main():
    # Setup TEPLogger
    log_dir = PROJECT_ROOT / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    logger = TEPLogger("step_030", str(log_dir / "step_030_geometric_elongation.log"))
    set_step_logger(logger)
    
    # Load data
    DATA_PATH = PROJECT_ROOT / "data" / "processed" / "INPOP19a_all_stations_residuals.csv"
    
    if not DATA_PATH.exists():
        logger.error("Data file not found!")
        return {"status": "FAIL", "reason": "Data file not found"}
    
    df = pd.read_csv(DATA_PATH)
    
    logger.info("Computing True Geometric Vectors from Astropy (J2000 Ephemeris)... (Vectorized for M4 Pro)")

    # Vectorized approach: process all times at once instead of row-by-row
    jds = df['date_julian'].values
    t = Time(jds, format='jd')

    try:
        # Vectorized sun and moon calculations
        sun = get_sun(t)
        moon = get_body('moon', t)
        
        # Vectorized angular separation
        sep = sun.separation(moon).radian
        true_elongations = sep
    except Exception as e:
        true_elongations = np.full(len(df), np.nan)

    df['true_elong_rad'] = true_elongations
    df = df.dropna(subset=['true_elong_rad'])  # PERFORMANCE FIX: Removed unnecessary .copy()

    df['cos_true_elong'] = np.cos(df['true_elong_rad'])
    df['cos_mean_elong'] = np.cos(df['elongation_rad'])

    # OLS
    X_mean = sm.add_constant(df['cos_mean_elong'])
    model_mean = sm.OLS(df['residual_m'], X_mean).fit()
    eta_mean = model_mean.params['cos_mean_elong'] / 13.0
    p_mean = model_mean.pvalues['cos_mean_elong']

    X_true = sm.add_constant(df['cos_true_elong'])
    model_true = sm.OLS(df['residual_m'], X_true).fit()
    eta_true = model_true.params['cos_true_elong'] / 13.0
    p_true = model_true.pvalues['cos_true_elong']

    logger.info("="*80)
    logger.info("GEOMETRIC VS MEAN PHASE VULNERABILITY INVESTIGATION")
    logger.info("="*80)
    logger.info(f"Mean Approximation Eta:  {eta_mean:.8f} (p-val: {p_mean:.4e})")
    logger.info(f"True Geometric Eta:      {eta_true:.8f} (p-val: {p_true:.4e})")

    # Does the true geometric predict it better?
    logger.info("\nPartial Regression (Competing predictors):")
    X_both = sm.add_constant(df[['cos_mean_elong', 'cos_true_elong']])
    model_both = sm.OLS(df['residual_m'], X_both).fit()
    logger.info(f"   Coefficient for Mean: {model_both.params['cos_mean_elong']}")
    logger.info(f"   Coefficient for True: {model_both.params['cos_true_elong']}")
    logger.info(f"   P-val for Mean: {model_both.pvalues['cos_mean_elong']}")
    logger.info(f"   P-val for True: {model_both.pvalues['cos_true_elong']}")

    if model_both.pvalues['cos_true_elong'] < model_both.pvalues['cos_mean_elong']:
        conclusion = "The TRUE GEOMETRIC boundary drives the signal, validating the physical claim and rejecting the mathematical proxy vulnerability!"
    else:
        conclusion = "The signal purely correlates to mathematical mean time instead of the actual physical geometry."
    logger.info(f"\nCONCLUSION: {conclusion}")
    
    # Save results
    output_data = {
        "step_id": "step_030",
        "status": "PASS",
        "mean_approximation": {
            "eta": float(eta_mean),
            "p_value": float(p_mean)
        },
        "true_geometric": {
            "eta": float(eta_true),
            "p_value": float(p_true)
        },
        "partial_regression": {
            "coefficient_mean": float(model_both.params['cos_mean_elong']),
            "coefficient_true": float(model_both.params['cos_true_elong']),
            "p_value_mean": float(model_both.pvalues['cos_mean_elong']),
            "p_value_true": float(model_both.pvalues['cos_true_elong'])
        },
        "conclusion": conclusion
    }
    
    logger.save_step_results(output_data, PROJECT_ROOT, "step_030_geometric_elongation")
    
    return output_data

if __name__ == "__main__":
    main()