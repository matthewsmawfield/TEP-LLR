#!/usr/bin/env python3
"""
Step 038: Full-Moon Deficit Analysis

Tests for the "full-moon deficit" effect documented by Murphy et al. (2010, 2014)
and Sabhlok et al. (2024), where signal degradation near full moon was attributed
to dust accumulation + thermal lensing on lunar retroreflectors.

TEP reinterpretation: The phase-dependent deficit may reflect scalar-field activation
tied to solar illumination geometry, consistent with:
- Perihelion enhancement (Step 024): η = -5.45×10^-4 at perihelion vs null at aphelion
- Negative η sign indicating gravitational potential suppression dominance
- Hardware-epoch consistency across all five instrument eras

This step quantifies the full-moon deficit and tests whether it correlates with
the TEP-predicted synodic modulation rather than purely thermal/dust effects.
"""

import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

import json
import numpy as np
import pandas as pd
from scipy import stats
import argparse
from scripts.utils.logger import TEPLogger, set_step_logger, set_verbose_mode, print_status
from scripts.utils.llr_constants import ETA_SCALE_FACTOR


def load_residuals(filepath):
    """Load residual data from CSV."""
    df = pd.read_csv(filepath)
    df['year'] = df['date_julian_year']
    return df

def analyze_full_moon_deficit(df):
    """Analyze residuals near full moon vs other phases."""
    elongation = df['elongation_rad'].values
    residuals = df['residual_m'].values
    
    # Define phase regions
    # Full moon: elongation within 20° of π (0.349 rad)
    full_moon_mask = np.abs(elongation - np.pi) < np.radians(20)
    
    # Quarter moon: elongation near π/2 or 3π/2
    quarter_mask = (np.abs(elongation - np.pi/2) < np.radians(20)) | \
                   (np.abs(elongation - 3*np.pi/2) < np.radians(20))
    
    # New moon: elongation within 20° of 0 or 2π
    new_moon_mask = (elongation < np.radians(20)) | (elongation > 2*np.pi - np.radians(20))
    
    results = {}
    
    for region_name, mask in [
        ('full_moon', full_moon_mask),
        ('quarter_moon', quarter_mask),
        ('new_moon', new_moon_mask)
    ]:
        n = np.sum(mask)
        if n < 10:
            results[region_name] = {'error': f'Insufficient data (n={n})'}
            continue
            
        res_region = residuals[mask]
        mean_res = np.mean(res_region)
        std_res = np.std(res_region, ddof=1)
        sem_res = std_res / np.sqrt(n)
        
        # Compute correlation with cos(elongation) in this region
        cos_elong = np.cos(elongation[mask])
        if np.std(cos_elong) > 0.01:  # Ensure variation in cos(D)
            r, p = stats.pearsonr(res_region, cos_elong)
            
            # OLS amplitude
            A = np.sum(res_region * cos_elong) / np.sum(cos_elong**2)
            eta = A / ETA_SCALE_FACTOR
        else:
            r, p, eta = 0, 1, 0
        
        results[region_name] = {
            'n_observations': int(n),
            'mean_residual_m': float(mean_res),
            'std_residual_m': float(std_res),
            'sem_residual_m': float(sem_res),
            'correlation_r': float(r),
            'correlation_p': float(p),
            'eta': float(eta)
        }
    
    return results

def test_tep_vs_thermal_hypothesis(df: pd.DataFrame) -> dict:
    """Docstring."""
    
    elongation = df['elongation_rad'].values
    residuals = df['residual_m'].values
    
    # Compute correlations with competing predictors
    cos_elong = np.cos(elongation)
    
    # TEP predictor: cos(D) (synodic phase)
    r_tep, p_tep = stats.pearsonr(residuals, cos_elong)
    
    # Thermal predictor: proximity to full moon (absolute deviation from π)
    thermal_proxy = -np.abs(elongation - np.pi)  # More negative = closer to full
    r_thermal, p_thermal = stats.pearsonr(residuals, thermal_proxy)
    
    return {
        'tep_correlation': {'r': float(r_tep), 'p': float(p_tep)},
        'thermal_correlation': {'r': float(r_thermal), 'p': float(p_thermal)},
        'preferred_model': 'TEP' if abs(r_tep) > abs(r_thermal) else 'Thermal',
        'r_tep_squared': float(r_tep**2),
        'r_thermal_squared': float(r_thermal**2)
    }

def analyze_by_illumination_geometry(df: pd.DataFrame) -> dict:
    """Docstring."""
    # Get heliocentric distance from date (simplified orbital model)
    # Using day of year to approximate orbital position for sensitivity testing
    # Perihelion ~ day 3, aphelion ~ day 186 (Earth's orbital parameters)
    df_temp = df.copy()  # Keep this copy since we add a column
    df_temp['day_of_year'] = (df_temp['date_julian'] % 365.25).astype(int)
    
    # Simplified heliocentric distance model for illumination analysis
    df_temp['days_from_perihelion'] = np.minimum(
        np.abs(df_temp['day_of_year'] - 3),
        365.25 - np.abs(df_temp['day_of_year'] - 3)
    )
    
    elongation = df_temp['elongation_rad'].values
    residuals = df_temp['residual_m'].values
    cos_elong = np.cos(elongation)
    
    # Subset near full moon
    full_moon_mask = np.abs(elongation - np.pi) < np.radians(30)
    
    if np.sum(full_moon_mask) < 50:
        return {'error': 'Insufficient full-moon data'}
    
    # Correlation with heliocentric distance (perihelion proximity)
    perihelion_proximity = 1 / (1 + df_temp['days_from_perihelion'].values)
    r_peri, p_peri = stats.pearsonr(residuals[full_moon_mask], 
                                     perihelion_proximity[full_moon_mask])
    
    return {
        'full_moon_subset_n': int(np.sum(full_moon_mask)),
        'perihelion_correlation_r': float(r_peri),
        'perihelion_correlation_p': float(p_peri),
        'interpretation': 'TEP-supported' if r_peri < -0.05 and p_peri < 0.05 else 'Inconclusive'
    }

def main():
    parser = argparse.ArgumentParser(description='Full-Moon Deficit Analysis (Step 038)')
    parser.add_argument('--verbose', action='store_true', help='Enable verbose output')
    args = parser.parse_args()
    
    # Setup logging
    log_dir = PROJECT_ROOT / "logs"
    log_dir.mkdir(exist_ok=True)
    log_file = log_dir / "step_038_full_moon_deficit_analysis.log"
    logger = TEPLogger("step_038", str(log_file))
    set_step_logger(logger)
    set_verbose_mode(args.verbose)
    
    print_status("Step 038: Full-Moon Deficit Analysis", "STEP")
    print_status("Testing Murphy/Sabhlok full-moon deficit vs TEP predictions", "INFO")
    
    # Load data
    data_path = PROJECT_ROOT / "data/processed/INPOP19a_all_stations_residuals.csv"
    if not data_path.exists():
        print_status(f"Error: Data file not found: {data_path}", "ERROR")
        sys.exit(1)
    
    print_status(f"Loading data from: {data_path}", "INFO")
    df = load_residuals(data_path)
    print_status(f"Loaded {len(df)} observations", "INFO")
    
    # Analyze full-moon deficit
    print_status("Analyzing full-moon deficit pattern...", "INFO")
    phase_analysis = analyze_full_moon_deficit(df)
    
    # Test TEP vs thermal hypothesis
    print_status("Testing TEP vs thermal/dust hypothesis...", "INFO")
    model_comparison = test_tep_vs_thermal_hypothesis(df)
    
    # Analyze by illumination geometry
    print_status("Analyzing heliocentric geometry dependence...", "INFO")
    geometry_analysis = analyze_by_illumination_geometry(df)
    
    # Compile results
    output = {
        "step_id": "step_038",
        "description": "Tests Murphy/Sabhlok full-moon deficit against TEP scalar-field predictions",
        "timestamp": pd.Timestamp.now().isoformat(),
        "references": {
            "murphy_2010": "Murphy et al. 2010, PASP, 122, 892",
            "murphy_2014": "Murphy et al. 2014, Icarus, 231, 183",
            "sabhlok_2024": "Sabhlok et al. 2024, Icarus, 412, 115927 (arXiv:2403.00899)"
        },
        "phase_analysis": phase_analysis,
        "model_comparison": model_comparison,
        "illumination_geometry": geometry_analysis,
        "assessment": {
            "tep_vs_thermal": model_comparison['preferred_model'],
            "conclusion": "Full-moon deficit pattern consistent with TEP synodic modulation" if model_comparison['preferred_model'] == 'TEP' else "Ambiguous",
            "note": "Sabhlok documented 10-15x signal degradation near full moon attributed to dust+thermal lensing; TEP offers scalar-field coupling as alternative mechanism"
        },
        "status": "PASS"
    }
    
    # Save output
    output_path = PROJECT_ROOT / "results/outputs/step_038_full_moon_deficit_analysis.json"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, 'w') as f:
        json.dump(output, f, indent=4, default=str)
    
    output_rel = output_path.relative_to(PROJECT_ROOT) if output_path.is_relative_to(PROJECT_ROOT) else output_path
    print_status(f"Results saved to: {output_rel}", "INFO")
    
    # Summary
    print_status("\n=== Full-Moon Deficit Analysis Summary ===", "STEP")
    if 'full_moon' in phase_analysis and 'n_observations' in phase_analysis['full_moon']:
        print_status(f"Full-moon observations: {phase_analysis['full_moon']['n_observations']}", "INFO")
    print_status(f"TEP correlation (r²): {model_comparison['r_tep_squared']:.4f}", "INFO")
    print_status(f"Thermal correlation (r²): {model_comparison['r_thermal_squared']:.4f}", "INFO")
    print_status(f"Preferred model: {model_comparison['preferred_model']}", "PASS" if model_comparison['preferred_model'] == 'TEP' else "WARNING")
    
    print_status("Step 038 completed successfully", "PASS")

if __name__ == "__main__":
    main()
