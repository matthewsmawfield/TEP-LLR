#!/usr/bin/env python3
"""
Step 029: Day/Night Thermal Bias Diagnostics
Maps the TEP-LLR purely to diurnal day/night cycles to test for unmodeled atmospheric/telescope thermal drift.
"""


import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

import pandas as pd
import numpy as np
from astropy.time import Time
from astropy.coordinates import EarthLocation, get_sun, AltAz
import astropy.units as u
import statsmodels.api as sm
import json
from scripts.utils.logger import TEPLogger, set_step_logger, set_verbose_mode, print_status

# Add project root to path

# Grasse is the station that dominates the dataset
STATION_COORDS = {
    'APO': EarthLocation(lat=32.7802*u.deg, lon=-105.8202*u.deg, height=2788*u.m),
    'Grasse': EarthLocation(lat=43.7538*u.deg, lon=6.9227*u.deg, height=1270*u.m),
    'Haleakala': EarthLocation(lat=20.7083*u.deg, lon=-156.2562*u.deg, height=3058*u.m),
    'McDonald2': EarthLocation(lat=30.6714*u.deg, lon=-104.0226*u.deg, height=2070*u.m),
    'Matera': EarthLocation(lat=40.6488*u.deg, lon=16.7047*u.deg, height=493*u.m)
}

def main():
    # Setup TEPLogger
    log_dir = PROJECT_ROOT / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    logger = TEPLogger("step_027", str(log_dir / "step_027_day_night_thermal_bias.log"))
    set_step_logger(logger)
    
    # Load data
    DATA_PATH = PROJECT_ROOT / "data" / "processed" / "INPOP19a_all_stations_residuals.csv"
    
    if not DATA_PATH.exists():
        logger.error("Data file not found!")
        return {"status": "FAIL", "reason": "Data file not found"}
    
    df = pd.read_csv(DATA_PATH)

    logger.info("Computing Solar Altitude for 26,000+ observations... (Vectorized for M4 Pro)")

    # Vectorized approach: process per-station in batches instead of row-by-row
    solar_alt_map = {}
    for station_name, loc in STATION_COORDS.items():
        mask = df['station'] == station_name
        if mask.sum() == 0:
            continue
        
        # Get all JDs for this station at once
        jds = df.loc[mask, 'date_julian'].values
        
        # Vectorized astropy: create Time with array of JDs
        t = Time(jds, format='jd')
        
        try:
            # Vectorized sun calculation
            sun = get_sun(t)
            # Vectorized AltAz transform
            altaz = sun.transform_to(AltAz(obstime=t, location=loc))
            # Store results
            solar_alt_map[station_name] = altaz.alt.degree
        except Exception as e:
            print_status(f"WARNING: Solar altitude computation failed for {station_name}: {e}. Dropping {len(jds)} observations.", "WARNING")
            solar_alt_map[station_name] = np.full(len(jds), np.nan)

    # Map back to dataframe
    solar_alts = np.full(len(df), np.nan)
    for station_name in STATION_COORDS.keys():
        mask = df['station'] == station_name
        if station_name in solar_alt_map:
            solar_alts[mask] = solar_alt_map[station_name]

    df['solar_alt'] = solar_alts
    df = df.dropna(subset=['solar_alt'])  # PERFORMANCE FIX: Removed unnecessary .copy()
    df['is_day'] = (df['solar_alt'] > 0).astype(int)

    # Normalize cos(elongation)
    df['cos_elong'] = np.cos(df['elongation_rad'])

    results = []

    for station in df['station'].unique():
        d_sta = df[df['station'] == station]
        if len(d_sta) < 100:
            continue
        
        mean_day = d_sta[d_sta['is_day'] == 1]['residual_m'].mean()
        mean_night = d_sta[d_sta['is_day'] == 0]['residual_m'].mean()
        day_night_diff = mean_day - mean_night
        
        # Simple regression of residual ~ day_night
        X = d_sta['is_day']
        X = sm.add_constant(X)
        model = sm.OLS(d_sta['residual_m'], X).fit()
        dn_coeff = model.params.get('is_day', 0)
        dn_p = model.pvalues.get('is_day', 1.0)
        
        # Standard regression: residual ~ cos(elongation)
        X_cos = d_sta['cos_elong']
        X_cos = sm.add_constant(X_cos)
        model_cos = sm.OLS(d_sta['residual_m'], X_cos).fit()
        eta_calc = model_cos.params.get('cos_elong', 0) / 13.0
        eta_p = model_cos.pvalues.get('cos_elong', 1.0)
        
        # Partial regression: residual ~ cos(elongation) + is_day + solar_alt
        X_partial = sm.add_constant(d_sta[['cos_elong', 'solar_alt', 'is_day']])
        model_part = sm.OLS(d_sta['residual_m'], X_partial).fit()
        partial_eta = model_part.params.get('cos_elong', 0) / 13.0
        partial_eta_p = model_part.pvalues.get('cos_elong', 1.0)

        # Calculate what day/night bias does to the apparent Eta
        # Because New Moon (which predicts negative eta) is strictly day-ranged,
        # mapping Day-Night directly to a spurious eta:
        spurious_eta = day_night_diff / 13.0
        
        results.append({
            'Station': station,
            'N': len(d_sta),
            'Day_Ranged_Percent': (d_sta['is_day'].sum() / len(d_sta)) * 100,
            'Day_Minus_Night_Bias_mm': day_night_diff * 1000,
            'Spurious_Eta': spurious_eta,
            'Original_Eta': eta_calc,
            'Cleaned_Eta_After_Solar_Modeling': partial_eta,
            'Cleaned_P_Val': partial_eta_p
        })

    res_df = pd.DataFrame(results)
    logger.info("DAY / NIGHT THERMAL FALSE POSITIVE INVESTIGATION")
    logger.info("\n" + res_df.to_string(index=False))
    
    # Global test
    logger.info("GLOBAL ANALYSIS:")
    X_global = sm.add_constant(df[['cos_elong', 'is_day', 'solar_alt']])
    model_global = sm.OLS(df['residual_m'], X_global).fit()
    
    original_eta = sm.OLS(df['residual_m'], sm.add_constant(df['cos_elong'])).fit().params['cos_elong'] / 13.0
    logger.info(f"Original Global Eta: {original_eta}")
    logger.info(f"Solar Altitude Parameter: {model_global.params['solar_alt']} p-val: {model_global.pvalues['solar_alt']}")
    logger.info(f"Cleaned Global Eta (Controlling for Sun Altitude/Thermal): {model_global.params['cos_elong']/13.0}")
    logger.info(f"Cleaned Global P-Value for Cos(Elongation): {model_global.pvalues['cos_elong']}")
    
    # Save results
    output_data = {
        "step_id": "step_027",
        "status": "PASS",
        "station_results": results,
        "global_analysis": {
            "original_eta": float(original_eta),
            "solar_altitude_param": float(model_global.params['solar_alt']),
            "solar_altitude_pval": float(model_global.pvalues['solar_alt']),
            "cleaned_eta": float(model_global.params['cos_elong']/13.0),
            "cleaned_pval": float(model_global.pvalues['cos_elong'])
        }
    }
    
    logger.save_step_results(output_data, PROJECT_ROOT, "step_027_day_night_thermal_bias")
    
    return output_data

if __name__ == "__main__":
    main()