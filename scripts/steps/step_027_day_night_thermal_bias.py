#!/usr/bin/env python3
"""
Step 027: Day/Night Thermal Bias Diagnostics
Maps the TEP-LLR purely to diurnal day/night cycles to test for unmodeled atmospheric/telescope thermal drift.
"""


import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

import pandas as pd
import numpy as np
from scripts.utils.numerics import stable_lstsq
from astropy.time import Time
from astropy.coordinates import EarthLocation, get_sun, AltAz
import astropy.units as u
from scipy import stats
import json
from scripts.utils.logger import TEPLogger, set_step_logger, set_verbose_mode, print_status
from scripts.utils.statistical_utils import linear_regression
from scripts.utils.llr_constants import ETA_SCALE_FACTOR

# Add project root to path

# Grasse is the station that dominates the dataset
STATION_COORDS = {
    'APO': EarthLocation(lat=32.7802*u.deg, lon=-105.8202*u.deg, height=2788*u.m),
    'Grasse': EarthLocation(lat=43.7538*u.deg, lon=6.9227*u.deg, height=1270*u.m),
    'Haleakala': EarthLocation(lat=20.7083*u.deg, lon=-156.2562*u.deg, height=3058*u.m),
    'McDonald2': EarthLocation(lat=30.6714*u.deg, lon=-104.0226*u.deg, height=2070*u.m),
    'Matera': EarthLocation(lat=40.6488*u.deg, lon=16.7047*u.deg, height=493*u.m)
}

def _multiple_regression(y, X, param_names):
    """Fit y = X @ beta using np.linalg.lstsq with proper diagnostics.

    Returns dict with coefficients, standard errors, and p-values,
    matching the diagnostics available in linear_regression.
    """
    n, k = X.shape
    beta, residuals, rank, _ = stable_lstsq(X, y)
    if rank < k:
        return {name: {'coeff': np.nan, 'se': np.nan, 'p': 1.0} for name in param_names}
    with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
        resid = y - X @ beta
    mse = np.sum(resid**2) / (n - k)
    XtX_inv = np.linalg.pinv(X.T @ X, rcond=1e-10, hermitian=True)
    cov = mse * XtX_inv
    se = np.sqrt(np.diag(cov))
    t_stats = beta / se
    pvals = 2 * (1 - stats.t.cdf(np.abs(t_stats), n - k))
    return {
        name: {'coeff': float(beta[i]), 'se': float(se[i]), 'p': float(pvals[i])}
        for i, name in enumerate(param_names)
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
        reg_dn = linear_regression(d_sta['residual_m'].values, d_sta['is_day'].values)
        dn_coeff = reg_dn['amplitude']
        dn_p = 2 * (1 - stats.t.cdf(abs(dn_coeff) / reg_dn['amplitude_error'], len(d_sta) - 2)) if reg_dn['amplitude_error'] > 0 else 1.0

        # Standard regression: residual ~ cos(elongation)
        reg_cos = linear_regression(d_sta['residual_m'].values, d_sta['cos_elong'].values)
        eta_calc = reg_cos['eta']
        eta_p = 2 * (1 - stats.t.cdf(abs(reg_cos['amplitude']) / reg_cos['amplitude_error'], len(d_sta) - 2)) if reg_cos['amplitude_error'] > 0 else 1.0

        # Partial regression: residual ~ cos(elongation) + is_day + solar_alt
        X_partial = np.column_stack([
            np.ones(len(d_sta)),
            d_sta['cos_elong'].values,
            d_sta['solar_alt'].values,
            d_sta['is_day'].values
        ])
        part_res = _multiple_regression(d_sta['residual_m'].values, X_partial,
                                        ['intercept', 'cos_elong', 'solar_alt', 'is_day'])
        partial_eta = part_res['cos_elong']['coeff'] / ETA_SCALE_FACTOR
        partial_eta_p = part_res['cos_elong']['p']

        # Calculate what day/night bias does to the apparent Eta
        # Because New Moon (which predicts negative eta) is strictly day-ranged,
        # mapping Day-Night directly to a spurious eta:
        spurious_eta = day_night_diff / ETA_SCALE_FACTOR

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
    X_global = np.column_stack([
        np.ones(len(df)),
        df['cos_elong'].values,
        df['solar_alt'].values,
        df['is_day'].values
    ])
    global_res = _multiple_regression(df['residual_m'].values, X_global,
                                      ['intercept', 'cos_elong', 'solar_alt', 'is_day'])

    reg_orig_global = linear_regression(df['residual_m'].values, df['cos_elong'].values)
    original_eta = reg_orig_global['eta']
    logger.info(f"Original Global Eta: {original_eta}")
    logger.info(f"Solar Altitude Parameter: {global_res['solar_alt']['coeff']} p-val: {global_res['solar_alt']['p']}")
    logger.info(f"Cleaned Global Eta (Controlling for Sun Altitude/Thermal): {global_res['cos_elong']['coeff']/ETA_SCALE_FACTOR}")
    logger.info(f"Cleaned Global P-Value for Cos(Elongation): {global_res['cos_elong']['p']}")

    # Save results
    output_data = {
        "step_id": "step_027",
        "status": "PASS",
        "station_results": results,
        "global_analysis": {
            "original_eta": float(original_eta),
            "solar_altitude_param": float(global_res['solar_alt']['coeff']),
            "solar_altitude_pval": float(global_res['solar_alt']['p']),
            "cleaned_eta": float(global_res['cos_elong']['coeff']/ETA_SCALE_FACTOR),
            "cleaned_pval": float(global_res['cos_elong']['p'])
        }
    }

    logger.save_step_results(output_data, PROJECT_ROOT, "step_027_day_night_thermal_bias")

    return output_data

if __name__ == "__main__":
    main()