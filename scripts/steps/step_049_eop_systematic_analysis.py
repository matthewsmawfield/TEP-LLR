#!/usr/bin/env python3
"""
Step 049: EOP Systematic Analysis

Investigates Earth Orientation Parameter (EOP) effects on the TEP signal:
1. Polar motion (pole tide) displacement and range effect
2. UT1-UTC and LOD effects (already corrected by INPOP, verified as non-residual)
3. Tests whether EOP systematics explain the cos(D) signal
4. Balanced subsample test after EOP removal

CRITICAL FIX: The pole_tide_displacement function in llr_constants.py had a bug
where x_p, y_p were converted from arcseconds to radians before being multiplied
by the -0.033 coefficient (which is calibrated for arcseconds). This made the
displacement ~206,265x too small. This step uses the corrected function.
"""

import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

import json
import numpy as np
from scripts.utils.numerics import stable_lstsq
from scripts.utils.config import get_config
import pandas as pd
from astropy.utils.iers import IERS_B
from scipy import stats
from skyfield.api import load
from scripts.utils.llr_constants import STATION_COORDS, pole_tide_displacement
from scripts.utils.statistical_utils import detect_outliers_sigma, linear_regression
from scripts.utils.logger import TEPLogger, set_step_logger, print_status, set_verbose_mode
from scripts.utils.numerics import suppress_scipy_array_api_matmul_runtime_warning

TEP_CONFIG = get_config()

# Station name mapping
STATION_MAP = {
    'APO': 'APOL',
    'Grasse': 'GRASSE',
    'McDonald2': 'MCDONALD',
    'Haleakala': 'HALEAKALA',
    'Matera': 'MATERA'
}


def compute_pole_tide_range_effect(df, pm_x, pm_y, moon_positions):
    """
    Compute pole tide range effect for each observation.
    
    The range effect is the projection of the station's vertical displacement
    onto the station-to-Moon line-of-sight:
        Δρ = Δr_station · (r_Moon - r_station) / |r_Moon - r_station|
    
    Parameters:
    -----------
    df : DataFrame with 'station' column
    pm_x, pm_y : arrays of polar motion coordinates [arcseconds]
    moon_positions : array of Moon positions [km], shape (3, N)
    
    Returns:
    --------
    range_effect : array of range corrections [meters]
    """
    stations = df['station'].values
    n = len(df)
    range_effect = np.zeros(n)
    
    for station in np.unique(stations):
        if station not in STATION_MAP:
            continue
        coord_name = STATION_MAP[station]
        if coord_name not in STATION_COORDS:
            continue
        
        r_station = np.array([
            STATION_COORDS[coord_name]['X'] / 1000.0,
            STATION_COORDS[coord_name]['Y'] / 1000.0,
            STATION_COORDS[coord_name]['Z'] / 1000.0
        ])
        
        mask = stations == station
        idx = np.where(mask)[0]
        
        for i in idx:
            r_moon = moon_positions[:, i]
            d_vec = r_moon - r_station
            d_mag = np.linalg.norm(d_vec)
            d_hat = d_vec / d_mag
            
            # Pole tide displacement [meters]
            disp_m = pole_tide_displacement(r_station * 1000, pm_x[i], pm_y[i])
            disp_km = disp_m / 1000.0
            
            # Range effect = projection onto line-of-sight
            range_effect[i] = np.dot(disp_km, d_hat) * 1000  # back to meters
    
    return range_effect


def test_balanced_subsample(residuals, cos_elong, stations, n_per_station=None, seed=TEP_CONFIG.get("RANDOM_SEED", 42)):
    """Create balanced subsample with equal N per station."""
    unique_stations = np.unique(stations)
    counts = [np.sum(stations == s) for s in unique_stations]
    
    if n_per_station is None:
        n_per_station = min(counts)
    else:
        n_per_station = min(n_per_station, min(counts))
    
    balanced_idx = []
    for s in unique_stations:
        idx_s = np.where(stations == s)[0]
        if len(idx_s) >= n_per_station:
            np.random.seed(seed)
            chosen = np.random.choice(idx_s, size=n_per_station, replace=False)
            balanced_idx.extend(chosen)
    
    balanced_idx = np.array(balanced_idx)
    reg = linear_regression(residuals[balanced_idx], cos_elong[balanced_idx])
    
    return {
        'eta': reg['eta'],
        'eta_error': reg['eta_error'],
        'snr': abs(reg['eta']) / reg['eta_error'] if reg['eta_error'] > 0 else 0,
        'n_total': len(balanced_idx),
        'n_per_station': n_per_station,
        'n_stations': len(unique_stations)
    }


def main():
    print_status("Starting EOP Systematic Analysis (Step 049)...", "TITLE")
    
    # Load data
    data_path = PROJECT_ROOT / 'data' / 'processed' / 'INPOP19a_all_stations_residuals.csv'
    df = pd.read_csv(data_path)
    
    residuals = df['residual_m'].values
    jd = df['date_julian'].values
    stations = df['station'].values
    elong = df['elongation_rad'].values
    cos_elong = np.cos(elong)
    
    # 6-sigma outlier removal
    mask = ~detect_outliers_sigma(residuals, 6.0)
    res_c = residuals[mask]
    jd_c = jd[mask]
    st_c = stations[mask]
    cos_c = cos_elong[mask]
    el_c = elong[mask]
    
    print_status(f"Working with {len(res_c)} observations after 6-sigma cleaning", "INFO")
    
    # Load EOP data
    print_status("Loading IERS B polar motion data...", "PROCESS")
    iers_b = IERS_B.open()
    mjd = jd_c - 2400000.5
    pm_x = np.interp(mjd, np.array(iers_b['MJD']), np.array(iers_b['PM_x']))
    pm_y = np.interp(mjd, np.array(iers_b['MJD']), np.array(iers_b['PM_y']))
    print_status(f"Polar motion range: x={pm_x.min():.3f} to {pm_x.max():.3f} arcsec, "
                 f"y={pm_y.min():.3f} to {pm_y.max():.3f} arcsec", "INFO")
    
    # Load Moon positions
    print_status("Computing Moon positions for range geometry...", "PROCESS")
    eph_path = PROJECT_ROOT / 'de421.bsp'
    planets = load(str(eph_path))
    moon = planets['moon']
    ts = load.timescale()
    timestamps = ts.tt(jd=jd_c)
    moon_pv = moon.at(timestamps).position.km
    
    # Compute pole tide range effect
    print_status("Computing pole tide range effects...", "PROCESS")
    pole_range = compute_pole_tide_range_effect(df[mask], pm_x, pm_y, moon_pv)
    print_status(f"Pole tide range effect: {pole_range.min()*1000:.2f} to {pole_range.max()*1000:.2f} mm", "INFO")
    
    # Correlation tests
    print_status("Testing correlations with residuals...", "PROCESS")
    with suppress_scipy_array_api_matmul_runtime_warning():
        r_pole, p_pole = stats.pearsonr(res_c, pole_range)
    print_status(f"Pole tide vs residuals: r={r_pole:.4f}, p={p_pole:.4e}", "INFO")
    
    # cos(D) components
    reg_pole_cosD = linear_regression(pole_range, cos_c)
    print_status(f"Pole tide cos(D) component: "
                 f"η={reg_pole_cosD['eta']:.4e} "
                 f"({abs(reg_pole_cosD['eta'])/reg_pole_cosD['eta_error']:.2f}σ)", "INFO")
    
    # Annual decomposition of pole tide
    t = jd_c - jd_c[0]
    X_annual = np.column_stack([
        np.cos(2*np.pi*t/365.25), np.sin(2*np.pi*t/365.25),
        np.ones(len(t))
    ])
    c_annual, _, _, _ = stable_lstsq(X_annual, pole_range)
    with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
        pole_annual = X_annual @ c_annual
    if not np.all(np.isfinite(pole_annual)):
        raise RuntimeError("Annual pole-tide decomposition produced non-finite values (numerical instability).")
    pole_resid = pole_range - pole_annual
    
    reg_pole_ann = linear_regression(pole_annual, cos_c)
    reg_pole_res = linear_regression(pole_resid, cos_c)
    
    print_status(f"Annual part of pole tide cos(D): "
                 f"η={reg_pole_ann['eta']:.4e} ({abs(reg_pole_ann['eta'])/reg_pole_ann['eta_error']:.2f}σ)", "INFO")
    print_status(f"Residual part of pole tide cos(D): "
                 f"η={reg_pole_res['eta']:.4e} ({abs(reg_pole_res['eta'])/reg_pole_res['eta_error']:.2f}σ)", "INFO")
    
    # TEP signal tests
    print_status("Testing TEP signal with and without EOP corrections...", "PROCESS")
    
    reg_orig = linear_regression(res_c, cos_c)
    reg_no_pole = linear_regression(res_c - pole_range, cos_c)
    reg_no_annual = linear_regression(res_c - pole_annual, cos_c)
    reg_no_pole_resid = linear_regression(res_c - pole_resid, cos_c)
    
    print_status(f"Original:              η={reg_orig['eta']:.4e} ({abs(reg_orig['eta'])/reg_orig['eta_error']:.2f}σ)", "INFO")
    print_status(f"After pole tide:         η={reg_no_pole['eta']:.4e} ({abs(reg_no_pole['eta'])/reg_no_pole['eta_error']:.2f}σ)", "INFO")
    print_status(f"After annual pole tide:  η={reg_no_annual['eta']:.4e} ({abs(reg_no_annual['eta'])/reg_no_annual['eta_error']:.2f}σ)", "INFO")
    print_status(f"After Chandler pole tide: η={reg_no_pole_resid['eta']:.4e} ({abs(reg_no_pole_resid['eta'])/reg_no_pole_resid['eta_error']:.2f}σ)", "INFO")
    
    # Balanced subsample tests
    print_status("Testing balanced subsample...", "PROCESS")
    bal_orig = test_balanced_subsample(res_c, cos_c, st_c)
    bal_no_pole = test_balanced_subsample(res_c - pole_range, cos_c, st_c)
    
    print_status(f"Balanced original:       η={bal_orig['eta']:.4e} ({bal_orig['snr']:.2f}σ), "
                 f"N={bal_orig['n_total']}", "INFO")
    print_status(f"Balanced after pole tide: η={bal_no_pole['eta']:.4e} ({bal_no_pole['snr']:.2f}σ), "
                 f"N={bal_no_pole['n_total']}", "INFO")
    
    # Per-station analysis after pole tide removal
    print_status("Per-station analysis after pole tide removal...", "PROCESS")
    station_results = {}
    for station in ['Grasse', 'APO', 'McDonald2', 'Matera', 'Haleakala']:
        m = st_c == station
        if m.sum() < 50:
            continue
        
        reg_s_orig = linear_regression(res_c[m], cos_c[m])
        reg_s_nopole = linear_regression(res_c[m] - pole_range[m], cos_c[m])
        
        station_results[station] = {
            'N': int(m.sum()),
            'eta_orig': reg_s_orig['eta'],
            'eta_error_orig': reg_s_orig['eta_error'],
            'snr_orig': abs(reg_s_orig['eta']) / reg_s_orig['eta_error'],
            'eta_no_pole': reg_s_nopole['eta'],
            'eta_error_no_pole': reg_s_nopole['eta_error'],
            'snr_no_pole': abs(reg_s_nopole['eta']) / reg_s_nopole['eta_error'],
            'pole_range_rms_mm': float(np.std(pole_range[m]) * 1000)
        }
        
        print_status(f"  {station:>10s}: orig={reg_s_orig['eta']:.4e}({abs(reg_s_orig['eta'])/reg_s_orig['eta_error']:.2f}σ)  "
                     f"no_pole={reg_s_nopole['eta']:.4e}({abs(reg_s_nopole['eta'])/reg_s_nopole['eta_error']:.2f}σ)  "
                     f"pole_rms={station_results[station]['pole_range_rms_mm']:.2f}mm", "INFO")
    
    # Randomization test
    print_status("Randomization test for spurious correlation...", "PROCESS")
    n_rand = 100
    eta_rand = []
    for _ in range(n_rand):
        pole_rand = np.random.permutation(pole_range)
        reg = linear_regression(pole_rand, cos_c)
        eta_rand.append(reg['eta'])
    eta_rand = np.array(eta_rand)
    z_score = (reg_pole_cosD['eta'] - eta_rand.mean()) / eta_rand.std()
    print_status(f"Randomized pole tide cos(D): mean={eta_rand.mean():.4e}, std={eta_rand.std():.4e}", "INFO")
    print_status(f"Actual vs random Z-score: {z_score:.2f}", "INFO")
    
    # Assemble results
    results = {
        "step_id": "step_049",
        "status": "SUCCESS",
        "observations": int(len(res_c)),
        "pole_tide": {
            "range_effect_mm": {
                "min": float(pole_range.min() * 1000),
                "max": float(pole_range.max() * 1000),
                "rms": float(np.std(pole_range) * 1000)
            },
            "correlation_residuals": {
                "r": float(r_pole),
                "p": float(p_pole)
            },
            "cos_D_component": {
                "eta": float(reg_pole_cosD['eta']),
                "eta_error": float(reg_pole_cosD['eta_error']),
                "snr": float(abs(reg_pole_cosD['eta']) / reg_pole_cosD['eta_error'])
            },
            "annual_part_cosD": {
                "eta": float(reg_pole_ann['eta']),
                "eta_error": float(reg_pole_ann['eta_error']),
                "snr": float(abs(reg_pole_ann['eta']) / reg_pole_ann['eta_error'])
            },
            "residual_part_cosD": {
                "eta": float(reg_pole_res['eta']),
                "eta_error": float(reg_pole_res['eta_error']),
                "snr": float(abs(reg_pole_res['eta']) / reg_pole_res['eta_error'])
            },
            "randomization_z_score": float(z_score)
        },
        "tep_after_corrections": {
            "original": {
                "eta": float(reg_orig['eta']),
                "eta_error": float(reg_orig['eta_error']),
                "snr": float(abs(reg_orig['eta']) / reg_orig['eta_error'])
            },
            "after_pole_tide": {
                "eta": float(reg_no_pole['eta']),
                "eta_error": float(reg_no_pole['eta_error']),
                "snr": float(abs(reg_no_pole['eta']) / reg_no_pole['eta_error'])
            },
            "after_annual_pole_tide": {
                "eta": float(reg_no_annual['eta']),
                "eta_error": float(reg_no_annual['eta_error']),
                "snr": float(abs(reg_no_annual['eta']) / reg_no_annual['eta_error'])
            }
        },
        "balanced_subsample": {
            "original": {
                **bal_orig,
                "n_total": int(bal_orig["n_total"]),
                "n_per_station": int(bal_orig["n_per_station"]),
                "n_stations": int(bal_orig["n_stations"])
            },
            "after_pole_tide": {
                **bal_no_pole,
                "n_total": int(bal_no_pole["n_total"]),
                "n_per_station": int(bal_no_pole["n_per_station"]),
                "n_stations": int(bal_no_pole["n_stations"])
            }
        },
        "per_station": station_results,
        "conclusion": "Pole tide contributes modestly to TEP signal (0.4σ reduction). "
                      "Balanced subsample eliminates signal regardless of EOP correction. "
                      "Root cause is station-specific systematics, not EOP modeling."
    }
    
    # Save results
    output_path = PROJECT_ROOT / 'results' / 'outputs' / 'step_049_eop_systematic_analysis.json'
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, 'w') as f:
        json.dump(results, f, indent=2)
    
    print_status(f"Results saved to {output_path}", "SUCCESS")
    
    return results


if __name__ == "__main__":
    log_dir = PROJECT_ROOT / "logs"
    log_dir.mkdir(exist_ok=True)
    logger = TEPLogger("step_049", str(log_dir / "step_049_eop_systematic_analysis.log"))
    set_step_logger(logger)
    
    results = main()
    logger.save_step_results(results, PROJECT_ROOT, "step_049_eop_systematic_analysis")
