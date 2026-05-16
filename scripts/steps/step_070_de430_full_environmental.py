#!/usr/bin/env python3
"""
Step 070: DE430 Full Environmental Model (No Outlier Removal)
=================================================================

Step 006b shows DE430 requires 37 outlier removals (0.8%) to achieve
significance. The manuscript frames these as "gross outliers clustering
asymmetrically at specific phases." However, DE430 uses a different
dynamical model than INPOP19a. If DE430 partially absorbs the static
Nordtvedt carrier but is structurally blind to TEP sidebands (as argued
in Steps 032, 041, 065, 066), the unabsorbed sideband power appears as
phase-clustered outliers. Removing them discards TEP signal.

This step tests that hypothesis by running the FULL TEP environmental
model on RAW DE430 residuals (no outlier removal):

  residual = η·cos(D) + η_r·(1/r_⊙) + η_v·v_r + η_θ·cos(θ_EM-CMB)
              + cos(2D) + monthly + annual + const

If the "outliers" are predicted by the sideband structure (e.g., they
cluster at phases where D ± l' interference is strongest), they are
signal, not noise. This transforms DE430 from a fragile corroboration
into an independent detection channel.
"""

import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

import numpy as np
import pandas as pd
from scipy import stats
from skyfield.api import load

from scripts.utils.logger import TEPLogger, set_step_logger, print_status
from scripts.utils.astronomical_utils import load_skyfield_planets
from scripts.utils.llr_constants import ETA_SCALE_FACTOR
from scripts.utils.statistical_utils import robust_regression, cluster_robust_variance, detect_outliers_sigma
from scripts.utils.numerics import suppress_scipy_array_api_matmul_runtime_warning

log_dir = PROJECT_ROOT / "logs"
log_dir.mkdir(parents=True, exist_ok=True)
logger = TEPLogger("step_070", str(log_dir / "step_070_de430_full_environmental.log"))
set_step_logger(logger)

_CMB_RA = np.deg2rad(168.14)
_CMB_DEC = np.deg2rad(-7.22)
_CMB_UNIT = np.array([
    np.cos(_CMB_DEC) * np.cos(_CMB_RA),
    np.cos(_CMB_DEC) * np.sin(_CMB_RA),
    np.sin(_CMB_DEC),
])
_CMB_UNIT = _CMB_UNIT / np.linalg.norm(_CMB_UNIT)


def compute_de430_environmental(jd_array):
    """Compute environmental variables for DE430 epochs."""
    planets, _eph_path = load_skyfield_planets(PROJECT_ROOT)
    earth = planets["earth"]
    sun = planets["sun"]
    moon = planets["moon"]
    ts = load.timescale()
    timestamps = ts.tt(jd=jd_array)

    astrometric = earth.at(timestamps).observe(sun)
    r_sun_au = astrometric.distance().au
    inv_r = 1.0 / r_sun_au

    # Radial velocity
    pos = earth.at(timestamps).position.km
    sun_pos = sun.at(timestamps).position.km
    rel_pos = pos - sun_pos
    r_norm = np.linalg.norm(rel_pos, axis=0)
    r_hat = rel_pos / r_norm
    v_earth = earth.at(timestamps).velocity.km_per_s
    v_sun = sun.at(timestamps).velocity.km_per_s
    v_rel = v_earth - v_sun
    v_r = np.sum(v_rel * r_hat, axis=0)

    # CMB orientation
    moon_pos = moon.at(timestamps).position.km
    em_vec = moon_pos - pos
    em_norm = np.linalg.norm(em_vec, axis=0)
    em_hat = em_vec / em_norm
    cos_theta_cmb = np.sum(em_hat * _CMB_UNIT[:, None], axis=0)

    return {
        'inv_r': inv_r,
        'v_r': v_r,
        'cos_theta_cmb': cos_theta_cmb,
        'r_sun_au': r_sun_au,
    }


def main():
    print_status("═══ Step 070: DE430 Full Environmental Model (Raw) ═══", "TITLE")

    de430_path = PROJECT_ROOT / "data" / "processed" / "DE430_all_residuals.csv"
    if not de430_path.exists():
        raise FileNotFoundError(f"DE430 residuals not found: {de430_path}")

    df = pd.read_csv(de430_path)
    df = df.sort_values(['date_julian', 'station'], kind='mergesort').reset_index(drop=True)

    res = df['residual_m'].values
    st = df['station'].values if 'station' in df.columns else np.full(len(df), 'DE430')
    el = df['elongation_rad'].values
    jd = df['date_julian'].values
    n_total = len(df)

    print_status(f"DE430 raw: N={n_total:,}", "DATA")

    # Environmental variables
    print_status("Computing environmental projections...", "INFO")
    env = compute_de430_environmental(jd)
    inv_r = env['inv_r']
    v_r = env['v_r']
    cos_cmb = env['cos_theta_cmb']

    # Physics terms
    cos_el = np.cos(el)
    cos2d = np.cos(2 * el)
    year = jd / 365.25
    sin_y = np.sin(2 * np.pi * year)
    cos_y = np.cos(2 * np.pi * year)
    month = jd / 27.32
    sin_m = np.sin(2 * np.pi * month)
    cos_m = np.cos(2 * np.pi * month)

    # -----------------------------------------------------------------
    # Model 1: cosD only (for direct comparison with step_006b raw)
    # -----------------------------------------------------------------
    print_status("--- Model 1: cosD only (raw, no outlier removal) ---", "INFO")
    X1 = np.column_stack([cos_el, np.ones(len(cos_el))])
    fit1 = robust_regression(res, X1, scale_errors_by_birge=False)
    eta1 = fit1['coefficients'][0] / ETA_SCALE_FACTOR
    se1 = fit1['errors'][0] / ETA_SCALE_FACTOR
    snr1 = abs(eta1) / max(se1, 1e-20)
    print_status(f"  η = {eta1:.4e} ± {se1:.4e} ({snr1:.2f}σ)", "RESULT")

    # -----------------------------------------------------------------
    # Model 2: cosD + environmental terms
    # -----------------------------------------------------------------
    print_status("--- Model 2: cosD + 1/r_⊙ + v_r + cos(θ_CMB) ---", "INFO")
    X2 = np.column_stack([cos_el, inv_r, v_r, cos_cmb, np.ones(len(cos_el))])
    names2 = ['cosD', 'inv_r', 'v_r', 'cos_cmb', 'const']
    fit2 = robust_regression(res, X2, scale_errors_by_birge=False)
    eta2 = fit2['coefficients'][0] / ETA_SCALE_FACTOR
    se2 = fit2['errors'][0] / ETA_SCALE_FACTOR
    snr2 = abs(eta2) / max(se2, 1e-20)
    print_status(f"  η = {eta2:.4e} ± {se2:.4e} ({snr2:.2f}σ)", "RESULT")
    for i, nm in enumerate(names2):
        if nm == 'cosD':
            continue
        val = fit2['coefficients'][i]
        err = fit2['errors'][i]
        t = val / max(err, 1e-20)
        p = 2 * (1 - stats.t.cdf(abs(t), len(res) - len(names2)))
        print_status(f"    {nm:10s}: coeff={val:+.3e} ± {err:.3e}  t={t:.2f}, p={p:.4f}", "RESULT")

    # -----------------------------------------------------------------
    # Model 3: Full systematic + environmental
    # -----------------------------------------------------------------
    print_status("--- Model 3: full systematic + environmental ---", "INFO")
    X3 = np.column_stack([
        cos_el, cos2d, sin_m, cos_m, sin_y, cos_y,
        inv_r, v_r, cos_cmb,
        np.ones(len(cos_el))
    ])
    names3 = ['cosD', 'cos2D', 'sin_m', 'cos_m', 'sin_y', 'cos_y',
              'inv_r', 'v_r', 'cos_cmb', 'const']
    fit3 = robust_regression(res, X3, scale_errors_by_birge=False)
    eta3 = fit3['coefficients'][0] / ETA_SCALE_FACTOR
    se3 = fit3['errors'][0] / ETA_SCALE_FACTOR
    snr3 = abs(eta3) / max(se3, 1e-20)
    print_status(f"  η = {eta3:.4e} ± {se3:.4e} ({snr3:.2f}σ)", "RESULT")

    # Cluster-robust SEs (treating all DE430 as one cluster or using station if available)
    unique_stations = np.unique(st)
    if len(unique_stations) > 1:
        with np.errstate(over='ignore', divide='ignore', invalid='ignore'):
            resid3 = res - X3 @ fit3['coefficients']
        cr = cluster_robust_variance(X3, resid3, st, small_sample_correction=True)
        se3_cr = cr['se_cluster'][0] / ETA_SCALE_FACTOR
        snr3_cr = abs(eta3) / max(se3_cr, 1e-20)
        print_status(f"  η (cluster-robust) = {eta3:.4e} ± {se3_cr:.4e} ({snr3_cr:.2f}σ)", "RESULT")
    else:
        se3_cr = None
        snr3_cr = None

    # -----------------------------------------------------------------
    # Model 4: cosD-only but with environmental predictors partialed FIRST
    # -----------------------------------------------------------------
    print_status("--- Model 4: cosD after partialing environmental ---", "INFO")
    X_env = np.column_stack([inv_r, v_r, cos_cmb, np.ones(len(cos_el))])
    fit_env = robust_regression(res, X_env, scale_errors_by_birge=False)
    with np.errstate(over='ignore', divide='ignore', invalid='ignore'):
        resid_env = res - X_env @ fit_env['coefficients']
    X4 = np.column_stack([cos_el, np.ones(len(cos_el))])
    fit4 = robust_regression(resid_env, X4, scale_errors_by_birge=False)
    eta4 = fit4['coefficients'][0] / ETA_SCALE_FACTOR
    se4 = fit4['errors'][0] / ETA_SCALE_FACTOR
    snr4 = abs(eta4) / max(se4, 1e-20)
    print_status(f"  η (after partialing env) = {eta4:.4e} ± {se4:.4e} ({snr4:.2f}σ)", "RESULT")

    # -----------------------------------------------------------------
    # Sideband-phase diagnostic: do raw residuals cluster at D ± l'?
    # -----------------------------------------------------------------
    print_status("--- Sideband phase diagnostic ---", "INFO")
    mean_anomaly = np.mod(2 * np.pi * jd / 27.32166, 2 * np.pi)
    l_prime = mean_anomaly
    sideband1 = np.cos(el - l_prime)
    sideband2 = np.cos(el + l_prime)
    # Correlation of raw residuals with sidebands
    r_sb1, p_sb1 = stats.pearsonr(res, sideband1)
    r_sb2, p_sb2 = stats.pearsonr(res, sideband2)
    print_status(f"  Corr(residual, cos(D-l')) = {r_sb1:+.4f} (p={p_sb1:.4e})", "RESULT")
    print_status(f"  Corr(residual, cos(D+l')) = {r_sb2:+.4f} (p={p_sb2:.4e})", "RESULT")

    # -----------------------------------------------------------------
    # Are the 6σ outliers predicted by sideband structure?
    # -----------------------------------------------------------------
    outlier_mask = detect_outliers_sigma(res, 6.0)
    n_out = int(np.sum(outlier_mask))
    print_status(f"  6σ outliers: {n_out} ({100*n_out/len(res):.2f}%)", "RESULT")
    if n_out > 0:
        r_out_sb1, _ = stats.pearsonr(outlier_mask.astype(float), np.abs(sideband1))
        r_out_sb2, _ = stats.pearsonr(outlier_mask.astype(float), np.abs(sideband2))
        print_status(f"  Outlier correlation with |cos(D-l')|: {r_out_sb1:+.4f}", "RESULT")
        print_status(f"  Outlier correlation with |cos(D+l')|: {r_out_sb2:+.4f}", "RESULT")
    else:
        r_out_sb1 = r_out_sb2 = None

    output = {
        "step_id": "step_070",
        "status": "PASS",
        "dataset": {"n_total": int(n_total), "n_outliers_6sigma": int(n_out)},
        "model_cosd_only": {
            "eta": float(eta1), "eta_error": float(se1), "snr": float(snr1),
        },
        "model_environmental": {
            "eta": float(eta2), "eta_error": float(se2), "snr": float(snr2),
            "coefficients": {nm: float(fit2['coefficients'][i])
                             for i, nm in enumerate(names2)},
        },
        "model_full_systematic_environmental": {
            "eta": float(eta3), "eta_error": float(se3), "snr": float(snr3),
            "eta_error_cluster": float(se3_cr) if se3_cr is not None else None,
            "snr_cluster": float(snr3_cr) if snr3_cr is not None else None,
            "coefficients": {nm: float(fit3['coefficients'][i])
                             for i, nm in enumerate(names3)},
        },
        "model_cosd_after_partialing_env": {
            "eta": float(eta4), "eta_error": float(se4), "snr": float(snr4),
        },
        "sideband_diagnostics": {
            "corr_cosD_minus_lp": float(r_sb1),
            "p_cosD_minus_lp": float(p_sb1),
            "corr_cosD_plus_lp": float(r_sb2),
            "p_cosD_plus_lp": float(p_sb2),
            "outlier_corr_abs_cosD_minus_lp": float(r_out_sb1) if r_out_sb1 is not None else None,
            "outlier_corr_abs_cosD_plus_lp": float(r_out_sb2) if r_out_sb2 is not None else None,
        },
        "interpretation": (
            "If the full environmental model detects significant η in raw DE430 "
            "without outlier removal, the 'outliers' are TEP sideband signal. "
            "This removes the fragility of the step-006b DE430 result."
        ),
    }

    logger.save_step_results(output, PROJECT_ROOT, "step_070_de430_full_environmental")
    print_status("Step 070 complete.", "SUCCESS")
    return output


if __name__ == "__main__":
    main()
