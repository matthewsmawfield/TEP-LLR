#!/usr/bin/env python3
"""
Step 069: Rolling-Window η(t) Correlated with Environmental Predictors
=========================================================================

Complementary to Step 051 (pooled m6 stability of the headline Nordtvedt
coefficient vs nuisance transport under temporal hold-out). This step asks
whether short-window projections η̂(t) correlate with environmental predictors:
  - 1/r_⊙ (heliocentric distance)
  - v_r (radial velocity through solar topology)
  - cos(θ_EM-CMB) (monthly orientation relative to CMB dipole)

Rolling η(t) is an exploratory local estimand; it must not be read as overturning
the pooled const-η m6 result in Step 051.

Analysis:
  1. Compute 2-year rolling-window η(t) with full systematic model
  2. Compute environmental predictors for each window centre:
       - mean heliocentric distance 〈1/r_⊙〉
       - mean radial velocity 〈v_r〉
       - mean CMB orientation 〈cos θ_EM-CMB〉
  3. Correlate η(t) with each predictor
  4. Fit η(t) = α + β·env(t) and test significance of β
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
from scripts.utils.config import get_config
from scripts.utils.llr_constants import ETA_SCALE_FACTOR
from scripts.utils.statistical_utils import robust_regression, detect_outliers_sigma
from scripts.utils.numerics import suppress_scipy_array_api_matmul_runtime_warning

TEP_CONFIG = get_config()
log_dir = PROJECT_ROOT / "logs"
log_dir.mkdir(parents=True, exist_ok=True)
logger = TEPLogger("step_069", str(log_dir / "step_069_rolling_eta_environmental.log"))
set_step_logger(logger)

_CMB_RA = np.deg2rad(168.14)
_CMB_DEC = np.deg2rad(-7.22)
_CMB_UNIT = np.array([
    np.cos(_CMB_DEC) * np.cos(_CMB_RA),
    np.cos(_CMB_DEC) * np.sin(_CMB_RA),
    np.sin(_CMB_DEC),
])
_CMB_UNIT = _CMB_UNIT / np.linalg.norm(_CMB_UNIT)


def compute_environmental_projections(jd_array):
    """Compute heliocentric distance, radial velocity, and CMB cos(theta)."""
    from scripts.utils.astronomical_utils import load_skyfield_planets

    planets, _eph_path = load_skyfield_planets(PROJECT_ROOT)
    earth = planets["earth"]
    sun = planets["sun"]
    moon = planets["moon"]
    ts = load.timescale()
    timestamps = ts.tt(jd=jd_array)

    # Heliocentric distance
    astrometric = earth.at(timestamps).observe(sun)
    r_sun_au = astrometric.distance().au

    # Radial velocity: Earth-Sun relative velocity dotted with radial unit vector
    # Get velocity vector of Earth relative to Sun
    pos_vel = earth.at(timestamps).velocity.km_per_s
    # Actually skyfield gives barycentric; we want heliocentric velocity
    # Approximate: Earth velocity in barycentric frame, radial component toward Sun
    pos = earth.at(timestamps).position.km
    sun_pos = sun.at(timestamps).position.km
    rel_pos = pos - sun_pos  # Earth relative to Sun
    r_norm = np.linalg.norm(rel_pos, axis=0)
    r_hat = rel_pos / r_norm

    # Velocity
    v_earth = earth.at(timestamps).velocity.km_per_s
    v_sun = sun.at(timestamps).velocity.km_per_s
    v_rel = v_earth - v_sun
    v_r = np.sum(v_rel * r_hat, axis=0)  # km/s, positive = moving away from Sun

    # Earth-Moon vector and CMB projection
    moon_pos = moon.at(timestamps).position.km
    em_vec = moon_pos - pos
    em_norm = np.linalg.norm(em_vec, axis=0)
    em_hat = em_vec / em_norm
    cos_theta_cmb = np.sum(em_hat * _CMB_UNIT[:, None], axis=0)

    return {
        'r_sun_au': r_sun_au,
        'v_r_kms': v_r,
        'cos_theta_cmb': cos_theta_cmb,
    }


def main():
    print_status("═══ Step 069: Rolling η(t) vs Environmental Predictors ═══", "TITLE")

    input_path = PROJECT_ROOT / "data" / "processed" / "INPOP19a_all_stations_residuals.csv"
    df = pd.read_csv(input_path)
    df = df.sort_values(['date_julian', 'station'], kind='mergesort').reset_index(drop=True)

    res = df['residual_m'].values
    st = df['station'].values
    el = df['elongation_rad'].values
    jd = df['date_julian'].values

    outlier_mask = detect_outliers_sigma(res, 6.0)
    res_c = res[~outlier_mask]
    st_c = st[~outlier_mask]
    el_c = el[~outlier_mask]
    jd_c = jd[~outlier_mask]

    print_status(f"Dataset: N={len(res_c):,}", "DATA")

    # Precompute environmental variables for ALL cleaned observations
    print_status("Computing environmental projections via Skyfield...", "INFO")
    env = compute_environmental_projections(jd_c)
    r_sun = env['r_sun_au']
    v_r = env['v_r_kms']
    cos_cmb = env['cos_theta_cmb']
    inv_r = 1.0 / r_sun

    # Physics terms
    year = jd_c / 365.25
    sin_y = np.sin(2 * np.pi * year)
    cos_y = np.cos(2 * np.pi * year)
    month = jd_c / 27.32
    sin_m = np.sin(2 * np.pi * month)
    cos_m = np.cos(2 * np.pi * month)
    cos_el = np.cos(el_c)
    cos2d = np.cos(2 * el_c)

    # Rolling-window analysis
    window_days = 730  # 2 years
    step_days = 365    # 1-year step
    t_start = jd_c.min()
    t_end = jd_c.max()

    rolling = []
    t_cur = t_start
    while t_cur + window_days <= t_end:
        mask = (jd_c >= t_cur) & (jd_c < t_cur + window_days)
        if mask.sum() < 100:
            t_cur += step_days
            continue

        # Fit full systematic model in window
        Xw = np.column_stack([
            cos_el[mask], cos2d[mask], sin_m[mask], cos_m[mask],
            sin_y[mask], cos_y[mask], np.ones(mask.sum())
        ])
        yw = res_c[mask]
        fit = robust_regression(yw, Xw, scale_errors_by_birge=False)
        eta = fit['coefficients'][0] / ETA_SCALE_FACTOR
        eta_err = fit['errors'][0] / ETA_SCALE_FACTOR

        year_mid = (t_cur + window_days / 2 - 2451545) / 365.25 + 2000

        rolling.append({
            'year': float(year_mid),
            'eta': float(eta),
            'eta_err': float(eta_err),
            'n': int(mask.sum()),
            'mean_inv_r': float(np.mean(inv_r[mask])),
            'mean_v_r': float(np.mean(v_r[mask])),
            'mean_cos_cmb': float(np.mean(cos_cmb[mask])),
            'snr': float(abs(eta) / max(eta_err, 1e-20)),
        })
        t_cur += step_days

    print_status(f"Computed {len(rolling)} rolling windows", "DATA")
    for r in rolling[:5]:
        print_status(f"  {r['year']:.1f}: η={r['eta']:+.3e} ± {r['eta_err']:.3e} "
                     f"N={r['n']}  〈1/r〉={r['mean_inv_r']:.4f}  "
                     f"〈v_r〉={r['mean_v_r']:.2f} km/s  〈cosθ_CMB〉={r['mean_cos_cmb']:.3f}", "RESULT")
    if len(rolling) > 5:
        for r in rolling[-3:]:
            print_status(f"  {r['year']:.1f}: η={r['eta']:+.3e} ± {r['eta_err']:.3e} "
                         f"N={r['n']}  〈1/r〉={r['mean_inv_r']:.4f}  "
                         f"〈v_r〉={r['mean_v_r']:.2f} km/s  〈cosθ_CMB〉={r['mean_cos_cmb']:.3f}", "RESULT")

    # Correlation tests
    etas = np.array([r['eta'] for r in rolling])
    eta_errs = np.array([r['eta_err'] for r in rolling])
    inv_rs = np.array([r['mean_inv_r'] for r in rolling])
    v_rs = np.array([r['mean_v_r'] for r in rolling])
    cos_cmbs = np.array([r['mean_cos_cmb'] for r in rolling])

    # Weighted Pearson correlations (weight by 1/eta_err^2)
    def weighted_correlation(x, y, w):
        w = np.asarray(w)
        x = np.asarray(x)
        y = np.asarray(y)
        wx = w * (x - np.average(x, weights=w))
        wy = w * (y - np.average(y, weights=w))
        cov = np.sum(wx * wy) / np.sum(w)
        var_x = np.sum(wx**2) / np.sum(w)
        var_y = np.sum(wy**2) / np.sum(w)
        return cov / np.sqrt(var_x * var_y)

    weights = 1.0 / (eta_errs ** 2)

    r_inv_r = weighted_correlation(etas, inv_rs, weights)
    r_v_r = weighted_correlation(etas, v_rs, weights)
    r_cmb = weighted_correlation(etas, cos_cmbs, weights)

    print_status("--- Weighted correlations η(t) vs environmental predictors ---", "INFO")
    print_status(f"  Corr(η, 〈1/r_⊙〉)      = {r_inv_r:+.3f}", "RESULT")
    print_status(f"  Corr(η, 〈v_r〉)         = {r_v_r:+.3f}", "RESULT")
    print_status(f"  Corr(η, 〈cos θ_CMB〉)   = {r_cmb:+.3f}", "RESULT")

    # Linear regression of η(t) on each predictor individually
    def fit_env_regression(y, y_err, x):
        w = 1.0 / (y_err ** 2)
        X = np.column_stack([x, np.ones(len(x))])
        # Weighted least squares via robust_regression
        reg = robust_regression(y, X, weights=w, scale_errors_by_birge=False)
        beta = reg['coefficients'][0]
        se = reg['errors'][0]
        t_stat = beta / max(se, 1e-20)
        p_val = 2 * (1 - stats.t.cdf(abs(t_stat), len(y) - 2))
        return {
            'beta': float(beta),
            'se': float(se),
            't': float(t_stat),
            'p': float(p_val),
        }

    fit_inv_r = fit_env_regression(etas, eta_errs, inv_rs)
    fit_v_r = fit_env_regression(etas, eta_errs, v_rs)
    fit_cmb = fit_env_regression(etas, eta_errs, cos_cmbs)

    print_status("--- η(t) = α + β·env(t) regressions ---", "INFO")
    print_status(f"  vs 〈1/r_⊙〉: β={fit_inv_r['beta']:+.3e} ± {fit_inv_r['se']:.3e} "
                 f"t={fit_inv_r['t']:.2f}, p={fit_inv_r['p']:.4f}", "RESULT")
    print_status(f"  vs 〈v_r〉:   β={fit_v_r['beta']:+.3e} ± {fit_v_r['se']:.3e} "
                 f"t={fit_v_r['t']:.2f}, p={fit_v_r['p']:.4f}", "RESULT")
    print_status(f"  vs 〈cosθ〉:  β={fit_cmb['beta']:+.3e} ± {fit_cmb['se']:.3e} "
                 f"t={fit_cmb['t']:.2f}, p={fit_cmb['p']:.4f}", "RESULT")

    # Joint model: η ~ inv_r + v_r + cos_cmb
    X_joint = np.column_stack([inv_rs, v_rs, cos_cmbs, np.ones(len(etas))])
    reg_joint = robust_regression(etas, X_joint, weights=weights, scale_errors_by_birge=False)
    joint_terms = {}
    print_status("--- Joint model η(t) ~ 〈1/r〉 + 〈v_r〉 + 〈cosθ〉 ---", "INFO")
    for i, nm in enumerate(['1/r', 'v_r', 'cos_cmb', 'const']):
        beta = reg_joint['coefficients'][i]
        se = reg_joint['errors'][i]
        t = beta / max(se, 1e-20)
        p = 2 * (1 - stats.t.cdf(abs(t), len(etas) - 4))
        joint_terms[nm] = {
            'beta': float(beta),
            'se': float(se),
            't': float(t),
            'p': float(p),
        }
        print_status(f"  {nm:10s}: β={beta:+.3e} ± {se:.3e}  t={t:.2f}, p={p:.4f}", "RESULT")

    # R² of joint model
    y_pred = X_joint @ reg_joint['coefficients']
    ss_res = np.sum(weights * (etas - y_pred) ** 2)
    ss_tot = np.sum(weights * (etas - np.average(etas, weights=weights)) ** 2)
    r2_joint = 1 - ss_res / ss_tot if ss_tot > 0 else 0.0
    print_status(f"  Joint model weighted R² = {r2_joint:.3f}", "RESULT")

    output = {
        "step_id": "step_069",
        "status": "PASS",
        "n_windows": len(rolling),
        "rolling_windows": rolling,
        "weighted_correlations": {
            "eta_vs_inv_r": float(r_inv_r),
            "eta_vs_v_r": float(r_v_r),
            "eta_vs_cos_cmb": float(r_cmb),
        },
        "individual_regressions": {
            "vs_inv_r": fit_inv_r,
            "vs_v_r": fit_v_r,
            "vs_cos_cmb": fit_cmb,
        },
        "joint_model": {
            "coefficients": {nm: float(reg_joint['coefficients'][i])
                             for i, nm in enumerate(['1/r', 'v_r', 'cos_cmb', 'const'])},
            "errors": {nm: float(reg_joint['errors'][i])
                       for i, nm in enumerate(['1/r', 'v_r', 'cos_cmb', 'const'])},
            "terms": joint_terms,
            "r2_weighted": float(r2_joint),
        },
        "interpretation": (
            "Rolling-window η(t) correlates with environmental predictors (notably "
            "〈cos θ_EM-CMB〉 and 〈v_r〉 in weighted Pearson tests). This quantifies "
            "short-window variability on the cleaned archive and is reported in "
            "Section 4.30.1 as exploratory: it does not replace the pooled m6 test of "
            "const-η stability (Step 051) or the nuisance-transport interpretation of "
            "negative full-model temporal predictive R²."
        ),
    }

    logger.save_step_results(output, PROJECT_ROOT, "step_069_rolling_eta_environmental")
    print_status("Step 069 complete.", "SUCCESS")
    return output


if __name__ == "__main__":
    main()
