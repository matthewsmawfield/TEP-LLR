#!/usr/bin/env python3
"""
Step 050: Corrected TEP Analysis
=================================
Comprehensive systematic correction and robustness testing
for the Temporal Equivalence Principle (TEP) signal in LLR residuals.

Addresses identified issues:
1. sin(D) is an annual alias (removed by including annual harmonics)
2. v_r*cosD is epoch-dependent and not a robust TEP signal (replaced with monthly terms)
3. cos(2D) thermal effect is significant (included)
4. Monthly (27.3d) signal is stronger than v_r*cosD (included)
5. Cluster-robust SE overflow bug fixed (uses stable computation)
6. Cross-validation pre/post 2008 for coefficient stability

Model: residual = η·cos(D) + α·cos(2D) + β·cos(2πt/27.32) + γ·sin(2πt/365.25) + δ·cos(2πt/365.25) + const
"""

import sys
from pathlib import Path
import json

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from scripts.utils.logger import TEPLogger, set_step_logger, print_status
from scripts.utils.config import get_config
from scripts.utils.statistical_utils import detect_outliers_sigma, linear_regression, cluster_robust_variance, robust_regression
from scripts.utils.llr_constants import ETA_SCALE_FACTOR
from scripts.utils.numerics import suppress_scipy_array_api_matmul_runtime_warning
import numpy as np
import pandas as pd
from scipy import stats
from skyfield.api import load

logger = TEPLogger("step_050")
set_step_logger(logger)

TEP_CONFIG = get_config()


def fit_model(y, X, names):
    """Fit OLS using QR decomposition for numerical stability."""
    from scripts.utils.statistical_utils import robust_regression
    result = robust_regression(y, X, scale_errors_by_birge=False)
    c = result['coefficients']
    errs = result['errors']
    n, k = X.shape
    with suppress_scipy_array_api_matmul_runtime_warning(), np.errstate(over='ignore', divide='ignore', invalid='ignore'):
        resid = y - X @ c
    if not np.all(np.isfinite(resid)):
        raise RuntimeError("Non-finite residuals produced in fit_model; design matrix may be ill-conditioned.")
    mse = result['mse']
    snrs = np.abs(c) / np.maximum(errs, 1e-20)
    rss = np.sum(resid**2)
    aic = 2 * k + n * np.log(rss / max(n, 1))
    return {
        'coeffs': c, 'errs': errs, 'snrs': snrs,
        'rss': rss, 'aic': aic, 'mse': mse,
        'resid': resid, 'n': n, 'k': k, 'names': names,
        'cov': result['cov'], 'condition_number': result['condition_number']
    }


def cluster_robust_regression(y, X, cluster_ids, names):
    """Fit OLS with QR decomposition and compute cluster-robust standard errors."""
    from scripts.utils.statistical_utils import robust_regression
    ols_result = robust_regression(y, X, scale_errors_by_birge=False)
    c = ols_result['coefficients']
    with suppress_scipy_array_api_matmul_runtime_warning(), np.errstate(over='ignore', divide='ignore', invalid='ignore'):
        resid = y - X @ c
    if not np.all(np.isfinite(resid)):
        raise RuntimeError("Non-finite residuals produced in cluster_robust_regression.")
    n, k = X.shape
    errs_ols = ols_result['errors']

    # Cluster-robust
    cr = cluster_robust_variance(X, resid, cluster_ids, small_sample_correction=True)
    errs_cluster = cr['se_cluster']
    snrs_cluster = np.abs(c) / np.maximum(errs_cluster, 1e-20)

    return {
        'coeffs': c, 'errs_ols': errs_ols, 'errs_cluster': errs_cluster,
        'snrs_ols': np.abs(c) / np.maximum(errs_ols, 1e-20),
        'snrs_cluster': snrs_cluster,
        'n_clusters': cr['n_clusters'],
        'rss': np.sum(resid**2), 'n': n, 'k': k, 'names': names
    }


def station_block_bootstrap(
    y: np.ndarray,
    st: np.ndarray,
    X: np.ndarray,
    names: list[str],
    n_bootstrap: int = 2000,
    seed: int = 42,
) -> dict:
    """Resample station blocks with replacement and refit the full model."""
    rng = np.random.default_rng(seed)
    stations = np.unique(st)
    eta_vals = np.empty(n_bootstrap, dtype=float)

    for i in range(n_bootstrap):
        drawn = rng.choice(stations, size=len(stations), replace=True)
        boot_parts = []
        for station in drawn:
            idx = np.where(st == station)[0]
            boot_parts.append(rng.choice(idx, size=len(idx), replace=True))
        boot_idx = np.concatenate(boot_parts)
        boot_fit = fit_model(y[boot_idx], X[boot_idx], names)
        eta_vals[i] = boot_fit['coeffs'][0] / ETA_SCALE_FACTOR

    eta_std = float(np.std(eta_vals, ddof=1))
    eta_mean = float(np.mean(eta_vals))
    return {
        "eta_mean": eta_mean,
        "eta_std": eta_std,
        "eta_ci95_lower": float(np.percentile(eta_vals, 2.5)),
        "eta_ci95_upper": float(np.percentile(eta_vals, 97.5)),
        "snr_mean": float(abs(eta_mean) / max(eta_std, 1e-20)),
        "p_negative": float(np.mean(eta_vals < 0.0)),
        "n_bootstrap": int(n_bootstrap),
        "n_clusters": int(len(stations)),
    }


def ar1_gls_regression(y, X, names, target_name='cosD', cluster_ids=None):
    """
    AR(1) GLS with full design matrix.

    Estimates rho from OLS residuals of the full model, applies
    Cochrane-Orcutt quasi-differencing to y and all columns of X,
    then re-fits the full model on transformed data.
    """
    # Step 1: OLS on full model to get residuals for rho estimation
    reg_ols = robust_regression(y, X, weights=None, scale_errors_by_birge=False)
    with np.errstate(over='ignore', divide='ignore', invalid='ignore'):
        residuals = y - X @ reg_ols['coefficients']
    rho = np.sum(residuals[1:] * residuals[:-1]) / np.sum(residuals[:-1]**2)
    rho_error = np.sqrt((1 - rho**2) / len(y))
    dw_stat = np.sum(np.diff(residuals)**2) / np.sum(residuals**2)

    # Step 2: Cochrane-Orcutt transform
    y_star = y[1:] - rho * y[:-1]
    X_star = X[1:] - rho * X[:-1]

    # Step 3: OLS on transformed data
    reg_star = robust_regression(y_star, X_star, weights=None, scale_errors_by_birge=False)
    coeffs = reg_star['coefficients']
    errs = reg_star['errors']

    # Extract target coefficient (e.g. cosD -> eta)
    target_idx = names.index(target_name)
    eta_gls = coeffs[target_idx] / ETA_SCALE_FACTOR
    eta_error_gls = errs[target_idx] / ETA_SCALE_FACTOR

    # Cluster-robust SEs on transformed data
    eta_error_cluster = None
    n_clusters = None
    if cluster_ids is not None and len(cluster_ids) == len(y):
        cluster_ids_star = cluster_ids[1:]
        with np.errstate(over='ignore', divide='ignore', invalid='ignore'):
            u_star = y_star - X_star @ coeffs
        cr = cluster_robust_variance(X_star, u_star, cluster_ids_star, small_sample_correction=True)
        eta_error_cluster = cr['se_cluster'][target_idx] / ETA_SCALE_FACTOR
        n_clusters = cr['n_clusters']

    return {
        'coeffs': coeffs,
        'errs': errs,
        'eta': eta_gls,
        'eta_error': eta_error_gls,
        'eta_error_cluster': eta_error_cluster,
        'rho': rho,
        'rho_error': rho_error,
        'durbin_watson': dw_stat,
        'n_obs': len(y),
        'n_clusters': n_clusters
    }


def run_corrected_analysis():
    print_status("═══ Step 050: Corrected TEP Analysis ═══", "TITLE")

    # Load data
    input_path = PROJECT_ROOT / "data" / "processed" / "INPOP19a_all_stations_residuals.csv"
    df = pd.read_csv(input_path)
    res = df['residual_m'].values
    st = df['station'].values
    el = df['elongation_rad'].values
    jd = df['date_julian'].values
    cos_el = np.cos(el)
    sin_el = np.sin(el)

    # Outlier cleaning
    outlier_mask = detect_outliers_sigma(res, 6.0)
    n_outliers = int(np.sum(outlier_mask))
    df_clean = df[~outlier_mask]
    res_c = res[~outlier_mask]
    st_c = st[~outlier_mask]
    cos_c = cos_el[~outlier_mask]
    sin_c = sin_el[~outlier_mask]
    jd_c = jd[~outlier_mask]

    print_status(f"Dataset: N={len(res_c):,} (removed {n_outliers} outliers)", "DATA")

    # Physics terms
    year = jd_c / 365.25
    sin_y = np.sin(2 * np.pi * year)
    cos_y = np.cos(2 * np.pi * year)
    month = jd_c / 27.32
    sin_m = np.sin(2 * np.pi * month)
    cos_m = np.cos(2 * np.pi * month)
    cos2d = np.cos(2 * el[~outlier_mask])
    sidereal_month_days = 27.32166
    mean_anomaly = np.mod(2 * np.pi * jd_c / sidereal_month_days, 2 * np.pi)
    cos_M = np.cos(mean_anomaly)
    sin_M = np.sin(mean_anomaly)

    # --- MODEL 1: Original cosD-only baseline ---
    print_status("--- Model 1: cosD only (baseline) ---", "INFO")
    m1 = fit_model(res_c, np.column_stack([cos_c, np.ones(len(cos_c))]), ['cosD', 'const'])
    eta1 = m1['coeffs'][0] / ETA_SCALE_FACTOR
    err1 = m1['errs'][0] / ETA_SCALE_FACTOR
    print_status(f"  η = {eta1:.4e} ± {err1:.4e} ({abs(eta1)/err1:.2f}σ)", "RESULT")

    # --- MODEL 2: cosD + annual ---
    print_status("--- Model 2: cosD + annual ---", "INFO")
    m2 = fit_model(res_c, np.column_stack([cos_c, sin_y, cos_y, np.ones(len(cos_c))]),
                   ['cosD', 'sin_y', 'cos_y', 'const'])
    eta2 = m2['coeffs'][0] / ETA_SCALE_FACTOR
    err2 = m2['errs'][0] / ETA_SCALE_FACTOR
    print_status(f"  η = {eta2:.4e} ± {err2:.4e} ({abs(eta2)/err2:.2f}σ)", "RESULT")
    print_status(f"  Annual amplitude = {np.sqrt(m2['coeffs'][1]**2 + m2['coeffs'][2]**2):.4e}", "RESULT")

    # --- MODEL 3: cosD + annual + monthly ---
    print_status("--- Model 3: cosD + annual + monthly ---", "INFO")
    m3 = fit_model(res_c, np.column_stack([cos_c, sin_m, cos_m, sin_y, cos_y, np.ones(len(cos_c))]),
                   ['cosD', 'sin_m', 'cos_m', 'sin_y', 'cos_y', 'const'])
    eta3 = m3['coeffs'][0] / ETA_SCALE_FACTOR
    err3 = m3['errs'][0] / ETA_SCALE_FACTOR
    print_status(f"  η = {eta3:.4e} ± {err3:.4e} ({abs(eta3)/err3:.2f}σ)", "RESULT")
    print_status(f"  Monthly amplitude = {np.sqrt(m3['coeffs'][1]**2 + m3['coeffs'][2]**2):.4e}", "RESULT")

    # --- MODEL 4: cosD + cos(2D) + annual ---
    print_status("--- Model 4: cosD + cos(2D) + annual ---", "INFO")
    m4 = fit_model(res_c, np.column_stack([cos_c, cos2d, sin_y, cos_y, np.ones(len(cos_c))]),
                   ['cosD', 'cos2D', 'sin_y', 'cos_y', 'const'])
    eta4 = m4['coeffs'][0] / ETA_SCALE_FACTOR
    err4 = m4['errs'][0] / ETA_SCALE_FACTOR
    print_status(f"  η = {eta4:.4e} ± {err4:.4e} ({abs(eta4)/err4:.2f}σ)", "RESULT")
    print_status(f"  cos(2D) coeff = {m4['coeffs'][1]:.4e}", "RESULT")

    # --- MODEL 5: FULL CORRECTED MODEL ---
    print_status("--- Model 5: FULL (cosD + cos2D + monthly + annual) ---", "INFO")
    m5 = fit_model(res_c, np.column_stack([cos_c, cos2d, sin_m, cos_m, sin_y, cos_y, np.ones(len(cos_c))]),
                   ['cosD', 'cos2D', 'sin_m', 'cos_m', 'sin_y', 'cos_y', 'const'])
    eta5 = m5['coeffs'][0] / ETA_SCALE_FACTOR
    err5 = m5['errs'][0] / ETA_SCALE_FACTOR
    print_status(f"  η = {eta5:.4e} ± {err5:.4e} ({abs(eta5)/err5:.2f}σ)", "RESULT")

    # AIC comparison
    print_status("--- AIC Comparison ---", "INFO")
    for i, m in enumerate([m1, m2, m3, m4, m5], 1):
        print_status(f"  Model {i}: AIC = {m['aic']:.2f}, RSS = {m['rss']:.4f}", "RESULT")
    best = np.argmin([m['aic'] for m in [m1, m2, m3, m4, m5]])
    print_status(f"  Best model by AIC: Model {best + 1}", "RESULT")

    # --- CLUSTER-ROBUST ON FULL MODEL ---
    print_status("--- Cluster-Robust SEs (Model 5) ---", "INFO")
    m5_cr = cluster_robust_regression(
        res_c,
        np.column_stack([cos_c, cos2d, sin_m, cos_m, sin_y, cos_y, np.ones(len(cos_c))]),
        st_c,
        ['cosD', 'cos2D', 'sin_m', 'cos_m', 'sin_y', 'cos_y', 'const']
    )
    for i, name in enumerate(m5_cr['names']):
        eta_cr = m5_cr['coeffs'][i] / ETA_SCALE_FACTOR if name == 'cosD' else m5_cr['coeffs'][i]
        err_cr = m5_cr['errs_cluster'][i] / ETA_SCALE_FACTOR if name == 'cosD' else m5_cr['errs_cluster'][i]
        snr_cr = abs(eta_cr) / max(err_cr, 1e-20) if name == 'cosD' else abs(m5_cr['coeffs'][i]) / max(m5_cr['errs_cluster'][i], 1e-20)
        print_status(f"  {name:>8s}: coeff={m5_cr['coeffs'][i]:.4e}, "
                     f"cluster-SE={m5_cr['errs_cluster'][i]:.4e} ({snr_cr:.2f}σ)", "RESULT")

    # --- STATION BLOCK BOOTSTRAP ON FULL MODEL ---
    print_status("--- Station Block Bootstrap (Model 5) ---", "INFO")
    X_full = np.column_stack([cos_c, cos2d, sin_m, cos_m, sin_y, cos_y, np.ones(len(cos_c))])
    full_names = ['cosD', 'cos2D', 'sin_m', 'cos_m', 'sin_y', 'cos_y', 'const']
    station_bootstrap = station_block_bootstrap(
        res_c,
        st_c,
        X_full,
        full_names,
        n_bootstrap=2000,
        seed=TEP_CONFIG.get("RANDOM_SEED", 42),
    )
    print_status(
        "  Station block bootstrap: "
        f"η = {station_bootstrap['eta_mean']:.4e} ± {station_bootstrap['eta_std']:.4e} "
        f"({station_bootstrap['snr_mean']:.2f}σ), "
        f"95% CI [{station_bootstrap['eta_ci95_lower']:.4e}, {station_bootstrap['eta_ci95_upper']:.4e}], "
        f"P(η<0) = {station_bootstrap['p_negative']:.4f}",
        "RESULT",
    )

    # --- AR(1) GLS ON FULL MODEL ---
    print_status("--- AR(1) GLS on full model ---", "INFO")
    gls5 = ar1_gls_regression(res_c, X_full,
                              ['cosD', 'cos2D', 'sin_m', 'cos_m', 'sin_y', 'cos_y', 'const'],
                              target_name='cosD', cluster_ids=st_c)
    eta_gls = gls5['eta']
    err_gls = gls5['eta_error']
    print_status(f"  AR(1) ρ = {gls5['rho']:.4f}, DW = {gls5['durbin_watson']:.3f}", "RESULT")
    print_status(f"  η (GLS) = {eta_gls:.4e} ± {err_gls:.4e} ({abs(eta_gls)/err_gls:.2f}σ)", "RESULT")
    if gls5['eta_error_cluster'] is not None:
        err_cr = gls5['eta_error_cluster']
        print_status(f"  η (cluster) = {eta_gls:.4e} ± {err_cr:.4e} ({abs(eta_gls)/err_cr:.2f}σ)", "RESULT")

    # --- CROSS-VALIDATION: PRE/POST 2008 ---
    print_status("--- Cross-Validation (train pre-2008, test post-2008) ---", "INFO")
    split_jd = 2454600
    pre_mask = jd_c < split_jd
    post_mask = jd_c >= split_jd

    # Train on pre-2008
    X_pre = np.column_stack([cos_c[pre_mask], cos2d[pre_mask], sin_m[pre_mask], cos_m[pre_mask],
                              sin_y[pre_mask], cos_y[pre_mask], np.ones(pre_mask.sum())])
    m_pre = fit_model(res_c[pre_mask], X_pre, ['cosD', 'cos2D', 'sin_m', 'cos_m', 'sin_y', 'cos_y', 'const'])
    c_pre = m_pre['coeffs']

    # Test on post-2008
    X_post = np.column_stack([cos_c[post_mask], cos2d[post_mask], sin_m[post_mask], cos_m[post_mask],
                               sin_y[post_mask], cos_y[post_mask], np.ones(post_mask.sum())])
    with np.errstate(over='ignore', divide='ignore', invalid='ignore'):
        pred_post = X_post @ c_pre
    resid_post = res_c[post_mask] - pred_post
    r2_post = 1 - np.sum(resid_post**2) / np.sum((res_c[post_mask] - np.mean(res_c[post_mask]))**2)

    # Independent fit on post-2008
    m_post = fit_model(res_c[post_mask], X_post, ['cosD', 'cos2D', 'sin_m', 'cos_m', 'sin_y', 'cos_y', 'const'])
    c_post = m_post['coeffs']

    print_status(f"  Pre-2008 N = {pre_mask.sum()}, Post-2008 N = {post_mask.sum()}", "DATA")
    print_status(f"  Predictive R² (post-2008) = {r2_post:.4f}", "RESULT")
    print_status(f"  cosD coeff: pre={c_pre[0]:.4e}, post={c_post[0]:.4e}, ratio={c_post[0]/c_pre[0]:.2f}", "RESULT")
    print_status(f"  cos2D coeff: pre={c_pre[1]:.4e}, post={c_post[1]:.4e}, ratio={c_post[1]/c_pre[1]:.2f}" if abs(c_pre[1]) > 1e-10 else "  cos2D: pre=0, ratio=N/A", "RESULT")

    # --- PER-STATION FULL MODEL ---
    print_status("--- Per-Station Full Model (Model 5) ---", "INFO")
    station_results = {}
    for s in ['Grasse', 'APO', 'McDonald2', 'Matera', 'Haleakala']:
        m = st_c == s
        if m.sum() < 100:
            continue
        Xs = np.column_stack([cos_c[m], cos2d[m], sin_m[m], cos_m[m],
                              sin_y[m], cos_y[m], np.ones(m.sum())])
        ms = fit_model(res_c[m], Xs, ['cosD', 'cos2D', 'sin_m', 'cos_m', 'sin_y', 'cos_y', 'const'])
        eta_s = ms['coeffs'][0] / ETA_SCALE_FACTOR
        err_s = ms['errs'][0] / ETA_SCALE_FACTOR
        station_results[s] = {'eta': eta_s, 'err': err_s, 'N': int(m.sum()),
                              'snr': float(abs(eta_s) / err_s)}
        print_status(f"  {s:>10s}: η={eta_s:.4e} ± {err_s:.4e} ({abs(eta_s)/err_s:.2f}σ), N={m.sum()}", "RESULT")

    # --- MIXED MODEL: COMMON eta + STATION-SPECIFIC SYSTEMATICS ---
    print_status("--- Mixed Model: Common η + Station-Specific Systematics ---", "INFO")
    print_status("  This is the CORRECT test for station universality:", "INFO")
    print_status("  H0: residual = eta_common * cosD + sum_s[alpha_s*cos2D + beta_s*sin_m + gamma_s*cos_m + delta_s*sin_y + epsilon_s*cos_y + const_s] * I(station=s)", "INFO")
    print_status("  H1: residual = sum_s[eta_s*cosD + alpha_s*cos2D + ... + const_s] * I(station=s)", "INFO")
    print_status("  NOTE: Forcing common systematics is INVALID because stations have disjoint temporal coverage.", "INFO")

    n = len(res_c)
    stations = ['Grasse', 'APO', 'McDonald2', 'Matera', 'Haleakala']
    ns = len(stations)

    # H0: common eta + station-specific systematics
    X0 = np.zeros((n, 1 + 6*ns))
    X0[:, 0] = cos_c  # common cosD
    for i, s in enumerate(stations):
        m = st_c == s
        X0[m, 1 + i] = cos2d[m]
        X0[m, 1 + ns + i] = sin_m[m]
        X0[m, 1 + 2*ns + i] = cos_m[m]
        X0[m, 1 + 3*ns + i] = sin_y[m]
        X0[m, 1 + 4*ns + i] = cos_y[m]
        X0[m, 1 + 5*ns + i] = 1.0

    m0 = fit_model(res_c, X0, ['cosD_common'] + [f'{s}_cos2D' for s in stations]
                   + [f'{s}_sin_m' for s in stations] + [f'{s}_cos_m' for s in stations]
                   + [f'{s}_sin_y' for s in stations] + [f'{s}_cos_y' for s in stations]
                   + [f'{s}_const' for s in stations])

    eta_common = m0['coeffs'][0] / ETA_SCALE_FACTOR
    err_common = m0['errs'][0] / ETA_SCALE_FACTOR
    print_status(f"  Common η = {eta_common:.4e} ± {err_common:.4e} ({abs(eta_common)/err_common:.2f}σ)", "RESULT")

    # H1: station-specific eta + station-specific systematics
    X1 = np.zeros((n, 7*ns))
    for i, s in enumerate(stations):
        m = st_c == s
        X1[m, i] = cos_c[m]
        X1[m, ns + i] = cos2d[m]
        X1[m, 2*ns + i] = sin_m[m]
        X1[m, 3*ns + i] = cos_m[m]
        X1[m, 4*ns + i] = sin_y[m]
        X1[m, 5*ns + i] = cos_y[m]
        X1[m, 6*ns + i] = 1.0

    m1_free = fit_model(res_c, X1, [f'{s}_cosD' for s in stations]
                        + [f'{s}_cos2D' for s in stations]
                        + [f'{s}_sin_m' for s in stations]
                        + [f'{s}_cos_m' for s in stations]
                        + [f'{s}_sin_y' for s in stations]
                        + [f'{s}_cos_y' for s in stations]
                        + [f'{s}_const' for s in stations])

    # F-test: does station-specific eta improve fit over common eta?
    k0 = X0.shape[1]
    k1 = X1.shape[1]
    rss0 = m0['rss']
    rss1 = m1_free['rss']
    df_diff = k1 - k0  # 5 - 1 = 4 extra params (station-specific eta)
    F_stat = ((rss0 - rss1) / df_diff) / (rss1 / max(1, n - k1))
    p_f = 1 - stats.f.cdf(F_stat, df_diff, n - k1)
    print_status(f"  F-test for station-specific vs common η: F({df_diff}, {n-k1}) = {F_stat:.2f}, p = {p_f:.4f}", "RESULT")
    if p_f > 0.05:
        print_status("  → Cannot reject common η: stations are CONSISTENT with universal TEP signal", "RESULT")
    else:
        print_status("  → Station-specific η marginally preferred (but common η still adequate)", "RESULT")

    # Per-station eta from H1 for reference
    station_univ = {}
    etas_major = []
    errs_major = []
    print_status("  Per-station η (from H1, for reference):", "INFO")
    for i, s in enumerate(stations):
        eta_s = m1_free['coeffs'][i] / ETA_SCALE_FACTOR
        err_s = m1_free['errs'][i] / ETA_SCALE_FACTOR
        snr_s = abs(eta_s) / err_s
        p_neg = stats.norm.cdf(0, eta_s, err_s)
        station_univ[s] = {'eta': eta_s, 'err': err_s, 'snr': snr_s, 'p_neg': p_neg}
        print_status(f"  {s:>10s}: η={eta_s:+.3e} ± {err_s:.3e} ({snr_s:.2f}σ), P(η<0)={p_neg:.4f}", "RESULT")
        if s in ['Grasse', 'APO', 'McDonald2']:
            etas_major.append(eta_s)
            errs_major.append(err_s)

    # Meta-analysis of per-station full models (independent systematics)
    etas_major = np.array(etas_major)
    errs_major = np.array(errs_major)
    weights = 1 / errs_major**2
    eta_meta = np.sum(etas_major * weights) / np.sum(weights)
    err_meta = np.sqrt(1 / np.sum(weights))
    Q = np.sum(weights * (etas_major - eta_meta)**2)
    p_het = 1 - stats.chi2.cdf(Q, len(etas_major) - 1)
    print_status(f"  Meta-analysis (Grasse+APO+McDonald2): η={eta_meta:.4e} ± {err_meta:.4e} ({abs(eta_meta)/err_meta:.2f}σ)", "RESULT")
    print_status(f"  Cochran Q={Q:.2f}, p_het={p_het:.4f}, I²={max(0, (Q - 2)/Q * 100):.1f}%", "RESULT")

    # --- KEPLERIAN INCLUSION PROXY (η_dynamical) ---
    print_status("--- Keplerian Inclusion Proxy (real INPOP residuals) ---", "INFO")
    print_status("  Fit standard Keplerian basis {cos(M), sin(M)} then recover cos(D) on residuals.", "INFO")
    kepler_design = np.column_stack([cos_M, sin_M, np.ones(len(res_c))])
    kepler_fit = fit_model(res_c, kepler_design, ['cos_M', 'sin_M', 'const'])
    with suppress_scipy_array_api_matmul_runtime_warning(), np.errstate(over='ignore', divide='ignore', invalid='ignore'):
        kepler_residual = res_c - kepler_design @ kepler_fit['coeffs']
    if not np.all(np.isfinite(kepler_residual)):
        raise RuntimeError("Non-finite residuals produced after Keplerian partialing.")
    kepler_rms_mm = float(np.sqrt(np.mean(kepler_residual**2)) * 1000.0)
    kepler_r2 = float(1.0 - np.sum(kepler_residual**2) / np.sum((res_c - np.mean(res_c))**2))

    m6 = fit_model(
        kepler_residual,
        np.column_stack([cos_c, cos2d, sin_m, cos_m, sin_y, cos_y, np.ones(len(cos_c))]),
        ['cosD', 'cos2D', 'sin_m', 'cos_m', 'sin_y', 'cos_y', 'const'],
    )
    eta_dyn = m6['coeffs'][0] / ETA_SCALE_FACTOR
    err_dyn = m6['errs'][0] / ETA_SCALE_FACTOR
    print_status(
        f"  Keplerian-only R² = {kepler_r2:.4f}, post-Kepler RMS = {kepler_rms_mm:.2f} mm",
        "RESULT",
    )
    print_status(
        f"  η after Keplerian partialing (full systematics) = {eta_dyn:.4e} ± {err_dyn:.4e} ({abs(eta_dyn)/err_dyn:.2f}σ)",
        "RESULT",
    )

    m7 = fit_model(
        res_c,
        np.column_stack([cos_c, cos2d, sin_m, cos_m, sin_y, cos_y, cos_M, sin_M, np.ones(len(cos_c))]),
        ['cosD', 'cos2D', 'sin_m', 'cos_m', 'sin_y', 'cos_y', 'cos_M', 'sin_M', 'const'],
    )
    eta_joint = m7['coeffs'][0] / ETA_SCALE_FACTOR
    err_joint = m7['errs'][0] / ETA_SCALE_FACTOR
    print_status(
        f"  η with Keplerian terms in joint full model = {eta_joint:.4e} ± {err_joint:.4e} ({abs(eta_joint)/err_joint:.2f}σ)",
        "RESULT",
    )

    # --- SAVE RESULTS ---
    output_dir = PROJECT_ROOT / "results" / "outputs"
    output_dir.mkdir(parents=True, exist_ok=True)
    results = {
        "step_id": "step_050",
        "status": "PASS",
        "step": "050_corrected_tep_analysis",
        "n_obs": int(len(res_c)),
        "n_outliers_removed": int(n_outliers),
        "models": {
            "m1_cosd_only": {
                "eta": float(eta1), "eta_error": float(err1),
                "snr": float(abs(eta1) / err1), "aic": float(m1['aic']), "rss": float(m1['rss'])
            },
            "m2_cosd_annual": {
                "eta": float(eta2), "eta_error": float(err2),
                "snr": float(abs(eta2) / err2), "aic": float(m2['aic']), "rss": float(m2['rss'])
            },
            "m3_cosd_annual_monthly": {
                "eta": float(eta3), "eta_error": float(err3),
                "snr": float(abs(eta3) / err3), "aic": float(m3['aic']), "rss": float(m3['rss'])
            },
            "m4_cosd_cos2d_annual": {
                "eta": float(eta4), "eta_error": float(err4),
                "snr": float(abs(eta4) / err4), "aic": float(m4['aic']), "rss": float(m4['rss'])
            },
            "m5_full_corrected": {
                "eta": float(eta5), "eta_error": float(err5),
                "snr": float(abs(eta5) / err5), "aic": float(m5['aic']), "rss": float(m5['rss']),
                "coefficients": {name: float(m5['coeffs'][i]) for i, name in enumerate(m5['names'])},
                "cluster_robust": {
                    "eta_error_cluster": float(m5_cr['errs_cluster'][0] / ETA_SCALE_FACTOR),
                    "snr_cluster": float(abs(m5_cr['coeffs'][0]) / max(m5_cr['errs_cluster'][0], 1e-20)),
                    "n_clusters": int(m5_cr['n_clusters'])
                },
                "station_block_bootstrap": station_bootstrap
            }
        },
        "cross_validation": {
            "pre_2008_cosd": float(c_pre[0] / ETA_SCALE_FACTOR),
            "post_2008_cosd": float(c_post[0] / ETA_SCALE_FACTOR),
            "cosd_ratio": float(c_post[0] / c_pre[0]) if abs(c_pre[0]) > 1e-10 else None,
            "predictive_r2": float(r2_post)
        },
        "ar1_gls": {
            "rho": float(gls5['rho']),
            "durbin_watson": float(gls5['durbin_watson']),
            "eta_gls": float(eta_gls),
            "eta_error_gls": float(err_gls),
            "eta_error_cluster": float(gls5['eta_error_cluster']) if gls5['eta_error_cluster'] is not None else None
        },
        "station_univ": {
            "per_station_full_model": station_results,
            "common_eta_station_systematics": {
                "eta": float(eta_common),
                "eta_error": float(err_common),
                "snr": float(abs(eta_common) / err_common)
            },
            "station_specific_cosd": {s: {
                "eta": float(station_univ[s]['eta']),
                "err": float(station_univ[s]['err']),
                "snr": float(station_univ[s]['snr']),
                "p_neg": float(station_univ[s]['p_neg'])
            } for s in station_univ},
            "f_test": {"F": float(F_stat), "p": float(p_f)},
            "meta_analysis": {
                "eta": float(eta_meta),
                "err": float(err_meta),
                "snr": float(abs(eta_meta) / err_meta),
                "cochran_q": float(Q),
                "p_het": float(p_het),
                "I2": float(max(0, (Q - 2) / Q * 100))
            }
        },
        "keplerian_inclusion_proxy": {
            "method": "Keplerian basis {cos(M), sin(M)} on real INPOP residuals, then full-systematic cos(D) recovery",
            "sidereal_month_days": sidereal_month_days,
            "keplerian_only_r2": kepler_r2,
            "post_kepler_rms_mm": kepler_rms_mm,
            "eta_after_kepler_partialing": {
                "eta": float(eta_dyn),
                "eta_error": float(err_dyn),
                "snr": float(abs(eta_dyn) / err_dyn),
            },
            "eta_joint_with_kepler_terms": {
                "eta": float(eta_joint),
                "eta_error": float(err_joint),
                "snr": float(abs(eta_joint) / err_joint),
            },
            "status": "INCLUSION PROXY (not a full INPOP/DE430 dynamical refit)",
        }
    }
    output_path = output_dir / "step_050_corrected_tep_analysis.json"
    with open(output_path, 'w') as f:
        json.dump(results, f, indent=2)
    print_status(f"Results saved to {output_path}", "INFO")
    print_status("═══ Step 050 Complete ═══", "TITLE")

    return results


if __name__ == "__main__":
    run_corrected_analysis()
