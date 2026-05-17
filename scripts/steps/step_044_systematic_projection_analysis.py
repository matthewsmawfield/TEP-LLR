#!/usr/bin/env python3
"""
Step 044: Systematic Projection Analysis for TEP-LLR

Computes the cos(elongation)-projected systematic bias for each error source.

Critical insight: In a linear regression y = A*x + B, only the component of
a systematic source s that is correlated with x biases the slope A. The
orthogonal component adds noise (increasing error bars) but does not bias A.

For TEP: η = A / ETA_SCALE_FACTOR, so the systematic bias to η is:
    δη_sys = cov(s, cos_elong) / var(cos_elong) / ETA_SCALE_FACTOR

This step also performs a phase-locked differential analysis (new moon vs
full moon), which cancels all common-mode systematics by construction.

Sources:
1. Ephemeris modeling: residual difference between INPOP19a and DE430
2. Atmospheric delay: monthly mean variation in detrended residuals
3. Instrumental: station-to-station mean differences
4. Tidal modeling: cos(2*elongation) harmonic amplitude
5. Thermal expansion: diurnal (24-hr) sinusoidal amplitude
"""

import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

import argparse
import json
import numpy as np
from scripts.utils.numerics import stable_lstsq, suppress_scipy_array_api_matmul_runtime_warning
from scripts.utils.config import get_config
import pandas as pd
from scipy import stats

from scripts.utils.llr_constants import ETA_SCALE_FACTOR, ELONGATION_MASK_WIDTH
from scripts.utils.logger import TEPLogger, set_step_logger, set_verbose_mode, print_status
from scripts.utils.statistical_utils import linear_regression

TEP_CONFIG = get_config()


def compute_systematic_projection(residuals, cos_elong, systematic_component):
    """Compute the cos(elongation)-projected bias from a systematic component.

    Returns the bias to the slope A and the corresponding η bias.
    """
    cos_centered = cos_elong - np.mean(cos_elong)
    sys_centered = systematic_component - np.mean(systematic_component)
    var_cos = np.sum(cos_centered**2)
    if var_cos == 0:
        return 0.0, 0.0
    cov_sys_cos = np.sum(sys_centered * cos_centered)
    bias_A = cov_sys_cos / var_cos  # bias to slope A [meters]
    bias_eta = bias_A / ETA_SCALE_FACTOR
    # Also compute correlation for diagnostics
    with suppress_scipy_array_api_matmul_runtime_warning():
        r_sc, p_sc = stats.pearsonr(systematic_component, cos_elong)
    return float(bias_A), float(bias_eta), float(r_sc), float(p_sc)


def run_systematic_projection_analysis(df, verbose=False):
    print_status("═══ Starting Step 044: Systematic Projection Analysis...", "TITLE")
    print_status("═══ STEP PURPOSE: Compute cos(elongation)-projected systematic bias for each error source", "INFO")
    print_status("═══ METHOD: Project each systematic component onto cos(elongation); phase-locked differential", "INFO")

    n = len(df)
    residuals = df['residual_m'].values
    cos_elong = np.cos(df['elongation_rad'].values)
    elongation = df['elongation_rad'].values
    jd = df['date_julian'].values

    print_status(f"═══ DATA SUMMARY", "INFO")
    print_status(f"    Dataset: N = {n:,} observations", "DATA")

    # ------------------------------------------------------------------
    # Baseline TEP fit (weighted, consistent with step_002)
    # ------------------------------------------------------------------
    weights = None
    if 'sigma_m' in df.columns:
        sigma = df['sigma_m'].values
        weights = 1.0 / sigma**2
        weights = np.where(np.isfinite(weights) & (weights > 0), weights, 1.0)
        print_status("    Using per-observation weights (1/σ²) from 'sigma_m'", "DATA")
    elif 'uncertainty_m' in df.columns:
        sigma = df['uncertainty_m'].values
        weights = 1.0 / sigma**2
        weights = np.where(np.isfinite(weights) & (weights > 0), weights, 1.0)
        print_status("    Using per-observation weights (1/σ²) from 'uncertainty_m'", "DATA")
    else:
        print_status("    No uncertainty column found; using unweighted fit", "DATA")

    reg_base = linear_regression(residuals, cos_elong, weights=weights)
    eta_base = reg_base['eta']
    eta_err_base = reg_base['eta_error']
    snr_base = abs(eta_base) / eta_err_base if eta_err_base > 0 else 0.0

    print_status(f"    Baseline η = {eta_base:.3e} ± {eta_err_base:.3e} ({snr_base:.2f}σ)", "DATA")

    # ------------------------------------------------------------------
    # Detrend residuals: remove best-fit TEP signal
    # ------------------------------------------------------------------
    X = np.column_stack([cos_elong, np.ones(n)])
    coeffs_tep, _, _, _ = stable_lstsq(X, residuals)
    detrended = residuals - coeffs_tep[0] * cos_elong  # m

    # ------------------------------------------------------------------
    # 1. Ephemeris systematic: INPOP19a vs DE430 residual difference
    # ------------------------------------------------------------------
    ephem_bias_A = 0.0
    ephem_bias_eta = 0.0
    ephem_r = 0.0
    ephem_p = 1.0
    ephem_available = False

    inpop_path = PROJECT_ROOT / 'data' / 'processed' / 'INPOP19a_all_stations_residuals.csv'
    de430_path = PROJECT_ROOT / 'data' / 'processed' / 'DE430_all_residuals.csv'

    # Load pre-computed cross-ephemeris comparison from step_006.
    # The ephemeris systematic is the scatter in η across ephemerides.
    # For N_ephem ephemerides, std(η_ephem) = sqrt(sum(η_i - mean)^2 / (N-1)).
    # This is the proper estimate because each η is computed independently on
    # its own residual set, avoiding the confounding of fitting difference
    # residuals (which mixes TEP signal with ephemeris modelling differences).
    step006_path = PROJECT_ROOT / 'results' / 'outputs' / 'step_006_multi_ephemeris_comparison.json'
    if step006_path.exists():
        with open(step006_path, 'r') as f:
            step006_data = json.load(f)
        etas = [v['eta'] for v in step006_data.get('comparisons', {}).values() if 'eta' in v]
        if len(etas) >= 2:
            ephem_eta_std = float(np.std(etas, ddof=1))
            ephem_bias_eta = ephem_eta_std
            ephem_bias_A = ephem_eta_std * ETA_SCALE_FACTOR
            ephem_available = True
            print_status(f"    Ephemeris scatter: N={len(etas)}, std(η) = {ephem_eta_std:.3e}", "CALC")
        elif len(etas) == 1:
            # Only one ephemeris available: use half the absolute value as conservative bound
            ephem_bias_eta = abs(etas[0]) * 0.5
            ephem_bias_A = ephem_bias_eta * ETA_SCALE_FACTOR
            ephem_available = True
            print_status(f"    Ephemeris scatter: N=1, conservative bound = {ephem_bias_eta:.3e}", "CALC")
        else:
            print_status("    WARNING: No ephemeris eta values in step_006 output", "WARNING")
    else:
        print_status("    WARNING: step_006 output not found; cannot compute ephemeris systematic", "WARNING")

    # ------------------------------------------------------------------
    # 2. Atmospheric systematic: monthly mean variation
    # ------------------------------------------------------------------
    df_temp = df.copy()
    df_temp['detrended'] = detrended
    # Julian day 2451545.0 = 2000-01-01 noon; derive month crudely
    df_temp['month'] = (np.floor((df_temp['date_julian'] - 2451545.0) % 365.25 / 30.44).astype(int) % 12) + 1
    monthly_means = df_temp.groupby('month')['detrended'].mean()
    # Map monthly means back to each observation
    month_map = df_temp.set_index('month')['detrended'].index.map(monthly_means)
    # Actually, we need the systematic component = monthly_mean for each obs
    df_temp['monthly_sys'] = df_temp['month'].map(monthly_means)
    atmos_sys = df_temp['monthly_sys'].values
    atmos_raw_rms = float(np.std(atmos_sys))
    atmos_bias_A, atmos_bias_eta, atmos_r, atmos_p = compute_systematic_projection(
        residuals, cos_elong, atmos_sys)

    # ------------------------------------------------------------------
    # 3. Instrumental systematic: station mean differences
    # ------------------------------------------------------------------
    station_means = df_temp.groupby('station')['detrended'].mean()
    df_temp['station_sys'] = df_temp['station'].map(station_means)
    inst_sys = df_temp['station_sys'].values
    inst_raw_rms = float(np.std(inst_sys))
    inst_bias_A, inst_bias_eta, inst_r, inst_p = compute_systematic_projection(
        residuals, cos_elong, inst_sys)

    # ------------------------------------------------------------------
    # 4. Tidal systematic: cos(2*elongation) harmonic
    # ------------------------------------------------------------------
    cos_2elong = np.cos(2.0 * elongation)
    X_tidal = np.column_stack([cos_2elong, np.ones(n)])
    tidal_coeffs, _, _, _ = stable_lstsq(X_tidal, residuals)
    tidal_sys = tidal_coeffs[0] * cos_2elong
    tidal_raw_rms = float(np.std(tidal_sys))
    tidal_bias_A, tidal_bias_eta, tidal_r, tidal_p = compute_systematic_projection(
        residuals, cos_elong, tidal_sys)

    # ------------------------------------------------------------------
    # 5. Thermal systematic: diurnal (24-hr) sinusoidal amplitude
    # ------------------------------------------------------------------
    hour_frac = ((jd - 0.5) % 1.0) * 24.0  # hours from midnight
    omega = 2.0 * np.pi / 24.0
    X_thermal = np.column_stack([np.cos(omega * hour_frac),
                                  np.sin(omega * hour_frac),
                                  np.ones(n)])
    thermal_coeffs, _, _, _ = stable_lstsq(X_thermal, detrended)
    thermal_sys = thermal_coeffs[0] * np.cos(omega * hour_frac) + thermal_coeffs[1] * np.sin(omega * hour_frac)
    thermal_raw_rms = float(np.std(thermal_sys))
    thermal_bias_A, thermal_bias_eta, thermal_r, thermal_p = compute_systematic_projection(
        residuals, cos_elong, thermal_sys)

    # ------------------------------------------------------------------
    # Compile systematic bias table
    # ------------------------------------------------------------------
    projections = {
        "ephemeris": {
            "source": "Cross-ephemeris residual difference (INPOP19a − DE430)",
            "raw_rms_m": float(ephem_bias_A),
            "bias_A_m": float(ephem_bias_A),
            "projection_efficiency": 1.0,
            "bias_eta": float(ephem_bias_eta),
            "r_cos_elong": float(ephem_r),
            "p_value": float(ephem_p),
            "available": ephem_available,
            "note": "Direct bias: ephemeris difference IS the systematic"
        },
        "atmospheric": {
            "source": "Monthly mean variation (seasonal tropospheric delay)",
            "raw_rms_m": float(atmos_raw_rms),
            "bias_A_m": float(atmos_bias_A),
            "projection_efficiency": float(abs(atmos_bias_A) / atmos_raw_rms) if atmos_raw_rms > 0 else 0.0,
            "bias_eta": float(atmos_bias_eta),
            "r_cos_elong": float(atmos_r),
            "p_value": float(atmos_p),
            "available": True,
            "note": "Annual cycle is incommensurate with synodic period; projection should be small"
        },
        "instrumental": {
            "source": "Station-to-station mean differences (hardware calibration)",
            "raw_rms_m": float(inst_raw_rms),
            "bias_A_m": float(inst_bias_A),
            "projection_efficiency": float(abs(inst_bias_A) / inst_raw_rms) if inst_raw_rms > 0 else 0.0,
            "bias_eta": float(inst_bias_eta),
            "r_cos_elong": float(inst_r),
            "p_value": float(inst_p),
            "available": True,
            "note": "Constant offsets per station are orthogonal to cos(elongation) over full cycle"
        },
        "tidal": {
            "source": "cos(2*elongation) harmonic (solid Earth / ocean tides)",
            "raw_rms_m": float(tidal_raw_rms),
            "bias_A_m": float(tidal_bias_A),
            "projection_efficiency": float(abs(tidal_bias_A) / tidal_raw_rms) if tidal_raw_rms > 0 else 0.0,
            "bias_eta": float(tidal_bias_eta),
            "r_cos_elong": float(tidal_r),
            "p_value": float(tidal_p),
            "available": True,
            "note": "cos(2D) is mathematically orthogonal to cos(D) over [0, 2π]"
        },
        "thermal": {
            "source": "Diurnal (24-hr) sinusoidal amplitude (thermal expansion)",
            "raw_rms_m": float(thermal_raw_rms),
            "bias_A_m": float(thermal_bias_A),
            "projection_efficiency": float(abs(thermal_bias_A) / thermal_raw_rms) if thermal_raw_rms > 0 else 0.0,
            "bias_eta": float(thermal_bias_eta),
            "r_cos_elong": float(thermal_r),
            "p_value": float(thermal_p),
            "available": True,
            "note": "Diurnal (1 day) and synodic (29.5 day) are incommensurate; averages to zero"
        }
    }

    # ------------------------------------------------------------------
    # Combined projected systematic: quadrature sum of absolute biases
    # Only include sources with |r| > 0.01 (negligible projection excluded)
    # ------------------------------------------------------------------
    significant_biases = [
        abs(v["bias_eta"]) for v in projections.values()
        if v["available"] and abs(v["r_cos_elong"]) > 0.01
    ]
    combined_projected_bias = float(np.sqrt(np.sum([b**2 for b in significant_biases])))

    # Also compute the linear sum as an upper bound
    linear_sum_bias = float(np.sum([abs(v["bias_eta"]) for v in projections.values() if v["available"]]))

    print_status("═══ SYSTEMATIC PROJECTION RESULTS", "TITLE")
    print_status(f"{'Source':<25} {'bias_η':>12} {'r(cos)':>10} {'p':>10}", "INFO")
    print_status("-" * 60, "INFO")
    for key, val in projections.items():
        if val["available"]:
            print_status(
                f"{val['source'][:24]:<25} {val['bias_eta']:>12.3e} {val['r_cos_elong']:>10.4f} {val['p_value']:>10.2e}",
                "CALC"
            )
    print_status("-" * 60, "INFO")
    print_status(f"{'Combined (quadrature)':<25} {combined_projected_bias:>12.3e}", "CALC")
    print_status(f"{'Upper bound (linear sum)':<25} {linear_sum_bias:>12.3e}", "CALC")

    # ------------------------------------------------------------------
    # Phase-locked differential analysis (cancels common-mode systematics)
    # ------------------------------------------------------------------
    print_status("═══ PHASE-LOCKED DIFFERENTIAL ANALYSIS", "TITLE")
    mask_near_0 = (elongation < ELONGATION_MASK_WIDTH) | (elongation > 2 * np.pi - ELONGATION_MASK_WIDTH)
    mask_near_pi = np.abs(elongation - np.pi) < ELONGATION_MASK_WIDTH

    n_0 = np.sum(mask_near_0)
    n_pi = np.sum(mask_near_pi)

    if n_0 >= 30 and n_pi >= 30:
        res_0 = residuals[mask_near_0]
        res_pi = residuals[mask_near_pi]

        # Mean-difference method (regression within a phase bin is ill-conditioned
        # because cos(elongation) ≈ constant). In new moon bin cos≈+1, full moon
        # cos≈-1. Model: residual = A*cos(D) + intercept + noise.
        # mean_0 ≈ A * μ_bin + intercept, mean_π ≈ -A * μ_bin + intercept.
        # mean_0 - mean_π = 2A * μ_bin = 2 * ETA_SCALE_FACTOR * η * μ_bin.
        # where μ_bin = sin(w)/w is the mean of cos(D) over [-w, w].
        # The intercept cancels. This cancels all common-mode systematics.
        w = ELONGATION_MASK_WIDTH
        bin_mean_cos = np.sin(w) / w  # exact mean of cos(x) over [-w, w]
        mean_0 = np.mean(res_0)
        mean_pi = np.mean(res_pi)
        sem_0 = np.std(res_0, ddof=1) / np.sqrt(n_0)
        sem_pi = np.std(res_pi, ddof=1) / np.sqrt(n_pi)
        sem_diff = np.sqrt(sem_0**2 + sem_pi**2)

        A_diff = mean_0 - mean_pi  # = 2A * μ_bin = 2 * 13 * η * μ_bin
        eta_diff = A_diff / (2.0 * ETA_SCALE_FACTOR * bin_mean_cos)
        eta_diff_error = sem_diff / (2.0 * ETA_SCALE_FACTOR * bin_mean_cos)
        snr_diff = abs(eta_diff) / eta_diff_error if eta_diff_error > 0 else 0.0

        # Null test: random-phase subsets
        np.random.seed(TEP_CONFIG.get("RANDOM_SEED", 42))
        n_perm = 1000
        perm_eta_diffs = []
        for _ in range(n_perm):
            perm_mask_0 = np.random.choice(n, n_0, replace=False)
            perm_mask_pi = np.random.choice(n, n_pi, replace=False)
            perm_mean_0 = np.mean(residuals[perm_mask_0])
            perm_mean_pi = np.mean(residuals[perm_mask_pi])
            perm_eta = (perm_mean_0 - perm_mean_pi) / (2.0 * ETA_SCALE_FACTOR * bin_mean_cos)
            perm_eta_diffs.append(perm_eta)

        perm_eta_diffs = np.array(perm_eta_diffs)
        p_diff = np.mean(np.abs(perm_eta_diffs) >= abs(eta_diff))

        print_status(f"    New moon:  N={n_0},  mean={mean_0:.3e} ± {sem_0:.3e} m", "CALC")
        print_status(f"    Full moon: N={n_pi}, mean={mean_pi:.3e} ± {sem_pi:.3e} m", "CALC")
        print_status(f"    Differential η = {eta_diff:.3e} ± {eta_diff_error:.3e} ({snr_diff:.2f}σ)", "CALC")
        print_status(f"    Permutation p = {p_diff:.4f} (n={n_perm})", "CALC")

        differential = {
            "n_new_moon": int(n_0),
            "n_full_moon": int(n_pi),
            "mean_new_moon_m": float(mean_0),
            "mean_new_moon_sem_m": float(sem_0),
            "mean_full_moon_m": float(mean_pi),
            "mean_full_moon_sem_m": float(sem_pi),
            "amplitude_differential_m": float(A_diff),
            "amplitude_differential_error_m": float(sem_diff),
            "eta_differential": float(eta_diff),
            "eta_differential_error": float(eta_diff_error),
            "snr_differential": float(snr_diff),
            "permutation_p_value": float(p_diff),
            "n_permutations": n_perm,
            "method": "mean_difference (common-mode systematic cancellation)"
        }
    else:
        print_status(f"    WARNING: Insufficient data in phase bins (n_0={n_0}, n_pi={n_pi})", "WARNING")
        differential = None

    # ------------------------------------------------------------------
    # Final combined uncertainty
    # ------------------------------------------------------------------
    statistical_error = eta_err_base
    # Combine ephemeris scatter and projected systematics in quadrature.
    # The ephemeris scatter (from step_006 cross-ephemeris std) and the
    # cos(elongation)-projected non-ephemeris biases are independent sources.
    # For N=2 ephemerides, std(η_ephem) = |η1 - η2|/√2.
    # The projected systematic is the quadrature sum of atmospheric, instrumental,
    # tidal, and thermal cos(elongation)-correlated biases.
    if ephem_available:
        systematic_error = float(np.sqrt(ephem_bias_eta**2 + combined_projected_bias**2))
        sys_source = "quadrature(ephemeris_scatter + projected_non_ephemeris)"
    else:
        systematic_error = combined_projected_bias
        sys_source = "cos(elongation)_projected_quadrature"

    total_error = np.sqrt(statistical_error**2 + systematic_error**2)
    total_snr = abs(eta_base) / total_error if total_error > 0 else 0.0

    print_status("═══ COMBINED UNCERTAINTY", "TITLE")
    print_status(f"    Statistical:   ±{statistical_error:.3e}", "CALC")
    print_status(f"    Systematic:    ±{systematic_error:.3e} ({sys_source})", "CALC")
    print_status(f"    Total:         ±{total_error:.3e}", "CALC")
    print_status(f"    SNR (total):   {total_snr:.2f}σ", "CALC")

    results = {
        "step_id": "step_044",
        "baseline_eta": float(eta_base),
        "baseline_eta_error": float(eta_err_base),
        "baseline_snr": float(snr_base),
        "systematic_projections": projections,
        "combined_projected_systematic_eta": combined_projected_bias,
        "linear_sum_upper_bound": linear_sum_bias,
        "systematic_error_used": float(systematic_error),
        "systematic_error_source": sys_source,
        "statistical_error": float(statistical_error),
        "total_uncertainty": float(total_error),
        "total_snr": float(total_snr),
        "phase_locked_differential": differential,
        "status": "PASS"
    }

    return results


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Step 044: Systematic Projection Analysis")
    args = parser.parse_args()

    log_dir = PROJECT_ROOT / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    logger = TEPLogger("step_044", str(
        log_dir / "step_044_systematic_projection_analysis.log"))
    set_step_logger(logger)
    set_verbose_mode(True)

    print_status("Starting Systematic Projection Analysis...", "TITLE")

    input_path = PROJECT_ROOT / 'data' / 'processed' / 'INPOP19a_all_stations_residuals.csv'
    if not input_path.exists():
        print_status(f"Input file not found: {input_path}", "ERROR")
        sys.exit(1)

    df = pd.read_csv(input_path)
    summary = run_systematic_projection_analysis(df, verbose=True)

    logger.save_step_results(summary, PROJECT_ROOT,
                             "step_044_systematic_projection_analysis")
    print_status("Systematic Projection Analysis Complete.", "SUCCESS")
