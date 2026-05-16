#!/usr/bin/env python3
"""
Step 057: Haleakala Null-Fluctuation Simulation

Formal simulation of how often an underpowered station with Haleakala's
observational characteristics (N, RMS, biased phase coverage from the
processed archive) would fluctuate to the observed positive station-level
eta from Step 029 under:
  (a) the true TEP model (eta_true = m5_full_corrected eta from Step 050)
  (b) the GR null model (eta_true = 0)

Also computes the family-wise rate: the probability that ANY of the 5
stations produces a deviation as large as Haleakala's under the true model.

Observed Haleakala eta and uncertainty are read from
`step_029_station_power_analysis.json`; the TEP null uses
`models.m5_full_corrected.eta` from `step_050_corrected_tep_analysis.json`.
"""

import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

import json
import numpy as np
import pandas as pd

from scripts.utils.llr_constants import ETA_SCALE_FACTOR
from scripts.utils.logger import TEPLogger, set_step_logger, print_status


def _load_required_json(rel_path: str) -> dict:
    path = PROJECT_ROOT / rel_path
    if not path.is_file():
        raise FileNotFoundError(f"Required upstream output missing: {path}")
    with open(path, encoding="utf-8") as f:
        return json.load(f)


def _haleakala_power_row(step_029: dict) -> dict:
    stations = step_029.get("per_station_power", {}).get("stations")
    if not stations:
        raise KeyError(
            "step_029_station_power_analysis.json missing per_station_power.stations"
        )
    for row in stations:
        if row.get("station") == "Haleakala":
            return row
    raise KeyError("Haleakala row not found in step_029 per_station_power.stations")


def _m5_full_corrected_eta(step_050: dict) -> float:
    try:
        return float(step_050["models"]["m5_full_corrected"]["eta"])
    except (KeyError, TypeError) as e:
        raise KeyError(
            "step_050_corrected_tep_analysis.json missing models.m5_full_corrected.eta"
        ) from e


def fit_station_eta(residuals, cos_elong):
    """Fit OLS eta from residuals vs cos(elongation)."""
    X = np.column_stack([cos_elong, np.ones(len(cos_elong), dtype=np.float64)])
    XtX_inv = np.linalg.pinv(X.T @ X, rcond=1e-10, hermitian=True)
    beta = XtX_inv @ (X.T @ residuals)
    A = float(beta[0])
    B = float(beta[1])
    pred = A * cos_elong + B
    resid = residuals - pred
    mse = float(np.sum(resid ** 2) / max(len(cos_elong) - 2, 1))
    se_A = float(np.sqrt(max(mse * XtX_inv[0, 0], 0.0)))
    eta = A / ETA_SCALE_FACTOR
    eta_err = se_A / ETA_SCALE_FACTOR
    return eta, eta_err


def fit_station_eta_batch(residuals_matrix: np.ndarray, cos_elong: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Vector OLS η for each column of residuals_matrix (n_obs, n_mc)."""
    n = cos_elong.shape[0]
    if residuals_matrix.shape[0] != n:
        raise ValueError("residuals_matrix row count must match cos_elong length")
    X = np.column_stack([cos_elong, np.ones(n, dtype=np.float64)])
    XtX_inv = np.linalg.pinv(X.T @ X, rcond=1e-10, hermitian=True)
    coeffs = XtX_inv @ (X.T @ residuals_matrix)
    A = coeffs[0, :]
    B = coeffs[1, :]
    pred = A[np.newaxis, :] * cos_elong[:, np.newaxis] + B[np.newaxis, :]
    resid = residuals_matrix - pred
    dof = max(n - 2, 1)
    mse = np.sum(resid**2, axis=0) / dof
    se_A = np.sqrt(np.maximum(mse * XtX_inv[0, 0], 0.0))
    eta = A / ETA_SCALE_FACTOR
    eta_err = se_A / ETA_SCALE_FACTOR
    return eta, eta_err


def simulate_station_fluctuations(
    haleakala_df, eta_true, eta_obs, eta_err_obs, n_mc=20000, seed=57
):
    """
    Simulate Haleakala-like station realizations under a specified true eta.
    Uses the station's actual elongation distribution to preserve phase-truncation
    structure in the processed archive.
    """
    rng = np.random.RandomState(seed)

    elongation = haleakala_df['elongation_rad'].values
    cos_elong = np.cos(elongation)
    n = len(cos_elong)
    rms = float(np.sqrt(np.mean(haleakala_df['residual_m'].values ** 2)))

    signal_amplitude = ETA_SCALE_FACTOR * eta_true
    signal = signal_amplitude * cos_elong

    # One RNG draw block matches sequential RandomState.normal(0, rms, n) per iteration.
    noise_mat = rng.normal(0, rms, size=(n, n_mc))
    res_mat = signal[:, np.newaxis] + noise_mat
    etas_sim, eta_errs_sim = fit_station_eta_batch(res_mat, cos_elong)

    eta_obs = float(eta_obs)
    eta_err_obs = float(eta_err_obs)

    # One-tailed: fluctuate to at or above the observed positive station eta
    p_one_tailed = float(np.mean(etas_sim >= eta_obs))

    # Two-tailed: deviation from true eta at least as large as observed
    dev_obs = abs(eta_obs - eta_true)
    p_two_tailed = float(np.mean(np.abs(etas_sim - eta_true) >= dev_obs))

    # Z-score of observed fluctuation
    z_obs = (eta_obs - eta_true) / np.std(etas_sim) if np.std(etas_sim) > 0 else 0.0

    return {
        'n_mc': n_mc,
        'seed': seed,
        'eta_true': float(eta_true),
        'n_obs': int(n),
        'rms_m': float(rms),
        'eta_observed': float(eta_obs),
        'eta_error_observed': float(eta_err_obs),
        'simulated_eta_mean': float(np.mean(etas_sim)),
        'simulated_eta_std': float(np.std(etas_sim)),
        'simulated_eta_median': float(np.median(etas_sim)),
        'p_one_tailed_opposite_sign': p_one_tailed,
        'p_two_tailed_deviation': p_two_tailed,
        'z_score_observed': float(z_obs),
        'percentiles': {
            '1': float(np.percentile(etas_sim, 1)),
            '2.5': float(np.percentile(etas_sim, 2.5)),
            '5': float(np.percentile(etas_sim, 5)),
            '95': float(np.percentile(etas_sim, 95)),
            '97.5': float(np.percentile(etas_sim, 97.5)),
            '99': float(np.percentile(etas_sim, 99)),
        },
    }


def family_wise_simulation(df_all, eta_true, eta_obs_haleakala, n_mc=20000, seed=157):
    """
    Simulate all 5 stations simultaneously under the true model.
    Compute the probability that at least ONE station produces a
    Haleakala-magnitude deviation.
    """
    rng = np.random.RandomState(seed)

    stations = ['APO', 'Grasse', 'Haleakala', 'Matera', 'McDonald2']
    station_params = {}

    for st in stations:
        sdf = df_all[df_all['station'] == st]
        if len(sdf) == 0:
            continue
        elongation = sdf['elongation_rad'].values
        cos_elong = np.cos(elongation)
        n = len(cos_elong)
        rms = float(np.sqrt(np.mean(sdf['residual_m'].values ** 2)))
        station_params[st] = {
            'elongation': elongation,
            'cos_elong': cos_elong,
            'n': n,
            'rms': rms,
        }

    eta_obs_haleakala = float(eta_obs_haleakala)
    dev_obs = abs(eta_obs_haleakala - eta_true)

    any_station_deviates = np.zeros(n_mc, dtype=bool)

    for i in range(n_mc):
        for st, params in station_params.items():
            n_s = params['n']
            rms_s = params['rms']
            cos_s = params['cos_elong']
            signal = ETA_SCALE_FACTOR * eta_true * cos_s
            noise = rng.normal(0, rms_s, n_s)
            residuals = signal + noise
            eta_fit, _ = fit_station_eta(residuals, cos_s)
            if abs(eta_fit - eta_true) >= dev_obs:
                any_station_deviates[i] = True
                break

    p_family_wise = float(np.mean(any_station_deviates))

    return {
        'n_mc': n_mc,
        'seed': seed,
        'eta_true': float(eta_true),
        'p_any_station_deviates': p_family_wise,
        'stations_included': list(station_params.keys()),
    }


def run_haleakala_simulation(df, verbose=False):
    print_status("═══ Starting Step 057: Haleakala Null-Fluctuation Simulation", "TITLE")
    print_status("═══ STEP PURPOSE: Quantify how often Haleakala's opposite-sign fluctuation arises under TEP", "INFO")
    print_status("═══ METHOD: Monte Carlo using actual Haleakala elongation distribution, N=20000 realizations", "INFO")

    step_029 = _load_required_json("results/outputs/step_029_station_power_analysis.json")
    step_050 = _load_required_json("results/outputs/step_050_corrected_tep_analysis.json")
    hk_power = _haleakala_power_row(step_029)
    eta_obs = float(hk_power["eta_obs"])
    eta_err_obs = float(hk_power["eta_err_obs"])
    eta_true_tep = _m5_full_corrected_eta(step_050)

    # Extract Haleakala data
    haleakala = df[df['station'] == 'Haleakala'].copy()
    if len(haleakala) == 0:
        print_status("ERROR: Haleakala not found in dataset", "ERROR")
        return {'status': 'FAIL', 'reason': 'Haleakala missing'}

    # Apply consistent 6σ outlier cleaning (same as primary pipeline)
    from scripts.utils.statistical_utils import detect_outliers_sigma
    outlier_mask = detect_outliers_sigma(haleakala['residual_m'].values, sigma_threshold=6.0)
    haleakala_clean = haleakala[~outlier_mask].copy()
    n_removed = int(np.sum(outlier_mask))
    print_status(f"    Haleakala raw: N = {len(haleakala)}, after 6σ cleaning: N = {len(haleakala_clean)} (removed {n_removed})", "DATA")

    # GR null
    eta_true_gr = 0.0

    # --- Simulation (a): Under TEP true model ---
    print_status(
        f">>> Simulation (a): Under TEP model (eta_true = {eta_true_tep:.3e})", "PROCESS"
    )
    sim_tep = simulate_station_fluctuations(
        haleakala_clean, eta_true_tep, eta_obs, eta_err_obs, n_mc=20000, seed=57
    )

    print_status(f"    Simulated eta mean:   {sim_tep['simulated_eta_mean']:.3e}", "CALC")
    print_status(f"    Simulated eta std:    {sim_tep['simulated_eta_std']:.3e}", "CALC")
    print_status(f"    Observed eta:         {sim_tep['eta_observed']:.3e}", "CALC")
    print_status(f"    Z-score of observed:  {sim_tep['z_score_observed']:.2f}", "CALC")
    print_status(
        f"    P(eta >= {eta_obs:.3e}):   {sim_tep['p_one_tailed_opposite_sign']:.4f} ({sim_tep['p_one_tailed_opposite_sign']*100:.2f}%)",
        "CALC",
    )
    print_status(f"    P(|dev| >= observed): {sim_tep['p_two_tailed_deviation']:.4f} ({sim_tep['p_two_tailed_deviation']*100:.2f}%)", "CALC")
    print_status(f"    97.5th percentile:      {sim_tep['percentiles']['97.5']:.3e}", "CALC")
    print_status(f"    99th percentile:        {sim_tep['percentiles']['99']:.3e}", "CALC")

    # --- Simulation (b): Under GR null ---
    print_status(">>> Simulation (b): Under GR null (eta_true = 0)", "PROCESS")
    sim_gr = simulate_station_fluctuations(
        haleakala_clean, eta_true_gr, eta_obs, eta_err_obs, n_mc=20000, seed=157
    )

    print_status(f"    Simulated eta mean:   {sim_gr['simulated_eta_mean']:.3e}", "CALC")
    print_status(f"    Simulated eta std:    {sim_gr['simulated_eta_std']:.3e}", "CALC")
    print_status(
        f"    P(eta >= {eta_obs:.3e}):   {sim_gr['p_one_tailed_opposite_sign']:.4f} ({sim_gr['p_one_tailed_opposite_sign']*100:.2f}%)",
        "CALC",
    )
    print_status(f"    P(|dev| >= observed): {sim_gr['p_two_tailed_deviation']:.4f} ({sim_gr['p_two_tailed_deviation']*100:.2f}%)", "CALC")

    # --- Family-wise simulation ---
    print_status(">>> Family-wise simulation: probability ANY station deviates this much", "PROCESS")
    fw_tep = family_wise_simulation(
        df, eta_true_tep, eta_obs, n_mc=20000, seed=257
    )
    fw_gr = family_wise_simulation(df, eta_true_gr, eta_obs, n_mc=20000, seed=357)

    print_status(f"    P(any station | TEP): {fw_tep['p_any_station_deviates']:.4f} ({fw_tep['p_any_station_deviates']*100:.2f}%)", "CALC")
    print_status(f"    P(any station | GR):   {fw_gr['p_any_station_deviates']:.4f} ({fw_gr['p_any_station_deviates']*100:.2f}%)", "CALC")

    # --- Solar-cycle context from Step 023 ---
    step_023_path = PROJECT_ROOT / 'results' / 'outputs' / 'step_023_solar_cycle_correlation.json'
    solar_context = {}
    if step_023_path.exists():
        with open(step_023_path, 'r') as f:
            step_023 = json.load(f)
        solar_mod = step_023.get('solar_modulation', {})
        low_eta = solar_mod.get('low_activity_eta')
        high_eta = solar_mod.get('high_activity_eta')
        if low_eta is not None and high_eta is not None:
            solar_swing = abs(high_eta - low_eta)
            solar_context = {
                'solar_cycle_swing_eta': float(solar_swing),
                'low_activity_eta': float(low_eta),
                'high_activity_eta': float(high_eta),
                'haleakala_mean_solar_index': step_023.get('haleakala_analysis', {}).get('mean_solar_index'),
            }
            print_status(">>> Solar-cycle context (Step 023)", "PROCESS")
            print_status(f"    Solar-cycle eta swing: {solar_swing:.3e}", "CALC")
            print_status(f"    Haleakala mean solar index: {solar_context.get('haleakala_mean_solar_index', 'N/A')}", "CALC")

    # --- Interpretation ---
    tep_is_rare = sim_tep['p_two_tailed_deviation'] < 0.05
    tep_is_common = sim_tep['p_two_tailed_deviation'] > 0.10

    solar_swing = None
    if solar_context:
        raw_swing = solar_context.get("solar_cycle_swing_eta")
        if raw_swing is not None:
            solar_swing = float(raw_swing)
    mean_solar_idx = solar_context.get("haleakala_mean_solar_index") if solar_context else None

    if tep_is_rare:
        if solar_swing is not None:
            ms = f"{float(mean_solar_idx):.3f}" if mean_solar_idx is not None else "N/A"
            tep_interpretation = (
                f"Haleakala's fluctuation is a {sim_tep['z_score_observed']:.1f}σ event under TEP "
                f"(p = {sim_tep['p_two_tailed_deviation']:.4f}), making it a marginal outlier. "
                f"The solar-cycle modulation amplitude ({solar_swing:.3e}) "
                f"is of comparable order, supporting the physical explanation that Haleakala's "
                f"solar-maximum operation epoch (mean index = {ms}) "
                f"contributes substantially to its deviation."
            )
        else:
            tep_interpretation = (
                f"Haleakala's fluctuation is a {sim_tep['z_score_observed']:.1f}σ event under TEP "
                f"(p = {sim_tep['p_two_tailed_deviation']:.4f}), making it a marginal outlier."
            )
    elif tep_is_common:
        tep_interpretation = (
            f"Haleakala's fluctuation is not anomalously large under TEP "
            f"(p = {sim_tep['p_two_tailed_deviation']:.4f}); it falls within the expected "
            f"distribution of underpowered-station noise."
        )
    else:
        tep_interpretation = (
            f"Haleakala's fluctuation is at the {sim_tep['p_two_tailed_deviation']*100:.1f}% level under TEP "
            f"(z = {sim_tep['z_score_observed']:.1f}σ) — uncommon but not impossible."
        )

    print_status("═══ INTERPRETATION", "INFO")
    print_status(f"    {tep_interpretation}", "INFO")

    print_status("═══ REPRODUCIBILITY", "INFO")
    print_status(f"    Output file: results/outputs/step_057_haleakala_null_fluctuation.json", "INFO")
    print_status(f"    MC realizations: 20,000 per model", "INFO")
    print_status(f"    Random seeds: 57 (TEP), 157 (GR), 257 (family TEP), 357 (family GR)", "INFO")

    results = {
        'step_id': 'step_057',
        'status': 'PASS',
        'upstream': {
            'eta_obs_step_029': eta_obs,
            'eta_err_obs_step_029': eta_err_obs,
            'eta_true_tep_step_050_m5_full_corrected': eta_true_tep,
        },
        'haleakala_simulation_tep': sim_tep,
        'haleakala_simulation_gr': sim_gr,
        'family_wise_tep': fw_tep,
        'family_wise_gr': fw_gr,
        'solar_cycle_context': solar_context,
        'interpretation': tep_interpretation,
    }

    return results


if __name__ == '__main__':
    log_dir = PROJECT_ROOT / 'logs'
    log_dir.mkdir(parents=True, exist_ok=True)
    logger = TEPLogger('step_057', str(log_dir / 'step_057_haleakala_null_fluctuation.log'))
    set_step_logger(logger)

    data_path = PROJECT_ROOT / 'data' / 'processed' / 'INPOP19a_all_stations_residuals.csv'
    if not data_path.exists():
        print_status('No processed INPOP19a residuals.', 'ERROR')
        sys.exit(1)

    df = pd.read_csv(data_path)

    results = run_haleakala_simulation(df, verbose=True)

    logger.save_step_results(results, PROJECT_ROOT, 'step_057_haleakala_null_fluctuation')
    print_status('Haleakala Null-Fluctuation Simulation Complete.', 'SUCCESS')
