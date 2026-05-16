#!/usr/bin/env python3
"""
Step 008: Systematic Error Analysis for TEP Nordtvedt Signal Detection
Enhanced with comprehensive systematic error budget table (DATA-DRIVEN)
"""

import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

import argparse
import json
import pandas as pd
import numpy as np
from scripts.utils.numerics import stable_lstsq

from scripts.utils.llr_constants import ETA_SCALE_FACTOR
from scripts.utils.logger import TEPLogger, set_step_logger, set_verbose_mode, print_status

def analyze_systematics(df, verbose=False, logger=None):
    stations = df['station'].unique()
    systematics = {}

    print_status(f"Analyzing {len(stations)} stations for systematic errors...", "INFO")
    print_status(f"Total observations: {len(df):,}", "INFO")
    print_status(f"Residual threshold for 'clean' status: 0.15 m", "INFO")
    print_status("", "INFO")

    clean_count = 0
    noisy_count = 0

    for s in stations:
        station_data = df[df['station'] == s]
        rms = station_data['residual_m'].std()
        is_clean = bool(rms < 0.15)
        systematics[s] = {
            "n_obs": len(station_data),
            "residual_rms": float(rms),
            "clean_profile": is_clean
        }

        # Always log station analysis (not just verbose)
        status_icon = "✓" if is_clean else "⚠"
        status_label = "CLEAN" if is_clean else "NOISY"
        print_status(f"{status_icon} Station {s}:", "PASS" if is_clean else "WARNING")
        print_status(f"    N obs     = {len(station_data):,}", "INFO")
        print_status(f"    RMS       = {rms:.4f} m", "INFO")
        print_status(f"    Status    = {status_label}", "PASS" if is_clean else "WARNING")

        if is_clean:
            clean_count += 1
        else:
            noisy_count += 1

    print_status("", "INFO")
    print_status(f"Systematic Analysis Summary:", "TITLE")
    print_status(f"  Clean stations:  {clean_count}/{len(stations)}", "PASS" if clean_count == len(stations) else "WARNING")
    print_status(f"  Noisy stations:  {noisy_count}/{len(stations)}", "INFO" if noisy_count == 0 else "WARNING")
    print_status(f"  Overall status:  {'PASS' if all(v.get('clean_profile', False) for v in systematics.values()) else 'WARNING'}", "PASS" if all(v.get('clean_profile', False) for v in systematics.values()) else "WARNING")

    return systematics

def correlated_systematic_amplitude(systematic_signal, cos_elong):
    """
    Compute the component of a systematic signal that is correlated with
    cos(elongation). Only this component can bias the TEP eta estimate;
    the orthogonal component contributes noise (already in statistical error).
    Returns the amplitude in metres of the cos(elongation)-correlated bias.
    """
    X = np.column_stack([cos_elong, np.ones(len(cos_elong))])
    coeffs, _, _, _ = stable_lstsq(X, systematic_signal)
    return float(abs(coeffs[0]))


def generate_systematic_error_budget(df, systematics, verbose=False, logger=None):
    """
    Generate comprehensive systematic error budget table.
    CRITICAL FIX v2: Quantifies only the ELONGATION-CORRELATED component
    of each systematic.  The total amplitude of a systematic (e.g. diurnal
    thermal expansion) is irrelevant if it is orthogonal to cos(elongation).
    Only the projection onto cos(elongation) can bias eta.
    """
    print_status("", "INFO")
    print_status("Generating Data-Driven Systematic Error Budget (v2: elongation-correlated)...", "TITLE")
    print_status("", "INFO")


    global_rms = df['residual_m'].std()
    global_mean = df['residual_m'].mean()
    n_total = len(df)

    # ------------------------------------------------------------------
    # Detrend residuals: remove the best-fit TEP cos(elongation) signal
    # so that the remaining variance isolates systematic sources.
    # ------------------------------------------------------------------
    cos_elong = np.cos(df['elongation_rad'].values)
    residuals = df['residual_m'].values
    X = np.column_stack([cos_elong, np.ones(len(cos_elong))])
    coeffs_tep, _, _, _ = stable_lstsq(X, residuals)
    detrended = residuals - coeffs_tep[0] * cos_elong  # m
    eta_fit = coeffs_tep[0] / ETA_SCALE_FACTOR

    # ------------------------------------------------------------------
    # 1. Ephemeris modeling uncertainty
    #    Estimated from cross-ephemeris scatter in eta (Step 006 output).
    #    With only two ephemerides available, we take half the absolute
    #    difference as a conservative upper bound.
    # ------------------------------------------------------------------
    ephem_json = PROJECT_ROOT / 'results/outputs/step_006_multi_ephemeris_comparison.json'
    ephem_uncertainty_cm = None
    n_ephem = 0
    if ephem_json.exists():
        with open(ephem_json, 'r') as f:
            ephem_data = json.load(f)
        etas = [v['eta'] for v in ephem_data.get('comparisons', {}).values() if 'eta' in v]
        n_ephem = len(etas)
        if n_ephem >= 2:
            ephem_eta_std = float(np.std(etas))
            ephem_uncertainty_cm = ephem_eta_std * ETA_SCALE_FACTOR * 100.0
        elif n_ephem == 1:
            ephem_uncertainty_cm = abs(etas[0]) * ETA_SCALE_FACTOR * 100.0 * 0.5
        else:
            raise RuntimeError("step_006_multi_ephemeris_comparison.json contains no ephemeris eta values; cannot compute data-driven ephemeris uncertainty.")
    else:
        raise FileNotFoundError(
            f"Required upstream data not found: {ephem_json}. "
            "Run step_006_multi_ephemeris_comparison.py first."
        )

    # ------------------------------------------------------------------
    # 2. Atmospheric delay modeling uncertainty
    #    v2 FIX: The total seasonal variation in detrended residuals includes
    #    genuine TEP heliocentric modulation (Jan-Jul sign flip confirmed in
    #    step_045).  We model the seasonal component and use only its
    #    ELONGATION-CORRELATED projection as the systematic bias.  The
    #    heliocentric TEP component is signal, not error.
    # ------------------------------------------------------------------
    df_temp = df.copy()
    df_temp['detrended'] = detrended
    df_temp['month'] = (np.floor((df_temp['date_julian'] - 2451545.0) % 365.25 / 30.44).astype(int) % 12) + 1
    monthly_means = df_temp.groupby('month')['detrended'].mean()

    # Build a seasonal model signal (monthly mean assigned to each obs)
    seasonal_model = df_temp['month'].map(monthly_means).values
    atmos_bias_m = correlated_systematic_amplitude(seasonal_model, cos_elong)
    atmos_uncertainty_cm = float(atmos_bias_m * 100.0)

    # ------------------------------------------------------------------
    # 3. Instrumental systematic uncertainty
    #    v2 FIX: Different hardware produces station-dependent biases, but
    #    only the component correlated with cos(elongation) can bias eta.
    #    We model each station's mean bias and project onto cos(elongation).
    # ------------------------------------------------------------------
    station_means = df_temp.groupby('station')['detrended'].mean()
    powered_stations = [s for s in station_means.index
                        if systematics.get(s, {}).get('n_obs', 0) >= 1000]
    if len(powered_stations) >= 2:
        inst_model = df_temp['station'].map(station_means[powered_stations]).fillna(0).values
    else:
        inst_model = df_temp['station'].map(station_means).fillna(0).values
    inst_bias_m = correlated_systematic_amplitude(inst_model, cos_elong)
    inst_uncertainty_cm = float(inst_bias_m * 100.0)

    # ------------------------------------------------------------------
    # 4. Tidal modeling uncertainty
    #    v2 FIX: Tides produce perturbations at TWICE the synodic frequency
    #    (cos(2*elongation)), which is ORTHOGONAL to cos(elongation) over a
    #    full period.  Only the small leakage onto cos(elongation) can bias
    #    eta.  We fit the tidal harmonic and project onto cos(elongation).
    # ------------------------------------------------------------------
    cos_2elong = np.cos(2.0 * df['elongation_rad'].values)
    X_tidal = np.column_stack([cos_2elong, np.ones(len(cos_2elong))])
    tidal_coeffs, _, _, _ = stable_lstsq(X_tidal, residuals)
    tidal_model = tidal_coeffs[0] * cos_2elong
    tidal_bias_m = correlated_systematic_amplitude(tidal_model, cos_elong)
    tidal_uncertainty_cm = float(tidal_bias_m * 100.0)

    # ------------------------------------------------------------------
    # 5. Thermal expansion uncertainty
    #    v2 FIX: Diurnal (24-hr) thermal expansion is orthogonal to
    #    cos(elongation) if observation times are uniformly distributed
    #    across lunar phases.  Only the elongation-correlated projection
    #    can bias eta.  We fit the diurnal model and project.
    # ------------------------------------------------------------------
    hour_frac = ((df['date_julian'].values - 0.5) % 1.0) * 24.0
    omega = 2.0 * np.pi / 24.0
    X_thermal = np.column_stack([np.cos(omega * hour_frac),
                                  np.sin(omega * hour_frac),
                                  np.ones(len(hour_frac))])
    thermal_coeffs, _, _, _ = stable_lstsq(X_thermal, detrended)
    thermal_model = (thermal_coeffs[0] * np.cos(omega * hour_frac) +
                     thermal_coeffs[1] * np.sin(omega * hour_frac))
    thermal_bias_m = correlated_systematic_amplitude(thermal_model, cos_elong)
    thermal_uncertainty_cm = float(thermal_bias_m * 100.0)

    # ------------------------------------------------------------------
    # Assemble budget
    # ------------------------------------------------------------------
    error_budget = {
        "ephemeris_modeling": {
            "source": "Ephemeris modeling (cross-ephemeris scatter)",
            "magnitude_cm": round(ephem_uncertainty_cm, 2),
            "description": "Uncertainty in lunar and planetary ephemeris fitting, derived from scatter of eta across INPOP19a and DE430",
            "reference": f"Data-driven from {n_ephem} ephemerides (step_006_multi_ephemeris_comparison)",
            "method": "std(eta_ephem) * 13.0 * 100",
            "data_driven": True
        },
        "atmospheric_delay": {
            "source": "Atmospheric delay modeling",
            "magnitude_cm": round(atmos_uncertainty_cm, 2),
            "description": "Tropospheric delay correction uncertainties: seasonal model projected onto cos(elongation). Only correlated component counts.",
            "reference": "Data-driven from monthly mean residual variation, elongation-correlated projection",
            "method": "correlated_systematic_amplitude(monthly_model, cos_elong) * 100",
            "data_driven": True
        },
        "instrumental": {
            "source": "Instrumental systematic",
            "magnitude_cm": round(inst_uncertainty_cm, 2),
            "description": "Detector calibration and timing electronics: station-bias model projected onto cos(elongation).",
            "reference": "Data-driven from powered-station mean residual variation, elongation-correlated projection",
            "method": "correlated_systematic_amplitude(station_bias_model, cos_elong) * 100",
            "data_driven": True
        },
        "tidal_modeling": {
            "source": "Tidal modeling",
            "magnitude_cm": round(tidal_uncertainty_cm, 2),
            "description": "Solid Earth and ocean tide model uncertainties: cos(2*elongation) harmonic projected onto cos(elongation). Orthogonal component excluded.",
            "reference": "Data-driven from 2nd-synodic-harmonic residual amplitude, elongation-correlated projection",
            "method": "correlated_systematic_amplitude(tidal_model, cos_elong) * 100",
            "data_driven": True
        },
        "thermal_expansion": {
            "source": "Thermal expansion",
            "magnitude_cm": round(thermal_uncertainty_cm, 2),
            "description": "Telescope and retroreflector thermal expansion: diurnal model projected onto cos(elongation). Orthogonal component excluded.",
            "reference": "Data-driven from diurnal residual variation, elongation-correlated projection",
            "method": "correlated_systematic_amplitude(diurnal_model, cos_elong) * 100",
            "data_driven": True
        }
    }

    total_budget_cm = sum(v["magnitude_cm"] for v in error_budget.values())
    for key in error_budget:
        error_budget[key]["percentage"] = round((error_budget[key]["magnitude_cm"] / total_budget_cm) * 100, 1)
        error_budget[key]["variance_contribution"] = (error_budget[key]["magnitude_cm"] / 100.0) ** 2

    observed_variance = (global_rms / 100.0) ** 2
    budget_variance = sum(v["variance_contribution"] for v in error_budget.values())
    unexplained_variance = max(0.0, observed_variance - budget_variance)
    unexplained_rms_cm = float(np.sqrt(unexplained_variance) * 100) if unexplained_variance > 0 else 0.0

    error_budget["unexplained"] = {
        "source": "Unexplained residual variance",
        "magnitude_cm": round(unexplained_rms_cm, 2),
        "description": "Residual variance not accounted for by known systematics",
        "percentage": round((unexplained_rms_cm / total_budget_cm) * 100, 1) if total_budget_cm > 0 else 0.0,
        "variance_contribution": unexplained_variance,
        "data_driven": True
    }

    # Summary statistics
    total_systematic_cm = float(np.sqrt(budget_variance) * 100)

    # Print error budget table
    print_status("Systematic Error Budget Table (DATA-DRIVEN):", "TITLE")
    print_status(f"{'Source':<40} {'Magnitude (cm)':<18} {'% Contribution':<15}", "INFO")
    print_status("-" * 75, "INFO")
    for key, value in error_budget.items():
        print_status(f"{value['source']:<40} {value['magnitude_cm']:<18.2f} {value['percentage']:<15.1f}", "INFO")
    print_status("-" * 75, "INFO")
    print_status(f"{'Total Budget (quadrature)':<40} {total_systematic_cm:<18.2f} {'-':<15}", "INFO")
    print_status(f"{'Total Budget (linear sum)':<40} {total_budget_cm:<18.2f} {'-':<15}", "INFO")
    print_status(f"{'Observed RMS':<40} {global_rms * 100:<18.2f} {'-':<15}", "INFO")
    print_status("", "INFO")

    return error_budget, total_systematic_cm, eta_fit

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Step 008: Systematic Error Analysis with Error Budget")
    args = parser.parse_args()

    log_dir = PROJECT_ROOT / "logs"
    logger = TEPLogger("step_008", str(
        log_dir / "step_008_systematic_error_analysis.log"))
    set_step_logger(logger)

    print_status("Starting Systematic Error Analysis...", "TITLE")

    input_path = PROJECT_ROOT / 'data/processed/INPOP19a_all_stations_residuals.csv'
    if not input_path.exists():
        print_status(f"Input file not found: {input_path}", "ERROR")
        sys.exit(1)

    df = pd.read_csv(input_path)
    print_status(f"Loaded {len(df):,} observations from {input_path.relative_to(PROJECT_ROOT)}", "INFO")
    print_status("", "INFO")

    sys_results = analyze_systematics(df)
    error_budget, total_systematic_cm, eta_fit = generate_systematic_error_budget(df, sys_results, logger=logger)

    all_clean = all(v.get("clean_profile", False)
                    for v in sys_results.values())

    results = {
        "step_id": "step_008",
        "systematics": sys_results,
        "error_budget": error_budget,
        "total_systematic_cm": total_systematic_cm,
        "eta_fit": float(eta_fit),
        "n_observations": len(df),
        "data_driven": True,
        "status": "PASS" if all_clean else "WARNING"
    }

    logger.save_step_results(results, PROJECT_ROOT,
                             "step_008_systematic_error_analysis")
    print_status("Systematic Error Analysis Complete.", "SUCCESS")
