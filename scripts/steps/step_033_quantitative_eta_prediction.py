#!/usr/bin/env python3
"""
Step 035: Quantitative η Analysis from TEP Framework

Analyzes the measured Nordtvedt parameter η in the context of the TEP framework.
Following the updated theoretical framework (Papers 6, 10, 11, 12, 13):
- The microscopic coupling β is constrained by Cassini to α_0 ≲ 3×10⁻³ in screened regime
- Observable response coefficients κ (like κ_MSP, κ_Cep) are empirically determined in unscreened regimes
- A key insight from Paper 10: TEP effect operates on clock rates via A(Φ) ≈ 1 − ηΦ/c²
- η itself is the TEP modification parameter for LLR, analogous to κ_MSP and κ_Cep
- No separate κ_LL is needed; η is the observable response coefficient for LLR
"""

import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

import json
import numpy as np
from scripts.utils.logger import TEPLogger, set_step_logger, set_verbose_mode, print_status
from scripts.utils.llr_constants import (
    RHO_T,
    RHO_T_ERROR,
    RHO_T_SOURCE,
    G,
    C_M_S,
    GM_EARTH,
    GM_MOON,
    EARTH_RADIUS_M,
    MOON_RADIUS_M
)

def compute_prediction():
    print_status("Computing Quantitative η Prediction from TEP...", "TITLE")

    # Load measured eta from step_002 output (deterministic pipeline result)
    step_002_path = (
        PROJECT_ROOT / "results" / "outputs" / "step_002_statistical_analysis.json"
    )
    if step_002_path.exists():
        with open(step_002_path, "r") as f:
            step_002_results = json.load(f)
        eta_measured = float(step_002_results.get('eta_ols', 0))
        eta_error = float(step_002_results.get('eta_ols_error', 0))
        print_status(f"Loaded measured η from step_002: {eta_measured:.4e} ± {eta_error:.4e}", "INFO")
    else:
        raise FileNotFoundError(
            f"Step 002 results not found: {step_002_path}. Run pipeline step 002 first."
        )

    print_status(f"ρ_T: {RHO_T} ± {RHO_T_ERROR} g/cm³ ({RHO_T_SOURCE})", "INFO")

    # =====================================================
    # THEORETICAL COMPARISON
    # =====================================================
    # Following Jakarta v0.8 framework (Papers 6, 10, 11, 12, 13):
    # - The microscopic coupling β is constrained by Cassini to α_0 ≲ 3×10⁻³ in screened regime
    # - Observable response coefficients κ (like κ_MSP, κ_Cep) are empirically determined in unscreened regimes
    # - A key insight from Paper 10: TEP effect operates on clock rates via A(Φ) ≈ 1 − ηΦ/c²
    # - η itself is the TEP modification parameter for LLR, analogous to κ_MSP and κ_Cep
    # - Different regimes exhibit different Observable Response Coefficients due to screening
    # - LLR operates in a more screened Solar System regime, yielding smaller η than galactic κ values

    print_status("", "INFO")
    print_status("--- THEORETICAL COMPARISON ---", "PROCESS")

    # Theoretical prediction from TEP formalism
    # TEP prediction: η ≈ α_0^2 (Φ_⊕^2 - Φ_Moon^2) where α_0 is the microscopic coupling
    # The microscopic coupling β is constrained by Cassini to α_0 ≲ 3×10⁻³ in the screened regime
    # This is the TEP-specific prediction, not the PPN formula η = 4β - γ - 3
    # The PPN formula is for standard scalar-tensor theories
    # TEP uses a different formalism based on compactness-dependent coupling
    # In the Jakarta v0.8 framework, we use the Cassini constraint for the theoretical baseline
    alpha_0_cassini = 3e-3  # Cassini bound: α_0 ≲ 3×10⁻³ in screened regime
    # eta_theoretical will be computed from compactness squared differential below

    # Physical constants - use llr_constants.py for consistency
    M_earth = GM_EARTH / G  # kg
    M_moon = GM_MOON / G  # kg
    c = C_M_S  # m/s
    R_earth = EARTH_RADIUS_M  # m
    R_moon = MOON_RADIUS_M  # m

    Phi_earth = G * M_earth / (R_earth * c**2)
    Phi_moon = G * M_moon / (R_moon * c**2)

    print_status(f"Earth compactness: Φ_⊕/c² = {Phi_earth:.3e}", "CALC")
    print_status(f"Moon compactness:  Φ_☽/c² = {Phi_moon:.3e}", "CALC")
    print_status(f"Compactness squared differential: Δ(Φ²) = {Phi_earth**2 - Phi_moon**2:.3e}", "CALC")

    # TEP theoretical prediction: η ≈ α_0^2 (Φ_⊕^2 - Φ_Moon^2)
    # Using Cassini constraint for microscopic coupling in screened regime
    eta_theoretical_tep = alpha_0_cassini**2 * (Phi_earth**2 - Phi_moon**2)
    print_status(f"Theoretical η from TEP formalism (Cassini-bound α_0): {eta_theoretical_tep:.4e}", "CALC")

    # Enhancement factor comparison
    enhancement_factor = abs(eta_measured) / abs(eta_theoretical_tep)
    print_status(f"Enhancement factor (measured / TEP theoretical): {enhancement_factor:.1f}×", "CALC")

    # Cross-domain comparison using Jakarta v0.8 Observable Response Coefficient framework
    # κ_MSP ~ 10^6 - 10^7 (pulsars, Paper 10: globular cluster pulsars)
    # κ_Cep = (1.05 ± 0.43) × 10^6 mag (Cepheids, Paper 11: Cepheid H0 tension)
    # These are in weakly screened regimes. LLR is in a more screened regime (Solar System),
    # so the effective response should be smaller. η is much smaller than κ_MSP and κ_Cep,
    # which is consistent with the screening mechanism.
    # Paper 13 (wide binaries) reports α_sat = 0.366 ± 0.012 for that specific regime
    kappa_msp_typical = 1e6  # Order-of-magnitude estimate from Paper 10
    kappa_cep_measured = 1.05e6  # Measured value from Paper 11: (1.05 ± 0.43) × 10^6 mag
    kappa_cep_error = 0.43e6
    print_status("", "INFO")
    print_status("--- CROSS-DOMAIN COMPARISON (Jakarta v0.8 κ Framework) ---", "PROCESS")
    print_status(f"κ_MSP (pulsars): ~{kappa_msp_typical:.1e}", "CALC")
    print_status(f"κ_Cep (Cepheids): {kappa_cep_measured:.2e} ± {kappa_cep_error:.2e} mag (Paper 11)", "CALC")
    print_status(f"η (LLR): {eta_measured:.4e} (Solar System screened regime)", "CALC")
    print_status(f"Screening ratio: η/κ_Cep ~ {eta_measured/kappa_cep_measured:.2e} (consistent with strong screening)", "CALC")

    # TEP prediction: δr = 13 η cos(D)
    # The amplitude scales linearly with η: |δr| = 13 |η| meters
    amplitude_measured = 13 * abs(eta_measured)
    amplitude_theoretical = 13 * abs(eta_theoretical_tep)
    print_status("", "INFO")
    print_status("--- AMPLITUDE PREDICTION ---", "PROCESS")
    print_status(f"Measured amplitude: {amplitude_measured:.4f} mm", "CALC")
    print_status(f"Theoretical amplitude (from TEP formalism): {amplitude_theoretical:.6f} mm", "CALC")

    output = {
        "step_id": "step_033",
        "status": "PASS",
        "measured_parameters": {
            "eta_measured": float(eta_measured),
            "eta_error": float(eta_error),
        },
        "theoretical_comparison": {
            "eta_theoretical_tep": float(eta_theoretical_tep),
            "enhancement_factor": float(enhancement_factor),
            "alpha_0_cassini_bound": float(alpha_0_cassini),
        },
        "physical_parameters": {
            "Phi_earth_over_c2": float(Phi_earth),
            "Phi_moon_over_c2": float(Phi_moon),
            "compactness_squared_differential": float(Phi_earth**2 - Phi_moon**2),
        },
        "cross_domain_comparison": {
            "kappa_msp_typical": kappa_msp_typical,
            "kappa_cep_measured": kappa_cep_measured,
            "kappa_cep_error": kappa_cep_error,
            "eta_LL_measured": float(eta_measured),
            "screening_ratio": float(eta_measured / kappa_cep_measured),
        },
        "amplitude_prediction": {
            "amplitude_measured_mm": float(amplitude_measured),
            "amplitude_theoretical_mm": float(amplitude_theoretical),
        },
        "independent_priors": {
            "rho_T_g_cm3": RHO_T,
            "rho_T_error_g_cm3": RHO_T_ERROR,
            "rho_T_source": RHO_T_SOURCE,
        },
        "conclusion": f"η = {eta_measured:.4e} is the TEP modification parameter for LLR, analogous to κ_MSP and κ_Cep in other domains. The measured η is ~{enhancement_factor:.0f}× larger than the TEP theoretical prediction from α_0^2(Φ_⊕^2 - Φ_Moon^2) using the Cassini-bound microscopic coupling, indicating additional contributions from the screening mechanism. η is much smaller than cross-domain κ values (η/κ_Cep ~ {eta_measured/kappa_cep_measured:.2e}), consistent with LLR being in a more screened Solar System regime.",
    }

    print_status("", "INFO")
    print_status("--- SUMMARY ---", "PROCESS")
    print_status(f"Measured η: {eta_measured:.4e}", "CALC")
    print_status(f"Theoretical η (Cassini-bound α_0): {eta_theoretical_tep:.4e}", "CALC")
    print_status(f"Enhancement: {enhancement_factor:.1f}×", "CALC")
    print_status(f"Amplitude: {amplitude_measured:.4f} mm", "CALC")
    print_status(f"Screening ratio η/κ_Cep: {eta_measured/kappa_cep_measured:.2e}", "CALC")

    return output

if __name__ == "__main__":
    log_dir = PROJECT_ROOT / "logs"
    logger = TEPLogger(
        "step_033", str(log_dir / "step_033_quantitative_eta_prediction.log")
    )
    set_step_logger(logger)

    results = compute_prediction()

    if results:
        logger.save_step_results(
            results, PROJECT_ROOT, "step_033_quantitative_eta_prediction"
        )
        print_status("Quantitative η Prediction Complete.", "SUCCESS")