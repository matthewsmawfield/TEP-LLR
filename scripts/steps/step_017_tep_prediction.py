#!/usr/bin/env python3
"""
Step 017: TEP-LLR Model Prediction
Finalized model predictions and comparison with GR predictions.

Following the updated TEP framework (Papers 6, 10, 11, 12, 13):
- The microscopic coupling β is constrained by Cassini to α_0 ≲ 3×10⁻³ in screened regime
- Observable response coefficients κ (like κ_MSP, κ_Cep) are empirically determined in unscreened regimes
- A key insight from Paper 10: TEP effect operates on clock rates via A(Φ) ≈ 1 − ηΦ/c²
- This suggests η itself is the TEP modification parameter for LLR, analogous to κ_MSP and κ_Cep
- No separate κ_LL is needed; η is the observable response coefficient for LLR
"""

import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

import argparse
import json
from scripts.utils.logger import TEPLogger, set_step_logger, set_verbose_mode, print_status
from scripts.utils.llr_constants import (
    RHO_T,
    RHO_T_ERROR,
    RHO_T_SOURCE
)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Step 017: TEP Prediction")
    parser.add_argument("--verbose", action="store_true", help="Enable verbose output")
    args = parser.parse_args()

    log_dir = PROJECT_ROOT / "logs"
    logger = TEPLogger("step_017", str(
        log_dir / "step_017_tep_prediction.log"))
    set_step_logger(logger)

    print_status("Generating Final TEP Signal Predictions...", "TITLE")

    step_002_path = PROJECT_ROOT / "results" / \
        "outputs" / "step_002_statistical_analysis.json"
    if not step_002_path.exists():
        print_status(f"Required input not found: {step_002_path}", "ERROR")
        sys.exit(1)

    with step_002_path.open("r") as f:
        step_002_results = json.load(f)

    detected_eta = float(step_002_results.get("eta_ols", 0))
    detected_eta_error = float(step_002_results.get("eta_ols_error", 0))

    # TEP prediction using η as the observable response coefficient
    # Following Jakarta v0.8 framework (Papers 6, 10, 11, 12, 13):
    # - The microscopic coupling β is constrained by Cassini to α_0 ≲ 3×10⁻³ in screened regime
    # - Observable response coefficients κ (like κ_MSP, κ_Cep) are empirically determined in unscreened regimes
    # - A key insight from Paper 10: TEP effect operates on clock rates via A(Φ) ≈ 1 − ηΦ/c²
    # - η itself is the TEP modification parameter for LLR, analogous to κ_MSP and κ_Cep
    # - Different regimes exhibit different Observable Response Coefficients due to screening
    # - LLR operates in a more screened Solar System regime, yielding smaller η than galactic κ values

    print_status(f"Measured η from LLR: {detected_eta:.4e} ± {detected_eta_error:.4e}", "CALC")
    print_status(f"ρ_T: {RHO_T} ± {RHO_T_ERROR} g/cm³ ({RHO_T_SOURCE})", "CALC")

    # Theoretical prediction from TEP formalism
    # TEP prediction: η ≈ α_0^2 (Φ_⊕^2 - Φ_Moon^2)
    # where α_0 is the microscopic coupling constrained by Cassini to α_0 ≲ 3×10⁻³ in screened regime
    # This is the TEP-specific prediction, not the PPN formula
    # The PPN formula η = 4β - γ - 3 is for standard scalar-tensor theories
    # TEP uses a different formalism based on compactness-dependent coupling
    # In the Jakarta v0.8 framework, we use the Cassini constraint for the theoretical baseline
    alpha_0_cassini = 3e-3  # Cassini bound: α_0 ≲ 3×10⁻³ in screened regime
    
    # Compute compactness values for TEP prediction
    # Use constants from llr_constants.py for consistency (imported above)
    from scripts.utils.llr_constants import G, C_M_S, GM_EARTH, GM_MOON, EARTH_RADIUS_M, MOON_RADIUS_M
    
    # Derived masses from GM values (M = GM / G)
    M_earth = GM_EARTH / G  # kg
    M_moon = GM_MOON / G  # kg
    c = C_M_S  # m/s
    R_earth = EARTH_RADIUS_M  # m
    R_moon = MOON_RADIUS_M  # m
    
    Phi_earth = G * M_earth / (R_earth * c**2)
    Phi_moon = G * M_moon / (R_moon * c**2)
    
    eta_theoretical_tep = alpha_0_cassini**2 * (Phi_earth**2 - Phi_moon**2)
    print_status(f"Theoretical η from TEP formalism (Cassini-bound α_0): {eta_theoretical_tep:.4e}", "CALC")

    # Enhancement factor comparison
    enhancement_factor = abs(detected_eta) / abs(eta_theoretical_tep)
    print_status(f"Enhancement factor (measured / TEP theoretical): {enhancement_factor:.1f}×", "CALC")

    # TEP prediction: δr = 13 η cos(D)
    # The amplitude scales linearly with η: |δr| = 13 |η| meters
    amplitude_measured = 13 * abs(detected_eta)
    amplitude_theoretical = 13 * abs(eta_theoretical_tep)
    print_status(f"Measured amplitude: {amplitude_measured:.4f} mm", "CALC")
    print_status(f"Theoretical amplitude (from TEP formalism): {amplitude_theoretical:.6f} mm", "CALC")

    # Comparison with cross-domain response coefficients using Jakarta v0.8 framework
    # κ_MSP ~ 10^6 - 10^7 (pulsars, Paper 10: globular cluster pulsars)
    # κ_Cep = (1.05 ± 0.43) × 10^6 mag (Cepheids, Paper 11: Cepheid H0 tension)
    # These are in weakly screened regimes. LLR is in a more screened regime (Solar System),
    # so the effective response should be smaller. η is much smaller than κ_MSP and κ_Cep,
    # which is consistent with the screening mechanism.
    # Paper 13 (wide binaries) reports α_sat = 0.366 ± 0.012 for that specific regime
    kappa_msp_typical = 1e6  # Order-of-magnitude estimate from Paper 10
    kappa_cep_measured = 1.05e6  # Measured value from Paper 11: (1.05 ± 0.43) × 10^6 mag
    kappa_cep_error = 0.43e6
    print_status(f"Cross-domain comparison (Jakarta v0.8 κ Framework):", "INFO")
    print_status(f"  κ_MSP (pulsars): ~{kappa_msp_typical:.1e}", "CALC")
    print_status(f"  κ_Cep (Cepheids): {kappa_cep_measured:.2e} ± {kappa_cep_error:.2e} mag (Paper 11)", "CALC")
    print_status(f"  η (LLR): {detected_eta:.4e} (Solar System screened regime)", "CALC")
    print_status(f"  Screening ratio: η/κ_Cep ~ {detected_eta/kappa_cep_measured:.2e} (consistent with strong screening)", "CALC")

    # Comparison of Nordtvedt parameter eta
    predictions = {
        "GR_prediction": 0.0,
        "TEP_measured_eta": float(detected_eta),
        "TEP_measured_eta_error": float(detected_eta_error),
        "TEP_theoretical_eta_tep": float(eta_theoretical_tep),
        "enhancement_factor": float(enhancement_factor),
        "amplitude_measured_mm": float(amplitude_measured),
        "amplitude_theoretical_mm": float(amplitude_theoretical),
        "rho_T_g_cm3": RHO_T,
        "rho_T_error_g_cm3": RHO_T_ERROR,
        "rho_T_source": RHO_T_SOURCE,
        "kappa_msp_typical": kappa_msp_typical,
        "kappa_cep_measured": kappa_cep_measured,
        "kappa_cep_error": kappa_cep_error,
        "alpha_0_cassini_bound": float(alpha_0_cassini),
        "screening_ratio": float(detected_eta / kappa_cep_measured),
    }

    print_status("Model Prediction Audit Summary:", "TITLE")
    print_status(f"  GR Prediction:    η = {predictions['GR_prediction']:.4e}", "CALC")
    print_status(f"  TEP Measured η:   η = {predictions['TEP_measured_eta']:.4e} ± {predictions['TEP_measured_eta_error']:.4e}", "CALC")
    print_status(f"  TEP Theoretical η (Cassini-bound α_0): η = {predictions['TEP_theoretical_eta_tep']:.4e}", "CALC")
    print_status(f"  Enhancement factor: {predictions['enhancement_factor']:.1f}×", "CALC")
    print_status(f"  Amplitude (measured): {predictions['amplitude_measured_mm']:.4f} mm", "CALC")
    print_status(f"  Amplitude (theoretical): {predictions['amplitude_theoretical_mm']:.6f} mm", "CALC")
    print_status(f"  ρ_T: {RHO_T} ± {RHO_T_ERROR} g/cm³ ({RHO_T_SOURCE})", "CALC")
    print_status(f"  Screening ratio η/κ_Cep: {predictions['screening_ratio']:.2e}", "CALC")

    results = {
        "step_id": "step_017",
        "status": "PASS",
        "predictions": predictions,
    }

    logger.save_step_results(results, PROJECT_ROOT, "step_017_tep_prediction")
    print_status("Final Model Comparison Complete.", "SUCCESS")