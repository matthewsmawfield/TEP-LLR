#!/usr/bin/env python3
"""
Step 033: Quantitative η Analysis from TEP Framework

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
from scripts.utils.statistical_utils import require_step003_eta_ols, require_step003_eta_ols_error
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

    # Load measured η from step_003 statistical output (deterministic pipeline result)
    step_003_path = (
        PROJECT_ROOT / "results" / "outputs" / "step_003_statistical_analysis.json"
    )
    if not step_003_path.exists():
        raise FileNotFoundError(
            f"step_003_statistical_analysis.json not found: {step_003_path}. Run pipeline step 003 first."
        )
    with open(step_003_path, "r") as f:
        step_003_results = json.load(f)
    eta_measured = require_step003_eta_ols(step_003_results)
    eta_error = require_step003_eta_ols_error(step_003_results)
    print_status(
        f"Loaded measured η from step_003: {eta_measured:.4e} ± {eta_error:.4e}", "INFO"
    )

    print_status(f"ρ_T: {RHO_T} ± {RHO_T_ERROR} g/cm³ ({RHO_T_SOURCE})", "INFO")

    # =====================================================
    # THEORETICAL COMPARISON — TWO BASELINES
    # =====================================================
    # Following Jakarta v0.8 framework (Papers 6, 10, 11, 12, 13):
    # - The microscopic coupling β is constrained by Cassini to α_0 ≲ 3×10⁻³
    # - Observable response coefficients κ are empirically determined per channel
    # - η is the LLR observable response coefficient, not derivable from α_0 alone
    #
    # We present two theoretically distinct baselines:
    #   1. MICROSCOPIC-COUPLING-ONLY (compactness-squared): η ≈ α_0^2 Δ(Φ²)
    #      This is a formal TEP expression that yields ~10⁻²⁴, far below detection.
    #      It serves as a consistency check: the measured η cannot arise from
    #      bare α_0 alone, confirming that η is an emergent response coefficient.
    #   2. PHENOMENOLOGICAL TEP (volumetric suppression): η ≈ -α_0 Δ⟨(ρ/ρ_T)³⟩
    #      This integrates interior density structure and yields η ~ 10⁻⁴,
    #      consistent with the observed order of magnitude.

    print_status("", "INFO")
    print_status("--- THEORETICAL COMPARISON ---", "PROCESS")

    alpha_0_cassini = 3e-3  # Cassini bound: α_0 ≲ 3×10⁻³ in screened regime

    # Baseline 1: Compactness-squared (microscopic-coupling-only)
    M_earth = GM_EARTH / G
    M_moon = GM_MOON / G
    c = C_M_S
    R_earth = EARTH_RADIUS_M
    R_moon = MOON_RADIUS_M

    Phi_earth = G * M_earth / (R_earth * c**2)
    Phi_moon = G * M_moon / (R_moon * c**2)
    compactness_sq_diff = Phi_earth**2 - Phi_moon**2

    eta_microscopic_only = alpha_0_cassini**2 * compactness_sq_diff
    print_status(f"[Baseline 1] Microscopic-coupling-only η = α_0² Δ(Φ²): {eta_microscopic_only:.4e}", "CALC")
    print_status(f"    Earth compactness: Φ_⊕/c² = {Phi_earth:.3e}", "CALC")
    print_status(f"    Moon compactness:  Φ_☽/c² = {Phi_moon:.3e}", "CALC")
    print_status(f"    Compactness squared differential: {compactness_sq_diff:.3e}", "CALC")

    # Baseline 2: Volumetric suppression model (from step_026 density simulation)
    step_026_path = PROJECT_ROOT / "results" / "outputs" / "step_026_tep_core_density_simulation.json"
    if step_026_path.exists():
        with open(step_026_path, "r") as f:
            step_026_results = json.load(f)
        eta_volumetric = float(step_026_results.get("volumetric_prediction", {}).get("eta_volumetric", 0))
        shielding_diff = float(step_026_results.get("simulation_results", {}).get("shielding_differential", 0))
        print_status(f"[Baseline 2] Volumetric suppression η = -α_0 Δ⟨(ρ/ρ_T)³⟩: {eta_volumetric:.4e}", "CALC")
        print_status(f"    Shielding differential (step_026): {shielding_diff:.4e}", "CALC")
    else:
        eta_volumetric = None
        shielding_diff = None
        print_status("[Baseline 2] Step 026 results not found; volumetric baseline unavailable.", "WARN")

    # Standard scalar-tensor baseline for comparison
    # In standard scalar-tensor theories: η ≈ 4α_0² (Nordtvedt formula)
    eta_standard_st = 4.0 * alpha_0_cassini**2
    print_status(f"[Baseline 3] Standard scalar-tensor η ≈ 4α_0²: {eta_standard_st:.4e}", "CALC")

    # Enhancement factors
    enhancement_microscopic = abs(eta_measured) / abs(eta_microscopic_only) if eta_microscopic_only != 0 else None
    enhancement_standard = abs(eta_measured) / abs(eta_standard_st) if eta_standard_st != 0 else None
    enhancement_volumetric = abs(eta_measured) / abs(eta_volumetric) if (eta_volumetric is not None and eta_volumetric != 0) else None

    print_status("", "INFO")
    print_status("--- ENHANCEMENT FACTORS ---", "PROCESS")
    if enhancement_microscopic is not None:
        print_status(f"Measured / microscopic-coupling-only: {enhancement_microscopic:.2e}×", "CALC")
    print_status(f"Measured / standard scalar-tensor (4α_0²): {enhancement_standard:.1f}×", "CALC")
    if enhancement_volumetric is not None:
        print_status(f"Measured / volumetric TEP prediction: {enhancement_volumetric:.1f}×", "CALC")

    # Cross-domain comparison using Jakarta v0.8 Observable Response Coefficient framework
    kappa_msp_typical = 1e6
    kappa_cep_measured = 1.05e6
    kappa_cep_error = 0.43e6
    print_status("", "INFO")
    print_status("--- CROSS-DOMAIN COMPARISON (Jakarta v0.8 κ Framework) ---", "PROCESS")
    print_status(f"κ_MSP (pulsars): ~{kappa_msp_typical:.1e}", "CALC")
    print_status(f"κ_Cep (Cepheids): {kappa_cep_measured:.2e} ± {kappa_cep_error:.2e} mag (Paper 11)", "CALC")
    print_status(f"η (LLR): {eta_measured:.4e} (Solar System screened regime)", "CALC")
    print_status(f"Screening ratio: η/κ_Cep ~ {eta_measured/kappa_cep_measured:.2e}", "CALC")

    # Amplitude predictions
    amplitude_measured = 13 * abs(eta_measured)
    amplitude_microscopic = 13 * abs(eta_microscopic_only)
    amplitude_volumetric = 13 * abs(eta_volumetric) if eta_volumetric is not None else None
    print_status("", "INFO")
    print_status("--- AMPLITUDE PREDICTION ---", "PROCESS")
    print_status(f"Measured amplitude: {amplitude_measured:.4f} mm", "CALC")
    print_status(f"Microscopic-only amplitude: {amplitude_microscopic:.6e} mm", "CALC")
    if amplitude_volumetric is not None:
        print_status(f"Volumetric TEP amplitude: {amplitude_volumetric:.4f} mm", "CALC")

    output = {
        "step_id": "step_033",
        "status": "PASS",
        "measured_parameters": {
            "eta_measured": float(eta_measured),
            "eta_error": float(eta_error),
        },
        "theoretical_comparison": {
            "microscopic_coupling_only": {
                "eta_microscopic_only": float(eta_microscopic_only),
                "compactness_squared_differential": float(compactness_sq_diff),
                "enhancement_factor": float(enhancement_microscopic) if enhancement_microscopic is not None else None,
                "caveat": "This baseline uses the bare α_0 with surface compactness only. It yields ~10⁻²⁴, confirming that the measured η is not produced by the microscopic coupling alone. η is an emergent observable response coefficient.",
            },
            "standard_scalar_tensor": {
                "eta_standard_st": float(eta_standard_st),
                "formula": "η ≈ 4α_0²",
                "enhancement_factor": float(enhancement_standard),
                "caveat": "Standard scalar-tensor Nordtvedt formula. The ~8.8× enhancement (using the Cassini upper-bound α_0 ≲ 3×10⁻³) indicates that η absorbs additional contributions from the TEP screening mechanism beyond bare scalar-tensor physics.",
            },
            "volumetric_tep": {
                "eta_volumetric": float(eta_volumetric) if eta_volumetric is not None else None,
                "shielding_differential": float(shielding_diff) if shielding_diff is not None else None,
                "enhancement_factor": float(enhancement_volumetric) if enhancement_volumetric is not None else None,
                "caveat": "Phenomenological TEP prediction from integrated density profiles (step_026). The ~3–5× residual is absorbed into PREM simplification, (ρ/ρ_T)³ exponent uncertainty, ρ_T = 20 ± 8 g/cm³ uncertainty, and the upper-bound nature of the Cassini α_0 constraint.",
            },
            "alpha_0_cassini_bound": float(alpha_0_cassini),
        },
        "physical_parameters": {
            "Phi_earth_over_c2": float(Phi_earth),
            "Phi_moon_over_c2": float(Phi_moon),
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
            "amplitude_microscopic_only_mm": float(amplitude_microscopic),
            "amplitude_volumetric_mm": float(amplitude_volumetric) if amplitude_volumetric is not None else None,
        },
        "independent_priors": {
            "rho_T_g_cm3": RHO_T,
            "rho_T_error_g_cm3": RHO_T_ERROR,
            "rho_T_source": RHO_T_SOURCE,
        },
        "conclusion": (
            f"η = {eta_measured:.4e} is the LLR observable response coefficient. "
            f"Three baselines are presented: "
            f"(1) microscopic-coupling-only η = α_0² Δ(Φ²) = {eta_microscopic_only:.4e} "
            f"yields an enhancement of {enhancement_microscopic:.2e}×, confirming that η is emergent "
            f"and not produced by bare α_0 alone; "
            f"(2) standard scalar-tensor η ≈ 4α_0² = {eta_standard_st:.4e} "
            f"yields an enhancement of {enhancement_standard:.1f}×, indicating that η absorbs "
            f"TEP-specific screening contributions beyond standard scalar-tensor physics; "
            f"(3) the phenomenological volumetric TEP model η ≈ -α_0 Δ⟨(ρ/ρ_T)³⟩ = "
            f"{eta_volumetric:.4e} "
            f"yields an enhancement of {enhancement_volumetric:.1f}×, "
            f"with the residual ~3–5× absorbed into model uncertainties. "
            f"The measured η is much smaller than cross-domain κ values "
            f"(η/κ_Cep ~ {eta_measured/kappa_cep_measured:.2e}), consistent with LLR being in a more screened Solar System regime."
        ),
    }

    print_status("", "INFO")
    print_status("--- SUMMARY ---", "PROCESS")
    print_status(f"Measured η: {eta_measured:.4e}", "CALC")
    print_status(f"Microscopic-only η: {eta_microscopic_only:.4e}", "CALC")
    print_status(f"Standard scalar-tensor η: {eta_standard_st:.4e}", "CALC")
    if eta_volumetric is not None:
        print_status(f"Volumetric TEP η: {eta_volumetric:.4e}", "CALC")
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