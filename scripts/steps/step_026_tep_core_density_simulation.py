#!/usr/bin/env python3
"""
Step 026: TEP Core Density Simulation

Simulates Earth's and Moon's internal density profiles to understand the
degree of Temporal Topology flattening in each planetary interior.

Following the updated TEP framework (Papers 6, 10, 11, 12, 13):
- ρ_T is the Temporal Topology saturation density (fundamental scale, not a local density switch)
- Cassini bounds the screened solar source charge |α_0 S_Σ^(⊙)| ≲ 8.2×10⁻⁶;
  in the standard unscreened limit it bounds the bare coupling α_0 ≲ 3.4×10⁻³
- Observable response coefficients κ (like κ_MSP, κ_Cep) are empirically determined in unscreened regimes
- A key insight from Paper 10: TEP effect operates on clock rates via A(Φ) ≈ 1 − ηΦ/c²
- η itself is the TEP modification parameter for LLR, analogous to κ_MSP and κ_Cep
- No separate κ_LL is needed; η is the observable response coefficient for LLR
- This step uses ρ_T from GNSS calibration as an independent prior (Paper 6)
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
    RHO_T_SOURCE
)

def run_core_simulation(logger, rho_T_input=None):
    """
    Run core density simulation with optional ρ_T input.

    CRITICAL NOTE: ρ_T is the Temporal Topology saturation density calibrated from
    GNSS atomic clock correlations (Paper 6, UCD). This is treated as an independent
    prior from a different observational domain (terrestrial geodesy), not a circular
    dependency on LLR results.
    """
    logger.info(">>> Starting TEP Core Density Simulation...")

    # Temporal Topology Saturation Density (ρ_T)
    # If not provided, use the value from GNSS calibration (Paper 6, UCD)
    # This is an independent prior from terrestrial geodesy, not LLR
    if rho_T_input is not None:
        RHO_T_G_CM3 = float(rho_T_input)
        logger.info(f"Using externally provided ρ_T = {RHO_T_G_CM3} g/cm³")
    else:
        RHO_T_G_CM3 = RHO_T
        logger.info(f"Using ρ_T = {RHO_T_G_CM3} g/cm³ from {RHO_T_SOURCE} (independent prior from GNSS)")

    # Earth Density Profile (Simplified PREM)
    # Inner Core, Outer Core, Lower Mantle, Upper Mantle, Crust
    # Earth radius ~ 6371 km
    earth_shells = [
        {"r_inner": 0, "r_outer": 1221, "rho": 13.0, "name": "Inner Core"},
        {"r_inner": 1221, "r_outer": 3480, "rho": 10.9, "name": "Outer Core"},
        {"r_inner": 3480, "r_outer": 5701, "rho": 5.0, "name": "Lower Mantle"},
        {"r_inner": 5701, "r_outer": 6371, "rho": 3.3, "name": "Upper Mantle/Crust"},
    ]

    # Moon Density Profile
    # Core is very small, predominantly homogeneous mantle ~ 3.3 g/cm^3
    moon_shells = [
        {"r_inner": 0, "r_outer": 330, "rho": 8.0, "name": "Small Core"},
        {"r_inner": 330, "r_outer": 1737, "rho": 3.3, "name": "Mantle/Crust"},
    ]

    # Calculating volumetric suppression integrals (phenomenological TEP evaluation)
    # Following updated framework (Papers 6, 10, 11, 12, 13):
    # - Cassini bounds the screened solar source charge |α_0 S_Σ^(⊙)| ≲ 8.2×10⁻⁶
    # - Observable response coefficients κ (like κ_MSP, κ_Cep) are empirically determined in unscreened regimes
    # - A key insight from Paper 10: TEP effect operates on clock rates via A(Φ) ≈ 1 − ηΦ/c²
    # - η itself is the TEP modification parameter for LLR, analogous to κ_MSP and κ_Cep
    # - No separate κ_LL is needed; η is the observable response coefficient for LLR
    # - This step calculates shielding factors from density profiles to understand the physics

    logger.info("Calculating Earth and Moon shielding factors from density profiles...")

    def calculate_shielding_parameter(shells):
        total_vol = 0
        shield_integral = 0
        for shell in shells:
            ri = shell["r_inner"]
            ro = shell["r_outer"]
            vol = (4.0 / 3.0) * np.pi * (ro**3 - ri**3)
            # Volumetric weighted suppression factor (rho / ρ_T)^3
            suppression_factor = (shell["rho"] / RHO_T_G_CM3) ** 3
            shield_integral += suppression_factor * vol
            total_vol += vol

        return shield_integral / total_vol

    earth_shielding = calculate_shielding_parameter(earth_shells)
    moon_shielding = calculate_shielding_parameter(moon_shells)

    logger.info(f"    Earth shielding factor: {earth_shielding:.4e}")
    logger.info(f"    Moon shielding factor: {moon_shielding:.4e}")

    # Load measured η from step_003 statistical output
    step_003_path = PROJECT_ROOT / "results" / "outputs" / "step_003_statistical_analysis.json"
    if not step_003_path.exists():
        raise FileNotFoundError(
            f"step_003_statistical_analysis.json not found: {step_003_path}. "
            "Run the empirical pipeline before the core-density comparison."
        )
    with open(step_003_path, "r") as f:
        step_003_results = json.load(f)
    measured_eta = require_step003_eta_ols(step_003_results)
    measured_eta_error = require_step003_eta_ols_error(step_003_results)
    logger.info(f"    Measured η from LLR: {measured_eta:.4e} ± {measured_eta_error:.4e}")

    # The shielding differential (Earth - Moon) quantifies the differential degree of
    # Temporal Topology flattening between Earth and Moon interiors. This differential
    # is the physical origin of the Nordtvedt effect in TEP: Earth's deeper potential
    # well leads to stronger screening, causing it to fall slightly differently toward
    # the Sun than the Moon.
    shielding_differential = earth_shielding - moon_shielding

    # Volumetric suppression model estimate for η
    # η = α_eff^(Earth-Moon) * (Earth_shielding - Moon_shielding)
    # Source-charge bookkeeping (issue 17-9): in any mediated-force
    # decomposition the differential-fall observable carries BOTH vertices —
    # the screened solar source charge and the screened body coupling —
    # δa/a_0 = (α_0 S_Sigma^(Sun)) * Δα_body. Cassini bounds the source
    # factor |α_0 S_Sigma^(Sun)| ≲ sqrt(2)·5.8e-6 ≈ 8.2e-6 (linear
    # source-charge map, α_0 = sqrt(2)·β_A). The volumetric number below
    # therefore admits only two readings, both fatal to a force-channel
    # interpretation:
    #   (a) α_eff is a body coupling — then the source factor is missing and
    #       the product-consistent estimate is
    #       |α_0 S_Sigma|·|α_eff|·Δ⟨(ρ/ρ_T)^3⟩ ≈ 8e-10 (body-coupling
    #       reading) or |α_0 S_Sigma|·|α_eff| ≈ 2.4e-8 (α_eff as the whole
    #       body-charge differential) — five orders below 10^-4;
    #   (b) α_eff is the full response — then it is NOT distinct from the
    #       solar source charge it must contain, and the "distinct
    #       projection" exemption is incoherent.
    # The quoted band is thus an effective observable response coefficient
    # (all microscopic factors folded into one number), not a decomposable
    # mediated-force prediction; it is additionally superseded by the
    # canonical pairwise operator (Paper 0 v0.15, Appendix E R11), which
    # suppresses the Earth-Moon force channel to ~10^-14.
    alpha_eff_earth_moon = -3.0e-3  # Illustrative Earth-Moon effective coupling
    eta_volumetric = alpha_eff_earth_moon * shielding_differential
    alpha_eff_solar_cassini = (2.0 ** 0.5) * 5.8e-6  # ≈ 8.2e-6
    eta_mediated_full = alpha_eff_solar_cassini * alpha_eff_earth_moon * shielding_differential
    eta_mediated_compact = alpha_eff_solar_cassini * alpha_eff_earth_moon
    logger.info(f"    Volumetric η estimate (effective ORC): {eta_volumetric:.4e}")
    logger.info(f"    (α_eff(Earth-Moon) = {alpha_eff_earth_moon}, shielding differential = {shielding_differential:.4e})")
    logger.info(f"    Mediated-force ceiling (α_⊙·α_eff·Δshielding): {eta_mediated_full:.4e}")
    logger.info(f"    Mediated-force ceiling (α_⊙·α_eff): {eta_mediated_compact:.4e}")

    # In the old framework, the relationship was η = -κ_LL * (Earth_Shielding - Moon_Shielding)
    # However, based on the new insight from Paper 10, η itself is the TEP modification parameter.
    # The shielding differential provides physical context for understanding the magnitude
    # of η, but η is not derived from it via a κ_LL coupling.
    #
    # The measured η can be compared with the shielding differential to understand the
    # efficiency of the screening mechanism in converting shielding differences into
    # the observable Nordtvedt parameter.

    if measured_eta is not None:
        # Calculate the effective coupling that would produce the measured η
        # This is for diagnostic purposes only, not a theoretical prediction
        effective_coupling = abs(measured_eta) / shielding_differential
        logger.info(f"    Effective coupling (diagnostic): η/Δshielding = {effective_coupling:.4e}")
        logger.info(f"    Shielding differential: {shielding_differential:.4e}")

        # Ratio of measured η to volumetric prediction
        volumetric_ratio = abs(measured_eta) / abs(eta_volumetric)
        logger.info(f"    Measured / volumetric prediction ratio: {volumetric_ratio:.2f}")

    results = {
        "step_id": "step_026",
        "data_type": "SYNTHETIC (THEORETICAL MODELING)",
        "status": "PASS",
        "independent_prior_note": f"ρ_T = {RHO_T_G_CM3} g/cm³ from {RHO_T_SOURCE} (independent prior from terrestrial geodesy)",
        "parameters": {
            "rho_T_g_cm3": float(RHO_T_G_CM3),
            "rho_T_source": RHO_T_SOURCE if rho_T_input is None else "external_input",
        },
        "simulation_results": {
            "earth_shielding_factor": float(earth_shielding),
            "moon_shielding_factor": float(moon_shielding),
            "shielding_differential": float(shielding_differential),
        },
        "comparison_with_measurement": {
            "measured_eta": float(measured_eta) if measured_eta is not None else None,
            "measured_eta_error": float(measured_eta_error) if measured_eta_error is not None else None,
            "effective_coupling_diagnostic": float(effective_coupling) if measured_eta is not None else None,
        } if measured_eta is not None else None,
        "volumetric_prediction": {
            "alpha_eff_earth_moon": float(alpha_eff_earth_moon),
            "eta_volumetric": float(eta_volumetric),
            "measured_to_volumetric_ratio": float(volumetric_ratio) if measured_eta is not None else None,
            "caveat": "Effective observable-response band only: the mediated-force decomposition is product-limited by the screened solar source charge and returns ≲10⁻⁹–10⁻⁸, so the quoted band cannot be a charge×charge differential-fall estimate. It is additionally superseded by the canonical pairwise operator (Paper 0 v0.15, Appendix E R11), which suppresses the Earth-Moon force channel to ~10⁻¹⁴.",
        },
        "source_charge_product_limit": {
            "alpha_eff_solar_cassini_ceiling": float(alpha_eff_solar_cassini),
            "eta_mediated_with_shielding_differential": float(eta_mediated_full),
            "eta_mediated_alpha_eff_as_full_differential": float(eta_mediated_compact),
            "verdict": "Under any mediated-force reading the differential-fall observable is bounded by ≲10⁻⁹–10⁻⁸ — five orders below the measured 10⁻⁴ band. The volumetric α_eff can be read only as an effective response coefficient (all factors folded in), not as a body coupling; the 'distinct projection' exemption from the source-charge bound is not available to a force channel.",
        },
        "conclusions": {
            "shielding_differential": f"{shielding_differential:.2e}",
            "explanation": "The simulation calculates Earth and Moon shielding factors from density profiles using ρ_T. The shielding differential quantifies the differential degree of Temporal Topology flattening between the two interiors. The volumetric suppression model yields an effective η band ≈ α_eff^(Earth-Moon) Δ⟨(ρ/ρ_T)^3⟩ ≈ -10⁻⁴ in the observed decade, but as a mediated-force estimate it is internally inconsistent: folding in the Cassini-bounded solar source charge gives ≲10⁻⁹–10⁻⁸, and the canonical pairwise operator independently suppresses the force channel to ~10⁻¹⁴. The surviving admissible channel for the measured residual is the clock-sector/ephemeris-reduction projection (Section 2.1), whose amplitude map is the open estimand question of the source-level integrator refit. η itself is the observable response coefficient for LLR, not derived from κ_LL.",
        },
    }

    logger.info(f"    Shielding differential: {shielding_differential:.2e}")
    logger.info(f"✓   TEP Core Density Simulation Complete.")

    return results

def main():
    # Setup TEPLogger
    log_dir = PROJECT_ROOT / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    logger = TEPLogger(
        "step_026", str(log_dir / "step_026_tep_core_density_simulation.log")
    )
    set_step_logger(logger)

    results = run_core_simulation(logger)

    output_path = (
        PROJECT_ROOT
        / "results"
        / "outputs"
        / "step_026_tep_core_density_simulation.json"
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w") as f:
        json.dump(results, f, indent=4)
    output_rel = (
        output_path.relative_to(PROJECT_ROOT)
        if output_path.is_relative_to(PROJECT_ROOT)
        else output_path
    )
    logger.info(f"    Saved results to {output_rel}")

if __name__ == "__main__":
    main()
