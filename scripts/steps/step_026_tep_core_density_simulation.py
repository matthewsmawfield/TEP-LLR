#!/usr/bin/env python3
"""
Step 028: TEP Core Density Simulation

Simulates Earth's and Moon's internal density profiles to understand the
degree of Temporal Topology flattening in each planetary interior.

Following the updated TEP framework (Papers 6, 10, 11, 12, 13):
- ρ_T is the Temporal Topology saturation density (fundamental scale, not a local density switch)
- The microscopic coupling β is constrained by Cassini to α_0 ≲ 3×10⁻³ in screened regime
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
    # - The microscopic coupling β is constrained by Cassini to α_0 ≲ 3×10⁻³ in screened regime
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

    # Load measured η from step_002 for comparison
    step_002_path = PROJECT_ROOT / "results" / "outputs" / "step_002_statistical_analysis.json"
    if step_002_path.exists():
        with open(step_002_path, "r") as f:
            step_002_results = json.load(f)
        measured_eta = float(step_002_results.get("eta_ols", 0))
        measured_eta_error = float(step_002_results.get("eta_ols_error", 0))
        logger.info(f"    Measured η from LLR: {measured_eta:.4e} ± {measured_eta_error:.4e}")
    else:
        measured_eta = None
        measured_eta_error = None
        logger.warning("    Step 002 results not found, cannot compare with measured η")

    # The shielding differential (Earth - Moon) quantifies the differential degree of
    # Temporal Topology flattening between Earth and Moon interiors. This differential
    # is the physical origin of the Nordtvedt effect in TEP: Earth's deeper potential
    # well leads to stronger screening, causing it to fall slightly differently toward
    # the Sun than the Moon.
    shielding_differential = earth_shielding - moon_shielding

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
        "conclusions": {
            "shielding_differential": f"{shielding_differential:.2e}",
            "explanation": "The simulation calculates Earth and Moon shielding factors from density profiles using ρ_T. The shielding differential quantifies the differential degree of Temporal Topology flattening, which is the physical origin of the Nordtvedt effect in TEP. η itself is the observable response coefficient for LLR, not derived from κ_LL.",
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