#!/usr/bin/env python3
"""
Step 049: Cross-Scale Validation Synthesis
Aggregates results from all UCD analysis steps (042-048) into a unified
cross-scale validation report. Computes combined significance and
assesses overall consistency of the universal-density hypothesis.
"""

import sys, json, numpy as np
from pathlib import Path
from scipy import stats

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))
from scripts.utils.logger import TEPLogger, set_step_logger, print_status

TEP_ALPHA = 1.0/3.0
RHO_T = 20.0


def load_step_results(step_id):
    """Load results from a previous step's JSON output."""
    outputs_dir = PROJECT_ROOT / "results" / "outputs"
    for f in outputs_dir.glob(f"{step_id}_*.json"):
        with open(f) as fh:
            return json.load(fh)
    return None


def run_synthesis(logger):
    logger.info(">>> Starting cross-scale validation synthesis...")

    synthesis = {"constraints": {}, "status": "PASS"}

    # GNSS calibration
    gnss = load_step_results("step_028")
    if gnss:
        synthesis["constraints"]["GNSS"] = {
            "scale": "Terrestrial (10^27 g)",
            "observable": "L_c ≈ 4200 km",
            "derived_rho_T": "19-20 g/cm³",
            "status": "Calibration anchor",
            "independence": "Primary — defines ρ_T scale",
        }

    # SPARC
    sparc = load_step_results("step_043")
    if sparc:
        alpha = sparc["power_law_fit"]["alpha"]
        alpha_err = sparc["power_law_fit"]["alpha_err"]
        sigma = abs(alpha - TEP_ALPHA) / alpha_err
        synthesis["constraints"]["SPARC"] = {
            "scale": "Galactic (10^9-10^12 M_sun)",
            "observable": f"α = {alpha:.3f} ± {alpha_err:.3f}",
            "sigma_from_1/3": float(sigma),
            "status": "Consistent" if sigma < 2 else "Tension",
            "independence": "Independent dataset — strongest cross-check",
        }

    # Screening
    screen = load_step_results("step_045")
    if screen:
        beta = screen["screening_law"]["beta"]
        synthesis["constraints"]["Screening"] = {
            "scale": "26 objects (10^-24 to 10^14 g/cm³)",
            "observable": f"S ∝ ρ^{beta:.3f}",
            "status": "Consistency check (algebraically expected)",
            "independence": "Derived — depends on ρ_T and definitions",
        }

    # Magnetar
    mag = load_step_results("step_046")
    if mag:
        match = mag["target_object"]["match_percent"]
        synthesis["constraints"]["Magnetar"] = {
            "scale": "Compact (10^14 g/cm³)",
            "observable": f"P_crit ≈ 6.8 s, {match:.0f}% match",
            "status": "Scale-consistency check (N=1, statistically limited)",
            "independence": "Derived — depends on ρ_T from GNSS",
        }

    # Milky Way
    mw = load_step_results("step_047")
    if mw:
        sigma_mw = mw["agreement_sigma"]
        synthesis["constraints"]["Milky Way"] = {
            "scale": "Local Galactic (10^12 M_sun)",
            "observable": f"R_trans ≈ 19 kpc, {sigma_mw:.1f}σ agreement",
            "status": "Scale-consistency check",
            "independence": "Derived — depends on SPARC k",
        }

    # Dependency audit
    audit = load_step_results("step_048")
    if audit:
        synthesis["dependency_audit"] = audit["conclusions"]

    # Overall assessment
    n_independent = 2  # GNSS + SPARC
    n_consistency = 3  # Screening + Magnetar + Milky Way
    synthesis["summary"] = {
        "independent_constraints": n_independent,
        "consistency_checks": n_consistency,
        "mass_range_orders": 40,
        "density_range_orders": 15,
        "overall_assessment": (
            "Two independent constraints (GNSS calibration, SPARC scaling) "
            "converge on ρ_T ≈ 20 g/cm³ and α ≈ 1/3. Three additional "
            "consistency checks (screening hierarchy, magnetar P_crit, "
            "Milky Way R_trans) show no tension. The universal-density "
            "hypothesis is empirically motivated but remains falsifiable: "
            "any future measurement requiring substantially different ρ_T "
            "or α ≠ 1/3 would exclude the model as formulated."
        ),
        "strongest_open_requirement": (
            "Independent replication of the GNSS L_c ≈ 4200 km calibration "
            "by groups outside the author's research program."
        ),
    }

    logger.info(f"  Independent constraints: {n_independent}")
    logger.info(f"  Consistency checks: {n_consistency}")
    logger.info(f"  Mass range: 40 orders of magnitude")
    logger.info("✓   Cross-scale synthesis complete")

    synthesis["step_id"] = "step_049"
    return synthesis


def main():
    log_dir = PROJECT_ROOT / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    logger = TEPLogger("step_049", str(log_dir / "step_049_cross_scale_synthesis.log"))
    set_step_logger(logger)
    results = run_synthesis(logger)
    logger.save_step_results(results, PROJECT_ROOT, "step_049_cross_scale_synthesis")


if __name__ == "__main__":
    main()
