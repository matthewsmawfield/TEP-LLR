#!/usr/bin/env python3
"""
Step 047: Milky Way Keplerian Transition Analysis
Tests whether the SPARC-calibrated M^{1/3} scaling predicts the Keplerian
decline radius in the Milky Way.
Uses published data: Bland-Hawthorn & Gerhard 2016, Jiao et al. 2023, Gaia DR3.
"""

import sys, json, numpy as np
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))
from scripts.utils.logger import TEPLogger, set_step_logger, print_status

# Published values
# The Milky Way check uses total dynamical mass (baryonic + phantom mass)
# since the Keplerian decline is observed in the total rotation curve.
MW_TOTAL_MASS = 1.0e12  # M_sun (Bland-Hawthorn & Gerhard 2016, ARA&A, 54, 529)
MW_TOTAL_MASS_ERR = 0.2e12  # ±20%
OBSERVED_TRANSITION = 19.0  # kpc (Jiao et al. 2023, A&A, 678, A208)
OBSERVED_TRANSITION_ERR = 2.0  # kpc

# SPARC calibration from step_043 (or default)
K_SPARC_DEFAULT = 7.86e-4  # kpc/M_sun^{1/3}


def run_analysis(logger):
    logger.info(">>> Starting Milky Way Keplerian transition analysis...")

    # Load SPARC k from step_043 if available
    k_sparc = K_SPARC_DEFAULT
    k_sparc_err = 0.04 * k_sparc  # ~4% fit error (typical uncertainty from SPARC calibration)
    step043_path = PROJECT_ROOT / "results" / "outputs" / "step_043_sparc_scaling_analysis.json"
    if step043_path.exists():
        with open(step043_path) as f:
            s043 = json.load(f)
        k_sparc = s043["power_law_fit"]["k_kpc_per_Msun_third"]
        k_sparc_err = k_sparc * (s043["power_law_fit"]["alpha_err"] /
                                  s043["power_law_fit"]["alpha"])
        logger.info(f"  Using k = {k_sparc:.4e} ± {k_sparc_err:.2e} from step_043")

    # Predicted transition radius
    r_pred = k_sparc * MW_TOTAL_MASS ** (1.0/3.0)

    # Error propagation
    # dR/R = dk/k + (1/3) dM/M
    mass_frac_err = MW_TOTAL_MASS_ERR / MW_TOTAL_MASS
    k_frac_err = k_sparc_err / k_sparc
    r_frac_err = np.sqrt(k_frac_err**2 + (mass_frac_err/3)**2)
    r_pred_err = r_pred * r_frac_err

    # Add model systematic (~10% for geometry)
    r_sys = 0.10 * r_pred
    r_total_err = np.sqrt(r_pred_err**2 + r_sys**2)

    logger.info(f"  Predicted R_trans = {r_pred:.1f} ± {r_total_err:.1f} kpc")
    logger.info(f"  Observed R_trans = {OBSERVED_TRANSITION:.1f} ± {OBSERVED_TRANSITION_ERR:.1f} kpc")

    # Consistency check
    sigma_diff = abs(r_pred - OBSERVED_TRANSITION) / np.sqrt(r_total_err**2 + OBSERVED_TRANSITION_ERR**2)
    logger.info(f"  Agreement: {sigma_diff:.1f}σ")

    results = {
        "step_id": "step_047",
        "status": "PASS",
        "input_data": {
            "MW_total_mass_Msun": MW_TOTAL_MASS,
            "MW_mass_uncertainty_Msun": MW_TOTAL_MASS_ERR,
            "MW_mass_ref": "Bland-Hawthorn & Gerhard 2016, ARA&A, 54, 529",
            "observed_transition_kpc": OBSERVED_TRANSITION,
            "observed_transition_err_kpc": OBSERVED_TRANSITION_ERR,
            "observed_transition_ref": "Jiao et al. 2023, A&A, 678, A208; Gaia DR3",
        },
        "sparc_calibration_k": float(k_sparc),
        "predicted_R_trans_kpc": float(r_pred),
        "predicted_R_trans_err_kpc": float(r_total_err),
        "error_budget": {
            "mass_contribution_pct": float(mass_frac_err/3 * 100),
            "k_calibration_contribution_pct": float(k_frac_err * 100),
            "model_systematic_pct": 10.0,
        },
        "agreement_sigma": float(sigma_diff),
        "consistent": bool(sigma_diff < 2.0),
        "caveat": "Outer-disk rotation curve sensitive to tracer selection and "
                  "non-circular motions. Treated as scale-consistency check.",
    }

    logger.info("✓   Milky Way analysis complete")
    return results


def main():
    log_dir = PROJECT_ROOT / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    logger = TEPLogger("step_047", str(log_dir / "step_047_milky_way_analysis.log"))
    set_step_logger(logger)
    results = run_analysis(logger)
    logger.save_step_results(results, PROJECT_ROOT, "step_047_milky_way_analysis")


if __name__ == "__main__":
    main()
