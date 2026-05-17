#!/usr/bin/env python3
"""
Step 075: Generic PPN Falsification Test
============================================

Tests whether the LLR data requires ANY scalar-tensor modification,
without invoking TEP.  Uses the standard PPN Nordtvedt parameter:

    η = 4β - γ - 3 - (10/3)ξ - α₁ + (2/3)α₂ - (2/3)ζ₁ - (1/3)ζ₂

In LLR, the observable is the synodic modulation amplitude, which in
standard PPN is a single static parameter.  With α₁=α₂=ξ=ζ₁=ζ₂=0 from
other experiments, this reduces to:

    η = 4β - γ - 3

Standard GR predicts η_GR = 0 (β = γ = 1).  We test:

  1. Is the measured η consistent with η = 0 at the 5σ level?
  2. Given external γ constraints (Cassini), what β is implied?
  3. Is the implied β consistent with β = 1?
  4. Joint (β, γ) confidence contours.

If the measured η deviates from zero and the implied β deviates from 1,
this establishes modified gravity independent of any TEP-specific
framework.
"""

import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

import json
import numpy as np
from scipy import stats

from scripts.utils.llr_constants import ETA_SCALE_FACTOR
from scripts.utils.logger import TEPLogger, set_step_logger, print_status


# External PPN constraints from Cassini and light deflection
# Cassini 2003: γ - 1 = (2.1 ± 2.3) × 10⁻⁵  (Bertotti et al. 2003)
# Perihelion precession (Mercury) constrains (2 + 2γ - β)/3, but less tightly
# than Cassini constrains γ.
CASSINI_GAMMA = 1.0
CASSINI_GAMMA_ERROR = 2.3e-5  # 1-sigma


def load_headline_eta():
    """Load the primary η estimate from Step 050 (corrected TEP analysis)."""
    path = PROJECT_ROOT / "results" / "outputs" / "step_050_corrected_tep_analysis.json"
    if path.exists():
        with open(path, "r") as f:
            data = json.load(f)
        # Try the precision-weighted full-systematic estimate (headline)
        pwr = data.get("precision_weighted_full_systematic", {})
        if pwr.get("eta") is not None:
            return float(pwr["eta"]), float(pwr["eta_error"])
        # Fallback to cooks-excised
        ols = data.get("cooks_excised_full_systematic", {})
        if ols.get("eta") is not None:
            return float(ols["eta"]), float(ols["eta_error"])
    # Fallback to Step 003
    path = PROJECT_ROOT / "results" / "outputs" / "step_003_statistical_analysis.json"
    if path.exists():
        with open(path, "r") as f:
            data = json.load(f)
        eta = data.get("eta", None)
        eta_err = data.get("eta_error", None)
        if eta is not None and eta_err is not None:
            return float(eta), float(eta_err)
    raise FileNotFoundError("No headline η found. Run Step 003 or 050 first.")


def compute_beta_from_eta(eta, eta_err, gamma, gamma_err):
    """
    From η = 4β - γ - 3, solve for β:
        β = (η + γ + 3) / 4
    """
    beta = (eta + gamma + 3) / 4.0
    var_beta = (eta_err ** 2 + gamma_err ** 2) / 16.0
    beta_err = np.sqrt(var_beta)
    return beta, beta_err


def test_gr_null(eta, eta_err):
    """Test H0: η = 0 (standard GR)."""
    z_gr = abs(eta) / eta_err if eta_err > 0 else 0.0
    p_gr = 2 * (1 - stats.norm.cdf(z_gr)) if z_gr > 0 else 1.0
    return {
        "gr_eta_expected": 0.0,
        "eta_measured": float(eta),
        "eta_error": float(eta_err),
        "z_score_vs_gr": float(z_gr),
        "p_value_two_tailed": float(p_gr),
        "gr_rejected_at_5sigma": bool(z_gr >= 5.0),
        "gr_rejected_at_3sigma": bool(z_gr >= 3.0),
        "conclusion": (
            f"Measured η = {eta:.4e} ± {eta_err:.4e} deviates from GR (η=0) "
            f"at {z_gr:.2f}σ (p = {p_gr:.2e}). "
            + ("GR is rejected at 5σ." if z_gr >= 5.0 else
               "GR is rejected at 3σ." if z_gr >= 3.0 else
               "GR is not rejected at 3σ.")
        ),
    }


def test_beta_given_gamma(eta, eta_err, gamma, gamma_err):
    """Test H0: β = 1 given external γ constraint."""
    beta, beta_err = compute_beta_from_eta(eta, eta_err, gamma, gamma_err)
    z_beta = abs(beta - 1.0) / beta_err if beta_err > 0 else 0.0
    p_beta = 2 * (1 - stats.norm.cdf(z_beta)) if z_beta > 0 else 1.0
    return {
        "beta_implied": float(beta),
        "beta_error": float(beta_err),
        "gamma_assumed": float(gamma),
        "gamma_error": float(gamma_err),
        "z_score_vs_beta_1": float(z_beta),
        "p_value_two_tailed": float(p_beta),
        "beta_1_rejected_at_5sigma": bool(z_beta >= 5.0),
        "beta_1_rejected_at_3sigma": bool(z_beta >= 3.0),
        "conclusion": (
            f"Given γ = {gamma:.6f} ± {gamma_err:.2e} (Cassini), "
            f"the LLR η implies β = {beta:.6f} ± {beta_err:.6f}. "
            f"This deviates from GR (β=1) at {z_beta:.2f}σ (p = {p_beta:.2e})."
        ),
    }


def joint_beta_gamma_contour(eta, eta_err, gamma_cassini, gamma_err, n_grid=200):
    """
    Compute joint (β, γ) confidence contours analytically.

    chi²(β,γ) = (4β - γ - 3 - η)² / η_err² + (γ - γ_Cassini)² / γ_err²

    The minimum is at β = 1 + η/4, γ = γ_Cassini (exactly, since the two
    constraints are compatible).  We evaluate contours analytically via
    the Hessian at the minimum rather than a coarse grid that misses the
    narrow valley.
    """
    # Analytical best fit
    beta_best = 1.0 + eta / 4.0
    gamma_best = float(gamma_cassini)
    chi2_min = 0.0  # exact at minimum

    # GR chi²
    gr_chi2 = float(
        ((4 * 1.0 - 1.0 - 3 - eta) ** 2) / (eta_err ** 2)
        + ((1.0 - gamma_cassini) ** 2) / (gamma_err ** 2)
    )

    # 2-DOF Δchi² thresholds
    delta_chi2 = {"68_percent": 2.30, "95_percent": 5.99, "99_percent": 9.21}
    levels = {k: chi2_min + v for k, v in delta_chi2.items()}
    gr_inside = {k: bool(gr_chi2 <= v) for k, v in levels.items()}

    # Covariance from Hessian at minimum
    # H = [[32, -8], [-8, 2]] / η_err²  +  [[0, 0], [0, 2]] / γ_err²
    H_11 = 32.0 / (eta_err ** 2)
    H_22 = 2.0 / (eta_err ** 2) + 2.0 / (gamma_err ** 2)
    H_12 = -8.0 / (eta_err ** 2)
    H = np.array([[H_11, H_12], [H_12, H_22]])
    cov = np.linalg.inv(H)
    beta_err = float(np.sqrt(cov[0, 0]))
    gamma_err_joint = float(np.sqrt(cov[1, 1]))

    # Also grid for visualization / sanity check (fine grid around minimum)
    half_width = max(8 * beta_err, 8 * gamma_err_joint, 0.001)
    beta_range = np.linspace(beta_best - half_width, beta_best + half_width, n_grid)
    gamma_range = np.linspace(gamma_best - half_width, gamma_best + half_width, n_grid)
    B, G = np.meshgrid(beta_range, gamma_range)
    chi2_llr = ((4 * B - G - 3 - eta) ** 2) / (eta_err ** 2)
    chi2_gamma = ((G - gamma_cassini) ** 2) / (gamma_err ** 2)
    chi2_total = chi2_llr + chi2_gamma
    grid_min = float(np.min(chi2_total))

    return {
        "beta_best_fit": float(beta_best),
        "beta_error_joint": beta_err,
        "gamma_best_fit": gamma_best,
        "gamma_error_joint": gamma_err_joint,
        "chi2_min": chi2_min,
        "grid_chi2_min_sanity_check": grid_min,
        "gr_chi2": gr_chi2,
        "gr_inside_68_percent": gr_inside["68_percent"],
        "gr_inside_95_percent": gr_inside["95_percent"],
        "gr_inside_99_percent": gr_inside["99_percent"],
        "conclusion": (
            f"Joint (β,γ) best fit: β={beta_best:.6f} ± {beta_err:.2e}, "
            f"γ={gamma_best:.6f} ± {gamma_err_joint:.2e}. "
            f"GR (β=γ=1) has Δchi²={gr_chi2:.1f}. "
            + (
                "GR lies inside 68% contour."
                if gr_inside["68_percent"]
                else (
                    "GR lies inside 95% contour — marginally excluded."
                    if gr_inside["95_percent"]
                    else "GR is outside 99% contour — rejected at >3σ."
                )
            )
        ),
    }


def run_generic_ppn_test() -> dict:
    print_status("═══ Step 075: Generic PPN Falsification Test ═══", "TITLE")

    eta, eta_err = load_headline_eta()
    print_status(f"Headline η = {eta:.4e} ± {eta_err:.4e}", "DATA")

    # 1. Test GR null: η = 0
    gr_test = test_gr_null(eta, eta_err)
    print_status(gr_test["conclusion"], "RESULT")

    # 2. Test β = 1 given Cassini γ
    beta_test = test_beta_given_gamma(eta, eta_err, CASSINI_GAMMA, CASSINI_GAMMA_ERROR)
    print_status(beta_test["conclusion"], "RESULT")

    # 3. Joint (β, γ) contour
    contour = joint_beta_gamma_contour(eta, eta_err, CASSINI_GAMMA, CASSINI_GAMMA_ERROR)
    print_status(contour["conclusion"], "RESULT")

    # Composite verdict
    gr_rejected = gr_test["gr_rejected_at_3sigma"]
    beta_rejected = beta_test["beta_1_rejected_at_3sigma"]

    if gr_rejected and beta_rejected:
        verdict = (
            "MODIFIED_GRAVITY_DETECTED: Both η ≠ 0 and β ≠ 1 are established "
            "at >3σ using standard PPN formalism, independent of TEP."
        )
        status = "PASS"
    elif gr_rejected:
        verdict = (
            "ETA_NONZERO_CONFIRMED: η deviates from GR at >3σ, but β = 1 "
            "is not independently rejected given γ uncertainty. "
            "Modified gravity is suggested but not conclusively established."
        )
        status = "PASS"
    else:
        verdict = (
            "GR_CONSISTENT: The measured η is consistent with GR (η=0) "
            "at the 3σ level. No PPN modification is required."
        )
        status = "FAIL"

    print_status(verdict, "RESULT")

    return {
        "step_id": "step_075",
        "status": status,
        "eta_headline": float(eta),
        "eta_headline_error": float(eta_err),
        "gr_null_test": gr_test,
        "beta_given_gamma_test": beta_test,
        "joint_beta_gamma_contour": contour,
        "ppn_parameter_constraints": {
            "gamma_cassini": CASSINI_GAMMA,
            "gamma_cassini_error": CASSINI_GAMMA_ERROR,
            "alpha_1": 0.0,
            "alpha_2": 0.0,
            "xi": 0.0,
            "zeta_1": 0.0,
            "zeta_2": 0.0,
            "constraint_source": "Cassini 2003 (Bertotti et al.), LLR Nordtvedt",
        },
        "interpretation": (
            "This step tests whether the LLR signal requires modified gravity "
            "in the standard PPN framework, completely independent of TEP. "
            f"The measured η = {eta:.3e} ± {eta_err:.3e} yields a PPN-interpreted "
            f"β = {beta_test['beta_implied']:.5f} ± {beta_test['beta_error']:.5f}. "
            "If β deviates from 1, scalar-tensor theories are implicated regardless "
            "of whether the underlying mechanism is TEP or another framework."
        ),
    }


def main() -> int:
    results = run_generic_ppn_test()
    output_path = PROJECT_ROOT / "results" / "outputs" / "step_075_generic_ppn_test.json"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8") as handle:
        json.dump(results, handle, indent=2)
    print_status(f"Saved results to {output_path.relative_to(PROJECT_ROOT)}", "SUCCESS")
    return 0


if __name__ == "__main__":
    log_dir = PROJECT_ROOT / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    logger = TEPLogger(
        "step_075", str(log_dir / "step_075_generic_ppn_test.log")
    )
    set_step_logger(logger)
    sys.exit(main())
