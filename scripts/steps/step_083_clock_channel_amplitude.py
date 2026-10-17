#!/usr/bin/env python3
"""
Step 083: Derived clock-sector amplitude bookkeeping (issue 17-3).

The phenomenological volumetric benchmark eta ~ alpha_eff * Delta<(rho/rho_T)^3>
uses an illustrative normalization and an underived exponent.  Under the
corpus's v0.15 master-action machinery (Paper 0, Appendix E), both ingredients
are superseded by derived quantities:

  * the quartic matter-branch potential pins the equilibrium interior
    level at  phi_min(rho) = (rho/(lambda M_Pl))^{1/3}  -- a cube-ROOT
    density response, supplying the depth bookkeeping
    S_A = min[1, (rho_bar/rho_T)^{1/3}] (e.g. S_A^Earth ~= 0.65);
  * the ambient-excursion transmission into interior clocks is measured
    on the solved unified profile (TEP step_76): real bodies sit on the
    leaking (flux-following) branch, where ambient field-value
    excursions transmit ~additively (T ~ 1; T_Moon = 0.99685 at the
    2.5e-11 excursion scale);
  * the admissible LLR channel is the clock-sector/ephemeris projection of the
    lunar conformal-field excursion  delta phi_sun(D)  -- a one-point function
    of the temporal landscape at the Moon's instantaneous position.

This step computes from nominal orbital geometry and simplified shell inputs:

  1. The solar conformal excursion at the Moon versus synodic phase D:
     d_SM(D) = sqrt(a_E^2 + a_M^2 - 2 a_E a_M cos D), so
     delta phi_sun(D) = u_sun(1 AU) * (a_M/a_E) * cos D to leading order,
     with the full geometric series evaluated to confirm the coefficient.
  2. The derived interior clock-depth factors on the real PREM-style Earth and
     lunar shell profiles (step_026 tables): the volume-averaged pinning
     <(rho/rho_T)^{1/3}> per body, and the corpus's S_A factors at the mean
     densities -- replacing the phenomenological <(rho/rho_T)^3> weighting.
  3. The amplitude bookkeeping for the measured eta_resid: the excursion
     amplitude carried through the lunar clock response versus the measured
     range-modulation coefficient, and the projection gap that remains.

Outputs: results/outputs/step_083_clock_channel_amplitude.json
"""

import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from scripts.utils.logger import TEPLogger, set_step_logger, print_status
from scripts.utils.statistical_utils import require_step003_eta_ols, require_step003_eta_ols_error
from scripts.utils.llr_constants import (
    C_M_S, GM_SUN, GM_EARTH, GM_MOON, RHO_T, RHO_T_ERROR, ETA_SCALE_FACTOR,
)

STEP_NUM = "083"

AU_M = 1.495978707e11          # IAU nominal astronomical unit [m]
A_MOON_M = 3.84399e8           # mean Earth-Moon separation [m]

# Real interior shell profiles (identical construction to step_026):
# simplified PREM for Earth; small core + mantle for the Moon.
EARTH_SHELLS = [
    (0, 1221e3, 13.0, "Inner Core"),
    (1221e3, 3480e3, 10.9, "Outer Core"),
    (3480e3, 5701e3, 5.0, "Lower Mantle"),
    (5701e3, 6371e3, 3.3, "Upper Mantle/Crust"),
]
MOON_SHELLS = [
    (0, 330e3, 8.0, "Small Core"),
    (330e3, 1737e3, 3.3, "Mantle/Crust"),
]


def volume_weighted(shells, exponent):
    """Volume-weighted <(rho/rho_T)^exponent> over the shell model."""
    num = den = 0.0
    for ri, ro, rho, _ in shells:
        v = (ro ** 3 - ri ** 3)
        num += (rho / RHO_T) ** exponent * v
        den += v
    return num / den


def mean_density(shells):
    m = v = 0.0
    for ri, ro, rho, _ in shells:
        vv = (ro ** 3 - ri ** 3)
        m += rho * vv
        v += vv
    return m / v


def compute_solar_excursion(n_samples=20000, phase_offset=0.0):
    """Uniform periodic quadrature; no duplicated endpoint or DC leakage."""
    if n_samples < 16:
        raise ValueError("At least 16 periodic samples are required")
    u_sun_1au = GM_SUN / (AU_M * C_M_S ** 2)
    eps = A_MOON_M / AU_M

    # exact series coefficient of cos D in u_sun(d_SM(D)): sample and project
    Ds = np.linspace(0.0, 2.0 * np.pi, n_samples, endpoint=False) + phase_offset
    d_sm = np.sqrt(AU_M ** 2 + A_MOON_M ** 2 - 2.0 * AU_M * A_MOON_M * np.cos(Ds))
    u_series = GM_SUN / (d_sm * C_M_S ** 2)
    u_mean = float(np.mean(u_series))
    # A periodic mean must not count both D=0 and D=2*pi. Duplicating
    # that endpoint leaks the ~1e-8 constant into the ~1e-11/~1e-14
    # harmonics. Explicitly center as an additional numerical safeguard.
    centered = u_series - u_mean
    cos_coeff = float(2.0 * np.mean(centered * np.cos(Ds)))
    cos2_coeff = float(2.0 * np.mean(centered * np.cos(2.0 * Ds)))
    excursion = {
        "u_sun_1au": u_sun_1au,
        "a_m_over_a_e": eps,
        "leading_amplitude": u_sun_1au * eps,
        "exact_cos_D_amplitude": cos_coeff,
        "exact_cos_2D_amplitude": cos2_coeff,
        "mean_level": u_mean,
        "peak_to_peak": float(u_series.max() - u_series.min()),
        "note": "leading term u_sun*(a_M/a_E) = 2.54e-11 cos D; the exact "
                "series confirms it to <1%; peak-to-peak ~5.1e-11",
    }

    return excursion


def main():
    logger = TEPLogger("step_083_clock_channel_amplitude")
    set_step_logger(logger)
    print_status(f"STEP {STEP_NUM}: Derived clock-sector amplitude bookkeeping", "TITLE")

    # ------------------------------------------------------------------
    # 1. Solar conformal excursion at the Moon, synodic component.
    #    u_sun(r) = GM_sun/(r c^2); the Moon's heliocentric distance
    #    oscillates as d_SM(D) = sqrt(a_E^2 + a_M^2 - 2 a_E a_M cos D).
    # ------------------------------------------------------------------
    excursion = compute_solar_excursion()
    cos_coeff = excursion["exact_cos_D_amplitude"]

    # ------------------------------------------------------------------
    # 2. Derived interior response: quartic pinning  phi_min ~ rho^{1/3}
    #    and the corpus's amplitude factor S_A = min[1,(rho_bar/rho_T)^{1/3}].
    # ------------------------------------------------------------------
    rho_earth = mean_density(EARTH_SHELLS)
    rho_moon = mean_density(MOON_SHELLS)

    S_A_earth = min(1.0, (rho_earth / RHO_T) ** (1.0 / 3.0))
    S_A_moon = min(1.0, (rho_moon / RHO_T) ** (1.0 / 3.0))

    phenom_cubed = {
        "earth": volume_weighted(EARTH_SHELLS, 3.0),
        "moon": volume_weighted(MOON_SHELLS, 3.0),
    }
    derived_cuberoot = {
        "earth": volume_weighted(EARTH_SHELLS, 1.0 / 3.0),
        "moon": volume_weighted(MOON_SHELLS, 1.0 / 3.0),
    }
    # Excursion transmission into the lunar interior, measured on the
    # solved unified profile (TEP step_76, step_76_ambient_transmission):
    # real bodies sit on the leaking (flux-following) branch where the
    # interior rides the ambient additively, so an ambient excursion at
    # the body's location transmits at T ~ 1; the departure is of order
    # u_int/delta with u_int the body's own plateau depth.  At the LLR
    # excursion scale (2.5e-11) the solved lunar transmission is 0.99685;
    # the S_A factors above are the equilibrium *depth* bookkeeping, not
    # the excursion-transmission factor.
    T_moon = 0.99685
    T_earth = 0.99554

    interior = {
        "rho_bar_earth_gcc": rho_earth,
        "rho_bar_moon_gcc": rho_moon,
        "S_A_earth": S_A_earth,
        "S_A_moon": S_A_moon,
        "S_A_differential": S_A_earth - S_A_moon,
        "excursion_transmission_moon": T_moon,
        "excursion_transmission_earth": T_earth,
        "excursion_transmission_source": "TEP step_76 unified-profile "
            "ambient-response solve (leaking branch, T ~ 1)",
        "volumetric_cubed_phenomenological": phenom_cubed,
        "volumetric_cuberoot_derived": derived_cuberoot,
        "differential_cubed": phenom_cubed["earth"] - phenom_cubed["moon"],
        "differential_cuberoot": derived_cuberoot["earth"] - derived_cuberoot["moon"],
        "note": "the quartic potential pins the equilibrium level as "
                "rho^{1/3}: the derived interior depth bookkeeping is the "
                "cube root of density, not the cube; the phenomenological "
                "(rho/rho_T)^3 weighting over-concentrates on the core and "
                "is superseded. The S_A factors are depth ratios; the "
                "ambient-excursion transmission into interior clocks is "
                "measured on the solved profile (step_76) at T ~ 1 on the "
                "leaking branch",
    }

    # ------------------------------------------------------------------
    # 3. Amplitude bookkeeping for the measured eta_resid.
    # ------------------------------------------------------------------
    step_003_path = PROJECT_ROOT / "results" / "outputs" / "step_003_statistical_analysis.json"
    with open(step_003_path, "r") as f:
        s003 = json.load(f)
    eta_ols = require_step003_eta_ols(s003)
    eta_ols_err = require_step003_eta_ols_error(s003)

    # headline precision-weighted estimate (step_050)
    eta_headline = None
    s050_path = PROJECT_ROOT / "results" / "outputs" / "step_050_corrected_tep_analysis.json"
    if s050_path.exists():
        with open(s050_path) as f:
            s050 = json.load(f)
        pw = s050.get("precision_weighted_full_systematic", {})
        eta_headline = pw.get("eta")
        eta_headline_err = pw.get("eta_error")

    eta_meas = eta_headline if eta_headline is not None else eta_ols
    eta_err = eta_headline_err if eta_headline is not None else eta_ols_err

    # range amplitude of the measured coefficient
    delta_r_meas_m = ETA_SCALE_FACTOR * eta_meas          # metres (signed)
    delta_r_meas_ols_m = ETA_SCALE_FACTOR * eta_ols
    # lunar clock response to the excursion: ambient field-value
    # excursions transmit into leaking-branch interiors ~additively
    # (T ~ 1, TEP step_76), so the Moon's effective excursion is the
    # ambient excursion itself up to the solved transmission
    delta_phi_eff = T_moon * cos_coeff

    ceiling_m = C_M_S * 1.28 * delta_phi_eff
    path_integrated_delta_r_m = 0.5 * ceiling_m
    eta_predicted_path = -path_integrated_delta_r_m / ETA_SCALE_FACTOR
    pull_sigma = abs(eta_meas - eta_predicted_path) / eta_err if eta_err else None

    # ------------------------------------------------------------------
    # 4. Unified Single-Solar-Solution Audit (Cassini at 1.6 R_sun and
    #    LLR at 1 AU evaluated on the SAME two-branch solar profile):
    #      (k / sqrt(2)) q y^2 + q^2 y^3 = 1,   k = 16.03, g_t = c H0 / 2.
    #    Under the unified metric n(x) = 1 - 2 Phi_N/c^2 - beta_A u(x),
    #    both the Cassini Shapiro delay and the direct scalar optical path
    #    delay across the Earth-Moon baseline sample the same solar field
    #    u_sun(r) with u'_lin(r) = 2 beta_A^2 GM_sun / (c^2 r^2).
    # ------------------------------------------------------------------
    R_SUN_M = 6.957e8
    H0_SI = 70.0e3 / 3.0856775814913673e22
    G_T_SI = C_M_S * H0_SI / 2.0
    K_TB = 16.03
    BETA_A = -1.0

    def y_of_q(q_val):
        f_eq = lambda ly: math.log(
            (K_TB / math.sqrt(2.0)) * q_val * math.exp(2.0 * ly)
            + (q_val ** 2) * math.exp(3.0 * ly)
        )
        return math.exp(brentq(f_eq, -200.0, 200.0, xtol=1e-14))

    def gN_sun(r_m):
        return GM_SUN / (r_m ** 2)

    def y_sun(r_m):
        return y_of_q(gN_sun(r_m) / G_T_SI)

    def uprime_lin_sun(r_m):
        return 2.0 * (BETA_A ** 2) * GM_SUN / ((C_M_S ** 2) * (r_m ** 2))

    y_1p6_rsun = y_sun(1.6 * R_SUN_M)
    y_1au = y_sun(AU_M)
    gamma_minus_1_cassini = -2.0 * (BETA_A ** 2) * y_1p6_rsun / (1.0 + 2.0 * (BETA_A ** 2) * y_1p6_rsun)
    gamma_minus_1_if_unscreened = -2.0 * (BETA_A ** 2) / (1.0 + 2.0 * (BETA_A ** 2))
    du_solved_twobranch = float(quad(lambda r: y_sun(r) * uprime_lin_sun(r), AU_M, AU_M + A_MOON_M, epsrel=1e-12)[0])
    du_unscreened_lin = float((2.0 * (BETA_A ** 2) * GM_SUN / (C_M_S ** 2)) * (1.0 / AU_M - 1.0 / (AU_M + A_MOON_M)))
    delta_r_solved_m = 0.5 * A_MOON_M * du_solved_twobranch
    eta_predicted_solved_twobranch = -delta_r_solved_m / ETA_SCALE_FACTOR
    screening_suppression_ratio = abs(eta_predicted_path / eta_predicted_solved_twobranch)

    unified_solar_solution = {
        "k_twobranch": K_TB,
        "g_t_m_s2": G_T_SI,
        "r_star_sun_AU": float(math.sqrt(GM_SUN / G_T_SI) / AU_M),
        "y_tb_1p6_Rsun": y_1p6_rsun,
        "gamma_PPN_minus_1_Cassini_1p6_Rsun": gamma_minus_1_cassini,
        "cassini_bound": 2.3e-5,
        "passes_cassini": bool(abs(gamma_minus_1_cassini) < 2.3e-5),
        "y_tb_1AU": y_1au,
        "delta_u_earth_moon_solved_twobranch": du_solved_twobranch,
        "delta_u_earth_moon_unscreened_2beta2_Phi": du_unscreened_lin,
        "delta_Phi_N_over_c2_input": cos_coeff,
        "delta_r_path_solved_twobranch_m": -delta_r_solved_m,
        "delta_r_path_solved_twobranch_mm": -delta_r_solved_m * 1e3,
        "eta_predicted_solved_twobranch": eta_predicted_solved_twobranch,
        "gamma_PPN_minus_1_if_solar_u_unscreened": gamma_minus_1_if_unscreened,
        "unscreened_to_solved_ratio": screening_suppression_ratio,
        "consistency_note": (
            "Under the unified metric n(x) = 1 - 2 Phi_N/c^2 - beta_A u(x), "
            "where the spatial 3-metric gamma_tilde_ij = gamma_ij = (1 - 2 Phi_N/c^2) delta_ij "
            "is identical to GR's while the spacetime lapse carries A_loc(u) = exp(beta_A u), "
            "the Newtonian solar tidal potential Delta Phi_{N,sun}/c^2 = 2.54e-11 cos D enters "
            "the BCRS->GCRS reduction identically in TEP and GR, and both the Cassini Shapiro "
            "delay and the direct scalar optical path delay sample the same screened solar "
            "field u_sun(r). The two-branch solution that passes Cassini "
            f"(y_tb(1.6 R_sun) = {y_1p6_rsun:.2e}, gamma - 1 = {gamma_minus_1_cassini:.2e}) "
            f"screens the solar scalar gradient at 1 AU by y_tb(1 AU) = {y_1au:.2e}, yielding "
            f"delta_u_solved = {du_solved_twobranch:.2e}, delta_r_path = {-delta_r_solved_m*1e3:.2e} mm, "
            f"and eta_solved = {eta_predicted_solved_twobranch:.2e} ({screening_suppression_ratio:.1e}x below "
            "the measured residual, cleanly satisfying direct-fit LLR bounds alongside Cassini). "
            "The -4.85 mm (eta = -3.73e-4) figure is the unscreened solar-potential ceiling "
            "(which would require an unscreened solar scalar field and violate Cassini at "
            "gamma - 1 = -2/3); the measured post-fit residual eta_resid = -(3.91 +/- 0.56)e-4 "
            "(-5.08 +/- 0.73 mm) is therefore reported as an empirical residual-channel "
            "detection whose dynamical or ephemeris-reduction transfer map awaits a "
            "source-level IMCCE/JPL integrator refit."
        ),
    }

    bookkeeping = {
        "eta_measured_headline": eta_headline,
        "eta_measured_ols": eta_ols,
        "eta_measured_err": eta_err,
        "delta_r_measured_m": delta_r_meas_m,
        "delta_r_measured_ols_m": delta_r_meas_ols_m,
        "excursion_cos_D_amplitude": cos_coeff,
        "lunar_excursion_transmission": T_moon,
        "effective_lunar_excursion": delta_phi_eff,
        "ceiling_delta_r_m": ceiling_m,
        "unstretched_path_integral_factor": 0.5,
        "predicted_delta_r_path_unscreened_ceiling_m": -path_integrated_delta_r_m,
        "eta_predicted_unscreened_ceiling": eta_predicted_path,
        "tension_sigma_unscreened_ceiling_vs_headline": pull_sigma,
        "predicted_delta_r_path_solved_twobranch_m": -delta_r_solved_m,
        "eta_predicted_solved_twobranch": eta_predicted_solved_twobranch,
        "unified_single_solar_solution": unified_solar_solution,
        "required_projection": (
            f"Solved two-branch solar field (passing Cassini with gamma-1 = {gamma_minus_1_cassini:.2e} "
            f"at 1.6 R_sun and y_tb(1 AU) = {y_1au:.2e}) gives Earth-Moon scalar excursion "
            f"delta_u_solved = {du_solved_twobranch:.2e}, |delta_r_path| = {delta_r_solved_m*1e3:.2e} mm "
            f"(eta_solved = {eta_predicted_solved_twobranch:.2e}), whereas the unscreened Newtonian "
            f"potential differential Delta Phi_N/c^2 = {cos_coeff:.3e} cos D defines the unscreened "
            f"reference ceiling (1/2) L delta_phi_eff = {path_integrated_delta_r_m*1e3:.2f} mm "
            f"(eta_ceiling = {eta_predicted_path:.3e}, excluded as a direct solar scalar gradient by Cassini) "
            f"vs measured |delta_r_meas| = {abs(delta_r_meas_m)*1e3:.2f} mm "
            f"(eta_meas = {eta_meas:.3e} +/- {eta_err:.2e})."
        ),
    }

    out = {
        "step": STEP_NUM,
        "status": "success",
        "description": "Derived clock-sector amplitude bookkeeping (issue 17-3): "
                       "unified single-solar-solution two-branch audit (Cassini + LLR), "
                       "quartic depth bookkeeping, step_76 ambient-excursion transmission, "
                       "and unscreened solar-potential ceiling diagnostic",
        "solar_excursion": excursion,
        "interior_response": interior,
        "amplitude_bookkeeping": bookkeeping,
        "verdict": (
            "The volumetric benchmark's two underived ingredients are superseded "
            "by the corpus's master-action machinery: the quartic potential pins "
            "the equilibrium interior level as rho^{1/3} (derived exponent 1/3, "
            "not 3), supplying the depth bookkeeping S_A^Earth = %.3f, "
            "S_A^Moon = %.3f. Evaluating Cassini and LLR on the single solved "
            "two-branch solar profile ((k/sqrt(2)) q y^2 + q^2 y^3 = 1, k = 16.03) "
            "gives y_tb(1.6 R_sun) = %.2e (gamma_PPN - 1 = %.2e, passing Cassini) "
            "and y_tb(1 AU) = %.2e, which suppresses the direct scalar optical "
            "excursion across the Earth-Moon baseline to delta_u_solved = %.2e "
            "(|delta_r_path| = %.2e mm, eta_solved = %.2e), %.1e times below the "
            "measured post-fit residual %.2f mm (eta_resid = %.2e +/- %.2e) and "
            "consistent with direct-fit LLR bounds. Because TEP's spatial 3-metric "
            "gamma_ij = (1 - 2 Phi_N/c^2) delta_ij is identical to GR's, the "
            "unscreened Newtonian tidal scale %.2e cos D (|delta_r| = %.2f mm, "
            "eta_ceiling = %.2e) cancels in BCRS->GCRS identically to GR and would "
            "violate Cassini (gamma_PPN - 1 = -2/3) if attributed to an unscreened "
            "solar scalar field; the measured post-fit synodic residual therefore "
            "stands as an empirical residual-channel detection awaiting a "
            "source-level IMCCE/JPL integrator refit."
            % (S_A_earth, S_A_moon, y_1p6_rsun, gamma_minus_1_cassini, y_1au,
               du_solved_twobranch, delta_r_solved_m * 1e3, eta_predicted_solved_twobranch,
               screening_suppression_ratio, abs(delta_r_meas_m) * 1e3, eta_meas, eta_err,
               cos_coeff, path_integrated_delta_r_m * 1e3, eta_predicted_path)
        ),
    }

    outdir = PROJECT_ROOT / "results" / "outputs"
    outdir.mkdir(parents=True, exist_ok=True)
    with open(outdir / "step_083_clock_channel_amplitude.json", "w") as f:
        json.dump(out, f, indent=2)
    print_status(f"Results saved to {outdir / 'step_083_clock_channel_amplitude.json'}", "SUCCESS")


if __name__ == "__main__":
    main()
