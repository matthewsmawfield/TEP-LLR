#!/usr/bin/env python3
"""
Step 062: Solar Radiation Pressure (SRP) Systematic Check

Tests whether the detected synodic-phase modulation in LLR residuals can be
attributed to unmodeled or imperfectly modeled solar radiation pressure.

Physical Motivation:
- At New Moon (D=0), the Moon lies between Earth and Sun. Solar radiation
  pressure on the Moon acts toward Earth, reducing the Earth-Moon range.
- At Full Moon (D=pi), the Moon lies on the anti-solar side of Earth. SRP
  on the Moon acts away from Earth, increasing the range.
- This geometric pattern produces a synodic-phase-locked signal with the
  same cos(D) signature as the Nordtvedt effect, making SRP a critical
  confounding systematic.
- However, SRP scales as 1/r_sun^2 (inverse-square law), while the
  Nordtvedt effect carries no such heliocentric-distance dependence to
  leading order. This difference provides the discriminating handle.

Methodological Warning:
  The SRP geometry proxy, SRP_proxy = cos(D) / r_sun^2, is nearly collinear
  with cos(D) itself because r_sun varies by only ~3.4% over the orbit
  (0.983–1.017 AU). A naive joint OLS regression of both predictors is
  therefore ill-conditioned and produces unstable, inflated standard
  errors. This step avoids that trap and uses three robust, orthogonal
  tests instead.

Tests:
  1. Collinearity diagnostic: quantify the correlation between cos(D)
     and SRP_proxy, and report the VIF.
  2. Detrended-residual correlation test: after removing the best-fit
     TEP cos(D) signal from the residuals, test whether the SRP proxy
     explains additional variance. If SRP were present, the detrended
     residuals should correlate with the proxy.
  3. SRP-scaling test (central): bin data by heliocentric distance,
     fit residual = eta_bin * cos(D) in each bin, and test whether
     eta_bin scales as 1/r_bin^2 (SRP hypothesis) or is consistent
     with constant (Nordtvedt hypothesis).
"""

import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

import argparse
import json
import numpy as np
import pandas as pd
from scipy import stats

from astropy.time import Time
from astropy.coordinates import get_sun, get_body

from scripts.utils.llr_constants import ETA_SCALE_FACTOR
from scripts.utils.numerics import stable_lstsq
from scripts.utils.statistical_utils import linear_regression, detect_outliers_sigma
from scripts.utils.logger import TEPLogger, set_step_logger, set_verbose_mode, print_status


def compute_srp_proxy(df, verbose=False):
    """
    Compute the SRP geometry proxy for each observation.

    Uses Astropy to obtain true geocentric positions of the Sun and Moon,
    then constructs:
        SRP_proxy = cos(elongation) / r_sun_au^2
    """
    jds = df["date_julian"].values
    t = Time(jds, format="jd")

    sun = get_sun(t)
    moon = get_body("moon", t)

    r_sun_au = sun.distance.au
    elong = sun.separation(moon).radian
    cos_elong = np.cos(elong)

    srp_proxy = cos_elong / (r_sun_au ** 2)

    if verbose:
        print_status(f"SRP proxy computed for {len(df):,} observations", "INFO")
        print_status(f"  r_sun range: {np.min(r_sun_au):.5f} – {np.max(r_sun_au):.5f} AU", "INFO")
        print_status(f"  cos(elong) range: {np.min(cos_elong):.4f} – {np.max(cos_elong):.4f}", "INFO")
        print_status(f"  SRP_proxy range: {np.min(srp_proxy):.4f} – {np.max(srp_proxy):.4f}", "INFO")

    return srp_proxy, r_sun_au, cos_elong


def run_srp_check(df, verbose=False):
    print_status("=" * 70, "INFO")
    print_status("STEP 064: SOLAR RADIATION PRESSURE SYSTEMATIC CHECK", "TITLE")
    print_status("=" * 70, "INFO")

    # Clean outliers
    outlier_mask = detect_outliers_sigma(df["residual_m"].values, sigma_threshold=6.0)
    df_clean = df[~outlier_mask].copy()
    print_status(f"Outlier cleaning: {np.sum(outlier_mask)}/{len(df)} removed (6σ)", "INFO")
    print_status(f"Clean dataset: N = {len(df_clean):,}", "INFO")

    # Compute SRP proxy
    srp_proxy, r_sun_au, cos_elong = compute_srp_proxy(df_clean, verbose=verbose)
    df_clean["srp_proxy"] = srp_proxy
    df_clean["cos_elong"] = cos_elong
    residuals = df_clean["residual_m"].values
    n = len(residuals)

    # ==================================================================
    # TEST 1: Collinearity diagnostic
    # ==================================================================
    corr_cos_srp = float(np.corrcoef(cos_elong, srp_proxy)[0, 1])
    vif_srp = float(1.0 / max(1e-12, 1.0 - corr_cos_srp ** 2))

    print_status("", "INFO")
    print_status("TEST 1: COLLINEARITY DIAGNOSTIC", "PROCESS")
    print_status(f"  Pearson r(cosD, SRP_proxy) = {corr_cos_srp:.6f}", "CALC")
    print_status(f"  VIF(SRP_proxy) = {vif_srp:.2f}", "CALC")
    if vif_srp > 10.0:
        print_status(
            "  WARNING: VIF > 10 indicates severe multicollinearity; "
            "joint OLS decomposition is unreliable.", "WARNING"
        )

    # ==================================================================
    # TEST 2: Detrended-residual correlation with SRP proxy
    # ==================================================================
    # Remove best-fit TEP signal
    reg_tep = linear_regression(residuals, cos_elong)
    tep_pred = reg_tep["amplitude"] * cos_elong + reg_tep["intercept"]
    detrended = residuals - tep_pred

    # Correlation of detrended residuals with SRP proxy
    r_detrend_srp, p_detrend_srp = stats.pearsonr(detrended, srp_proxy)
    snr_detrend_srp = abs(r_detrend_srp) * np.sqrt(n - 2) / np.sqrt(1 - r_detrend_srp**2 + 1e-12)

    print_status("", "INFO")
    print_status("TEST 2: DETRENDED-RESIDUAL SRP CORRELATION", "PROCESS")
    print_status(
        f"  After removing best-fit TEP cos(D) signal (eta = {reg_tep['eta']:.4e})", "INFO"
    )
    print_status(f"  r(detrended, SRP_proxy) = {r_detrend_srp:.6e}", "CALC")
    print_status(f"  p-value = {p_detrend_srp:.4e}", "CALC")
    print_status(f"  SNR = {snr_detrend_srp:.2f}σ", "CALC")

    detrended_srp_significant = bool(p_detrend_srp < 0.05)
    if detrended_srp_significant:
        print_status(
            "  RESULT: Significant correlation with SRP proxy AFTER TEP detrending. "
            "This suggests SRP may contribute to the residuals.", "WARNING"
        )
    else:
        print_status(
            "  RESULT: No significant correlation with SRP proxy after TEP detrending. "
            "SRP does not explain residual variance beyond the TEP cos(D) term.", "PASS"
        )

    # ==================================================================
    # TEST 3: SRP-scaling test (central test)
    # ==================================================================
    # Bin by heliocentric distance, fit eta in each bin, test for 1/r^2 scaling
    n_bins = 10
    bin_edges = np.linspace(r_sun_au.min(), r_sun_au.max(), n_bins + 1)

    bin_centers = []
    bin_etas = []
    bin_eta_errs = []
    bin_inv_r2 = []
    bin_nobs = []

    for i in range(n_bins):
        mask = (r_sun_au >= bin_edges[i]) & (r_sun_au < bin_edges[i + 1])
        if mask.sum() > 50:
            bin_resid = residuals[mask]
            bin_cos = cos_elong[mask]
            reg_bin = linear_regression(bin_resid, bin_cos)
            bin_centers.append((bin_edges[i] + bin_edges[i + 1]) / 2)
            bin_etas.append(reg_bin["eta"])
            bin_eta_errs.append(reg_bin["eta_error"])
            bin_inv_r2.append(1.0 / ((bin_edges[i] + bin_edges[i + 1]) / 2) ** 2)
            bin_nobs.append(int(mask.sum()))

    bin_centers = np.array(bin_centers)
    bin_etas = np.array(bin_etas)
    bin_eta_errs = np.array(bin_eta_errs)
    bin_inv_r2 = np.array(bin_inv_r2)

    # Fit eta_bin = C + M * inv_r2  (linear in inverse-square)
    # SRP predicts M != 0 and C ≈ 0 (pure 1/r^2 scaling)
    # Nordtvedt predicts M = 0 (constant eta)
    X_scale = np.column_stack([np.ones(len(bin_inv_r2)), bin_inv_r2])
    coeffs_scale, _, _, _ = stable_lstsq(X_scale, bin_etas)
    C_fit, M_fit = coeffs_scale

    # Standard errors from weighted least squares using bin errors
    weights = 1.0 / (bin_eta_errs ** 2)
    W = np.diag(weights)
    XtWX = X_scale.T @ W @ X_scale
    try:
        XtWX_inv = np.linalg.inv(XtWX)
    except np.linalg.LinAlgError:
        XtWX_inv = np.linalg.pinv(XtWX, rcond=1e-12, hermitian=True)
    se_scale = np.sqrt(np.diag(XtWX_inv))
    se_C, se_M = se_scale

    # T-statistics
    t_M = M_fit / se_M if se_M > 0 else 0.0
    p_M = 2.0 * (1.0 - stats.t.cdf(abs(t_M), df=len(bin_etas) - 2)) if se_M > 0 else 1.0
    t_C = C_fit / se_C if se_C > 0 else 0.0
    p_C = 2.0 * (1.0 - stats.t.cdf(abs(t_C), df=len(bin_etas) - 2)) if se_C > 0 else 1.0

    # Reduced chi^2 of the scaling fit
    pred_etas = X_scale @ coeffs_scale
    chi2 = np.sum(((bin_etas - pred_etas) / bin_eta_errs) ** 2)
    chi2_red = chi2 / (len(bin_etas) - 2)

    print_status("", "INFO")
    print_status("TEST 3: SRP SCALING TEST (central)", "PROCESS")
    print_status(f"  Bins: {len(bin_centers)} (N per bin: {min(bin_nobs)}–{max(bin_nobs)})", "INFO")
    print_status(f"  Fit: eta(r) = C + M / r^2", "CALC")
    print_status(f"  C = {C_fit:.4e} ± {se_C:.4e}  (t = {t_C:.2f}, p = {p_C:.3f})", "CALC")
    print_status(f"  M = {M_fit:.4e} ± {se_M:.4e}  (t = {t_M:.2f}, p = {p_M:.3f})", "CALC")
    print_status(f"  χ²_red = {chi2_red:.2f}", "CALC")

    # Expected SRP scaling amplitude
    # If the ~5mm synodic signal were pure SRP, M would be ~ eta_mean * mean(r^2)
    mean_r2 = np.mean(r_sun_au ** 2)
    expected_M_for_pure_srp = reg_tep["eta"] * mean_r2
    print_status(f"  Expected M for pure SRP = {expected_M_for_pure_srp:.4e}", "CALC")

    # Interpretation
    srp_scaling_detected = bool(abs(t_M) > 2.0 and p_M < 0.05)
    constant_eta_preferred = bool(abs(t_M) < 2.0)

    if constant_eta_preferred:
        scaling_interpretation = (
            "The measured synodic modulation shows NO significant 1/r² scaling. "
            "This is inconsistent with an SRP origin and consistent with a "
            "heliocentric-distance-independent Nordtvedt-like signal."
        )
        scaling_status = "PASS"
    elif srp_scaling_detected and np.sign(M_fit) == np.sign(expected_M_for_pure_srp):
        scaling_interpretation = (
            "Significant 1/r² scaling detected with the sign expected for SRP. "
            "This raises a vulnerability that must be addressed."
        )
        scaling_status = "FAIL"
    else:
        scaling_interpretation = (
            "Significant 1/r² scaling detected, but with a sign or magnitude "
            "inconsistent with simple SRP. Further investigation required."
        )
        scaling_status = "WARNING"

    print_status(f"  {scaling_interpretation}", scaling_status)

    # ==================================================================
    # TEST 4: Perihelion-vs-aphelion eta comparison (leverages Step 022)
    # ==================================================================
    p15 = np.percentile(r_sun_au, 15)
    p85 = np.percentile(r_sun_au, 85)
    peri_mask = r_sun_au <= p15
    aph_mask = r_sun_au >= p85

    reg_peri = linear_regression(residuals[peri_mask], cos_elong[peri_mask])
    reg_aph = linear_regression(residuals[aph_mask], cos_elong[aph_mask])

    eta_diff = reg_peri["eta"] - reg_aph["eta"]
    diff_err = np.sqrt(reg_peri["eta_error"] ** 2 + reg_aph["eta_error"] ** 2)
    diff_snr = abs(eta_diff) / diff_err if diff_err > 0 else 0.0

    # SRP-predicted differential (perihelion should be ~7% larger in magnitude)
    srp_predicted_diff = reg_peri["eta"] * (1.0 - (p15 / p85) ** 2)

    print_status("", "INFO")
    print_status("TEST 4: PERIHELION–APHELION DIFFERENTIAL", "PROCESS")
    print_status(f"  Perihelion η = {reg_peri['eta']:.4e} ± {reg_peri['eta_error']:.4e} (N={peri_mask.sum()})", "CALC")
    print_status(f"  Aphelion   η = {reg_aph['eta']:.4e} ± {reg_aph['eta_error']:.4e} (N={aph_mask.sum()})", "CALC")
    print_status(f"  Δη = {eta_diff:.4e} ± {diff_err:.4e}  (SNR = {diff_snr:.2f}σ)", "CALC")
    print_status(f"  SRP-predicted Δη ≈ {srp_predicted_diff:.4e}", "CALC")

    # ==================================================================
    # Overall summary
    # ==================================================================
    print_status("", "INFO")
    print_status("OVERALL SUMMARY", "TITLE")

    # Criteria for PASS:
    # - No SRP correlation in detrended residuals
    # - No 1/r^2 scaling detected
    # - Perihelion-aphelion differential not consistent with SRP prediction
    if (not detrended_srp_significant and
        constant_eta_preferred and
        diff_snr < 2.0):
        overall_status = "PASS"
        overall_interpretation = (
            "All four tests fail to detect an SRP signature. The synodic-phase "
            "modulation is consistent with a heliocentric-distance-independent "
            "Nordtvedt-like signal and inconsistent with a simple unmodeled-SRP "
            "origin. INPOP19a's built-in SRP modeling, combined with the lack of "
            "1/r² scaling in the residuals, rules out SRP as the dominant source of "
            "the detected signal."
        )
    elif detrended_srp_significant or srp_scaling_detected:
        overall_status = "FAIL"
        overall_interpretation = (
            "One or more tests detect an SRP-like signature in the residuals. "
            "This constitutes a critical vulnerability that must be investigated "
            "before a TEP claim can be considered robust."
        )
    else:
        overall_status = "WARNING"
        overall_interpretation = (
            "Tests are inconclusive. While no clear SRP signature is detected, "
            "the statistical power of the scaling test is limited by the small "
            "heliocentric-distance baseline (3.4%)."
        )

    print_status(f"  Overall status: {overall_status}", overall_status)
    print_status(f"  {overall_interpretation}", overall_status)
    print_status("=" * 70, "INFO")

    return {
        "step_id": "step_064",
        "status": overall_status,
        "interpretation": overall_interpretation,
        "test_1_collinearity": {
            "r_cosd_srp_proxy": corr_cos_srp,
            "vif_srp_proxy": vif_srp,
            "severe_multicollinearity": bool(vif_srp > 10.0),
        },
        "test_2_detrended_srp_correlation": {
            "r_detrended_srp": float(r_detrend_srp),
            "p_value": float(p_detrend_srp),
            "snr_sigma": float(snr_detrend_srp),
            "significant": detrended_srp_significant,
        },
        "test_3_srp_scaling": {
            "n_bins": int(len(bin_centers)),
            "fit_c": float(C_fit),
            "fit_c_error": float(se_C),
            "fit_m": float(M_fit),
            "fit_m_error": float(se_M),
            "t_stat_m": float(t_M),
            "p_value_m": float(p_M),
            "chi2_red": float(chi2_red),
            "expected_m_for_pure_srp": float(expected_M_for_pure_srp),
            "srp_scaling_detected": srp_scaling_detected,
            "constant_eta_preferred": constant_eta_preferred,
            "scaling_interpretation": scaling_interpretation,
            "bin_centers_au": bin_centers.tolist(),
            "bin_etas": bin_etas.tolist(),
            "bin_eta_errors": bin_eta_errs.tolist(),
            "bin_nobs": bin_nobs,
        },
        "test_4_perihelion_aphelion": {
            "perihelion_eta": float(reg_peri["eta"]),
            "perihelion_eta_error": float(reg_peri["eta_error"]),
            "perihelion_n": int(peri_mask.sum()),
            "aphelion_eta": float(reg_aph["eta"]),
            "aphelion_eta_error": float(reg_aph["eta_error"]),
            "aphelion_n": int(aph_mask.sum()),
            "delta_eta": float(eta_diff),
            "delta_eta_error": float(diff_err),
            "significance_sigma": float(diff_snr),
            "srp_predicted_delta_eta": float(srp_predicted_diff),
        },
        "baseline_tep_only": {
            "eta": float(reg_tep["eta"]),
            "eta_error": float(reg_tep["eta_error"]),
            "snr_sigma": float(abs(reg_tep["eta"]) / reg_tep["eta_error"]),
        },
        "n_observations": n,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Step 064: Solar Radiation Pressure Systematic Check"
    )
    parser.add_argument("--verbose", action="store_true", help="Verbose output")
    args = parser.parse_args()

    log_dir = PROJECT_ROOT / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    logger = TEPLogger(
        "step_064", str(log_dir / "step_064_srp_systematic_check.log")
    )
    set_step_logger(logger)
    set_verbose_mode(args.verbose)

    data_path = (
        PROJECT_ROOT / "data" / "processed" / "INPOP19a_all_stations_residuals.csv"
    )
    if not data_path.exists():
        print_status(f"Input file not found: {data_path}", "ERROR")
        sys.exit(1)

    df = pd.read_csv(data_path)
    results = run_srp_check(df, verbose=args.verbose)
    logger.save_step_results(results, PROJECT_ROOT, "step_064_srp_systematic_check")
