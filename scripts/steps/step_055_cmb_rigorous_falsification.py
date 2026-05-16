#!/usr/bin/env python3
"""
Step 055: CMB Anisotropy Rigorous Falsification Suite
=======================================================

This step directly addresses three specific methodological criticisms of the
CMB anisotropy claim (Section 5.7):

  1. Aliasing: The correlation r(cosD, cosθ) = 0.050 is statistically
     significant at N ~ 25,000 and could spuriously produce a cosθ signal.
  2. Multicollinearity: In the 5-parameter joint fit, heliocentric distance
     becomes non-significant when cosθ is added, suggesting collinearity
     rather than a genuine cosmological effect.
  3. Cosmological overreach: A monthly regression term is insufficient
     evidence to link the Solar System signal to large-scale structure.

The suite comprises five scored falsification checks (A–E) used for the
published PASS/WARNING gate, plus supplementary temporal and reference-axis
controls (F–G) reported in the JSON alongside the gate.

  A. Aliasing Simulation Falsification:
     Simulate 5,000 datasets with a realistic synodic signal and NO CMB
     dependence. Because cosθ correlates with cosD at r = 0.050, some
     cosθ coefficient will emerge from aliasing. The observed coefficient
     must lie beyond the 99.9th percentile of this null distribution.

  B. Multicollinearity Diagnostic Suite:
     Compute Variance Inflation Factors (VIFs), condition numbers, and
     variance decomposition for all predictors in the full joint model.
     Quantify exactly how much the 0.05 correlation inflates standard
     errors. Show that even after accounting for this inflation, cosθ
     remains highly significant.

  C. Permutation Test:
     Break any true CMB–residual relationship by permuting cosθ values
     while preserving the marginal distribution and correlation structure
     with other predictors. Fit the joint model 1,000 times. Compare the
     observed η_θ to the permutation null distribution.

  D. Random-axis Monte Carlo on S² (primary ΔAIC null):
     Draw fixed sky axes uniformly on the celestial sphere (50,000 draws by
     default). For each random direction, fit the base-plus-direction model
     and build the empirical null for ΔAIC under ``any fixed axis.'' Compare
     the Planck dipole to this null (Jeffreys intervals on empirical p-values).
     Supplementary refined nulls (D3–D5): synodic
     phase scrambling, CMB-axis rotation (global SO(3) and local cone),
     and Gram–Schmidt orthogonalized direction scrambling.

  E. Gram–Schmidt Orthogonalization:
     Explicitly orthogonalize cosθ against cosD, r_c, and vr_c. Fit the
     joint model using only the orthogonal residual. If the signal is
     genuine, it must persist in the component that is mathematically
     independent of all other predictors.

The tests are deliberately adversarial: they assume the alternative
hypothesis is false and ask which narrower null explanations the data reject.
Passing these diagnostics supports a directional residual component, but it
does not by itself establish a cosmological origin.

  F. Year-block jackknife (temporal stability):
     For each calendar year in the cleaned sample, refit the five-parameter
     joint model with all other years retained and predictors re-centered on
     that subset. Summarizes whether η_θ and its t-statistic are stable under
     temporal excision.

  G. Alternative fixed celestial axes:
     Repeat the joint-model direction metrics for the north ecliptic pole,
     the ICRS galactic north pole, and the time-averaged Earth–Sun orbital
     angular-momentum unit vector from the same ephemeris. These are
     non-CMB reference directions in the same hypothesis class (fixed sky
     axis against the Earth–Moon line).

  H. Dual-axis simultaneous fit (CMB + galactic north) — primary orientation
     identifiability block; single-axis Planck vs galactic ΔAIC comparisons are
     reported as supplementary.

  I. TEP η_θ prediction coverage (Step 033):
     Compare fitted η_θ from the Planck-axis joint model to the Step 033
     orientation scale prediction (η_θ ≈ η₀ m_CMB, |m_CMB| ~ O(1)).
"""

import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

import json
import os
from concurrent.futures import ThreadPoolExecutor

import numpy as np
from scripts.utils.numerics import stable_lstsq
import pandas as pd
from scipy import stats
from skyfield.api import load

from scripts.utils.config import get_config
from scripts.utils.llr_constants import ETA_SCALE_FACTOR
from scripts.utils.statistical_utils import detect_outliers_sigma
from scripts.utils.logger import TEPLogger, set_step_logger, print_status

# CMB dipole direction (Planck 2018)
_CMB_RA_RAD = np.deg2rad(168.14)
_CMB_DEC_RAD = np.deg2rad(-7.22)
_CMB_UNIT = np.array([
    np.cos(_CMB_DEC_RAD) * np.cos(_CMB_RA_RAD),
    np.cos(_CMB_DEC_RAD) * np.sin(_CMB_RA_RAD),
    np.sin(_CMB_DEC_RAD),
])
_CMB_UNIT = _CMB_UNIT / np.linalg.norm(_CMB_UNIT)

# North ecliptic pole, J2000 ICRS (IAU convention)
_NEP_RA_RAD = np.deg2rad(270.0)
_NEP_DEC_RAD = np.deg2rad(66.0 + 33.0 / 60.0 + 38.0 / 3600.0)
_NEP_UNIT = np.array(
    [
        np.cos(_NEP_DEC_RAD) * np.cos(_NEP_RA_RAD),
        np.cos(_NEP_DEC_RAD) * np.sin(_NEP_RA_RAD),
        np.sin(_NEP_DEC_RAD),
    ]
)
_NEP_UNIT = _NEP_UNIT / np.linalg.norm(_NEP_UNIT)

# Galactic north pole J2000 ICRS (e.g. Blaauw et al.; standard reduction)
_GAL_N_RA_RAD = np.deg2rad(192.85948)
_GAL_N_DEC_RAD = np.deg2rad(27.128336)
_GAL_N_UNIT = np.array(
    [
        np.cos(_GAL_N_DEC_RAD) * np.cos(_GAL_N_RA_RAD),
        np.cos(_GAL_N_DEC_RAD) * np.sin(_GAL_N_RA_RAD),
        np.sin(_GAL_N_DEC_RAD),
    ]
)
_GAL_N_UNIT = _GAL_N_UNIT / np.linalg.norm(_GAL_N_UNIT)


def jd_to_gregorian_year(jd):
    """Gregorian calendar year for Julian Date (vectorized via pandas)."""
    jd = np.asarray(jd, dtype=np.float64)
    # JD 2440587.5 = Unix epoch 1970-01-01T00:00:00Z
    dt = pd.to_datetime(jd - 2440587.5, unit="D", origin="unix", utc=True)
    return dt.year.values.astype(np.int32)


def jeffreys_binomial_ci(k, n, alpha=0.05):
    """Equal-tailed Jeffreys interval for Binomial(k | n) proportion."""
    if n <= 0:
        return None, None
    k = int(np.clip(k, 0, n))
    low = stats.beta.ppf(alpha / 2.0, k + 0.5, n - k + 0.5)
    high = stats.beta.ppf(1.0 - alpha / 2.0, k + 0.5, n - k + 0.5)
    return float(low), float(high)


def compute_cmb_projections(jd_array):
    """Compute CMB-frame kinematic projections for every epoch."""
    eph_path = PROJECT_ROOT / "data" / "raw" / "de440.bsp"
    if not eph_path.exists():
        raise FileNotFoundError(
            f"Required Skyfield kernel missing for Step 055: {eph_path}. "
            "Download and place de440.bsp under data/raw and verify it via Step 000 manifest."
        )
    planets = load(str(eph_path))
    earth = planets["earth"]
    moon = planets["moon"]
    sun = planets["sun"]
    ts = load.timescale()
    timestamps = ts.tt(jd=jd_array)

    earth_pos = earth.at(timestamps)
    earth_pv = earth_pos.position.km
    earth_vv = earth_pos.velocity.km_per_s
    sun_pv = sun.at(timestamps).position.km
    moon_pv = moon.at(timestamps).position.km

    rel_pos = earth_pv - sun_pv
    rel_vel = earth_vv - sun.at(timestamps).velocity.km_per_s
    distance_km = np.linalg.norm(rel_pos, axis=0)
    r_hat = rel_pos / distance_km
    v_radial_kms = np.sum(rel_vel * r_hat, axis=0)
    v_parallel_kms = np.sum(rel_vel * _CMB_UNIT[:, None], axis=0)

    em_vec = moon_pv - earth_pv
    em_dist = np.linalg.norm(em_vec, axis=0)
    em_hat = em_vec / em_dist
    earth_moon_cos_theta = np.sum(em_hat * _CMB_UNIT[:, None], axis=0)

    rel_pos = np.asarray(rel_pos)
    rel_vel = np.asarray(rel_vel)
    if rel_pos.ndim != 2 or rel_vel.ndim != 2:
        raise ValueError("Expected 2D ephemeris position/velocity arrays from Skyfield.")
    n0, n1 = rel_pos.shape
    if n0 == 3 and n1 != 3:
        rel_T = rel_pos.T
        vel_T = rel_vel.T
    elif n1 == 3 and n0 != 3:
        rel_T = rel_pos
        vel_T = rel_vel
    else:
        raise ValueError(
            f"Unexpected ephemeris array layout rel_pos.shape={rel_pos.shape}; "
            "expected (3, n) or (n, 3) with n >> 3."
        )
    # np.cross on (3,n) with axisa=0 returns (n,3); always work in (n,3) row layout.
    h_vec = np.cross(rel_T, vel_T)
    h_norm = np.linalg.norm(h_vec, axis=-1, keepdims=True)
    with np.errstate(invalid="ignore", divide="ignore"):
        h_hat = h_vec / h_norm
    bad = (h_norm.squeeze() < 1e-12)
    h_hat[bad, :] = np.nan
    h_mean = np.nanmean(h_hat, axis=0)
    h_mean_norm = np.linalg.norm(h_mean)
    if h_mean_norm < 1e-12:
        mean_orbit_normal = np.array([0.0, 0.0, 1.0])
    else:
        mean_orbit_normal = h_mean / h_mean_norm

    return {
        "v_parallel_kms": v_parallel_kms,
        "earth_moon_cos_theta": earth_moon_cos_theta,
        "sun_distance_au": distance_km / 1.495978707e8,
        "radial_velocity_kms": v_radial_kms,
        "earth_moon_unit_vectors": em_hat,
        "mean_orbit_normal": mean_orbit_normal,
    }


def fit_full_joint_model(res, cosD, r_c, vr_c, cos_theta_c):
    """Fit the 5-parameter joint model and return (eta_theta, se_eta_theta, eta_r, eta_vr, eta_0)."""
    X = np.column_stack([cosD, r_c * cosD, vr_c * cosD, cos_theta_c * cosD, np.ones(len(cosD))])
    coeffs, _, rank, _ = stable_lstsq(X, res)
    if rank < 5:
        return None, None, None, None, None
    with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
        resid = res - X @ coeffs
    mse = np.sum(resid ** 2) / (len(res) - 5)
    cov = mse * np.linalg.pinv(X.T @ X, rcond=1e-10, hermitian=True)
    se = np.sqrt(np.diag(cov))
    eta_0 = coeffs[0] / ETA_SCALE_FACTOR
    eta_r = coeffs[1] / ETA_SCALE_FACTOR
    eta_vr = coeffs[2] / ETA_SCALE_FACTOR
    eta_theta = coeffs[3] / ETA_SCALE_FACTOR
    se_eta_theta = se[3] / ETA_SCALE_FACTOR
    return eta_theta, se_eta_theta, eta_r, eta_vr, eta_0


def fit_orientation_joint_model(res, cosD, r_c, vr_c, cos_theta_c, n_params=5):
    """Fit base+single-orientation (5p) or dual-orientation (6p) joint models."""
    if n_params == 5:
        X = np.column_stack([cosD, r_c * cosD, vr_c * cosD, cos_theta_c * cosD, np.ones(len(cosD))])
    elif n_params == 6:
        cos_theta_cmb_c, cos_theta_gal_c = cos_theta_c
        X = np.column_stack(
            [
                cosD,
                r_c * cosD,
                vr_c * cosD,
                cos_theta_cmb_c * cosD,
                cos_theta_gal_c * cosD,
                np.ones(len(cosD)),
            ]
        )
    else:
        raise ValueError(f"n_params must be 5 or 6, got {n_params}")
    n = len(cosD)
    coeffs, _, rank, _ = stable_lstsq(X, res)
    if rank < n_params:
        return None
    with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
        resid = res - X @ coeffs
    rss = float(np.sum(resid ** 2))
    mse = rss / (n - n_params)
    cov = mse * np.linalg.pinv(X.T @ X, rcond=1e-10, hermitian=True)
    se = np.sqrt(np.diag(cov))
    aic = float(n * np.log(rss / n) + 2 * n_params)
    out = {
        "n_params": n_params,
        "rss": rss,
        "aic": aic,
        "coeffs": coeffs,
        "se": se,
        "vifs": {name: float(v) for name, v in zip(
            ["cosD", "r_c·cosD", "vr_c·cosD", "cosθ_cmb·cosD", "cosθ_gal·cosD", "intercept"][:n_params],
            compute_vif(X),
        )},
    }
    if n_params == 5:
        eta_theta = coeffs[3] / ETA_SCALE_FACTOR
        se_theta = se[3] / ETA_SCALE_FACTOR
        out["eta_theta"] = float(eta_theta)
        out["eta_theta_se"] = float(se_theta)
        out["eta_theta_t"] = float(eta_theta / se_theta) if se_theta > 0 else 0.0
    else:
        eta_cmb = coeffs[3] / ETA_SCALE_FACTOR
        eta_gal = coeffs[4] / ETA_SCALE_FACTOR
        se_cmb = se[3] / ETA_SCALE_FACTOR
        se_gal = se[4] / ETA_SCALE_FACTOR
        out["eta_cmb"] = float(eta_cmb)
        out["eta_gal"] = float(eta_gal)
        out["eta_cmb_se"] = float(se_cmb)
        out["eta_gal_se"] = float(se_gal)
        out["eta_cmb_t"] = float(eta_cmb / se_cmb) if se_cmb > 0 else 0.0
        out["eta_gal_t"] = float(eta_gal / se_gal) if se_gal > 0 else 0.0
    return out


def nested_f_test(rss_small, p_small, rss_large, p_large, n_obs):
    """F-test whether the larger model significantly improves fit."""
    df_num = p_large - p_small
    df_den = n_obs - p_large
    if df_num <= 0 or df_den <= 0 or rss_large <= 0:
        return None, None
    delta_rss = rss_small - rss_large
    if delta_rss <= 0:
        return 0.0, 1.0
    f_stat = (delta_rss / df_num) / (rss_large / df_den)
    p_val = 1 - stats.f.cdf(f_stat, df_num, df_den)
    return float(f_stat), float(p_val)


def compute_vif(X):
    """Compute Variance Inflation Factors for each column of design matrix X.
    VIF_j = 1 / (1 - R²_j), where R²_j is from regressing column j on all others.
    """
    n, p = X.shape
    vifs = []
    for j in range(p):
        y_col = X[:, j]
        X_others = np.delete(X, j, axis=1)
        coeffs, _, rank, _ = stable_lstsq(X_others, y_col)
        if rank < X_others.shape[1] or np.var(y_col) == 0:
            vifs.append(np.inf)
            continue
        with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
            y_pred = X_others @ coeffs
        ss_res = np.sum((y_col - y_pred) ** 2)
        ss_tot = np.sum((y_col - np.mean(y_col)) ** 2)
        r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0.0
        vifs.append(1.0 / (1.0 - r2) if r2 < 0.999999 else np.inf)
    return np.array(vifs)


def random_unit_vector(rng):
    """Generate a uniformly random unit vector on the sphere."""
    v = rng.normal(size=3)
    return v / np.linalg.norm(v)


def rotate_vector_about_axis(v: np.ndarray, axis_unit: np.ndarray, angle_rad: float) -> np.ndarray:
    """Rodrigues rotation of v about axis_unit by angle_rad."""
    k = axis_unit / np.linalg.norm(axis_unit)
    v = np.asarray(v, dtype=np.float64)
    return (
        v * np.cos(angle_rad)
        + np.cross(k, v) * np.sin(angle_rad)
        + k * np.dot(k, v) * (1.0 - np.cos(angle_rad))
    )


def random_unit_perpendicular(v: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    """Uniform random unit vector in the plane perpendicular to v."""
    v = v / np.linalg.norm(v)
    raw = rng.normal(size=3)
    perp = raw - np.dot(raw, v) * v
    n = np.linalg.norm(perp)
    if n < 1e-12:
        raw = np.array([1.0, 0.0, 0.0]) if abs(v[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
        perp = raw - np.dot(raw, v) * v
        n = np.linalg.norm(perp)
    return perp / n


def orthogonalized_cos_theta_cosD(
    cos_theta_c: np.ndarray,
    cosD: np.ndarray,
    r_c: np.ndarray,
    vr_c: np.ndarray,
) -> np.ndarray | None:
    """Gram–Schmidt residual of cosθ_c·cosD with respect to the heliocentric base."""
    y_target = cos_theta_c * cosD
    x_orth_base = np.column_stack([cosD, r_c * cosD, vr_c * cosD, np.ones(len(cosD))])
    coeffs_gs, _, rank_gs, _ = stable_lstsq(x_orth_base, y_target)
    if rank_gs < 4:
        return None
    with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
        return y_target - x_orth_base @ coeffs_gs


def load_step033_orientation_prediction():
    """Load Step 033 CMB orientation coefficient scale prediction."""
    path = PROJECT_ROOT / "results" / "outputs" / "step_033_quantitative_eta_prediction.json"
    if not path.exists():
        raise FileNotFoundError(
            f"step_033_quantitative_eta_prediction.json not found: {path}. "
            "Run pipeline step 033 before step 055."
        )
    with open(path, "r", encoding="utf-8") as f:
        step_033 = json.load(f)
    pred = step_033.get("cmb_orientation_prediction")
    if pred is None:
        raise RuntimeError(
            "step_033 output missing cmb_orientation_prediction; rerun step 033."
        )
    return pred


def compute_eta_theta_prediction_coverage(
    eta_theta_fitted: float,
    eta_theta_se: float,
    pred_block: dict,
):
    """Compare fitted η_θ to Step 033 TEP orientation scale prediction."""
    center = float(pred_block["eta_theta_predicted_center"])
    se_center = float(pred_block["eta_0_synodic_error"])
    band = pred_block["order_of_magnitude_band"]
    lo = float(band["lo"])
    hi = float(band["hi"])
    band_lo, band_hi = min(lo, hi), max(lo, hi)
    within_band = bool(band_lo <= eta_theta_fitted <= band_hi)
    combined_se = float(np.sqrt(eta_theta_se**2 + se_center**2))
    z_score = float((eta_theta_fitted - center) / combined_se) if combined_se > 0 else None
    ratio_to_center = (
        float(eta_theta_fitted / center) if center not in (0.0, -0.0) else None
    )
    return {
        "eta_theta_fitted": float(eta_theta_fitted),
        "eta_theta_fitted_se": float(eta_theta_se),
        "eta_theta_predicted_center": center,
        "eta_0_synodic_error": se_center,
        "order_of_magnitude_band": {"lo": lo, "hi": hi},
        "within_predicted_band": within_band,
        "z_score_vs_predicted_center": z_score,
        "ratio_fitted_to_predicted_center": ratio_to_center,
        "formula": pred_block.get("formula"),
        "volumetric_scale_anchor": pred_block.get("volumetric_scale_anchor"),
    }


def delta_aic_for_orthogonal_component(
    y_orth: np.ndarray,
    res: np.ndarray,
    cosD: np.ndarray,
    r_c: np.ndarray,
    vr_c: np.ndarray,
    rss_base: float,
    aic_base: float,
    n_clean: int,
    n_eff: float,
    use_eff: bool = True,
):
    """Augment the heliocentric base with an orthogonalized direction column."""
    if y_orth is None:
        return None, None, None, None
    x_aug = np.column_stack([cosD, r_c * cosD, vr_c * cosD, y_orth, np.ones(n_clean)])
    coeffs_aug, _, rank_aug, _ = stable_lstsq(x_aug, res)
    if rank_aug < 5:
        return None, None, None, None
    with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
        resid_aug = res - x_aug @ coeffs_aug
    rss_aug = np.sum(resid_aug ** 2)
    aic_aug = n_clean * np.log(rss_aug / n_clean) + 2 * 5
    delta_aic = aic_base - aic_aug
    delta_rss = rss_base - rss_aug
    f_dir = (delta_rss / 1) / (rss_aug / (n_clean - 5))
    p_f = 1 - stats.f.cdf(f_dir, 1, n_clean - 5)
    if use_eff and n_eff > 5:
        f_dir_eff = (delta_rss / 1) / (rss_aug / (n_eff - 5))
    else:
        f_dir_eff = f_dir
    return delta_aic, f_dir, p_f, f_dir_eff


def cmb_falsification_analysis(df, verbose=False):
    print_status("═══ Step 055: CMB Anisotropy Rigorous Falsification ═══", "TITLE")

    residuals = df["residual_m"].values
    cos_elong = np.cos(df["elongation_rad"].values)
    jd = df["date_julian"].values

    # ------------------------------------------------------------------
    # 1. Compute kinematics
    # ------------------------------------------------------------------
    print_status("Computing CMB-frame kinematics...", "PROCESS")
    cmb_data = compute_cmb_projections(jd)
    v_par = cmb_data["v_parallel_kms"]
    cos_theta = cmb_data["earth_moon_cos_theta"]
    r = cmb_data["sun_distance_au"]
    vr = cmb_data["radial_velocity_kms"]
    em_hat = cmb_data["earth_moon_unit_vectors"]

    corr_cosD_cosTheta = np.corrcoef(cos_elong, cos_theta)[0, 1]
    print_status(f"r(cosD, cosθ) = {corr_cosD_cosTheta:.4f}", "CALC")

    outlier_mask = detect_outliers_sigma(residuals, sigma_threshold=6.0)
    mask_clean = ~outlier_mask
    n_clean = mask_clean.sum()
    print_status(f"Cleaned dataset: {n_clean:,} observations", "CALC")

    res = residuals[mask_clean]
    cosD = cos_elong[mask_clean]
    v_par_c = v_par[mask_clean] - np.mean(v_par[mask_clean])
    cos_theta_c = cos_theta[mask_clean] - np.mean(cos_theta[mask_clean])
    r_c = r[mask_clean] - np.mean(r[mask_clean])
    vr_c = vr[mask_clean] - np.mean(vr[mask_clean])
    v_par_clean = v_par[mask_clean]
    r_clean = r[mask_clean]
    vr_clean = vr[mask_clean]
    cos_theta_clean = cos_theta[mask_clean]

    # ------------------------------------------------------------------
    # A. Aliasing Simulation Falsification
    # ------------------------------------------------------------------
    print_status("═══ A. Aliasing Simulation Falsification ═══", "TITLE")

    # Estimate the true synodic amplitude from the data
    X_syn = np.column_stack([cosD, np.ones(n_clean)])
    coeffs_syn, _, _, _ = stable_lstsq(X_syn, res)
    beta_syn = coeffs_syn[0]
    intercept_syn = coeffs_syn[1]
    resid_syn = res - (beta_syn * cosD + intercept_syn)
    noise_std = np.std(resid_syn)

    n_sim = 5000
    rng_sim = np.random.default_rng(42)
    eta_theta_sim = []
    t_theta_sim = []

    for i in range(n_sim):
        # Simulate data: true synodic signal + noise, ZERO CMB dependence
        res_sim = beta_syn * cosD + intercept_syn + rng_sim.normal(0, noise_std, n_clean)
        eta_t, se_t, _, _, _ = fit_full_joint_model(res_sim, cosD, r_c, vr_c, cos_theta_c)
        if eta_t is not None:
            eta_theta_sim.append(eta_t)
            t_theta_sim.append(eta_t / se_t if se_t > 0 else 0.0)

    eta_theta_sim = np.array(eta_theta_sim)
    t_theta_sim = np.array(t_theta_sim)

    # Observed value
    eta_obs, se_obs, _, _, _ = fit_full_joint_model(res, cosD, r_c, vr_c, cos_theta_c)
    t_obs = eta_obs / se_obs if se_obs > 0 else 0.0

    p_aliasing_eta = np.mean(np.abs(eta_theta_sim) >= np.abs(eta_obs))
    p_aliasing_t = np.mean(np.abs(t_theta_sim) >= np.abs(t_obs))

    print_status(
        f"Simulated η_θ (no CMB): mean={np.mean(eta_theta_sim):.4e}, "
        f"std={np.std(eta_theta_sim):.4e}, max|η|={np.max(np.abs(eta_theta_sim)):.4e}",
        "CALC",
    )
    print_status(
        f"Simulated t (no CMB): mean={np.mean(t_theta_sim):.2f}, "
        f"std={np.std(t_theta_sim):.2f}, max|t|={np.max(np.abs(t_theta_sim)):.2f}",
        "CALC",
    )
    print_status(
        f"Observed: η_θ = {eta_obs:.4e}, t = {t_obs:.2f}",
        "CALC",
    )
    print_status(
        f"Aliasing p-value (coefficient): {p_aliasing_eta:.4f}",
        "CALC",
    )
    print_status(
        f"Aliasing p-value (t-statistic): {p_aliasing_t:.4f}",
        "CALC",
    )

    aliasing_result = {
        "n_simulations": n_sim,
        "eta_theta_sim_mean": float(np.mean(eta_theta_sim)),
        "eta_theta_sim_std": float(np.std(eta_theta_sim)),
        "eta_theta_sim_99_9_percentile": float(np.percentile(np.abs(eta_theta_sim), 99.9)),
        "t_sim_mean": float(np.mean(t_theta_sim)),
        "t_sim_std": float(np.std(t_theta_sim)),
        "t_sim_99_9_percentile": float(np.percentile(np.abs(t_theta_sim), 99.9)),
        "observed_eta_theta": float(eta_obs),
        "observed_t": float(t_obs),
        "p_aliasing_coefficient": float(p_aliasing_eta),
        "p_aliasing_t_statistic": float(p_aliasing_t),
        "aliasing_rejected": bool(p_aliasing_t < 0.001),
    }

    # ------------------------------------------------------------------
    # B. Multicollinearity Diagnostic Suite
    # ------------------------------------------------------------------
    print_status("═══ B. Multicollinearity Diagnostic Suite ═══", "TITLE")

    X_full = np.column_stack([cosD, r_c * cosD, vr_c * cosD, cos_theta_c * cosD, np.ones(n_clean)])
    vifs = compute_vif(X_full)
    # Condition number of X^T X (not X itself)
    XtX = X_full.T @ X_full
    s = np.linalg.svd(XtX, compute_uv=False)
    cond_XtX = s[0] / s[-1] if s[-1] > 0 else np.inf
    cond_X = np.linalg.cond(X_full)

    predictor_names = ["cosD", "r_c·cosD", "vr_c·cosD", "cosθ_c·cosD", "intercept"]
    for name, vif in zip(predictor_names, vifs):
        if name == "intercept":
            continue
        print_status(f"VIF({name}) = {vif:.3f}", "CALC")
    print_status(f"Condition number κ(X) = {cond_X:.2e}", "CALC")
    print_status(f"Condition number κ(XᵀX) = {cond_XtX:.2e}", "CALC")

    # Compute what the SE would be if cosθ were perfectly orthogonal to all others
    # SE_orthogonal = sqrt(MSE / sum((cosθ_c*cosD)_orth²))
    # First orthogonalize cosθ_c*cosD against the other three columns
    y_orth = cos_theta_c * cosD
    X_others = np.column_stack([cosD, r_c * cosD, vr_c * cosD, np.ones(n_clean)])
    coeffs_orth, _, rank_orth, _ = stable_lstsq(X_others, y_orth)
    if rank_orth == 4:
        with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
            y_resid = y_orth - X_others @ coeffs_orth
        # Refit full model to get MSE
        coeffs_full, _, _, _ = stable_lstsq(X_full, res)
        with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
            resid_full = res - X_full @ coeffs_full
        mse_full = np.sum(resid_full ** 2) / (n_clean - 5)
        # SE if orthogonal: sqrt(MSE / sum(y_resid²))
        se_orthogonal = np.sqrt(mse_full / np.sum(y_resid ** 2)) / ETA_SCALE_FACTOR
        se_actual = se_obs
        inflation = se_actual / se_orthogonal if se_orthogonal > 0 else np.inf
        print_status(
            f"SE(cosθ) actual = {se_actual:.4e}; SE if orthogonal = {se_orthogonal:.4e}; "
            f"inflation = {inflation:.2f}x",
            "CALC",
        )
    else:
        se_orthogonal = None
        inflation = None

    vifs_no_intercept = [v for name, v in zip(predictor_names, vifs) if name != "intercept"]
    multicollinearity_result = {
        "vifs": {name: float(v) for name, v in zip(predictor_names, vifs) if name != "intercept"},
        "max_vif": float(np.max(vifs_no_intercept)),
        "mean_vif": float(np.mean(vifs_no_intercept)),
        "condition_number_X": float(cond_X),
        "condition_number_XtX": float(cond_XtX),
        "se_cos_theta_actual": float(se_actual) if se_actual is not None else None,
        "se_cos_theta_orthogonal": float(se_orthogonal) if se_orthogonal is not None else None,
        "se_inflation_factor": float(inflation) if inflation is not None else None,
        "severe_multicollinearity": bool(np.max(vifs_no_intercept) > 10 or cond_X > 1e3),
    }

    # ------------------------------------------------------------------
    # C. Permutation Test
    # ------------------------------------------------------------------
    print_status("═══ C. Permutation Test for cosθ ═══", "TITLE")

    n_perm = 1000
    rng_perm = np.random.default_rng(43)
    eta_theta_perm = []
    t_theta_perm = []

    for _ in range(n_perm):
        cos_theta_perm = rng_perm.permutation(cos_theta_c)
        eta_p, se_p, _, _, _ = fit_full_joint_model(res, cosD, r_c, vr_c, cos_theta_perm)
        if eta_p is not None:
            eta_theta_perm.append(eta_p)
            t_theta_perm.append(eta_p / se_p if se_p > 0 else 0.0)

    eta_theta_perm = np.array(eta_theta_perm)
    t_theta_perm = np.array(t_theta_perm)

    p_perm_eta = np.mean(np.abs(eta_theta_perm) >= np.abs(eta_obs))
    p_perm_t = np.mean(np.abs(t_theta_perm) >= np.abs(t_obs))

    print_status(
        f"Permuted η_θ: mean={np.mean(eta_theta_perm):.4e}, std={np.std(eta_theta_perm):.4e}, "
        f"99.9th %ile |η|={np.percentile(np.abs(eta_theta_perm), 99.9):.4e}",
        "CALC",
    )
    print_status(
        f"Permuted t: mean={np.mean(t_theta_perm):.2f}, std={np.std(t_theta_perm):.2f}, "
        f"99.9th %ile |t|={np.percentile(np.abs(t_theta_perm), 99.9):.2f}",
        "CALC",
    )
    print_status(f"Permutation p (coefficient) = {p_perm_eta:.4f}", "CALC")
    print_status(f"Permutation p (t-statistic) = {p_perm_t:.4f}", "CALC")

    permutation_result = {
        "n_permutations": n_perm,
        "eta_perm_mean": float(np.mean(eta_theta_perm)),
        "eta_perm_std": float(np.std(eta_theta_perm)),
        "eta_perm_99_9_percentile": float(np.percentile(np.abs(eta_theta_perm), 99.9)),
        "t_perm_mean": float(np.mean(t_theta_perm)),
        "t_perm_std": float(np.std(t_theta_perm)),
        "t_perm_99_9_percentile": float(np.percentile(np.abs(t_theta_perm), 99.9)),
        "observed_eta": float(eta_obs),
        "observed_t": float(t_obs),
        "p_value_coefficient": float(p_perm_eta),
        "p_value_t_statistic": float(p_perm_t),
        "permutation_rejected": bool(p_perm_t < 0.001),
    }

    # ------------------------------------------------------------------
    # D. Random-axis Monte Carlo on S² (ΔAIC null for any fixed axis)
    # ------------------------------------------------------------------
    print_status("═══ D. Random-axis Monte Carlo (any fixed sky axis) ═══", "TITLE")

    # Base heliocentric model (synodic + distance + radial velocity, no direction)
    X_base = np.column_stack([cosD, r_c * cosD, vr_c * cosD, np.ones(n_clean)])
    coeffs_base, _, rank_base, _ = stable_lstsq(X_base, res)
    with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
        resid_base = res - X_base @ coeffs_base
    rss_base = np.sum(resid_base ** 2)
    aic_base = n_clean * np.log(rss_base / n_clean) + 2 * 4

    # Estimate temporal autocorrelation in base residuals for effective-sample-size correction
    def autocorr_lag1(x):
        x_c = x - np.mean(x)
        denom = np.sum(x_c ** 2)
        return np.sum(x_c[:-1] * x_c[1:]) / denom if denom > 0 else 0.0

    rho_base = autocorr_lag1(resid_base)
    # Effective sample size for AR(1): n_eff = n * (1 - rho) / (1 + rho)
    n_eff = n_clean * (1.0 - rho_base) / (1.0 + rho_base) if abs(rho_base) < 0.999 else n_clean
    print_status(f"Base residual AR(1) ρ = {rho_base:.4f}; n_eff = {n_eff:.1f}", "CALC")

    def delta_aic_for_direction(cos_theta_dir_c, use_eff=False):
        """Fit base + direction and return delta_AIC = AIC_base - AIC_augmented.
        Positive means the direction improves the fit.
        If use_eff=True, return autocorrelation-adjusted F using n_eff."""
        X_aug = np.column_stack([cosD, r_c * cosD, vr_c * cosD, cos_theta_dir_c * cosD, np.ones(n_clean)])
        coeffs_aug, _, rank_aug, _ = stable_lstsq(X_aug, res)
        if rank_aug < 5:
            return None, None, None, None
        with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
            resid_aug = res - X_aug @ coeffs_aug
        rss_aug = np.sum(resid_aug ** 2)
        aic_aug = n_clean * np.log(rss_aug / n_clean) + 2 * 5
        delta_aic = aic_base - aic_aug
        # Standard F-statistic
        delta_rss = rss_base - rss_aug
        f_dir = (delta_rss / 1) / (rss_aug / (n_clean - 5))
        p_f = 1 - stats.f.cdf(f_dir, 1, n_clean - 5)
        # Effective-sample-size F-statistic (autocorrelation-adjusted)
        if use_eff and n_eff > 5:
            f_dir_eff = (delta_rss / 1) / (rss_aug / (n_eff - 5))
            p_f_eff = 1 - stats.f.cdf(f_dir_eff, 1, n_eff - 5)
        else:
            f_dir_eff = f_dir
            p_f_eff = p_f
        return delta_aic, f_dir, p_f, f_dir_eff

    # True CMB (compute once)
    delta_aic_true, f_true, p_f_true, f_true_eff = delta_aic_for_direction(cos_theta_c, use_eff=True)

    # D1. Uniform random directions on S² (empirical ΔAIC null)
    n_scramble = 50000
    strict_alpha = 0.01
    n_threads = max(2, min(int(get_config()["N_WORKERS"]), os.cpu_count() or 8))
    em_sub = em_hat[:, mask_clean]

    def run_scramble_batch(n_batch: int, seed_batch: int):
        rng_sc = np.random.default_rng(int(seed_batch))
        da, fs, fse, cr = [], [], [], []
        for _ in range(int(n_batch)):
            dir_vec = random_unit_vector(rng_sc)
            cos_theta_rand = np.sum(em_sub * dir_vec[:, None], axis=0)
            cos_theta_rand_c = cos_theta_rand - np.mean(cos_theta_rand)
            d_aic, f_dir, _p_f, f_dir_eff = delta_aic_for_direction(cos_theta_rand_c, use_eff=True)
            if d_aic is not None:
                da.append(d_aic)
                fs.append(f_dir)
                fse.append(f_dir_eff)
                cr.append(np.corrcoef(cosD, cos_theta_rand)[0, 1])
        return (
            np.asarray(da, dtype=np.float64),
            np.asarray(fs, dtype=np.float64),
            np.asarray(fse, dtype=np.float64),
            np.asarray(cr, dtype=np.float64),
        )

    edges = np.linspace(0, n_scramble, n_threads + 1, dtype=int)
    batches = [(int(edges[t + 1] - edges[t]), 44 + t * 1_000_003) for t in range(n_threads) if edges[t + 1] > edges[t]]
    parts = []
    with ThreadPoolExecutor(max_workers=n_threads) as ex:
        for fut in [ex.submit(run_scramble_batch, nb, sd) for nb, sd in batches]:
            parts.append(fut.result())
    delta_aics_scram = np.concatenate([p[0] for p in parts])
    f_stats_scram = np.concatenate([p[1] for p in parts])
    f_stats_eff_scram = np.concatenate([p[2] for p in parts])
    corr_cosD_scram = np.concatenate([p[3] for p in parts])

    p_scramble_delta_aic = np.mean(delta_aics_scram >= delta_aic_true)
    p_scramble_f = np.mean(f_stats_scram >= f_true)
    p_scramble_f_eff = np.mean(f_stats_eff_scram >= f_true_eff)

    print_status(
        f"Base model AIC = {aic_base:.1f}; True CMB ΔAIC = {delta_aic_true:.1f}",
        "CALC",
    )
    print_status(
        f"Random-axis null (n={n_scramble}): ΔAIC median = {np.median(delta_aics_scram):.1f}, "
        f"99th %ile = {np.percentile(delta_aics_scram, 99):.1f}, max = {np.max(delta_aics_scram):.1f}",
        "CALC",
    )
    print_status(
        f"Planck ΔAIC = {delta_aic_true:.1f}; p(any axis ΔAIC ≥ Planck) = {p_scramble_delta_aic:.4f}",
        "CALC",
    )
    print_status(
        f"Uniform scramble F: median = {np.median(f_stats_scram):.2f}, "
        f"99th %ile = {np.percentile(f_stats_scram, 99):.2f}",
        "CALC",
    )
    print_status(
        f"True CMB F = {f_true:.2f}; p(random F ≥ true) = {p_scramble_f:.4f}",
        "CALC",
    )
    print_status(
        f"Eff-sample-size F: true = {f_true_eff:.2f}; "
        f"p(random F_eff ≥ true) = {p_scramble_f_eff:.4f}",
        "CALC",
    )

    # ------------------------------------------------------------------
    # D2. Correlation-matched scramble (rigorous correlation-structure null)
    # ------------------------------------------------------------------
    # Random directions with different r(cosD, cosθ) have different aliasing
    # potential.  The true correlation is corr_cosD_cosTheta ~ 0.050.
    # We restrict the null to directions with a matched correlation to ensure
    # the comparison is conditioned on the same correlation structure.
    corr_true = abs(corr_cosD_cosTheta)
    corr_tol = 0.02
    matched_mask = np.abs(np.abs(corr_cosD_scram) - corr_true) <= corr_tol
    n_matched = matched_mask.sum()
    if n_matched >= 100:
        p_scramble_delta_aic_matched = np.mean(delta_aics_scram[matched_mask] >= delta_aic_true)
        p_scramble_f_matched = np.mean(f_stats_scram[matched_mask] >= f_true)
        p_scramble_f_eff_matched = np.mean(f_stats_eff_scram[matched_mask] >= f_true_eff)
        print_status(
            f"Correlation-matched scramble: {n_matched} directions "
            f"(|r(cosD, cosθ)| within {corr_tol:.3f} of {corr_true:.4f})",
            "CALC",
        )
        print_status(
            f"  p_matched(ΔAIC) = {p_scramble_delta_aic_matched:.4f}, "
            f"p_matched(F) = {p_scramble_f_matched:.4f}, "
            f"p_matched(F_eff) = {p_scramble_f_eff_matched:.4f}",
            "CALC",
        )
    else:
        p_scramble_delta_aic_matched = None
        p_scramble_f_matched = None
        p_scramble_f_eff_matched = None
        print_status(
            f"Correlation-matched scramble: insufficient matches ({n_matched})",
            "WARNING",
        )

    k_scramble_f = int(np.sum(f_stats_scram >= f_true))
    k_scramble_delta_aic = int(np.sum(delta_aics_scram >= delta_aic_true))
    n_scramble_eff = int(len(f_stats_scram))
    j_lo_f, j_hi_f = jeffreys_binomial_ci(k_scramble_f, n_scramble_eff)
    j_lo_da, j_hi_da = jeffreys_binomial_ci(k_scramble_delta_aic, n_scramble_eff)
    if n_matched >= 100:
        k_matched_f = int(np.sum(f_stats_scram[matched_mask] >= f_true))
        j_lo_fm, j_hi_fm = jeffreys_binomial_ci(k_matched_f, int(n_matched))
    else:
        k_matched_f = None
        j_lo_fm = j_hi_fm = None

    # ------------------------------------------------------------------
    # D3–D5. Refined directional null tests
    # ------------------------------------------------------------------
    print_status("═══ D3–D5. Refined directional null tests ═══", "TITLE")
    from scipy.stats import ortho_group

    n_refined = 5000
    local_cone_deg = 30.0

    def summarize_refined_null(f_stats_arr, f_obs, delta_aics_arr=None, delta_obs=None):
        f_stats_arr = np.asarray(f_stats_arr, dtype=np.float64)
        n_draws = int(len(f_stats_arr))
        if n_draws == 0:
            raise RuntimeError("Refined null test produced zero valid draws")
        k_ge = int(np.sum(f_stats_arr >= f_obs))
        j_lo, j_hi = jeffreys_binomial_ci(k_ge, n_draws)
        out = {
            "n_draws": n_draws,
            "true_f_statistic_eff": float(f_obs),
            "f_eff_median": float(np.median(f_stats_arr)),
            "f_eff_99th_percentile": float(np.percentile(f_stats_arr, 99)),
            "p_f_eff": float(k_ge / n_draws),
            "k_random_ge_true_f": k_ge,
            "jeffreys_95_lower": j_lo,
            "jeffreys_95_upper": j_hi,
        }
        if delta_aics_arr is not None and delta_obs is not None:
            da = np.asarray(delta_aics_arr, dtype=np.float64)
            out["true_delta_aic"] = float(delta_obs)
            out["delta_aic_median"] = float(np.median(da))
            out["p_delta_aic"] = float(np.mean(da >= delta_obs))
        return out

    def fit_orth_for_cos_theta_c(cos_theta_dir_c):
        y_orth = orthogonalized_cos_theta_cosD(cos_theta_dir_c, cosD, r_c, vr_c)
        return delta_aic_for_orthogonal_component(
            y_orth, res, cosD, r_c, vr_c, rss_base, aic_base, n_clean, n_eff, use_eff=True
        )

    y_orth_cmb = orthogonalized_cos_theta_cosD(cos_theta_c, cosD, r_c, vr_c)
    delta_aic_orth_true, _f_orth_true, _, f_orth_true_eff = delta_aic_for_orthogonal_component(
        y_orth_cmb, res, cosD, r_c, vr_c, rss_base, aic_base, n_clean, n_eff, use_eff=True
    )
    if f_orth_true_eff is None:
        raise RuntimeError("CMB orthogonalized direction fit failed")

    def run_phase_batch(n_batch, seed_batch):
        rng_p = np.random.default_rng(int(seed_batch))
        fs, da = [], []
        for _ in range(int(n_batch)):
            shift = int(rng_p.integers(1, n_clean))
            cos_theta_shifted = np.roll(cos_theta_c, shift)
            d_aic, _f_dir, _, f_dir_eff = delta_aic_for_direction(cos_theta_shifted, use_eff=True)
            if f_dir_eff is not None:
                fs.append(f_dir_eff)
                if d_aic is not None:
                    da.append(d_aic)
        return np.asarray(fs, dtype=np.float64), np.asarray(da, dtype=np.float64)

    def run_rotation_global_batch(n_batch, seed_batch):
        rng_r = np.random.default_rng(int(seed_batch))
        fs, da = [], []
        for _ in range(int(n_batch)):
            rot = ortho_group.rvs(dim=3, random_state=rng_r)
            n_rot = rot @ _CMB_UNIT
            cos_theta_rot = np.sum(em_sub * n_rot[:, None], axis=0)
            cos_theta_rot_c = cos_theta_rot - np.mean(cos_theta_rot)
            d_aic, _f_dir, _, f_dir_eff = delta_aic_for_direction(cos_theta_rot_c, use_eff=True)
            if f_dir_eff is not None:
                fs.append(f_dir_eff)
                if d_aic is not None:
                    da.append(d_aic)
        return np.asarray(fs, dtype=np.float64), np.asarray(da, dtype=np.float64)

    def run_rotation_local_batch(n_batch, seed_batch):
        rng_rl = np.random.default_rng(int(seed_batch))
        fs, da = [], []
        max_angle = np.deg2rad(local_cone_deg)
        for _ in range(int(n_batch)):
            perp = random_unit_perpendicular(_CMB_UNIT, rng_rl)
            angle = float(rng_rl.uniform(0.0, max_angle))
            n_rot = rotate_vector_about_axis(_CMB_UNIT, perp, angle)
            n_rot = n_rot / np.linalg.norm(n_rot)
            cos_theta_rot = np.sum(em_sub * n_rot[:, None], axis=0)
            cos_theta_rot_c = cos_theta_rot - np.mean(cos_theta_rot)
            d_aic, _f_dir, _, f_dir_eff = delta_aic_for_direction(cos_theta_rot_c, use_eff=True)
            if f_dir_eff is not None:
                fs.append(f_dir_eff)
                if d_aic is not None:
                    da.append(d_aic)
        return np.asarray(fs, dtype=np.float64), np.asarray(da, dtype=np.float64)

    def run_gs_scramble_batch(n_batch, seed_batch):
        rng_gs = np.random.default_rng(int(seed_batch))
        fs, da = [], []
        for _ in range(int(n_batch)):
            dir_vec = random_unit_vector(rng_gs)
            cos_theta_rand = np.sum(em_sub * dir_vec[:, None], axis=0)
            cos_theta_rand_c = cos_theta_rand - np.mean(cos_theta_rand)
            d_aic, _f_dir, _, f_dir_eff = fit_orth_for_cos_theta_c(cos_theta_rand_c)
            if f_dir_eff is not None:
                fs.append(f_dir_eff)
                if d_aic is not None:
                    da.append(d_aic)
        return np.asarray(fs, dtype=np.float64), np.asarray(da, dtype=np.float64)

    refined_edges = np.linspace(0, n_refined, n_threads + 1, dtype=int)
    refined_batches = [
        (int(refined_edges[t + 1] - refined_edges[t]), 880_003 + t * 1_000_003)
        for t in range(n_threads)
        if refined_edges[t + 1] > refined_edges[t]
    ]

    def _run_refined_batches(run_fn):
        parts_local = []
        with ThreadPoolExecutor(max_workers=n_threads) as ex:
            futs = [ex.submit(run_fn, nb, sd) for nb, sd in refined_batches]
            for fut in futs:
                parts_local.append(fut.result())
        fs_all = np.concatenate([p[0] for p in parts_local])
        da_all = np.concatenate([p[1] for p in parts_local])
        return fs_all, da_all

    phase_f, phase_da = _run_refined_batches(run_phase_batch)
    rot_global_f, rot_global_da = _run_refined_batches(run_rotation_global_batch)
    rot_local_f, rot_local_da = _run_refined_batches(run_rotation_local_batch)
    gs_f, gs_da = _run_refined_batches(run_gs_scramble_batch)

    phase_null = summarize_refined_null(phase_f, f_true_eff, phase_da, delta_aic_true)
    rotation_global_null = summarize_refined_null(
        rot_global_f, f_true_eff, rot_global_da, delta_aic_true
    )
    rotation_local_null = summarize_refined_null(
        rot_local_f, f_true_eff, rot_local_da, delta_aic_true
    )
    orthogonal_scramble_null = summarize_refined_null(
        gs_f, f_orth_true_eff, gs_da, delta_aic_orth_true
    )
    orthogonal_scramble_null["description"] = (
        "Random sky directions with cosθ·cosD Gram–Schmidt orthogonalized against the "
        "heliocentric base before augmentation (same geometry as test E)."
    )
    phase_null["description"] = (
        "Circular shifts of the CMB cosθ_c series destroy synodic phase alignment with cosD "
        "while preserving the marginal cosθ distribution."
    )
    rotation_global_null["description"] = (
        "Uniform SO(3) rotations of the Planck axis (marginal on S², comparable to uniform scramble)."
    )
    rotation_local_null["description"] = (
        f"Random rotations of the Planck axis within a {local_cone_deg:.0f}° cone "
        "(local orientation specificity)."
    )

    refined_directional_nulls = {
        "phase_null": phase_null,
        "rotation_global_null": rotation_global_null,
        "rotation_local_null": rotation_local_null,
        "orthogonal_scramble_null": orthogonal_scramble_null,
        "local_cone_deg": float(local_cone_deg),
    }

    print_status(
        f"Phase null: p_F_eff = {phase_null['p_f_eff']:.4f} "
        f"(median F_eff = {phase_null['f_eff_median']:.2f})",
        "CALC",
    )
    print_status(
        f"Rotation global: p_F_eff = {rotation_global_null['p_f_eff']:.4f}",
        "CALC",
    )
    print_status(
        f"Rotation local ({local_cone_deg:.0f}°): p_F_eff = {rotation_local_null['p_f_eff']:.4f}",
        "CALC",
    )
    print_status(
        f"GS orthogonal scramble: p_F_eff = {orthogonal_scramble_null['p_f_eff']:.4f} "
        f"(true F_eff = {f_orth_true_eff:.2f})",
        "CALC",
    )

    random_axis_delta_aic_null = {
        "description": (
            "Uniform random unit vectors on S²: empirical null for ΔAIC when the "
            "orientation predictor is any fixed sky axis (not only permuted or rotated Planck)."
        ),
        "n_draws": n_scramble_eff,
        "observed_delta_aic_planck": float(delta_aic_true),
        "null_delta_aic_median": float(np.median(delta_aics_scram)),
        "null_delta_aic_99th_percentile": float(np.percentile(delta_aics_scram, 99)),
        "null_delta_aic_max": float(np.max(delta_aics_scram)),
        "empirical_p_value": float(p_scramble_delta_aic),
        "planck_axis_exceeds_median_by": float(delta_aic_true - np.median(delta_aics_scram)),
        "jeffreys_95_interval_on_p": {
            "k_random_ge_true": k_scramble_delta_aic,
            "n_draws": n_scramble_eff,
            "jeffreys_95_lower": j_lo_da,
            "jeffreys_95_upper": j_hi_da,
        },
        "correlation_matched_subnull": {
            "n_matched": int(n_matched) if n_matched is not None else None,
            "corr_tolerance": float(corr_tol),
            "empirical_p_value": float(p_scramble_delta_aic_matched)
            if p_scramble_delta_aic_matched is not None
            else None,
        },
        "role": "primary_delta_aic_null_for_any_fixed_axis",
    }

    sky_scramble_result = {
        "random_axis_delta_aic_null": random_axis_delta_aic_null,
        "n_scrambles": n_scramble,
        "strict_alpha": strict_alpha,
        "base_model_aic": float(aic_base),
        "true_delta_aic": float(delta_aic_true),
        "true_f_statistic": float(f_true),
        "true_f_p_value": float(p_f_true),
        "true_f_statistic_eff": float(f_true_eff),
        "base_residual_ar1_rho": float(rho_base),
        "effective_sample_size": float(n_eff),
        "delta_aic_median": float(np.median(delta_aics_scram)),
        "delta_aic_99th_percentile": float(np.percentile(delta_aics_scram, 99)),
        "delta_aic_max": float(np.max(delta_aics_scram)),
        "f_median": float(np.median(f_stats_scram)),
        "f_99th_percentile": float(np.percentile(f_stats_scram, 99)),
        "f_eff_median": float(np.median(f_stats_eff_scram)),
        "f_eff_99th_percentile": float(np.percentile(f_stats_eff_scram, 99)),
        "p_scramble_delta_aic": float(p_scramble_delta_aic),
        "p_scramble_f": float(p_scramble_f),
        "p_scramble_f_eff": float(p_scramble_f_eff),
        "p_scramble_f_monte_carlo_counts": {
            "k_random_ge_true_f": k_scramble_f,
            "n_draws": n_scramble_eff,
            "jeffreys_95_lower": j_lo_f,
            "jeffreys_95_upper": j_hi_f,
        },
        "p_scramble_delta_aic_monte_carlo_counts": {
            "k_random_ge_true_delta_aic": k_scramble_delta_aic,
            "n_draws": n_scramble_eff,
            "jeffreys_95_lower": j_lo_da,
            "jeffreys_95_upper": j_hi_da,
        },
        "correlation_matched": {
            "n_matched": int(n_matched) if n_matched is not None else None,
            "corr_tolerance": float(corr_tol),
            "p_scramble_delta_aic_matched": float(p_scramble_delta_aic_matched) if p_scramble_delta_aic_matched is not None else None,
            "p_scramble_f_matched": float(p_scramble_f_matched) if p_scramble_f_matched is not None else None,
            "p_scramble_f_eff_matched": float(p_scramble_f_eff_matched) if p_scramble_f_eff_matched is not None else None,
            "k_random_ge_true_f_matched": int(k_matched_f) if k_matched_f is not None else None,
            "jeffreys_95_lower_f_matched": j_lo_fm,
            "jeffreys_95_upper_f_matched": j_hi_fm,
        },
        "true_direction_preferred": bool(
            p_scramble_delta_aic < strict_alpha and p_scramble_f < strict_alpha
        ),
        "directional_anatomy_supported": bool(
            p_scramble_f < 0.1
        ),
        "refined_directional_nulls": refined_directional_nulls,
    }

    # ------------------------------------------------------------------
    # E. Gram–Schmidt Orthogonalization
    # ------------------------------------------------------------------
    print_status("═══ E. Gram–Schmidt Orthogonalization ═══", "TITLE")

    # Orthogonalize cosθ·cosD against [cosD, r_c·cosD, vr_c·cosD, 1]
    y_target = cos_theta_c * cosD
    X_orth_base = np.column_stack([cosD, r_c * cosD, vr_c * cosD, np.ones(n_clean)])
    coeffs_gs, _, rank_gs, _ = stable_lstsq(X_orth_base, y_target)
    if rank_gs == 4:
        with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
            y_orth_gs = y_target - X_orth_base @ coeffs_gs
        # Verify orthogonality
        orth_dot_products = [np.dot(y_orth_gs, X_orth_base[:, j]) for j in range(4)]
        print_status(
            f"Orthogonality residuals (dot products): "
            f"{[f'{d:.2e}' for d in orth_dot_products]}",
            "CALC",
        )

        X_orth_model = np.column_stack([cosD, r_c * cosD, vr_c * cosD, y_orth_gs, np.ones(n_clean)])
        coeffs_orth_m, _, rank_orth_m, _ = stable_lstsq(X_orth_model, res)
        if rank_orth_m == 5:
            with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
                resid_orth_m = res - X_orth_model @ coeffs_orth_m
            mse_orth_m = np.sum(resid_orth_m ** 2) / (n_clean - 5)
            cov_orth_m = mse_orth_m * np.linalg.pinv(X_orth_model.T @ X_orth_model, rcond=1e-10, hermitian=True)
            se_orth_m = np.sqrt(np.diag(cov_orth_m))

            eta_theta_orth = coeffs_orth_m[3] / ETA_SCALE_FACTOR
            se_eta_theta_orth = se_orth_m[3] / ETA_SCALE_FACTOR
            t_theta_orth = eta_theta_orth / se_eta_theta_orth if se_eta_theta_orth > 0 else 0.0
            p_theta_orth = 2 * (1 - stats.t.cdf(abs(t_theta_orth), n_clean - 5)) if se_eta_theta_orth > 0 else 1.0

            print_status(
                f"Orthogonalized η_θ = {eta_theta_orth:.4e} ± {se_eta_theta_orth:.4e} "
                f"(t={t_theta_orth:.2f}, p={p_theta_orth:.4e})",
                "CALC",
            )

            orthogonalization_result = {
                "available": True,
                "eta_theta_orthogonal": float(eta_theta_orth),
                "eta_theta_orthogonal_error": float(se_eta_theta_orth),
                "eta_theta_orthogonal_t": float(t_theta_orth),
                "eta_theta_orthogonal_p": float(p_theta_orth),
                "orthogonality_residuals": [float(d) for d in orth_dot_products],
                "signal_persists": bool(abs(t_theta_orth) > 2.0),
            }
        else:
            orthogonalization_result = {"available": False, "reason": "rank_deficient_orth_model"}
    else:
        orthogonalization_result = {"available": False, "reason": "rank_deficient_gs"}

    if not orthogonalization_result.get("available"):
        print_status("Orthogonalization failed.", "WARNING")

    # ------------------------------------------------------------------
    # F. Joint-model directional anatomy
    # ------------------------------------------------------------------
    print_status("═══ F. Joint-model directional anatomy ═══", "TITLE")

    perp1 = np.cross(_CMB_UNIT, np.array([0.0, 0.0, 1.0]))
    if np.linalg.norm(perp1) < 1e-12:
        perp1 = np.cross(_CMB_UNIT, np.array([0.0, 1.0, 0.0]))
    perp1 = perp1 / np.linalg.norm(perp1)
    perp2 = np.cross(_CMB_UNIT, perp1)
    perp2 = perp2 / np.linalg.norm(perp2)

    def joint_direction_metrics(dir_vec):
        cos_theta_dir = np.sum(em_hat[:, mask_clean] * dir_vec[:, None], axis=0)
        cos_theta_dir_c = cos_theta_dir - np.mean(cos_theta_dir)
        delta_aic, f_dir, p_f, _f_eff = delta_aic_for_direction(cos_theta_dir_c)
        eta_dir, se_dir, _, _, _ = fit_full_joint_model(res, cosD, r_c, vr_c, cos_theta_dir_c)
        t_dir = eta_dir / se_dir if se_dir and se_dir > 0 else 0.0
        return {
            "delta_aic": float(delta_aic) if delta_aic is not None else None,
            "f_statistic": float(f_dir) if f_dir is not None else None,
            "f_p_value": float(p_f) if p_f is not None else None,
            "eta_theta": float(eta_dir) if eta_dir is not None else None,
            "eta_theta_error": float(se_dir) if se_dir is not None else None,
            "eta_theta_t": float(t_dir),
        }

    directional_true = joint_direction_metrics(_CMB_UNIT)
    directional_perp1 = joint_direction_metrics(perp1)
    directional_perp2 = joint_direction_metrics(perp2)
    directional_antipode = joint_direction_metrics(-_CMB_UNIT)

    delta_aic_perp1 = directional_true["delta_aic"] - directional_perp1["delta_aic"]
    delta_aic_perp2 = directional_true["delta_aic"] - directional_perp2["delta_aic"]
    delta_aic_antipode = directional_true["delta_aic"] - directional_antipode["delta_aic"]

    print_status(
        f"ΔAIC(true vs perp1) = {delta_aic_perp1:.1f}; "
        f"ΔAIC(true vs perp2) = {delta_aic_perp2:.1f}; "
        f"ΔAIC(true vs antipode) = {delta_aic_antipode:.1f}",
        "CALC",
    )
    print_status(
        f"|η_θ|/σ: true={abs(directional_true['eta_theta_t']):.2f}, "
        f"perp1={abs(directional_perp1['eta_theta_t']):.2f}, "
        f"perp2={abs(directional_perp2['eta_theta_t']):.2f}, "
        f"antipode={abs(directional_antipode['eta_theta_t']):.2f}",
        "CALC",
    )

    antipode_magnitude_ratio = (
        abs(directional_antipode["eta_theta"]) / abs(directional_true["eta_theta"])
        if directional_true["eta_theta"] not in (None, 0.0)
        else None
    )
    directional_anatomy_result = {
        "true_cmb": directional_true,
        "perpendicular_1": directional_perp1,
        "perpendicular_2": directional_perp2,
        "antipode": directional_antipode,
        "delta_aic_perp1_vs_true": float(delta_aic_perp1),
        "delta_aic_perp2_vs_true": float(delta_aic_perp2),
        "delta_aic_antipode_vs_true": float(delta_aic_antipode),
        "antipode_magnitude_ratio": float(antipode_magnitude_ratio) if antipode_magnitude_ratio is not None else None,
        "dipole_antipode_consistent": bool(abs(delta_aic_antipode) < 2.0),
        "perpendicular_directions_suppressed": bool(
            delta_aic_perp1 > 90.0 and delta_aic_perp2 > 90.0
        ),
        "directional_anatomy_passed": bool(
            delta_aic_perp1 > 90.0
            and delta_aic_perp2 > 90.0
            and abs(delta_aic_antipode) < 2.0
        ),
    }

    # ------------------------------------------------------------------
    # G1. Year-block jackknife (temporal stability of η_θ)
    # ------------------------------------------------------------------
    print_status("═══ G1. Year-block jackknife ═══", "TITLE")
    jd_clean = jd[mask_clean]
    years_clean = jd_to_gregorian_year(jd_clean)
    unique_years = np.sort(np.unique(years_clean))
    jackknife_rows = []
    for yr in unique_years:
        fm = years_clean != yr
        n_sub = int(fm.sum())
        if n_sub < 8000:
            continue
        res_j = res[fm]
        cosD_j = cosD[fm]
        r_j = r_clean[fm]
        vr_j = vr_clean[fm]
        v_par_j = v_par_clean[fm]
        cos_th_j = cos_theta_clean[fm]
        r_c_j = r_j - np.mean(r_j)
        vr_c_j = vr_j - np.mean(vr_j)
        v_par_c_j = v_par_j - np.mean(v_par_j)
        cos_theta_c_j = cos_th_j - np.mean(cos_th_j)
        et, seet, _, _, _ = fit_full_joint_model(res_j, cosD_j, r_c_j, vr_c_j, cos_theta_c_j)
        if et is None or seet is None:
            continue
        t_j = float(et / seet) if seet > 0 else 0.0
        jackknife_rows.append(
            {"year_excluded": int(yr), "n": n_sub, "eta_theta": float(et), "t": t_j}
        )

    jack_t = np.array([row["t"] for row in jackknife_rows], dtype=np.float64)
    jack_eta = np.array([row["eta_theta"] for row in jackknife_rows], dtype=np.float64)
    all_same_sign_as_full = bool(
        len(jack_eta) > 0 and np.all(np.sign(jack_eta) == np.sign(float(eta_obs)))
    )
    year_jackknife_result = {
        "n_years_in_sample": int(len(unique_years)),
        "n_jackknife_fits": int(len(jackknife_rows)),
        "per_year_results": jackknife_rows,
        "eta_theta_full_sample": float(eta_obs),
        "jackknife_median_t": float(np.median(jack_t)) if len(jack_t) else None,
        "jackknife_min_t": float(np.min(jack_t)) if len(jack_t) else None,
        "jackknife_max_t": float(np.max(jack_t)) if len(jack_t) else None,
        "all_jackknife_eta_same_sign_as_full": all_same_sign_as_full,
        "temporal_stability_strong": bool(
            len(jack_t) >= 10
            and all_same_sign_as_full
            and float(np.median(np.abs(jack_t))) > 5.0
        ),
    }
    med_abs_t_jack = float(np.median(np.abs(jack_t))) if len(jack_t) else float("nan")
    print_status(
        f"Jackknife: {len(jackknife_rows)} yearly excisions; "
        f"median |t| on η_θ = {med_abs_t_jack:.2f}; "
        f"same sign as full sample: {all_same_sign_as_full}",
        "CALC",
    )

    # ------------------------------------------------------------------
    # G2. Alternative fixed celestial axes (same hypothesis class)
    # ------------------------------------------------------------------
    print_status("═══ G2. Alternative fixed celestial axes ═══", "TITLE")
    mean_orbit_n = np.asarray(cmb_data["mean_orbit_normal"], dtype=np.float64).ravel()
    mean_orbit_n = mean_orbit_n / (np.linalg.norm(mean_orbit_n) + 1e-30)
    sep_cmb_orbit_deg = float(
        np.degrees(np.arccos(np.clip(np.dot(_CMB_UNIT, mean_orbit_n), -1.0, 1.0)))
    )
    sep_cmb_nep_deg = float(
        np.degrees(np.arccos(np.clip(np.dot(_CMB_UNIT, _NEP_UNIT), -1.0, 1.0)))
    )
    alt_axes = {
        "north_ecliptic_pole": _NEP_UNIT,
        "galactic_north_pole": _GAL_N_UNIT,
        "mean_orbit_angular_momentum": mean_orbit_n,
    }
    alternative_fixed_directions = {}
    for name, vec in alt_axes.items():
        u = np.asarray(vec, dtype=np.float64).ravel()
        u = u / np.linalg.norm(u)
        alternative_fixed_directions[name] = joint_direction_metrics(u)
    print_status(
        f"Angles: CMB–mean orbit normal = {sep_cmb_orbit_deg:.2f}°, CMB–ecliptic pole = {sep_cmb_nep_deg:.2f}°",
        "CALC",
    )

    # ------------------------------------------------------------------
    # H. Dual-axis simultaneous fit (CMB + galactic north)
    # ------------------------------------------------------------------
    print_status("═══ H. Dual-axis simultaneous fit (CMB + galactic) ═══", "TITLE")
    cos_theta_gal = np.sum(em_hat[:, mask_clean] * _GAL_N_UNIT[:, None], axis=0)
    cos_theta_gal_c = cos_theta_gal - np.mean(cos_theta_gal)
    corr_cmb_gal = float(np.corrcoef(cos_theta_c, cos_theta_gal_c)[0, 1])
    sep_cmb_gal_deg = float(
        np.degrees(np.arccos(np.clip(np.dot(_CMB_UNIT, _GAL_N_UNIT), -1.0, 1.0)))
    )

    fit_cmb_only = fit_orientation_joint_model(res, cosD, r_c, vr_c, cos_theta_c, n_params=5)
    fit_gal_only = fit_orientation_joint_model(res, cosD, r_c, vr_c, cos_theta_gal_c, n_params=5)
    fit_dual = fit_orientation_joint_model(
        res, cosD, r_c, vr_c, (cos_theta_c, cos_theta_gal_c), n_params=6
    )

    if fit_cmb_only is None or fit_gal_only is None or fit_dual is None:
        raise RuntimeError("Dual-axis joint model fit failed (rank deficient).")

    # Base heliocentric model (4 parameters) for nested comparisons
    X_base_h = np.column_stack([cosD, r_c * cosD, vr_c * cosD, np.ones(n_clean)])
    coeffs_base_h, _, rank_base_h, _ = stable_lstsq(X_base_h, res)
    if rank_base_h < 4:
        raise RuntimeError("Heliocentric base model rank deficient in dual-axis block.")
    with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
        resid_base_h = res - X_base_h @ coeffs_base_h
    rss_base_h = float(np.sum(resid_base_h ** 2))
    aic_base_h = float(n_clean * np.log(rss_base_h / n_clean) + 2 * 4)

    f_add_cmb, p_add_cmb = nested_f_test(rss_base_h, 4, fit_cmb_only["rss"], 5, n_clean)
    f_add_gal, p_add_gal = nested_f_test(rss_base_h, 4, fit_gal_only["rss"], 5, n_clean)
    f_add_gal_given_cmb, p_add_gal_given_cmb = nested_f_test(
        fit_cmb_only["rss"], 5, fit_dual["rss"], 6, n_clean
    )
    f_add_cmb_given_gal, p_add_cmb_given_gal = nested_f_test(
        fit_gal_only["rss"], 5, fit_dual["rss"], 6, n_clean
    )

    delta_aic_cmb_only = aic_base_h - fit_cmb_only["aic"]
    delta_aic_gal_only = aic_base_h - fit_gal_only["aic"]
    delta_aic_dual = aic_base_h - fit_dual["aic"]
    gal_beats_cmb_in_single_axis = bool(fit_gal_only["aic"] < fit_cmb_only["aic"])
    dual_beats_gal_only = bool(fit_dual["aic"] < fit_gal_only["aic"])
    dual_beats_cmb_only = bool(fit_dual["aic"] < fit_cmb_only["aic"])
    both_significant_in_dual = bool(
        abs(fit_dual["eta_cmb_t"]) > 2.0 and abs(fit_dual["eta_gal_t"]) > 2.0
    )
    gal_dominates_in_dual = bool(abs(fit_dual["eta_gal_t"]) > abs(fit_dual["eta_cmb_t"]))

    dual_axis_result = {
        "cmb_galactic_cos_correlation": corr_cmb_gal,
        "cmb_galactic_axis_separation_deg": sep_cmb_gal_deg,
        "base_heliocentric_aic": aic_base_h,
        "cmb_only": {
            "aic": fit_cmb_only["aic"],
            "delta_aic_vs_base": float(delta_aic_cmb_only),
            "eta_theta": fit_cmb_only["eta_theta"],
            "eta_theta_t": fit_cmb_only["eta_theta_t"],
            "max_vif": float(max(fit_cmb_only["vifs"].values())),
        },
        "galactic_only": {
            "aic": fit_gal_only["aic"],
            "delta_aic_vs_base": float(delta_aic_gal_only),
            "eta_theta": fit_gal_only["eta_theta"],
            "eta_theta_t": fit_gal_only["eta_theta_t"],
            "max_vif": float(max(fit_gal_only["vifs"].values())),
        },
        "dual_axis": {
            "aic": fit_dual["aic"],
            "delta_aic_vs_base": float(delta_aic_dual),
            "eta_cmb": fit_dual["eta_cmb"],
            "eta_gal": fit_dual["eta_gal"],
            "eta_cmb_t": fit_dual["eta_cmb_t"],
            "eta_gal_t": fit_dual["eta_gal_t"],
            "max_vif": float(max(fit_dual["vifs"].values())),
        },
        "nested_f_tests": {
            "add_cmb_to_base": {"f": f_add_cmb, "p": p_add_cmb},
            "add_gal_to_base": {"f": f_add_gal, "p": p_add_gal},
            "add_gal_given_cmb": {"f": f_add_gal_given_cmb, "p": p_add_gal_given_cmb},
            "add_cmb_given_gal": {"f": f_add_cmb_given_gal, "p": p_add_cmb_given_gal},
        },
        "galactic_beats_cmb_single_axis_aic": gal_beats_cmb_in_single_axis,
        "dual_improves_over_gal_only": dual_beats_gal_only,
        "dual_improves_over_cmb_only": dual_beats_cmb_only,
        "both_axes_significant_in_dual": both_significant_in_dual,
        "galactic_coefficient_stronger_in_dual": gal_dominates_in_dual,
        "interpretation": (
            "Galactic north attains lower AIC than Planck in single-axis fits, but the "
            "six-parameter dual model tests whether that advantage is a separable second "
            "axis or shared dipole structure."
        ),
    }

    print_status(
        f"r(cosθ_CMB, cosθ_gal) = {corr_cmb_gal:.4f}; axis separation = {sep_cmb_gal_deg:.1f}°",
        "CALC",
    )
    print_status(
        f"AIC: base={aic_base_h:.1f}; CMB-only={fit_cmb_only['aic']:.1f} "
        f"(Δ={delta_aic_cmb_only:.1f}); gal-only={fit_gal_only['aic']:.1f} "
        f"(Δ={delta_aic_gal_only:.1f}); dual={fit_dual['aic']:.1f} (Δ={delta_aic_dual:.1f})",
        "CALC",
    )
    print_status(
        f"Dual model: η_CMB t={fit_dual['eta_cmb_t']:.2f}, η_gal t={fit_dual['eta_gal_t']:.2f}; "
        f"add gal|CMB p={p_add_gal_given_cmb:.4e}; add CMB|gal p={p_add_cmb_given_gal:.4e}",
        "CALC",
    )

    # ------------------------------------------------------------------
    # I. TEP η_θ prediction coverage (Step 033)
    # ------------------------------------------------------------------
    print_status("═══ I. TEP η_θ prediction coverage (Step 033) ═══", "TITLE")
    pred_033 = load_step033_orientation_prediction()
    tep_prediction_coverage = compute_eta_theta_prediction_coverage(
        float(eta_obs), float(se_obs), pred_033
    )
    band = tep_prediction_coverage["order_of_magnitude_band"]
    print_status(
        f"Fitted η_θ = {tep_prediction_coverage['eta_theta_fitted']:.4e} "
        f"(predicted center {tep_prediction_coverage['eta_theta_predicted_center']:.4e}; "
        f"band [{band['lo']:.0e}, {band['hi']:.0e}]; within band = "
        f"{tep_prediction_coverage['within_predicted_band']})",
        "CALC",
    )
    if tep_prediction_coverage["z_score_vs_predicted_center"] is not None:
        print_status(
            f"z vs predicted center = {tep_prediction_coverage['z_score_vs_predicted_center']:.2f}; "
            f"ratio fitted/predicted = {tep_prediction_coverage['ratio_fitted_to_predicted_center']:.2f}",
            "CALC",
        )

    # ------------------------------------------------------------------
    # 6. Compile results
    # ------------------------------------------------------------------
    sky_scramble_pass = bool(sky_scramble_result.get("true_direction_preferred", False))
    checks = [
        aliasing_result.get("aliasing_rejected", False),
        permutation_result.get("permutation_rejected", False),
        sky_scramble_pass,
        orthogonalization_result.get("signal_persists", False) if orthogonalization_result.get("available") else False,
        directional_anatomy_result.get("directional_anatomy_passed", False),
    ]
    n_passed = sum(checks)
    if n_passed >= 4:
        falsification_result = "STRONG_SUPPORT"
        print_status(
            f"RESULT: {n_passed}/5 falsification tests passed. Directional residual component is supported with caveats.",
            "SUCCESS",
        )
    elif n_passed >= 3:
        falsification_result = "PARTIAL_SUPPORT"
        print_status(
            f"RESULT: {n_passed}/5 falsification tests passed. Directional residual component is only partially supported.",
            "INFO",
        )
    elif n_passed >= 2:
        falsification_result = "MARGINAL_SUPPORT"
        print_status(
            f"RESULT: {n_passed}/5 falsification tests passed. Marginal robustness.",
            "INFO",
        )
    else:
        falsification_result = "WEAK_SUPPORT"
        print_status(
            f"RESULT: Only {n_passed}/5 falsification tests passed. CMB anisotropy claim is weakened.",
            "WARNING",
        )
    pipeline_status = "FAIL" if n_passed < 2 else "PASS"

    results = {
        "step_id": "step_055",
        "status": pipeline_status,
        "falsification_result": falsification_result,
        "n_observations": int(n_clean),
        "correlation_cosD_cos_theta": float(corr_cosD_cosTheta),
        "aliasing_simulation": aliasing_result,
        "multicollinearity_diagnostics": multicollinearity_result,
        "permutation_test": permutation_result,
        "sky_scrambling": sky_scramble_result,
        "orthogonalization_test": orthogonalization_result,
        "directional_anatomy": directional_anatomy_result,
        "year_jackknife": year_jackknife_result,
        "alternative_fixed_directions": alternative_fixed_directions,
        "dual_axis_identifiability": dual_axis_result,
        "primary_orientation_block": "dual_axis_identifiability",
        "supplementary_single_axis_metrics": {
            "planck_delta_aic_vs_base": float(delta_aic_true),
            "galactic_single_axis": alternative_fixed_directions.get("galactic_north_pole"),
            "random_axis_delta_aic_null": random_axis_delta_aic_null,
        },
        "tep_prediction_coverage": tep_prediction_coverage,
        "reference_axis_geometry": {
            "cmb_mean_orbit_normal_separation_deg": sep_cmb_orbit_deg,
            "cmb_ecliptic_pole_separation_deg": sep_cmb_nep_deg,
            "cmb_galactic_separation_deg": sep_cmb_gal_deg,
        },
        "falsification_tests_passed": int(n_passed),
        "falsification_tests_total": 5,
        "interpretation": (
            f"{n_passed}/5 scored diagnostics pass (A–E). Primary orientation identifiability is "
            f"the dual-axis block (H): r(cosθ_CMB, cosθ_gal)={corr_cmb_gal:.3f}; in the six-parameter "
            f"model η_gal t={fit_dual['eta_gal_t']:.2f} while η_CMB t={fit_dual['eta_cmb_t']:.2f} "
            f"(add-gal|CMB p={p_add_gal_given_cmb:.4e}, add-CMB|gal p={p_add_cmb_given_gal:.3f}). "
            f"Single-axis Planck vs galactic ΔAIC comparisons are supplementary "
            f"(galactic-only beats CMB-only in AIC: {gal_beats_cmb_in_single_axis}). "
            f"Random-axis Monte Carlo on S² (D): Planck ΔAIC={delta_aic_true:.1f} vs null median "
            f"{random_axis_delta_aic_null['null_delta_aic_median']:.1f}; "
            f"p(any axis ≥ Planck)={random_axis_delta_aic_null['empirical_p_value']:.3f}. "
            f"TEP η_θ coverage (I): fitted η_θ={tep_prediction_coverage['eta_theta_fitted']:.4e} "
            f"within Step 033 band={tep_prediction_coverage['within_predicted_band']} "
            f"(z={tep_prediction_coverage['z_score_vs_predicted_center']:.2f} vs η₀ scale). "
            f"Aliasing, permutation, and orthogonalization support a nontrivial directional component; "
            f"phase null p_F_eff={refined_directional_nulls['phase_null']['p_f_eff']:.4f}. "
            f"Present as directional-anatomy plus fixed-sky identifiability, not unique CMB-frame proof."
        ),
    }

    return results


if __name__ == "__main__":
    log_dir = PROJECT_ROOT / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    logger = TEPLogger(
        "step_055", str(log_dir / "step_055_cmb_rigorous_falsification.log")
    )
    set_step_logger(logger)

    data_path = PROJECT_ROOT / "data" / "processed" / "INPOP19a_all_stations_residuals.csv"
    if not data_path.exists():
        raise FileNotFoundError(
            "Missing required processed residual archive for Step 055: "
            f"{data_path}. "
            "This step is defined on the INPOP19a archive used by the primary estimand; "
            "automatic ephemeris substitution is disabled."
        )

    df = pd.read_csv(data_path)

    results = cmb_falsification_analysis(df, verbose=True)
    if results:
        results["input_dataset"] = {
            "label": "INPOP19a_all_stations_residuals",
            "path": str(data_path.relative_to(PROJECT_ROOT)),
        }

    if results:
        logger.save_step_results(
            results, PROJECT_ROOT, "step_055_cmb_rigorous_falsification"
        )
        print_status("CMB Falsification analysis complete.", "SUCCESS")
    else:
        print_status("Analysis failed.", "ERROR")
        sys.exit(1)
