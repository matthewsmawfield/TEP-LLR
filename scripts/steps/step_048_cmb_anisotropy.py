#!/usr/bin/env python3
"""
Step 048: CMB Dipole Anisotropy Test
====================================

Tests whether residual-channel structure is consistent with anisotropic coupling
defined relative to the Planck 2018 CMB dipole as an operational fixed celestial
axis (velocity projection and Earth–Moon orientation), with Step 055 supplying
adversarial nulls on uniqueness and synodic phase coupling.

Physical Mechanism:
The CMB dipole defines a standard kinematic reference direction in cosmology:
(l, b) = (264.02°, 48.25°) in galactic coordinates, corresponding to
(α, δ) ≈ (168.14°, -7.22°) in J2000 equatorial coordinates, with amplitude
v_CMB ≈ 369 km/s (Planck 2018). In TEP-motivated embeddings that single out a
large-scale rest frame, motion and orientation relative to that axis should
modulate the effective coupling; this step fits those predictors without
asserting a unique cosmological origin (see Step 055).

Two distinct predictions arise:

  1. Annual velocity projection: Earth's orbital velocity (~30 km/s) projects
     onto the CMB dipole direction with amplitude varying sinusoidally over the
     year. The parallel component
         v_parallel(t) = v_orbital(t) · n_CMB
     should modulate η if the coupling is velocity-dependent in the CMB frame.

  2. Monthly orientation anisotropy: The Earth-Moon line sweeps across the
     celestial sphere with synodic period. If the scalar gradient has a
     preferred direction aligned with the CMB dipole, the effective η should
     depend on the cosine of the angle θ between the Earth-Moon vector and
     the CMB dipole:
         η ∝ cos(θ_EM-CMB)
     This creates a monthly-period anisotropy superimposed on the synodic
     signal, distinct from the heliocentric distance/velocity modulation.

Crucially, the CMB dipole direction in ecliptic coordinates is at longitude
≈ 173° and latitude ≈ 11°, meaning the annual velocity projection and the
heliocentric distance modulation (perihelion at longitude ≈ 103°) are offset
by approximately 70° in orbital phase. This phase offset makes them
statistically distinguishable.

This step uses DE421 ephemeris via Skyfield to compute:
  - Earth's orbital velocity vector in the Barycentric frame
  - The Earth-Moon unit vector for every observation
  - Dot products with the CMB dipole direction
  - Tests for correlation with η in both annual and monthly frequencies.
"""

import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

import numpy as np
from scripts.utils.numerics import stable_lstsq
import pandas as pd
from scipy import stats
from skyfield.api import load

from scripts.utils.llr_constants import ETA_SCALE_FACTOR
from scripts.utils.statistical_utils import linear_regression, detect_outliers_sigma
from scripts.utils.logger import TEPLogger, set_step_logger, print_status


# CMB dipole direction (Planck 2018)
# Galactic: l = 264.021°, b = 48.253°
# J2000 equatorial: α = 168.14°, δ = -7.22°
# Unit vector in ICRS/Barycentric frame (km basis, any consistent frame)
_CMB_RA_RAD = np.deg2rad(168.14)
_CMB_DEC_RAD = np.deg2rad(-7.22)
_CMB_UNIT = np.array([
    np.cos(_CMB_DEC_RAD) * np.cos(_CMB_RA_RAD),
    np.cos(_CMB_DEC_RAD) * np.sin(_CMB_RA_RAD),
    np.sin(_CMB_DEC_RAD),
])
# Normalise to machine precision
_CMB_UNIT = _CMB_UNIT / np.linalg.norm(_CMB_UNIT)


def compute_cmb_projections(jd_array):
    """Compute CMB-frame kinematic projections for every epoch.

    Parameters:
    -----------
    jd_array : array-like
        Julian dates (TT scale)

    Returns:
    --------
    dict with keys:
        v_parallel_kms : ndarray
            Earth orbital velocity projected onto CMB dipole [km/s]
        v_perp_kms : ndarray
            Perpendicular component magnitude [km/s]
        earth_moon_cos_theta : ndarray
            Cosine of angle between Earth-Moon vector and CMB dipole
        sun_distance_au : ndarray
            Heliocentric distance [AU]
        speed_kms : ndarray
            Total orbital speed [km/s]
    """
    from scripts.utils.astronomical_utils import load_skyfield_planets

    planets, _eph_path = load_skyfield_planets(PROJECT_ROOT)
    earth = planets["earth"]
    moon = planets["moon"]
    sun = planets["sun"]
    ts = load.timescale()

    timestamps = ts.tt(jd=jd_array)

    # Earth barycentric position and velocity
    earth_pos = earth.at(timestamps)
    earth_pv = earth_pos.position.km  # shape (3, N)
    earth_vv = earth_pos.velocity.km_per_s  # shape (3, N)

    # Sun barycentric position
    sun_pv = sun.at(timestamps).position.km  # shape (3, N)

    # Moon barycentric position
    moon_pv = moon.at(timestamps).position.km  # shape (3, N)

    # Earth-Sun relative vectors (heliocentric Earth)
    rel_pos = earth_pv - sun_pv
    rel_vel = earth_vv - sun.at(timestamps).velocity.km_per_s

    distance_km = np.linalg.norm(rel_pos, axis=0)
    speed_kms = np.linalg.norm(rel_vel, axis=0)

    # Radial and parallel velocity components
    r_hat = rel_pos / distance_km
    v_radial_kms = np.sum(rel_vel * r_hat, axis=0)

    # CMB dipole projection of orbital velocity
    v_parallel_kms = np.sum(rel_vel * _CMB_UNIT[:, None], axis=0)
    v_perp_vec = rel_vel - v_parallel_kms * _CMB_UNIT[:, None]
    v_perp_kms = np.linalg.norm(v_perp_vec, axis=0)

    # Earth-Moon vector (Earth -> Moon)
    em_vec = moon_pv - earth_pv
    em_dist = np.linalg.norm(em_vec, axis=0)
    em_hat = em_vec / em_dist

    # Cosine of angle between Earth-Moon line and CMB dipole
    earth_moon_cos_theta = np.sum(em_hat * _CMB_UNIT[:, None], axis=0)

    return {
        "v_parallel_kms": v_parallel_kms,
        "v_perp_kms": v_perp_kms,
        "earth_moon_cos_theta": earth_moon_cos_theta,
        "sun_distance_au": distance_km / 1.495978707e8,
        "speed_kms": speed_kms,
        "radial_velocity_kms": v_radial_kms,
        "earth_moon_unit_vectors": em_hat,
    }


def cmb_anisotropy_analysis(df, verbose=False):
    """Main analysis: test for CMB-frame anisotropic TEP modulation.

    Parameters:
    -----------
    df : pd.DataFrame
        DataFrame with columns: date_julian, residual_m, elongation_rad
    verbose : bool
        Print detailed output

    Returns:
    --------
    dict : Analysis results
    """
    print_status(
        "═══ Step 048: CMB Dipole Anisotropy Test ═══",
        "TITLE",
    )

    n = len(df)
    residuals = df["residual_m"].values
    cos_elong = np.cos(df["elongation_rad"].values)
    jd = df["date_julian"].values

    # ------------------------------------------------------------------
    # 1. Compute CMB-frame kinematics
    # ------------------------------------------------------------------
    print_status("Computing CMB-frame projections via DE421/Skyfield...", "PROCESS")
    cmb_data = compute_cmb_projections(jd)

    v_par = cmb_data["v_parallel_kms"]
    v_perp = cmb_data["v_perp_kms"]
    cos_theta = cmb_data["earth_moon_cos_theta"]
    r = cmb_data["sun_distance_au"]
    speed = cmb_data["speed_kms"]
    vr = cmb_data["radial_velocity_kms"]
    em_hat = cmb_data["earth_moon_unit_vectors"]

    print_status(
        f"v_parallel range: {v_par.min():.3f} – {v_par.max():.3f} km/s",
        "CALC",
    )
    print_status(
        f"v_perp range: {v_perp.min():.3f} – {v_perp.max():.3f} km/s",
        "CALC",
    )
    print_status(
        f"cos(θ_EM-CMB) range: {cos_theta.min():.4f} – {cos_theta.max():.4f}",
        "CALC",
    )

    # Correlation: v_parallel with heliocentric distance (should be weak
    # if CMB direction is not aligned with perihelion/aphelion axis)
    corr_vpar_r = np.corrcoef(v_par, r)[0, 1]
    print_status(
        f"Correlation r(v_parallel, r) = {corr_vpar_r:.4f}",
        "CALC",
    )

    # Correlation: cos_theta with cos(elongation) (should be weak —
    # monthly vs synodic are different frequencies)
    corr_cosD_cosTheta = np.corrcoef(cos_elong, cos_theta)[0, 1]
    print_status(
        f"Correlation r(cos(D), cosθ_EM-CMB) = {corr_cosD_cosTheta:.4f}",
        "CALC",
    )

    # Pre-filter outliers
    outlier_mask = detect_outliers_sigma(residuals, sigma_threshold=6.0)
    mask_clean = ~outlier_mask
    n_clean = mask_clean.sum()
    print_status(
        f"Cleaned dataset: {n_clean:,} observations "
        f"({n - n_clean} outliers removed)",
        "CALC",
    )

    res = residuals[mask_clean]
    cosD = cos_elong[mask_clean]
    v_par_c = v_par[mask_clean] - np.mean(v_par[mask_clean])
    cos_theta_c = cos_theta[mask_clean] - np.mean(cos_theta[mask_clean])
    r_c = r[mask_clean] - np.mean(r[mask_clean])
    vr_c = vr[mask_clean] - np.mean(vr[mask_clean])

    # ------------------------------------------------------------------
    # 2. Model A: Annual CMB velocity projection
    # ------------------------------------------------------------------
    print_status("═══ MODEL A: Annual CMB velocity projection ═══", "TITLE")

    # Split by v_parallel sign and magnitude
    vp_p15 = np.percentile(v_par_c, 15)
    vp_p85 = np.percentile(v_par_c, 85)
    mask_vp_high = v_par_c >= vp_p85
    mask_vp_low = v_par_c <= vp_p15

    reg_vp_high = linear_regression(res[mask_vp_high], cosD[mask_vp_high])
    reg_vp_low = linear_regression(res[mask_vp_low], cosD[mask_vp_low])

    diff_vp = reg_vp_high["eta"] - reg_vp_low["eta"]
    diff_vp_err = np.sqrt(
        reg_vp_high["eta_error"] ** 2 + reg_vp_low["eta_error"] ** 2
    )
    diff_vp_sig = abs(diff_vp) / diff_vp_err if diff_vp_err > 0 else 0.0

    print_status(
        f"High v_parallel (> {vp_p85:.3f}): η = {reg_vp_high['eta']:.4e} "
        f"± {reg_vp_high['eta_error']:.4e}",
        "CALC",
    )
    print_status(
        f"Low v_parallel (< {vp_p15:.3f}):  η = {reg_vp_low['eta']:.4e} "
        f"± {reg_vp_low['eta_error']:.4e}",
        "CALC",
    )
    print_status(
        f"Δη = {diff_vp:.4e} ± {diff_vp_err:.4e} ({diff_vp_sig:.2f}σ)",
        "CALC",
    )

    # ------------------------------------------------------------------
    # 3. Model B: Monthly Earth-Moon orientation anisotropy
    # ------------------------------------------------------------------
    print_status(
        "═══ MODEL B: Monthly Earth-Moon orientation anisotropy ═══",
        "TITLE",
    )

    # Test: does η correlate with cos(θ_EM-CMB)?
    # If the scalar field has a preferred direction, the amplitude should
    # depend on how aligned the Earth-Moon test mass dipole is with that
    # direction.
    ct_p85 = np.percentile(cos_theta_c, 85)
    ct_p15 = np.percentile(cos_theta_c, 15)
    mask_ct_high = cos_theta_c >= ct_p85
    mask_ct_low = cos_theta_c <= ct_p15

    reg_ct_high = linear_regression(res[mask_ct_high], cosD[mask_ct_high])
    reg_ct_low = linear_regression(res[mask_ct_low], cosD[mask_ct_low])

    diff_ct = reg_ct_high["eta"] - reg_ct_low["eta"]
    diff_ct_err = np.sqrt(
        reg_ct_high["eta_error"] ** 2 + reg_ct_low["eta_error"] ** 2
    )
    diff_ct_sig = abs(diff_ct) / diff_ct_err if diff_ct_err > 0 else 0.0

    print_status(
        f"High cosθ (> {ct_p85:.3f}): η = {reg_ct_high['eta']:.4e} "
        f"± {reg_ct_high['eta_error']:.4e}",
        "CALC",
    )
    print_status(
        f"Low cosθ (< {ct_p15:.3f}):  η = {reg_ct_low['eta']:.4e} "
        f"± {reg_ct_low['eta_error']:.4e}",
        "CALC",
    )
    print_status(
        f"Δη = {diff_ct:.4e} ± {diff_ct_err:.4e} ({diff_ct_sig:.2f}σ)",
        "CALC",
    )

    # ------------------------------------------------------------------
    # 4. Model C: Joint CMB anisotropy fit
    # ------------------------------------------------------------------
    print_status("═══ MODEL C: Joint CMB anisotropy fit ═══", "TITLE")

    # Base model: η = η0 + η_syn cos(D) + η_vpar v_parallel cos(D) + η_theta cos(θ) cos(D)
    # In the residual formulation: residual = β0 cos(D) + β1 v_parallel cos(D)
    #                                     + β2 cos(θ) cos(D) + β3 + ε
    X_joint = np.column_stack([
        cosD,
        v_par_c * cosD,
        cos_theta_c * cosD,
        np.ones(n_clean),
    ])
    with np.errstate(divide='ignore', invalid='ignore', over='ignore'):
        coeffs_joint, _, rank_joint, _ = stable_lstsq(X_joint, res)

    if rank_joint < 4:
        print_status("WARNING: Joint model rank-deficient", "WARNING")
        joint_result = {"available": False, "reason": "rank_deficient"}
    else:
        with np.errstate(divide='ignore', invalid='ignore', over='ignore'):
            resid_joint = res - X_joint @ coeffs_joint
        mse_joint = np.sum(resid_joint ** 2) / (n_clean - 4)
        XtX_inv = np.linalg.pinv(X_joint.T @ X_joint, rcond=1e-10, hermitian=True)
        cov_joint = mse_joint * XtX_inv
        se_joint = np.sqrt(np.diag(cov_joint))

        eta_0 = coeffs_joint[0] / ETA_SCALE_FACTOR
        eta_vpar = coeffs_joint[1] / ETA_SCALE_FACTOR
        eta_theta = coeffs_joint[2] / ETA_SCALE_FACTOR

        se_eta_0 = se_joint[0] / ETA_SCALE_FACTOR
        se_eta_vpar = se_joint[1] / ETA_SCALE_FACTOR
        se_eta_theta = se_joint[2] / ETA_SCALE_FACTOR

        t_vpar = eta_vpar / se_eta_vpar if se_eta_vpar > 0 else 0.0
        t_theta = eta_theta / se_eta_theta if se_eta_theta > 0 else 0.0
        p_vpar = 2 * (1 - stats.t.cdf(abs(t_vpar), n_clean - 4)) if se_eta_vpar > 0 else 1.0
        p_theta = 2 * (1 - stats.t.cdf(abs(t_theta), n_clean - 4)) if se_eta_theta > 0 else 1.0

        rss_joint = np.sum(resid_joint ** 2)
        k_joint = 4
        aic_joint = n_clean * np.log(rss_joint / n_clean) + 2 * k_joint
        bic_joint = n_clean * np.log(rss_joint / n_clean) + k_joint * np.log(n_clean)

        # Simple model (cosD only, 2 params)
        X_simple = np.column_stack([cosD, np.ones(n_clean)])
        with np.errstate(divide='ignore', invalid='ignore', over='ignore'):
            coeffs_simple, _, _, _ = stable_lstsq(X_simple, res)
            resid_simple = res - X_simple @ coeffs_simple
        rss_simple = np.sum(resid_simple ** 2)
        k_simple = 2
        aic_simple = n_clean * np.log(rss_simple / n_clean) + 2 * k_simple
        bic_simple = n_clean * np.log(rss_simple / n_clean) + k_simple * np.log(n_clean)

        # Distance-velocity model (from Step 047, 4 params)
        X_dv = np.column_stack([cosD, r_c * cosD, vr_c * cosD, np.ones(n_clean)])
        with np.errstate(divide='ignore', invalid='ignore', over='ignore'):
            coeffs_dv, _, _, _ = stable_lstsq(X_dv, res)
            resid_dv = res - X_dv @ coeffs_dv
        rss_dv = np.sum(resid_dv ** 2)
        k_dv = 4
        aic_dv = n_clean * np.log(rss_dv / n_clean) + 2 * k_dv
        bic_dv = n_clean * np.log(rss_dv / n_clean) + k_dv * np.log(n_clean)

        # CMB-only model (3 params)
        X_cmb = np.column_stack([cosD, v_par_c * cosD, np.ones(n_clean)])
        with np.errstate(divide='ignore', invalid='ignore', over='ignore'):
            coeffs_cmb, _, _, _ = stable_lstsq(X_cmb, res)
            resid_cmb = res - X_cmb @ coeffs_cmb
        rss_cmb = np.sum(resid_cmb ** 2)
        k_cmb = 3
        aic_cmb = n_clean * np.log(rss_cmb / n_clean) + 2 * k_cmb
        bic_cmb = n_clean * np.log(rss_cmb / n_clean) + k_cmb * np.log(n_clean)

        print_status(
            f"Joint CMB: η = {eta_0:.4e} + {eta_vpar:.4e}·v_parallel "
            f"+ {eta_theta:.4e}·cos(θ)",
            "CALC",
        )
        print_status(
            f"  η_vpar: {eta_vpar:.4e} ± {se_eta_vpar:.4e} "
            f"(t={t_vpar:.2f}, p={p_vpar:.4f})",
            "CALC",
        )
        print_status(
            f"  η_theta: {eta_theta:.4e} ± {se_eta_theta:.4e} "
            f"(t={t_theta:.2f}, p={p_theta:.4f})",
            "CALC",
        )
        print_status(
            f"  AIC simple={aic_simple:.1f}, CMB-only={aic_cmb:.1f}, "
            f"dist-vel={aic_dv:.1f}, joint-CMB={aic_joint:.1f}",
            "CALC",
        )

        delta_aic_cmb = aic_cmb - aic_simple
        delta_aic_joint = aic_joint - aic_simple
        delta_aic_dv = aic_dv - aic_simple

        print_status(
            f"  ΔAIC(CMB-only vs simple) = {delta_aic_cmb:.1f}",
            "CALC",
        )
        print_status(
            f"  ΔAIC(dist-vel vs simple) = {delta_aic_dv:.1f}",
            "CALC",
        )
        print_status(
            f"  ΔAIC(joint-CMB vs simple) = {delta_aic_joint:.1f}",
            "CALC",
        )

        if delta_aic_joint < -2:
            model_pref = "joint_cmb_preferred"
            print_status("Joint CMB model preferred (ΔAIC < -2)", "SUCCESS")
        elif delta_aic_cmb < -2:
            model_pref = "cmb_only_preferred"
            print_status("CMB-only model preferred (ΔAIC < -2)", "SUCCESS")
        elif delta_aic_dv < -2:
            model_pref = "dist_vel_preferred"
            print_status("Dist-vel model preferred; no CMB signal", "WARNING")
        else:
            model_pref = "simple_preferred"
            print_status("No CMB anisotropy detected", "WARNING")

        joint_result = {
            "available": True,
            "eta_0": float(eta_0),
            "eta_0_error": float(se_eta_0),
            "eta_vpar": float(eta_vpar),
            "eta_vpar_error": float(se_eta_vpar),
            "eta_vpar_t": float(t_vpar),
            "eta_vpar_p": float(p_vpar),
            "eta_theta": float(eta_theta),
            "eta_theta_error": float(se_eta_theta),
            "eta_theta_t": float(t_theta),
            "eta_theta_p": float(p_theta),
            "rss": float(rss_joint),
            "aic": float(aic_joint),
            "bic": float(bic_joint),
            "aic_simple": float(aic_simple),
            "aic_cmb_only": float(aic_cmb),
            "aic_dist_vel": float(aic_dv),
            "delta_aic_cmb_vs_simple": float(delta_aic_cmb),
            "delta_aic_joint_vs_simple": float(delta_aic_joint),
            "delta_aic_dist_vel_vs_simple": float(delta_aic_dv),
            "model_preference": model_pref,
        }

    # ------------------------------------------------------------------
    # 4a. Refinement A: Orthogonalize v_parallel against distance
    # ------------------------------------------------------------------
    print_status("═══ REFINEMENT A: Orthogonalized v_parallel ═══", "TITLE")

    # Regress v_parallel on distance, use the residual as the CMB-frame
    # velocity predictor that is genuinely independent from heliocentric
    # distance modulation.
    X_vr = np.column_stack([np.ones(n_clean), r_c])
    with np.errstate(divide='ignore', invalid='ignore', over='ignore'):
        beta_vr, _, rank_vr, _ = stable_lstsq(X_vr, v_par_c)
    if rank_vr == 2:
        with np.errstate(divide='ignore', invalid='ignore', over='ignore'):
            v_par_perp = v_par_c - X_vr @ beta_vr
        corr_vpar_perp_r = np.corrcoef(v_par_perp, r_c)[0, 1]
        print_status(
            f"Residual v_parallel (after removing r): "
            f"r(v_perp, r) = {corr_vpar_perp_r:.4f}",
            "CALC",
        )

        # Joint fit with orthogonalized velocity
        X_orth = np.column_stack([
            cosD,
            v_par_perp * cosD,
            cos_theta_c * cosD,
            np.ones(n_clean),
        ])
        with np.errstate(divide='ignore', invalid='ignore', over='ignore'):
            coeffs_orth, _, rank_orth, _ = stable_lstsq(X_orth, res)

        orth_result = {"available": False}
        if rank_orth == 4:
            with np.errstate(divide='ignore', invalid='ignore', over='ignore'):
                resid_orth = res - X_orth @ coeffs_orth
            mse_orth = np.sum(resid_orth ** 2) / (n_clean - 4)
            XtX_orth_inv = np.linalg.pinv(X_orth.T @ X_orth, rcond=1e-10, hermitian=True)
            cov_orth = mse_orth * XtX_orth_inv
            se_orth = np.sqrt(np.diag(cov_orth))

            eta_0_orth = coeffs_orth[0] / ETA_SCALE_FACTOR
            eta_vpar_perp = coeffs_orth[1] / ETA_SCALE_FACTOR
            eta_theta_orth = coeffs_orth[2] / ETA_SCALE_FACTOR

            se_eta_vpar_perp = se_orth[1] / ETA_SCALE_FACTOR
            se_eta_theta_orth = se_orth[2] / ETA_SCALE_FACTOR

            t_vpar_perp = (eta_vpar_perp / se_eta_vpar_perp
                           if se_eta_vpar_perp > 0 else 0.0)
            t_theta_orth = (eta_theta_orth / se_eta_theta_orth
                            if se_eta_theta_orth > 0 else 0.0)
            p_vpar_perp = (2 * (1 - stats.t.cdf(abs(t_vpar_perp), n_clean - 4))
                           if se_eta_vpar_perp > 0 else 1.0)
            p_theta_orth = (2 * (1 - stats.t.cdf(abs(t_theta_orth), n_clean - 4))
                            if se_eta_theta_orth > 0 else 1.0)

            rss_orth = np.sum(resid_orth ** 2)
            aic_orth = n_clean * np.log(rss_orth / n_clean) + 2 * 4

            print_status(
                f"Orthogonalized: η_vpar_perp = {eta_vpar_perp:.4e} "
                f"± {se_eta_vpar_perp:.4e} (t={t_vpar_perp:.2f}, p={p_vpar_perp:.4f})",
                "CALC",
            )
            print_status(
                f"Orthogonalized: η_theta     = {eta_theta_orth:.4e} "
                f"± {se_eta_theta_orth:.4e} (t={t_theta_orth:.2f}, p={p_theta_orth:.4f})",
                "CALC",
            )

            orth_result = {
                "available": True,
                "eta_0": float(eta_0_orth),
                "eta_0_error": float(se_orth[0] / ETA_SCALE_FACTOR),
                "eta_vpar_perp": float(eta_vpar_perp),
                "eta_vpar_perp_error": float(se_eta_vpar_perp),
                "eta_vpar_perp_t": float(t_vpar_perp),
                "eta_vpar_perp_p": float(p_vpar_perp),
                "eta_theta": float(eta_theta_orth),
                "eta_theta_error": float(se_eta_theta_orth),
                "eta_theta_t": float(t_theta_orth),
                "eta_theta_p": float(p_theta_orth),
                "corr_vpar_perp_vs_r": float(corr_vpar_perp_r),
                "rss": float(rss_orth),
                "aic": float(aic_orth),
                "delta_aic_vs_simple": float(aic_orth - aic_simple),
            }
    else:
        orth_result = {"available": False, "reason": "rank_deficient"}
        print_status("Orthogonalization failed (rank deficient)", "WARNING")

    # ------------------------------------------------------------------
    # 4b. Refinement B: Full joint regression (synodic + distance + v_r + cosθ)
    # ------------------------------------------------------------------
    print_status(
        "═══ REFINEMENT B: Full joint (synodic + distance + v_r + cosθ) ═══",
        "TITLE",
    )

    X_full = np.column_stack([
        cosD,
        r_c * cosD,
        vr_c * cosD,
        cos_theta_c * cosD,
        np.ones(n_clean),
    ])
    with np.errstate(divide='ignore', invalid='ignore', over='ignore'):
        coeffs_full, _, rank_full, _ = stable_lstsq(X_full, res)

    full_result = {"available": False}
    if rank_full == 5:
        with np.errstate(divide='ignore', invalid='ignore', over='ignore'):
            resid_full = res - X_full @ coeffs_full
        mse_full = np.sum(resid_full ** 2) / (n_clean - 5)
        XtX_full_inv = np.linalg.pinv(X_full.T @ X_full, rcond=1e-10, hermitian=True)
        cov_full = mse_full * XtX_full_inv
        se_full = np.sqrt(np.diag(cov_full))

        eta_0_full = coeffs_full[0] / ETA_SCALE_FACTOR
        eta_r_full = coeffs_full[1] / ETA_SCALE_FACTOR
        eta_vr_full = coeffs_full[2] / ETA_SCALE_FACTOR
        eta_theta_full = coeffs_full[3] / ETA_SCALE_FACTOR

        se_eta_r_full = se_full[1] / ETA_SCALE_FACTOR
        se_eta_vr_full = se_full[2] / ETA_SCALE_FACTOR
        se_eta_theta_full = se_full[3] / ETA_SCALE_FACTOR

        t_r_full = eta_r_full / se_eta_r_full if se_eta_r_full > 0 else 0.0
        t_vr_full = eta_vr_full / se_eta_vr_full if se_eta_vr_full > 0 else 0.0
        t_theta_full = (eta_theta_full / se_eta_theta_full
                        if se_eta_theta_full > 0 else 0.0)
        p_r_full = (2 * (1 - stats.t.cdf(abs(t_r_full), n_clean - 5))
                    if se_eta_r_full > 0 else 1.0)
        p_vr_full = (2 * (1 - stats.t.cdf(abs(t_vr_full), n_clean - 5))
                     if se_eta_vr_full > 0 else 1.0)
        p_theta_full = (2 * (1 - stats.t.cdf(abs(t_theta_full), n_clean - 5))
                        if se_eta_theta_full > 0 else 1.0)

        rss_full = np.sum(resid_full ** 2)
        aic_full = n_clean * np.log(rss_full / n_clean) + 2 * 5

        print_status(
            f"Full joint: η_r = {eta_r_full:.4e} ± {se_eta_r_full:.4e} "
            f"(t={t_r_full:.2f}, p={p_r_full:.4f})",
            "CALC",
        )
        print_status(
            f"Full joint: η_vr = {eta_vr_full:.4e} ± {se_eta_vr_full:.4e} "
            f"(t={t_vr_full:.2f}, p={p_vr_full:.4f})",
            "CALC",
        )
        print_status(
            f"Full joint: η_theta = {eta_theta_full:.4e} "
            f"± {se_eta_theta_full:.4e} (t={t_theta_full:.2f}, p={p_theta_full:.4f})",
            "CALC",
        )

        full_result = {
            "available": True,
            "eta_0": float(eta_0_full),
            "eta_0_error": float(se_full[0] / ETA_SCALE_FACTOR),
            "eta_r": float(eta_r_full),
            "eta_r_error": float(se_eta_r_full),
            "eta_r_t": float(t_r_full),
            "eta_r_p": float(p_r_full),
            "eta_vr": float(eta_vr_full),
            "eta_vr_error": float(se_eta_vr_full),
            "eta_vr_t": float(t_vr_full),
            "eta_vr_p": float(p_vr_full),
            "eta_theta": float(eta_theta_full),
            "eta_theta_error": float(se_eta_theta_full),
            "eta_theta_t": float(t_theta_full),
            "eta_theta_p": float(p_theta_full),
            "rss": float(rss_full),
            "aic": float(aic_full),
            "delta_aic_vs_simple": float(aic_full - aic_simple),
            "delta_aic_vs_dist_vel": float(aic_full - aic_dv),
        }
    else:
        print_status("Full joint model rank-deficient", "WARNING")

    # ------------------------------------------------------------------
    # 4c. Refinement C: Nested model comparison for cosθ marginal contribution
    # ------------------------------------------------------------------
    print_status(
        "═══ REFINEMENT C: cosθ marginal contribution (nested models) ═══",
        "TITLE",
    )

    # Model (a): synodic + distance + radial velocity (Step 047 model)
    X_base = np.column_stack([
        cosD,
        r_c * cosD,
        vr_c * cosD,
        np.ones(n_clean),
    ])
    with np.errstate(divide='ignore', invalid='ignore', over='ignore'):
        coeffs_base, _, rank_base, _ = stable_lstsq(X_base, res)

    nested_result = {"available": False}
    if rank_base == 4:
        with np.errstate(divide='ignore', invalid='ignore', over='ignore'):
            resid_base = res - X_base @ coeffs_base
        rss_base = np.sum(resid_base ** 2)
        k_base = 4
        aic_base = n_clean * np.log(rss_base / n_clean) + 2 * k_base

        # Model (b): add cosθ
        X_cos = np.column_stack([
            cosD,
            r_c * cosD,
            vr_c * cosD,
            cos_theta_c * cosD,
            np.ones(n_clean),
        ])
        with np.errstate(divide='ignore', invalid='ignore', over='ignore'):
            coeffs_cos, _, rank_cos, _ = stable_lstsq(X_cos, res)

        if rank_cos == 5:
            with np.errstate(divide='ignore', invalid='ignore', over='ignore'):
                resid_cos = res - X_cos @ coeffs_cos
            rss_cos = np.sum(resid_cos ** 2)
            k_cos = 5
            aic_cos = n_clean * np.log(rss_cos / n_clean) + 2 * k_cos

            delta_aic_cos_base = aic_cos - aic_base
            delta_rss = rss_base - rss_cos
            f_stat = (delta_rss / 1) / (rss_cos / (n_clean - 5))
            p_f = 1 - stats.f.cdf(f_stat, 1, n_clean - 5)

            # Partial eta for cosθ in the nested model
            eta_theta_nested = coeffs_cos[3] / ETA_SCALE_FACTOR
            se_theta_nested = np.sqrt(
                np.sum(resid_cos ** 2) / (n_clean - 5)
                * np.linalg.pinv(X_cos.T @ X_cos, rcond=1e-10, hermitian=True)[3, 3]
            ) / ETA_SCALE_FACTOR
            t_theta_nested = (eta_theta_nested / se_theta_nested
                              if se_theta_nested > 0 else 0.0)
            p_theta_nested = (2 * (1 - stats.t.cdf(
                abs(t_theta_nested), n_clean - 5))
                if se_theta_nested > 0 else 1.0)

            print_status(
                f"Base model (synodic + dist + v_r): AIC = {aic_base:.1f}",
                "CALC",
            )
            print_status(
                f"+cosθ model: AIC = {aic_cos:.1f}, "
                f"ΔAIC = {delta_aic_cos_base:.1f}",
                "CALC",
            )
            print_status(
                f"F(1, {n_clean - 5}) = {f_stat:.2f}, p = {p_f:.4f}",
                "CALC",
            )
            print_status(
                f"cosθ partial η = {eta_theta_nested:.4e} "
                f"± {se_theta_nested:.4e} (t={t_theta_nested:.2f}, "
                f"p={p_theta_nested:.4f})",
                "CALC",
            )

            nested_result = {
                "available": True,
                "aic_base": float(aic_base),
                "aic_with_cos_theta": float(aic_cos),
                "delta_aic_cos_vs_base": float(delta_aic_cos_base),
                "f_statistic": float(f_stat),
                "f_p_value": float(p_f),
                "eta_theta_partial": float(eta_theta_nested),
                "eta_theta_partial_error": float(se_theta_nested),
                "eta_theta_partial_t": float(t_theta_nested),
                "eta_theta_partial_p": float(p_theta_nested),
            }

    # ------------------------------------------------------------------
    # 5. Phase-locked test: compare orbital longitude of perihelion
    #    vs CMB dipole direction
    # ------------------------------------------------------------------
    print_status("═══ Phase offset test ═══", "TITLE")

    # Perihelion longitude ≈ 103° (ecliptic), CMB dipole longitude ≈ 173°
    # The 70° offset means v_parallel peaks ~70 days after perihelion
    # We test whether a 70° phase-shifted signal exists in the residuals

    # Mean anomaly M ≈ 2π(t - T_perihelion)/365.25
    # T_perihelion ≈ JD 2451545.0 + 4 days (Jan 4, 2000)
    # Test for annual power at the CMB dipole phase by including both
    # sin and cos at the annual frequency.  The joint F-test is
    # phase-independent and properly detects any annual signal regardless
    # of whether it is in-phase or quadrature with the 70° offset.
    t_days = jd - 2451545.0
    annual_sin_cmb = np.sin(2 * np.pi * t_days / 365.25 - np.deg2rad(70))
    annual_cos_cmb = np.cos(2 * np.pi * t_days / 365.25 - np.deg2rad(70))

    annual_sin_cmb_c = annual_sin_cmb[mask_clean] - np.mean(annual_sin_cmb[mask_clean])
    annual_cos_cmb_c = annual_cos_cmb[mask_clean] - np.mean(annual_cos_cmb[mask_clean])

    X_phase = np.column_stack([
        cosD,
        annual_sin_cmb_c * cosD,
        annual_cos_cmb_c * cosD,
        np.ones(n_clean),
    ])
    with np.errstate(divide='ignore', invalid='ignore', over='ignore'):
        coeffs_phase, _, rank_phase, _ = stable_lstsq(X_phase, res)

    phase_result = None
    if rank_phase == 4:
        with np.errstate(divide='ignore', invalid='ignore', over='ignore'):
            resid_phase = res - X_phase @ coeffs_phase
        rss_phase = np.sum(resid_phase ** 2)
        mse_phase = rss_phase / (n_clean - 4)
        cov_phase = mse_phase * np.linalg.pinv(X_phase.T @ X_phase, rcond=1e-10, hermitian=True)
        se_phase = np.sqrt(np.diag(cov_phase))

        # Individual coefficients (for reference; joint F-test is primary)
        eta_sin = coeffs_phase[1] / ETA_SCALE_FACTOR
        eta_cos = coeffs_phase[2] / ETA_SCALE_FACTOR
        se_eta_sin = se_phase[1] / ETA_SCALE_FACTOR
        se_eta_cos = se_phase[2] / ETA_SCALE_FACTOR
        t_sin = eta_sin / se_eta_sin if se_eta_sin > 0 else 0.0
        t_cos = eta_cos / se_eta_cos if se_eta_cos > 0 else 0.0
        p_sin = 2 * (1 - stats.t.cdf(abs(t_sin), n_clean - 4)) if se_eta_sin > 0 else 1.0
        p_cos = 2 * (1 - stats.t.cdf(abs(t_cos), n_clean - 4)) if se_eta_cos > 0 else 1.0

        # Joint F-test: does the sin+cos pair improve over synodic-only?
        X_phase_null = np.column_stack([cosD, np.ones(n_clean)])
        with np.errstate(divide='ignore', invalid='ignore', over='ignore'):
            coeffs_phase_null, _, _, _ = stable_lstsq(X_phase_null, res)
            resid_phase_null = res - X_phase_null @ coeffs_phase_null
        rss_phase_null = np.sum(resid_phase_null ** 2)
        delta_rss_phase = rss_phase_null - rss_phase
        f_phase = (delta_rss_phase / 2) / (rss_phase / (n_clean - 4))
        p_f_phase = 1 - stats.f.cdf(f_phase, 2, n_clean - 4)

        phase_result = {
            "eta_sin": float(eta_sin),
            "eta_sin_error": float(se_eta_sin),
            "eta_sin_t": float(t_sin),
            "eta_sin_p": float(p_sin),
            "eta_cos": float(eta_cos),
            "eta_cos_error": float(se_eta_cos),
            "eta_cos_t": float(t_cos),
            "eta_cos_p": float(p_cos),
            "f_statistic": float(f_phase),
            "f_p_value": float(p_f_phase),
            "phase_shift_deg": 70.0,
        }

        print_status(
            f"CMB-phase annual sin: η = {eta_sin:.4e} "
            f"± {se_eta_sin:.4e} (t={t_sin:.2f}, p={p_sin:.4f})",
            "CALC",
        )
        print_status(
            f"CMB-phase annual cos: η = {eta_cos:.4e} "
            f"± {se_eta_cos:.4e} (t={t_cos:.2f}, p={p_cos:.4f})",
            "CALC",
        )
        print_status(
            f"Joint F(2, {n_clean - 4}) = {f_phase:.2f}, p = {p_f_phase:.4f}",
            "CALC",
        )

    # ------------------------------------------------------------------
    # 6. Binned anisotropy analysis
    # ------------------------------------------------------------------
    print_status("═══ Binned anisotropy analysis ═══", "TITLE")

    n_bins = 8
    bin_edges = np.linspace(-1, 1, n_bins + 1)
    bin_centers = []
    bin_etas = []
    bin_eta_errs = []

    for i in range(n_bins):
        mask = (cos_theta_c >= bin_edges[i]) & (cos_theta_c < bin_edges[i + 1])
        if mask.sum() > 100:
            reg_bin = linear_regression(res[mask], cosD[mask])
            bin_centers.append((bin_edges[i] + bin_edges[i + 1]) / 2)
            bin_etas.append(reg_bin["eta"])
            bin_eta_errs.append(reg_bin["eta_error"])

    # Fit linear trend in η vs cos(θ)
    orient_trend = None
    if len(bin_centers) >= 4:
        X_trend = np.column_stack([np.ones(len(bin_centers)), bin_centers])
        coeffs_trend, _, rank_trend, _ = stable_lstsq(
            X_trend, bin_etas)
        if rank_trend == 2:
            resid_trend = np.array(bin_etas) - X_trend @ coeffs_trend
            mse_t = np.sum(resid_trend ** 2) / (len(bin_centers) - 2)
            cov_t = mse_t * np.linalg.pinv(X_trend.T @ X_trend, rcond=1e-10, hermitian=True)
            se_t = np.sqrt(np.diag(cov_t))
            slope = coeffs_trend[1]
            slope_err = se_t[1]
            t_slope = slope / slope_err if slope_err > 0 else 0.0
            p_slope = 2 * (1 - stats.t.cdf(abs(t_slope), len(bin_centers) - 2))

            orient_trend = {
                "slope_eta_per_cos_theta": float(slope),
                "slope_error": float(slope_err),
                "t_statistic": float(t_slope),
                "p_value": float(p_slope),
                "n_bins": len(bin_centers),
            }

            print_status(
                f"Binned trend: dη/dcos(θ) = {slope:.4e} ± {slope_err:.4e} "
                f"(t={t_slope:.2f}, p={p_slope:.4f})",
                "CALC",
            )

    # ------------------------------------------------------------------
    # 7. Directional specificity test: is the anisotropy genuinely
    #    aligned with the CMB dipole, or would any sky direction work?
    # ------------------------------------------------------------------
    print_status("═══ Directional specificity test ═══", "TITLE")

    def direction_unit_vector(ra_deg, dec_deg):
        ra = np.deg2rad(ra_deg)
        dec = np.deg2rad(dec_deg)
        vec = np.array([
            np.cos(dec) * np.cos(ra),
            np.cos(dec) * np.sin(ra),
            np.sin(dec),
        ])
        return vec / np.linalg.norm(vec)

    # Test three geometrically well-defined null directions:
    # 1. True antipode (flip both RA and Dec)
    # 2. Perpendicular in equatorial plane (cross with z-axis)
    # 3. Perpendicular in meridian plane (cross CMB with perp-1)
    z_axis = np.array([0.0, 0.0, 1.0])
    perp1 = np.cross(z_axis, _CMB_UNIT)
    perp1 = perp1 / np.linalg.norm(perp1)
    perp2 = np.cross(_CMB_UNIT, perp1)
    perp2 = perp2 / np.linalg.norm(perp2)

    null_directions = [
        ("Perpendicular-1", perp1),
        ("Anti-CMB", -_CMB_UNIT),
        ("Perpendicular-2", perp2),
    ]

    direction_results = []
    true_sig = diff_ct_sig

    for name, dir_vec in null_directions:
        cos_theta_null = np.sum(em_hat * dir_vec[:, None], axis=0)
        cos_theta_null_c = cos_theta_null[mask_clean] - np.mean(cos_theta_null[mask_clean])

        ct_p85_null = np.percentile(cos_theta_null_c, 85)
        ct_p15_null = np.percentile(cos_theta_null_c, 15)
        mask_high_null = cos_theta_null_c >= ct_p85_null
        mask_low_null = cos_theta_null_c <= ct_p15_null

        reg_high_null = linear_regression(res[mask_high_null], cosD[mask_high_null])
        reg_low_null = linear_regression(res[mask_low_null], cosD[mask_low_null])

        diff_null = reg_high_null["eta"] - reg_low_null["eta"]
        diff_err_null = np.sqrt(
            reg_high_null["eta_error"] ** 2 + reg_low_null["eta_error"] ** 2
        )
        sig_null = diff_null / diff_err_null if diff_err_null > 0 else 0.0

        direction_results.append({
            "direction": name,
            "delta_eta": float(diff_null),
            "delta_eta_error": float(diff_err_null),
            "significance_sigma": float(sig_null),
        })
        print_status(
            f"{name}: Δη = {diff_null:.4e} ± {diff_err_null:.4e} "
            f"({abs(sig_null):.2f}σ)",
            "CALC",
        )

    n_null_significant = sum(1 for d in direction_results if abs(d["significance_sigma"]) > 2.0)

    # Anti-parallel direction (180° in RA) should show reversed sign, same magnitude
    anti_cmb = direction_results[1]
    dipole_consistency = abs(anti_cmb["delta_eta"]) / abs(diff_ct)
    print_status(
        f"Anti-CMB |Δη| / True CMB |Δη| = {dipole_consistency:.3f} (expect ≈1.0 for dipole)",
        "CALC",
    )

    # Perpendicular directions should be suppressed relative to CMB
    perp_max_sig = max(abs(direction_results[0]["significance_sigma"]),
                       abs(direction_results[2]["significance_sigma"]))
    perp_suppression = perp_max_sig / true_sig
    print_status(
        f"Max perpendicular σ / True CMB σ = {perp_suppression:.3f} "
        f"(expect ≈0 for pure dipole)",
        "CALC",
    )
    print_status(
        f"True CMB: {true_sig:.2f}σ; Null directions >2σ: {n_null_significant}/3",
        "CALC",
    )

    # ------------------------------------------------------------------
    # 7b. Joint-model null-direction control
    # ------------------------------------------------------------------
    print_status("═══ Joint-model null-direction control ═══", "TITLE")

    def run_joint_for_direction(dir_vec, label):
        """Run the full joint regression for an arbitrary dipole direction.

        Model: residual = eta_0 cosD + eta_r r_c cosD + eta_vr vr_c cosD
                         + eta_theta cos_theta_dir cosD + intercept + epsilon
        """
        cos_theta_dir = np.sum(em_hat * dir_vec[:, None], axis=0)
        cos_theta_dir_c = cos_theta_dir[mask_clean] - np.mean(cos_theta_dir[mask_clean])

        X_dir = np.column_stack([
            cosD, r_c * cosD, vr_c * cosD,
            cos_theta_dir_c * cosD, np.ones(n_clean),
        ])
        with np.errstate(divide='ignore', invalid='ignore', over='ignore'):
            coeffs_dir, _, rank_dir, _ = stable_lstsq(X_dir, res)

        if rank_dir != 5:
            return {"available": False, "reason": "rank_deficient"}

        with np.errstate(divide='ignore', invalid='ignore', over='ignore'):
            resid_dir = res - X_dir @ coeffs_dir
        rss_dir = np.sum(resid_dir ** 2)
        mse_dir = rss_dir / (n_clean - 5)
        cov_dir = mse_dir * np.linalg.pinv(X_dir.T @ X_dir, rcond=1e-10, hermitian=True)
        se_dir = np.sqrt(np.diag(cov_dir))

        eta_0_dir = coeffs_dir[0] / ETA_SCALE_FACTOR
        eta_r_dir = coeffs_dir[1] / ETA_SCALE_FACTOR
        eta_vr_dir = coeffs_dir[2] / ETA_SCALE_FACTOR
        eta_theta_dir = coeffs_dir[3] / ETA_SCALE_FACTOR

        se_eta_r_dir = se_dir[1] / ETA_SCALE_FACTOR
        se_eta_vr_dir = se_dir[2] / ETA_SCALE_FACTOR
        se_eta_theta_dir = se_dir[3] / ETA_SCALE_FACTOR

        t_r_dir = eta_r_dir / se_eta_r_dir if se_eta_r_dir > 0 else 0.0
        t_vr_dir = eta_vr_dir / se_eta_vr_dir if se_eta_vr_dir > 0 else 0.0
        t_theta_dir = eta_theta_dir / se_eta_theta_dir if se_eta_theta_dir > 0 else 0.0
        p_r_dir = (2 * (1 - stats.t.cdf(abs(t_r_dir), n_clean - 5))
                    if se_eta_r_dir > 0 else 1.0)
        p_vr_dir = (2 * (1 - stats.t.cdf(abs(t_vr_dir), n_clean - 5))
                     if se_eta_vr_dir > 0 else 1.0)
        p_theta_dir = (2 * (1 - stats.t.cdf(abs(t_theta_dir), n_clean - 5))
                        if se_eta_theta_dir > 0 else 1.0)

        aic_dir = n_clean * np.log(rss_dir / n_clean) + 2 * 5
        bic_dir = n_clean * np.log(rss_dir / n_clean) + 5 * np.log(n_clean)

        print_status(
            f"{label}: η_θ = {eta_theta_dir:.4e} ± {se_eta_theta_dir:.4e} "
            f"(t={t_theta_dir:.2f}, p={p_theta_dir:.4f}); AIC={aic_dir:.1f}",
            "CALC",
        )

        return {
            "available": True,
            "label": label,
            "eta_0": float(eta_0_dir),
            "eta_0_error": float(se_dir[0] / ETA_SCALE_FACTOR),
            "eta_r": float(eta_r_dir),
            "eta_r_error": float(se_eta_r_dir),
            "eta_r_t": float(t_r_dir),
            "eta_r_p": float(p_r_dir),
            "eta_vr": float(eta_vr_dir),
            "eta_vr_error": float(se_eta_vr_dir),
            "eta_vr_t": float(t_vr_dir),
            "eta_vr_p": float(p_vr_dir),
            "eta_theta": float(eta_theta_dir),
            "eta_theta_error": float(se_eta_theta_dir),
            "eta_theta_t": float(t_theta_dir),
            "eta_theta_p": float(p_theta_dir),
            "rss": float(rss_dir),
            "aic": float(aic_dir),
            "bic": float(bic_dir),
        }

    # True CMB (already computed in Refinement B, re-run for consistent comparison)
    joint_true = run_joint_for_direction(_CMB_UNIT, "True CMB")

    # 90° null-1: perpendicular in equatorial plane (cross with celestial pole).
    # For a pure dipole, an exactly 90° direction should yield η_θ ≈ 0.
    joint_90a = run_joint_for_direction(perp1, "CMB +90° (perp1)")

    # 90° null-2: perpendicular in meridian plane (cross CMB × perp1).
    # Tests whether the null result is specific to one perpendicular plane.
    joint_90b = run_joint_for_direction(perp2, "CMB +90° (perp2)")

    # 180° control: true antipode (-CMB). For a dipole field, reversing the
    # axis preserves the magnitude and reverses the sign (η_θ → -η_θ).
    joint_180 = run_joint_for_direction(-_CMB_UNIT, "CMB +180° (antipode)")

    # Aggregate comparison
    null_control_results = {
        "true_cmb": joint_true,
        "perpendicular_1": joint_90a,
        "perpendicular_2": joint_90b,
        "antipode": joint_180,
    }

    if all(r.get("available") for r in null_control_results.values()):
        aic_true = joint_true["aic"]
        aic_90a = joint_90a["aic"]
        aic_90b = joint_90b["aic"]
        aic_180 = joint_180["aic"]

        delta_aic_90a = aic_90a - aic_true
        delta_aic_90b = aic_90b - aic_true
        delta_aic_180 = aic_180 - aic_true

        sig_true = abs(joint_true["eta_theta"] / joint_true["eta_theta_error"])
        sig_90a = abs(joint_90a["eta_theta"] / joint_90a["eta_theta_error"])
        sig_90b = abs(joint_90b["eta_theta"] / joint_90b["eta_theta_error"])
        sig_180 = abs(joint_180["eta_theta"] / joint_180["eta_theta_error"])

        print_status(
            f"ΔAIC(perp1 vs true) = {delta_aic_90a:.1f}; "
            f"ΔAIC(perp2 vs true) = {delta_aic_90b:.1f}; "
            f"ΔAIC(antipode vs true) = {delta_aic_180:.1f}",
            "CALC",
        )
        print_status(
            f"|η_θ|/σ: true={sig_true:.2f}, perp1={sig_90a:.2f}, "
            f"perp2={sig_90b:.2f}, antipode={sig_180:.2f}",
            "CALC",
        )

        null_control_summary = {
            "delta_aic_perp1_vs_true": float(delta_aic_90a),
            "delta_aic_perp2_vs_true": float(delta_aic_90b),
            "delta_aic_antipode_vs_true": float(delta_aic_180),
            "eta_theta_significance_true": float(sig_true),
            "eta_theta_significance_perp1": float(sig_90a),
            "eta_theta_significance_perp2": float(sig_90b),
            "eta_theta_significance_antipode": float(sig_180),
            "true_perp1_preferred": delta_aic_90a > 2,
            "true_perp2_preferred": delta_aic_90b > 2,
            "dipole_antipode_consistent": abs(delta_aic_180) < 2,
            "true_direction_preferred": delta_aic_90a > 2 and delta_aic_90b > 2,
        }
    else:
        null_control_summary = {
            "available": False,
            "reason": "one_or_more_models_failed",
        }

    # ------------------------------------------------------------------
    # 7c. Higher-order multipole test
    # ------------------------------------------------------------------
    print_status("═══ Higher-order multipole test ═══", "TITLE")

    # Legendre polynomials: P2(x) = (3x^2 - 1)/2, P3(x) = (5x^3 - 3x)/2
    cos_theta_full = cos_theta_c
    p2 = 0.5 * (3 * cos_theta_full ** 2 - 1)
    p3 = 0.5 * (5 * cos_theta_full ** 3 - 3 * cos_theta_full)

    # Joint dipole + quadrupole + octupole model
    X_multi = np.column_stack([
        cosD,
        r_c * cosD,
        vr_c * cosD,
        cos_theta_full * cosD,
        p2 * cosD,
        p3 * cosD,
        np.ones(n_clean),
    ])
    with np.errstate(divide='ignore', invalid='ignore', over='ignore'):
        coeffs_multi, _, rank_multi, _ = stable_lstsq(X_multi, res)

    multipole_result = None
    if rank_multi == 7:
        with np.errstate(divide='ignore', invalid='ignore', over='ignore'):
            resid_multi = res - X_multi @ coeffs_multi
        rss_multi = np.sum(resid_multi ** 2)
        mse_multi = rss_multi / (n_clean - 7)
        cov_multi = mse_multi * np.linalg.pinv(X_multi.T @ X_multi, rcond=1e-10, hermitian=True)
        se_multi = np.sqrt(np.diag(cov_multi))

        eta_multi = {
            "dipole": coeffs_multi[3] / ETA_SCALE_FACTOR,
            "quadrupole": coeffs_multi[4] / ETA_SCALE_FACTOR,
            "octupole": coeffs_multi[5] / ETA_SCALE_FACTOR,
        }
        se_multi_eta = {
            "dipole": se_multi[3] / ETA_SCALE_FACTOR,
            "quadrupole": se_multi[4] / ETA_SCALE_FACTOR,
            "octupole": se_multi[5] / ETA_SCALE_FACTOR,
        }

        # F-test: dipole+quadrupole+octupole vs. dipole-only (controlling for r, vr)
        X_multi_null = np.column_stack([
            cosD, r_c * cosD, vr_c * cosD, cos_theta_full * cosD, np.ones(n_clean),
        ])
        with np.errstate(divide='ignore', invalid='ignore', over='ignore'):
            coeffs_multi_null, _, _, _ = stable_lstsq(X_multi_null, res)
            resid_multi_null = res - X_multi_null @ coeffs_multi_null
        rss_multi_null = np.sum(resid_multi_null ** 2)
        delta_rss_multi = rss_multi_null - rss_multi
        f_multi = (delta_rss_multi / 2) / (rss_multi / (n_clean - 7))
        p_f_multi = 1 - stats.f.cdf(f_multi, 2, n_clean - 7)

        multipole_result = {
            "eta_dipole": float(eta_multi["dipole"]),
            "eta_dipole_error": float(se_multi_eta["dipole"]),
            "eta_dipole_t": float(eta_multi["dipole"] / se_multi_eta["dipole"]),
            "eta_quadrupole": float(eta_multi["quadrupole"]),
            "eta_quadrupole_error": float(se_multi_eta["quadrupole"]),
            "eta_quadrupole_t": float(eta_multi["quadrupole"] / se_multi_eta["quadrupole"]),
            "eta_octupole": float(eta_multi["octupole"]),
            "eta_octupole_error": float(se_multi_eta["octupole"]),
            "eta_octupole_t": float(eta_multi["octupole"] / se_multi_eta["octupole"]),
            "f_statistic": float(f_multi),
            "f_p_value": float(p_f_multi),
            "aic": float(n_clean * np.log(rss_multi / n_clean) + 2 * 7),
        }

        print_status(
            f"Dipole: η = {eta_multi['dipole']:.4e} ± {se_multi_eta['dipole']:.4e}",
            "CALC",
        )
        print_status(
            f"Quadrupole: η = {eta_multi['quadrupole']:.4e} ± {se_multi_eta['quadrupole']:.4e}",
            "CALC",
        )
        print_status(
            f"Octupole: η = {eta_multi['octupole']:.4e} ± {se_multi_eta['octupole']:.4e}",
            "CALC",
        )
        print_status(
            f"Joint F(2, {n_clean - 7}) = {f_multi:.2f}, p = {p_f_multi:.4f}",
            "CALC",
        )

    # ------------------------------------------------------------------
    # 8. Bootstrap robustness analysis for full joint regression
    # ------------------------------------------------------------------
    print_status("═══ Bootstrap robustness (n=200) ═══", "TITLE")

    rng = np.random.default_rng(42)
    n_boot = 200
    boot_coeffs = []

    for b in range(n_boot):
        idx = rng.integers(0, n_clean, size=n_clean)
        X_b = np.column_stack([
            cosD[idx],
            r_c[idx] * cosD[idx],
            vr_c[idx] * cosD[idx],
            cos_theta_c[idx] * cosD[idx],
            np.ones(n_clean),
        ])
        coeffs_b, _, rank_b, _ = stable_lstsq(X_b, res[idx])
        if rank_b == 5:
            boot_coeffs.append(coeffs_b)

    boot_result = {"available": False}
    if boot_coeffs:
        boot_arr = np.array(boot_coeffs)
        # Convert to eta-scale
        eta_theta_boot = boot_arr[:, 3] / ETA_SCALE_FACTOR
        eta_vr_boot = boot_arr[:, 2] / ETA_SCALE_FACTOR
        eta_r_boot = boot_arr[:, 1] / ETA_SCALE_FACTOR

        def boot_stats(arr):
            return {
                "mean": float(np.mean(arr)),
                "std": float(np.std(arr)),
                "ci95_lower": float(np.percentile(arr, 2.5)),
                "ci95_upper": float(np.percentile(arr, 97.5)),
            }

        boot_result = {
            "available": True,
            "n_successful": len(boot_coeffs),
            "eta_theta": boot_stats(eta_theta_boot),
            "eta_vr": boot_stats(eta_vr_boot),
            "eta_r": boot_stats(eta_r_boot),
        }

        print_status(
            f"η_θ bootstrap: {boot_result['eta_theta']['mean']:.4e} "
            f"± {boot_result['eta_theta']['std']:.4e} "
            f"[95% CI: {boot_result['eta_theta']['ci95_lower']:.4e}, "
            f"{boot_result['eta_theta']['ci95_upper']:.4e}]",
            "CALC",
        )
        print_status(
            f"η_vr bootstrap: {boot_result['eta_vr']['mean']:.4e} "
            f"± {boot_result['eta_vr']['std']:.4e} "
            f"[95% CI: {boot_result['eta_vr']['ci95_lower']:.4e}, "
            f"{boot_result['eta_vr']['ci95_upper']:.4e}]",
            "CALC",
        )
        print_status(
            f"η_r bootstrap: {boot_result['eta_r']['mean']:.4e} "
            f"± {boot_result['eta_r']['std']:.4e} "
            f"[95% CI: {boot_result['eta_r']['ci95_lower']:.4e}, "
            f"{boot_result['eta_r']['ci95_upper']:.4e}]",
            "CALC",
        )

    # ------------------------------------------------------------------
    # 9. Annual envelope of monthly anisotropy
    # ------------------------------------------------------------------
    print_status("═══ Annual envelope of monthly anisotropy ═══", "TITLE")

    # Annual phase (zero at JD 2451545.0, Jan 1 2000)
    t_days_c = jd[mask_clean] - 2451545.0
    annual_phase_sin = np.sin(2 * np.pi * t_days_c / 365.25)
    annual_phase_cos = np.cos(2 * np.pi * t_days_c / 365.25)

    # Model: synodic + cosθ + cosθ×annual_sin + cosθ×annual_cos
    X_env = np.column_stack([
        cosD,
        cos_theta_c * cosD,
        cos_theta_c * annual_phase_sin * cosD,
        cos_theta_c * annual_phase_cos * cosD,
        np.ones(n_clean),
    ])
    with np.errstate(divide='ignore', invalid='ignore', over='ignore'):
        coeffs_env, _, rank_env, _ = stable_lstsq(X_env, res)

    envelope_result = None
    if rank_env == 5:
        with np.errstate(divide='ignore', invalid='ignore', over='ignore'):
            resid_env = res - X_env @ coeffs_env
        mse_env = np.sum(resid_env ** 2) / (n_clean - 5)
        XtX_env_inv = np.linalg.pinv(X_env.T @ X_env, rcond=1e-10, hermitian=True)
        cov_env = mse_env * XtX_env_inv
        se_env = np.sqrt(np.diag(cov_env))

        eta_env_sin = coeffs_env[2] / ETA_SCALE_FACTOR
        eta_env_cos = coeffs_env[3] / ETA_SCALE_FACTOR
        se_env_sin = se_env[2] / ETA_SCALE_FACTOR
        se_env_cos = se_env[3] / ETA_SCALE_FACTOR

        t_env_sin = eta_env_sin / se_env_sin if se_env_sin > 0 else 0.0
        t_env_cos = eta_env_cos / se_env_cos if se_env_cos > 0 else 0.0
        p_env_sin = (2 * (1 - stats.t.cdf(abs(t_env_sin), n_clean - 5))
                     if se_env_sin > 0 else 1.0)
        p_env_cos = (2 * (1 - stats.t.cdf(abs(t_env_cos), n_clean - 5))
                     if se_env_cos > 0 else 1.0)

        # F-test for both envelope terms jointly
        rss_env = np.sum(resid_env ** 2)
        rss_no_env = rss_full if full_result.get("available") else rss_base
        if full_result.get("available"):
            delta_rss_env = rss_no_env - rss_env
            f_env = (delta_rss_env / 2) / (rss_env / (n_clean - 5))
            p_f_env = 1 - stats.f.cdf(f_env, 2, n_clean - 5)
        else:
            f_env = None
            p_f_env = None

        envelope_result = {
            "available": True,
            "eta_env_sin": float(eta_env_sin),
            "eta_env_sin_error": float(se_env_sin),
            "eta_env_sin_t": float(t_env_sin),
            "eta_env_sin_p": float(p_env_sin),
            "eta_env_cos": float(eta_env_cos),
            "eta_env_cos_error": float(se_env_cos),
            "eta_env_cos_t": float(t_env_cos),
            "eta_env_cos_p": float(p_env_cos),
            "f_statistic": float(f_env) if f_env is not None else None,
            "f_p_value": float(p_f_env) if p_f_env is not None else None,
            "rss": float(rss_env),
            "aic": float(n_clean * np.log(rss_env / n_clean) + 2 * 5),
        }

        print_status(
            f"Annual envelope sin: η = {eta_env_sin:.4e} ± {se_env_sin:.4e} "
            f"(t={t_env_sin:.2f}, p={p_env_sin:.4f})",
            "CALC",
        )
        print_status(
            f"Annual envelope cos: η = {eta_env_cos:.4e} ± {se_env_cos:.4e} "
            f"(t={t_env_cos:.2f}, p={p_env_cos:.4f})",
            "CALC",
        )
        if f_env is not None:
            print_status(
                f"Joint F(2, {n_clean - 5}) = {f_env:.2f}, p = {p_f_env:.4f}",
                "CALC",
            )

    # ------------------------------------------------------------------
    # 10. Cross-station consistency
    # ------------------------------------------------------------------
    print_status("═══ Cross-station consistency ═══", "TITLE")

    station_ids = df["station"].values[mask_clean]
    unique_stations = np.unique(station_ids)
    station_results = []

    for st in unique_stations:
        mask_st = station_ids == st
        n_st = mask_st.sum()
        if n_st < 500:
            continue

        cosD_st = cosD[mask_st]
        res_st = res[mask_st]
        cos_theta_st = cos_theta_c[mask_st]
        r_st = r_c[mask_st]
        vr_st = vr_c[mask_st]

        # --- Simple split test (for backward compatibility) ---
        ct_p85_st = np.percentile(cos_theta_st, 85)
        ct_p15_st = np.percentile(cos_theta_st, 15)
        mask_high_st = cos_theta_st >= ct_p85_st
        mask_low_st = cos_theta_st <= ct_p15_st

        if mask_high_st.sum() < 30 or mask_low_st.sum() < 30:
            continue

        reg_high_st = linear_regression(res_st[mask_high_st], cosD_st[mask_high_st])
        reg_low_st = linear_regression(res_st[mask_low_st], cosD_st[mask_low_st])

        diff_st = reg_high_st["eta"] - reg_low_st["eta"]
        diff_err_st = np.sqrt(
            reg_high_st["eta_error"] ** 2 + reg_low_st["eta_error"] ** 2
        )
        sig_st = diff_st / diff_err_st if diff_err_st > 0 else 0.0

        # --- Full joint regression per station ---
        X_st = np.column_stack([
            cosD_st, r_st * cosD_st, vr_st * cosD_st,
            cos_theta_st * cosD_st, np.ones(n_st),
        ])
        with np.errstate(divide='ignore', invalid='ignore', over='ignore'):
            coeffs_st, _, rank_st, _ = stable_lstsq(X_st, res_st)

        joint_st = {"eta_r": None, "eta_vr": None, "eta_theta": None}
        if rank_st == 5:
            with np.errstate(divide='ignore', invalid='ignore', over='ignore'):
                resid_st = res_st - X_st @ coeffs_st
            mse_st = np.sum(resid_st ** 2) / (n_st - 5)
            cov_st = mse_st * np.linalg.pinv(X_st.T @ X_st, rcond=1e-10, hermitian=True)
            se_st = np.sqrt(np.diag(cov_st))
            joint_st = {
                "eta_r": float(coeffs_st[1] / ETA_SCALE_FACTOR),
                "eta_r_error": float(se_st[1] / ETA_SCALE_FACTOR),
                "eta_r_t": float((coeffs_st[1] / se_st[1]) if se_st[1] > 0 else 0.0),
                "eta_vr": float(coeffs_st[2] / ETA_SCALE_FACTOR),
                "eta_vr_error": float(se_st[2] / ETA_SCALE_FACTOR),
                "eta_vr_t": float((coeffs_st[2] / se_st[2]) if se_st[2] > 0 else 0.0),
                "eta_theta": float(coeffs_st[3] / ETA_SCALE_FACTOR),
                "eta_theta_error": float(se_st[3] / ETA_SCALE_FACTOR),
                "eta_theta_t": float((coeffs_st[3] / se_st[3]) if se_st[3] > 0 else 0.0),
            }

        station_results.append({
            "station": str(st),
            "n_observations": int(n_st),
            "delta_eta": float(diff_st),
            "delta_eta_error": float(diff_err_st),
            "significance_sigma": float(sig_st),
            "joint_regression": joint_st,
        })
        print_status(
            f"{st}: split Δη = {diff_st:.4e} ({abs(sig_st):.2f}σ); "
            f"joint η_r={joint_st.get('eta_r', 0):.2e} "
            f"η_vr={joint_st.get('eta_vr', 0):.2e} "
            f"η_θ={joint_st.get('eta_theta', 0):.2e}",
            "CALC",
        )

    n_stations_significant = sum(
        1 for s in station_results if abs(s["significance_sigma"]) > 2.0
    )
    if station_results:
        print_status(
            f"Stations with significant cosθ effect: "
            f"{n_stations_significant}/{len(station_results)}",
            "CALC",
        )

    # ------------------------------------------------------------------
    # 11. Compile results
    # ------------------------------------------------------------------
    if joint_result.get("available"):
        if joint_result.get("eta_vpar_p", 1.0) < 0.05 or joint_result.get("eta_theta_p", 1.0) < 0.05:
            cmb_detection_result = "SIGNIFICANT"
            print_status(
                "RESULT: Significant CMB-frame anisotropy detected.",
                "SUCCESS",
            )
        elif diff_vp_sig > 2.0 or diff_ct_sig > 2.0:
            cmb_detection_result = "MARGINAL"
            print_status(
                "RESULT: Marginal CMB anisotropy detected.",
                "SUCCESS",
            )
        else:
            cmb_detection_result = "NOT_SIGNIFICANT"
            print_status(
                "RESULT: No significant CMB anisotropy detected.",
                "INFO",
            )
        pipeline_status = "PASS"
    else:
        cmb_detection_result = "MODEL_FAILED"
        pipeline_status = "FAIL"
        print_status("RESULT: Joint model failed to converge.", "ERROR")

    results = {
        "step_id": "step_048",
        "status": pipeline_status,
        "cmb_detection_result": cmb_detection_result,
        "n_observations": int(n_clean),
        "n_outliers_removed": int(n - n_clean),
        "cmb_direction": {
            "galactic_l_deg": 264.02,
            "galactic_b_deg": 48.25,
            "equatorial_ra_deg": 168.14,
            "equatorial_dec_deg": -7.22,
        },
        "kinematic_ranges": {
            "v_parallel_kms_min": float(v_par.min()),
            "v_parallel_kms_max": float(v_par.max()),
            "v_perp_kms_min": float(v_perp.min()),
            "v_perp_kms_max": float(v_perp.max()),
            "cos_theta_min": float(cos_theta.min()),
            "cos_theta_max": float(cos_theta.max()),
        },
        "correlations": {
            "v_parallel_vs_distance": float(corr_vpar_r),
            "cosD_vs_cos_theta": float(corr_cosD_cosTheta),
        },
        "model_A_annual_velocity_projection": {
            "high_vpar_eta": float(reg_vp_high["eta"]),
            "high_vpar_eta_error": float(reg_vp_high["eta_error"]),
            "low_vpar_eta": float(reg_vp_low["eta"]),
            "low_vpar_eta_error": float(reg_vp_low["eta_error"]),
            "delta_eta": float(diff_vp),
            "delta_eta_error": float(diff_vp_err),
            "significance_sigma": float(diff_vp_sig),
        },
        "model_B_monthly_orientation_anisotropy": {
            "high_cos_theta_eta": float(reg_ct_high["eta"]),
            "high_cos_theta_eta_error": float(reg_ct_high["eta_error"]),
            "low_cos_theta_eta": float(reg_ct_low["eta"]),
            "low_cos_theta_eta_error": float(reg_ct_low["eta_error"]),
            "delta_eta": float(diff_ct),
            "delta_eta_error": float(diff_ct_err),
            "significance_sigma": float(diff_ct_sig),
        },
        "model_C_joint_cmb_fit": joint_result,
        "refinement_A_orthogonalized_vpar": orth_result,
        "refinement_B_full_joint_regression": full_result,
        "refinement_C_nested_cos_theta": nested_result,
        "directional_specificity_test": {
            "true_cmb_significance_sigma": float(true_sig),
            "null_directions": direction_results,
            "n_null_significant_gt_2sigma": int(n_null_significant),
            "dipole_consistency_ratio": float(dipole_consistency),
            "perpendicular_suppression_ratio": float(perp_suppression),
        },
        "joint_model_null_direction_control": {
            "models": null_control_results,
            "summary": null_control_summary,
        },
        "higher_order_multipole_test": multipole_result,
        "bootstrap_robustness": boot_result,
        "annual_envelope_test": envelope_result,
        "cross_station_consistency": {
            "stations": station_results,
            "n_stations_significant_gt_2sigma": int(n_stations_significant),
            "n_stations_tested": len(station_results),
        },
        "phase_shifted_annual_test": phase_result,
        "binned_orientation_trend": orient_trend,
    }

    return results


if __name__ == "__main__":
    log_dir = PROJECT_ROOT / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    logger = TEPLogger(
        "step_048", str(log_dir / "step_048_cmb_anisotropy.log")
    )
    set_step_logger(logger)

    data_path = (
        PROJECT_ROOT / "data" / "processed" / "INPOP19a_all_stations_residuals.csv"
    )
    if not data_path.exists():
        print_status(f"No processed INPOP19a residuals at {data_path}", "ERROR")
        sys.exit(1)

    df = pd.read_csv(data_path)

    results = cmb_anisotropy_analysis(df, verbose=True)

    if results:
        logger.save_step_results(
            results, PROJECT_ROOT, "step_048_cmb_anisotropy"
        )
        print_status("CMB Anisotropy analysis complete.", "SUCCESS")
    else:
        print_status("Analysis failed.", "ERROR")
        sys.exit(1)
