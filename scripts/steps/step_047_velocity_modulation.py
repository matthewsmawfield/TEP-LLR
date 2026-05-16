#!/usr/bin/env python3
r"""
Step 047: Orbital Velocity Modulation of Temporal Shear
========================================================

Tests the TEP-specific prediction that temporal shear experienced by the
Earth-Moon system depends not only on heliocentric distance (scalar gradient
steepness) but also on the orbital velocity through the solar temporal topology.

Physical Mechanism:
In TEP, the scalar field $\phi$ possesses a spatial gradient $\nabla\phi$ that
defines the Temporal Topology. A body moving through this topology at velocity
$\mathbf{v}$ experiences an effective temporal shear rate:

    $d\phi/dt = \partial\phi/\partial t + \mathbf{v} \cdot \nabla\phi$

For a quasi-static solar field, $\partial\phi/\partial t \approx 0$, so the
dominant term is $\mathbf{v} \cdot \nabla\phi \approx v_r \, |\nabla\phi|$,
where $v_r$ is the heliocentric radial velocity.

Crucially, in a Kepler orbit with small eccentricity $e$:
    $r \approx a(1 - e\cos E)$          (distance, max at aphelion)
    $v_r \approx e\sin E \sqrt{GM/a}$   (radial velocity, zero at apsides)

Distance and radial velocity are approximately in **quadrature** (90° out of
phase). A distance-only model and a velocity-only model make orthogonal
predictions. A joint fit can distinguish both effects.

This step computes Earth's orbital velocity vector for every LLR epoch using
the DE440 ephemeris via Skyfield, then tests:
  1. Does $\eta$ correlate with orbital speed $|v|$?
  2. Does $\eta$ correlate with radial velocity $v_r$?
  3. Joint model: $\eta = \eta_0 + \eta_r \Delta r + \eta_v v_r + \epsilon$
  4. Compare AIC/BIC to determine whether velocity adds explanatory power.

This directly tests whether the TEP temporal topology is dynamical (velocity-
dependent) or purely static (distance-only).
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

from scripts.utils.astronomical_utils import load_skyfield_planets
from scripts.utils.llr_constants import ETA_SCALE_FACTOR
from scripts.utils.statistical_utils import linear_regression, detect_outliers_sigma
from scripts.utils.logger import TEPLogger, set_step_logger, print_status

# CMB dipole direction (Planck 2018, J2000)
_CMB_RA = np.deg2rad(168.14)
_CMB_DEC = np.deg2rad(-7.22)
_CMB_UNIT = np.array([
    np.cos(_CMB_DEC) * np.cos(_CMB_RA),
    np.cos(_CMB_DEC) * np.sin(_CMB_RA),
    np.sin(_CMB_DEC),
])


def compute_orbital_velocities(jd_array):
    """Compute Earth's heliocentric orbital velocity components.

    Uses Skyfield with DE440 ephemeris to compute the Earth-Sun relative
    velocity vector at each Julian date.

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
        radial_velocity_kms : ndarray
            Radial component v_r = (v · r_hat) [km/s]
            Positive = moving away from Sun, Negative = moving toward Sun
        tangential_velocity_kms : ndarray
            Tangential component |v_⊥| [km/s]
        distance_au : ndarray
            Heliocentric distance [AU]
        earth_moon_cos_theta : ndarray
            Cosine of angle between Earth-Moon vector and CMB dipole
    """
    planets, _eph_path = load_skyfield_planets(PROJECT_ROOT)
    earth = planets["earth"]
    moon = planets["moon"]
    sun = planets["sun"]
    ts = load.timescale()

    timestamps = ts.tt(jd=jd_array)

    # Earth position and velocity in Barycentric frame
    earth_pos = earth.at(timestamps)
    sun_pos = sun.at(timestamps)
    moon_pos = moon.at(timestamps)

    # Relative position: Earth - Sun
    rel_pos = earth_pos.position.km - sun_pos.position.km  # shape (3, N)
    rel_vel = earth_pos.velocity.km_per_s - sun_pos.velocity.km_per_s  # shape (3, N)

    # Heliocentric distance
    distance_km = np.linalg.norm(rel_pos, axis=0)
    distance_au = distance_km / 1.495978707e8  # 1 AU in km

    # Radial unit vector (Earth -> Sun direction, outward from Sun)
    r_hat = rel_pos / distance_km

    # Orbital speed
    speed_kms = np.linalg.norm(rel_vel, axis=0)

    # Radial velocity: positive = moving away from Sun
    radial_velocity_kms = np.sum(rel_vel * r_hat, axis=0)

    # Tangential velocity (perpendicular component)
    v_radial_vec = radial_velocity_kms * r_hat
    v_tangential_vec = rel_vel - v_radial_vec
    tangential_velocity_kms = np.linalg.norm(v_tangential_vec, axis=0)

    # Earth-Moon vector (Earth -> Moon)
    em_vec = moon_pos.position.km - earth_pos.position.km  # shape (3, N)
    em_dist = np.linalg.norm(em_vec, axis=0)
    em_hat = em_vec / em_dist

    # CMB dipole projection of Earth-Moon orientation
    earth_moon_cos_theta = np.sum(em_hat * _CMB_UNIT[:, None], axis=0)

    return {
        "speed_kms": speed_kms,
        "radial_velocity_kms": radial_velocity_kms,
        "tangential_velocity_kms": tangential_velocity_kms,
        "distance_au": distance_au,
        "earth_moon_cos_theta": earth_moon_cos_theta,
    }


def velocity_modulation_analysis(df, verbose=False):
    """Main analysis: test for velocity-dependent TEP modulation.

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
        "═══ Step 047: Orbital Velocity Modulation of Temporal Shear ═══",
        "TITLE",
    )

    n = len(df)
    residuals = df["residual_m"].values
    cos_elong = np.cos(df["elongation_rad"].values)
    jd = df["date_julian"].values

    # ------------------------------------------------------------------
    # 1. Compute orbital kinematics for every observation
    # ------------------------------------------------------------------
    print_status("Computing Earth orbital velocities via DE440/Skyfield...", "PROCESS")
    vel_data = compute_orbital_velocities(jd)

    df = df.copy()
    df["speed_kms"] = vel_data["speed_kms"]
    df["radial_velocity_kms"] = vel_data["radial_velocity_kms"]
    df["tangential_velocity_kms"] = vel_data["tangential_velocity_kms"]
    df["sun_distance_au"] = vel_data["distance_au"]
    df["cos_theta"] = vel_data["earth_moon_cos_theta"]

    print_status(
        f"Orbital speed range: {vel_data['speed_kms'].min():.3f} – "
        f"{vel_data['speed_kms'].max():.3f} km/s",
        "CALC",
    )
    print_status(
        f"Radial velocity range: {vel_data['radial_velocity_kms'].min():.3f} – "
        f"{vel_data['radial_velocity_kms'].max():.3f} km/s",
        "CALC",
    )
    print_status(
        f"Heliocentric distance range: {vel_data['distance_au'].min():.5f} – "
        f"{vel_data['distance_au'].max():.5f} AU",
        "CALC",
    )

    # Correlation between distance and radial velocity (should be ~0 for quad)
    corr_r_vr = np.corrcoef(df["sun_distance_au"], df["radial_velocity_kms"])[0, 1]
    print_status(
        f"Correlation r(distance, v_r) = {corr_r_vr:.4f} "
        f"(quadrature prediction: ~0)",
        "CALC",
    )

    # Pre-filter outliers
    outlier_mask = detect_outliers_sigma(df["residual_m"].values, sigma_threshold=6.0)
    df_clean = df[~outlier_mask]
    n_clean = len(df_clean)
    print_status(
        f"Cleaned dataset: {n_clean:,} observations "
        f"({n - n_clean} outliers removed)",
        "CALC",
    )

    res = df_clean["residual_m"].values
    cosD = np.cos(df_clean["elongation_rad"].values)
    r = df_clean["sun_distance_au"].values
    vr = df_clean["radial_velocity_kms"].values
    v = df_clean["speed_kms"].values
    cos_theta = df_clean["cos_theta"].values

    # Center predictors to reduce collinearity
    r_c = r - np.mean(r)
    vr_c = vr - np.mean(vr)
    v_c = v - np.mean(v)
    cos_theta_c = cos_theta - np.mean(cos_theta)

    # ------------------------------------------------------------------
    # 2. Model A: Distance-only (reproduction of Step 022 logic)
    # ------------------------------------------------------------------
    print_status("═══ MODEL A: Distance-only modulation ═══", "TITLE")

    # Split by distance: deep perihelion (closest 15%) vs deep aphelion (furthest 15%)
    p15 = np.percentile(r, 15)
    p85 = np.percentile(r, 85)
    mask_peri = r <= p15
    mask_aph = r >= p85

    reg_peri = linear_regression(res[mask_peri], cosD[mask_peri])
    reg_aph = linear_regression(res[mask_aph], cosD[mask_aph])

    eta_r_peri = reg_peri["eta"]
    eta_r_aph = reg_aph["eta"]
    diff_r = eta_r_peri - eta_r_aph
    diff_r_err = np.sqrt(reg_peri["eta_error"] ** 2 + reg_aph["eta_error"] ** 2)
    diff_r_sig = abs(diff_r) / diff_r_err if diff_r_err > 0 else 0.0

    print_status(
        f"Deep perihelion: η = {eta_r_peri:.4e} ± {reg_peri['eta_error']:.4e} "
        f"({abs(eta_r_peri)/reg_peri['eta_error']:.2f}σ)",
        "CALC",
    )
    print_status(
        f"Deep aphelion:    η = {eta_r_aph:.4e} ± {reg_aph['eta_error']:.4e} "
        f"({abs(eta_r_aph)/reg_aph['eta_error']:.2f}σ)",
        "CALC",
    )
    print_status(
        f"Δη (peri - aph) = {diff_r:.4e} ± {diff_r_err:.4e} ({diff_r_sig:.2f}σ)",
        "CALC",
    )

    # ------------------------------------------------------------------
    # 3. Model B: Velocity-only modulation
    # ------------------------------------------------------------------
    print_status("═══ MODEL B: Velocity-only modulation ═══", "TITLE")

    # Split by radial velocity: fastest approach vs fastest recession
    # For a Kepler orbit, max |v_r| occurs at quadrature (~90° from apsides)
    # Positive v_r = moving away from Sun (post-perihelion)
    # Negative v_r = moving toward Sun (pre-perihelion)
    vr_p15 = np.percentile(vr, 15)
    vr_p85 = np.percentile(vr, 85)
    mask_fast_away = vr >= vr_p85  # fastest recession from Sun
    mask_fast_in = vr <= vr_p15  # fastest approach toward Sun

    reg_away = linear_regression(res[mask_fast_away], cosD[mask_fast_away])
    reg_in = linear_regression(res[mask_fast_in], cosD[mask_fast_in])

    eta_v_away = reg_away["eta"]
    eta_v_in = reg_in["eta"]
    diff_v = eta_v_away - eta_v_in
    diff_v_err = np.sqrt(reg_away["eta_error"] ** 2 + reg_in["eta_error"] ** 2)
    diff_v_sig = abs(diff_v) / diff_v_err if diff_v_err > 0 else 0.0

    print_status(
        f"Fast recession (v_r > {vr_p85:.3f} km/s): "
        f"η = {eta_v_away:.4e} ± {reg_away['eta_error']:.4e} "
        f"({abs(eta_v_away)/reg_away['eta_error']:.2f}σ)",
        "CALC",
    )
    print_status(
        f"Fast approach (v_r < {vr_p15:.3f} km/s):   "
        f"η = {eta_v_in:.4e} ± {reg_in['eta_error']:.4e} "
        f"({abs(eta_v_in)/reg_in['eta_error']:.2f}σ)",
        "CALC",
    )
    print_status(
        f"Δη (away - in) = {diff_v:.4e} ± {diff_v_err:.4e} ({diff_v_sig:.2f}σ)",
        "CALC",
    )

    # ------------------------------------------------------------------
    # 4. Model C: Speed-only modulation (total speed |v|)
    # ------------------------------------------------------------------
    print_status("═══ MODEL C: Total speed modulation ═══", "TITLE")

    v_p15 = np.percentile(v, 15)
    v_p85 = np.percentile(v, 85)
    mask_fast = v >= v_p85
    mask_slow = v <= v_p15

    reg_fast = linear_regression(res[mask_fast], cosD[mask_fast])
    reg_slow = linear_regression(res[mask_slow], cosD[mask_slow])

    eta_speed_fast = reg_fast["eta"]
    eta_speed_slow = reg_slow["eta"]
    diff_speed = eta_speed_fast - eta_speed_slow
    diff_speed_err = np.sqrt(reg_fast["eta_error"] ** 2 + reg_slow["eta_error"] ** 2)
    diff_speed_sig = abs(diff_speed) / diff_speed_err if diff_speed_err > 0 else 0.0

    print_status(
        f"High speed (|v| > {v_p85:.3f} km/s): "
        f"η = {eta_speed_fast:.4e} ± {reg_fast['eta_error']:.4e} "
        f"({abs(eta_speed_fast)/reg_fast['eta_error']:.2f}σ)",
        "CALC",
    )
    print_status(
        f"Low speed (|v| < {v_p15:.3f} km/s):  "
        f"η = {eta_speed_slow:.4e} ± {reg_slow['eta_error']:.4e} "
        f"({abs(eta_speed_slow)/reg_slow['eta_error']:.2f}σ)",
        "CALC",
    )
    print_status(
        f"Δη (fast - slow) = {diff_speed:.4e} ± {diff_speed_err:.4e} "
        f"({diff_speed_sig:.2f}σ)",
        "CALC",
    )

    # ------------------------------------------------------------------
    # 5. Model D: Joint distance + radial velocity fit
    # ------------------------------------------------------------------
    print_status("═══ MODEL D: Joint distance + radial velocity ═══", "TITLE")

    # Joint model: residual = A*cos(D) + B*r*cos(D) + C*vr*cos(D) + intercept
    # Equivalently: η = η_0 + η_r * r_c + η_v * vr_c for the cos(D) amplitude
    # We fit: residual = (a0 + a1*r_c + a2*vr_c) * cos(D) + intercept
    # But this is nonlinear in parameters. Instead, use interaction formulation:
    # residual = β0*cos(D) + β1*r_c*cos(D) + β2*vr_c*cos(D) + β3 + ε

    X_joint = np.column_stack([cosD, r_c * cosD, vr_c * cosD, np.ones(n_clean)])
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

        # Convert to η parameterization
        # residual = [β0 + β1*r_c + β2*vr_c] * cos(D) + β3
        # η(r, vr) = (β0 + β1*r_c + β2*vr_c) / ETA_SCALE_FACTOR
        eta_0_joint = coeffs_joint[0] / ETA_SCALE_FACTOR
        eta_r_joint = coeffs_joint[1] / ETA_SCALE_FACTOR
        eta_vr_joint = coeffs_joint[2] / ETA_SCALE_FACTOR
        intercept_joint = coeffs_joint[3]

        se_eta_0 = se_joint[0] / ETA_SCALE_FACTOR
        se_eta_r = se_joint[1] / ETA_SCALE_FACTOR
        se_eta_vr = se_joint[2] / ETA_SCALE_FACTOR

        t_eta_r = eta_r_joint / se_eta_r if se_eta_r > 0 else 0.0
        t_eta_vr = eta_vr_joint / se_eta_vr if se_eta_vr > 0 else 0.0
        p_eta_r = 2 * (1 - stats.t.cdf(abs(t_eta_r), n_clean - 4)) if se_eta_r > 0 else 1.0
        p_eta_vr = 2 * (1 - stats.t.cdf(abs(t_eta_vr), n_clean - 4)) if se_eta_vr > 0 else 1.0

        # AIC and BIC for model comparison
        rss_joint = np.sum(resid_joint ** 2)
        k_joint = 4
        aic_joint = n_clean * np.log(rss_joint / n_clean) + 2 * k_joint
        bic_joint = n_clean * np.log(rss_joint / n_clean) + k_joint * np.log(n_clean)

        # Simple model (distance only, 2 params: cosD + intercept)
        X_simple = np.column_stack([cosD, np.ones(n_clean)])
        with np.errstate(divide='ignore', invalid='ignore', over='ignore'):
            coeffs_simple, _, _, _ = stable_lstsq(X_simple, res)
            resid_simple = res - X_simple @ coeffs_simple
        rss_simple = np.sum(resid_simple ** 2)
        k_simple = 2
        aic_simple = n_clean * np.log(rss_simple / n_clean) + 2 * k_simple
        bic_simple = n_clean * np.log(rss_simple / n_clean) + k_simple * np.log(n_clean)

        # Velocity-only model (2 params: cosD + vr*cosD + intercept... wait that's 3)
        # Actually: cosD, vr_c*cosD, intercept = 3 params
        X_vel = np.column_stack([cosD, vr_c * cosD, np.ones(n_clean)])
        with np.errstate(divide='ignore', invalid='ignore', over='ignore'):
            coeffs_vel, _, _, _ = stable_lstsq(X_vel, res)
            resid_vel = res - X_vel @ coeffs_vel
        rss_vel = np.sum(resid_vel ** 2)
        k_vel = 3
        aic_vel = n_clean * np.log(rss_vel / n_clean) + 2 * k_vel
        bic_vel = n_clean * np.log(rss_vel / n_clean) + k_vel * np.log(n_clean)

        print_status(
            f"Joint model: η(r, v_r) = "
            f"{eta_0_joint:.4e} + {eta_r_joint:.4e}·Δr + {eta_vr_joint:.4e}·Δv_r",
            "CALC",
        )
        print_status(
            f"  η_r (distance coeff):   {eta_r_joint:.4e} ± {se_eta_r:.4e} "
            f"(t={t_eta_r:.2f}, p={p_eta_r:.4f})",
            "CALC",
        )
        print_status(
            f"  η_vr (velocity coeff):  {eta_vr_joint:.4e} ± {se_eta_vr:.4e} "
            f"(t={t_eta_vr:.2f}, p={p_eta_vr:.4f})",
            "CALC",
        )
        print_status(
            f"  AIC simple={aic_simple:.1f}, vel-only={aic_vel:.1f}, "
            f"joint={aic_joint:.1f}",
            "CALC",
        )
        print_status(
            f"  BIC simple={bic_simple:.1f}, vel-only={bic_vel:.1f}, "
            f"joint={bic_joint:.1f}",
            "CALC",
        )

        delta_aic_vel = aic_vel - aic_simple
        delta_bic_vel = bic_vel - bic_simple
        delta_aic_joint = aic_joint - aic_simple
        delta_bic_joint = bic_joint - bic_simple

        print_status(
            f"  ΔAIC(vel-only vs simple) = {delta_aic_vel:.1f}",
            "CALC",
        )
        print_status(
            f"  ΔAIC(joint vs simple) = {delta_aic_joint:.1f}",
            "CALC",
        )

        # Model selection interpretation
        if delta_aic_joint < -2:
            model_pref = "joint_preferred"
            print_status(
                "Joint model preferred over distance-only (ΔAIC < -2)",
                "SUCCESS",
            )
        elif delta_aic_vel < -2 and delta_aic_joint >= -2:
            model_pref = "velocity_only_preferred"
            print_status(
                "Velocity-only model preferred over distance-only (ΔAIC < -2)",
                "SUCCESS",
            )
        else:
            model_pref = "simple_preferred"
            print_status(
                "Distance-only model sufficient; no velocity signal detected",
                "WARNING",
            )

        joint_result = {
            "available": True,
            "eta_0": float(eta_0_joint),
            "eta_0_error": float(se_eta_0),
            "eta_r": float(eta_r_joint),
            "eta_r_error": float(se_eta_r),
            "eta_r_t": float(t_eta_r),
            "eta_r_p": float(p_eta_r),
            "eta_vr": float(eta_vr_joint),
            "eta_vr_error": float(se_eta_vr),
            "eta_vr_t": float(t_eta_vr),
            "eta_vr_p": float(p_eta_vr),
            "intercept_m": float(intercept_joint),
            "rss": float(rss_joint),
            "aic": float(aic_joint),
            "bic": float(bic_joint),
            "aic_simple": float(aic_simple),
            "bic_simple": float(bic_simple),
            "aic_vel_only": float(aic_vel),
            "bic_vel_only": float(bic_vel),
            "delta_aic_vel_vs_simple": float(delta_aic_vel),
            "delta_aic_joint_vs_simple": float(delta_aic_joint),
            "delta_bic_vel_vs_simple": float(delta_bic_vel),
            "delta_bic_joint_vs_simple": float(delta_bic_joint),
            "model_preference": model_pref,
            "n_params": k_joint,
        }

    # ------------------------------------------------------------------
    # 6. Model E: CMB-controlled joint fit (cos_theta as control)
    # ------------------------------------------------------------------
    print_status("═══ MODEL E: CMB-controlled joint fit ═══", "TITLE")

    # Add CMB orientation as a control variable to isolate velocity effects
    # from cosmological anisotropy. This tests whether the distance/velocity
    # signals persist when the CMB dipole alignment is accounted for.
    X_cmb = np.column_stack([
        cosD, r_c * cosD, vr_c * cosD, cos_theta_c * cosD, np.ones(n_clean)
    ])
    with np.errstate(divide='ignore', invalid='ignore', over='ignore'):
        coeffs_cmb, _, rank_cmb, _ = stable_lstsq(X_cmb, res)

    if rank_cmb < 5:
        print_status("WARNING: CMB-controlled model rank-deficient", "WARNING")
        cmb_controlled_result = {"available": False, "reason": "rank_deficient"}
    else:
        with np.errstate(divide='ignore', invalid='ignore', over='ignore'):
            resid_cmb = res - X_cmb @ coeffs_cmb
        mse_cmb = np.sum(resid_cmb ** 2) / (n_clean - 5)
        XtX_cmb_inv = np.linalg.pinv(X_cmb.T @ X_cmb, rcond=1e-10, hermitian=True)
        cov_cmb = mse_cmb * XtX_cmb_inv
        se_cmb = np.sqrt(np.diag(cov_cmb))

        eta_0_cmb = coeffs_cmb[0] / ETA_SCALE_FACTOR
        eta_r_cmb = coeffs_cmb[1] / ETA_SCALE_FACTOR
        eta_vr_cmb = coeffs_cmb[2] / ETA_SCALE_FACTOR
        eta_theta_cmb = coeffs_cmb[3] / ETA_SCALE_FACTOR
        intercept_cmb = coeffs_cmb[4]

        se_eta_0_cmb = se_cmb[0] / ETA_SCALE_FACTOR
        se_eta_r_cmb = se_cmb[1] / ETA_SCALE_FACTOR
        se_eta_vr_cmb = se_cmb[2] / ETA_SCALE_FACTOR
        se_eta_theta_cmb = se_cmb[3] / ETA_SCALE_FACTOR

        t_eta_r_cmb = eta_r_cmb / se_eta_r_cmb if se_eta_r_cmb > 0 else 0.0
        t_eta_vr_cmb = eta_vr_cmb / se_eta_vr_cmb if se_eta_vr_cmb > 0 else 0.0
        t_eta_theta_cmb = eta_theta_cmb / se_eta_theta_cmb if se_eta_theta_cmb > 0 else 0.0

        p_eta_r_cmb = 2 * (1 - stats.t.cdf(abs(t_eta_r_cmb), n_clean - 5)) if se_eta_r_cmb > 0 else 1.0
        p_eta_vr_cmb = 2 * (1 - stats.t.cdf(abs(t_eta_vr_cmb), n_clean - 5)) if se_eta_vr_cmb > 0 else 1.0
        p_eta_theta_cmb = 2 * (1 - stats.t.cdf(abs(t_eta_theta_cmb), n_clean - 5)) if se_eta_theta_cmb > 0 else 1.0

        rss_cmb = np.sum(resid_cmb ** 2)
        k_cmb = 5
        aic_cmb = n_clean * np.log(rss_cmb / n_clean) + 2 * k_cmb
        bic_cmb = n_clean * np.log(rss_cmb / n_clean) + k_cmb * np.log(n_clean)

        print_status(
            f"CMB-controlled: η(r, v_r, θ) = "
            f"{eta_0_cmb:.4e} + {eta_r_cmb:.4e}·Δr + {eta_vr_cmb:.4e}·Δv_r + "
            f"{eta_theta_cmb:.4e}·cos(θ)",
            "CALC",
        )
        print_status(
            f"  η_r (distance):   {eta_r_cmb:.4e} ± {se_eta_r_cmb:.4e} "
            f"(t={t_eta_r_cmb:.2f}, p={p_eta_r_cmb:.4f})",
            "CALC",
        )
        print_status(
            f"  η_vr (velocity):  {eta_vr_cmb:.4e} ± {se_eta_vr_cmb:.4e} "
            f"(t={t_eta_vr_cmb:.2f}, p={p_eta_vr_cmb:.4e})",
            "CALC",
        )
        print_status(
            f"  η_θ (CMB orient): {eta_theta_cmb:.4e} ± {se_eta_theta_cmb:.4e} "
            f"(t={t_eta_theta_cmb:.2f}, p={p_eta_theta_cmb:.4e})",
            "CALC",
        )
        print_status(
            f"  AIC joint={aic_joint:.1f}, CMB-ctrl={aic_cmb:.1f}",
            "CALC",
        )

        cmb_controlled_result = {
            "available": True,
            "eta_0": float(eta_0_cmb),
            "eta_0_error": float(se_eta_0_cmb),
            "eta_r": float(eta_r_cmb),
            "eta_r_error": float(se_eta_r_cmb),
            "eta_r_t": float(t_eta_r_cmb),
            "eta_r_p": float(p_eta_r_cmb),
            "eta_vr": float(eta_vr_cmb),
            "eta_vr_error": float(se_eta_vr_cmb),
            "eta_vr_t": float(t_eta_vr_cmb),
            "eta_vr_p": float(p_eta_vr_cmb),
            "eta_theta": float(eta_theta_cmb),
            "eta_theta_error": float(se_eta_theta_cmb),
            "eta_theta_t": float(t_eta_theta_cmb),
            "eta_theta_p": float(p_eta_theta_cmb),
            "intercept_m": float(intercept_cmb),
            "rss": float(rss_cmb),
            "aic": float(aic_cmb),
            "bic": float(bic_cmb),
            "delta_aic_vs_joint": float(aic_cmb - aic_joint),
            "delta_bic_vs_joint": float(bic_cmb - bic_joint),
            "n_params": k_cmb,
        }

    # ------------------------------------------------------------------
    # 7. Binning analysis: 10 bins in radial velocity
    # ------------------------------------------------------------------
    print_status("═══ Binned radial-velocity analysis ═══", "TITLE")

    n_bins = 10
    bin_edges = np.linspace(vr.min(), vr.max(), n_bins + 1)
    bin_centers = []
    bin_etas = []
    bin_eta_errs = []

    for i in range(n_bins):
        mask = (vr >= bin_edges[i]) & (vr < bin_edges[i + 1])
        if mask.sum() > 30:
            reg_bin = linear_regression(res[mask], cosD[mask])
            bin_centers.append((bin_edges[i] + bin_edges[i + 1]) / 2)
            bin_etas.append(reg_bin["eta"])
            bin_eta_errs.append(reg_bin["eta_error"])

    # Fit linear trend in η vs v_r
    vel_trend = None
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

            vel_trend = {
                "slope_eta_per_kms": float(slope),
                "slope_error": float(slope_err),
                "t_statistic": float(t_slope),
                "p_value": float(p_slope),
                "n_bins": len(bin_centers),
            }

            print_status(
                f"Binned trend: dη/dv_r = {slope:.4e} ± {slope_err:.4e} "
                f"(t={t_slope:.2f}, p={p_slope:.4f})",
                "CALC",
            )

    # ------------------------------------------------------------------
    # 7. Spectral test: Lomb-Scargle for v_r-modulated signal
    # ------------------------------------------------------------------
    print_status("═══ Spectral sideband test (velocity-modulated) ═══", "TITLE")

    # TEP theory: if η modulates with v_r ~ sin(E), and the base signal is
    # cos(D), then the modulated signal is approximately cos(D) * (1 + m*sin(E)).
    # This creates sidebands at D ± E, where E is the eccentric anomaly frequency
    # (annual). But E ≈ M + e*sin(M) for small e, so the sideband is near D ± l'.

    # For a more direct test: construct a velocity-weighted residual and
    # check if it shows enhanced correlation with cos(D) at specific orbital phases.

    # Test: compute η in 4 orbital quadrants based on (r, v_r) phase
    # Quadrant I:  r < mean, vr > 0  (near perihelion, moving away)
    # Quadrant II: r > mean, vr > 0  (post-aphelion, moving away)
    # Quadrant III: r > mean, vr < 0 (near aphelion, moving in)
    # Quadrant IV: r < mean, vr < 0 (pre-perihelion, moving in)
    r_mean = np.mean(r)
    vr_mean = np.mean(vr)

    quadrants = {
        "QI_peri_away": (r < r_mean) & (vr > vr_mean),
        "QII_post_peri": (r > r_mean) & (vr > vr_mean),
        "QIII_aph_in": (r > r_mean) & (vr < vr_mean),
        "QIV_pre_peri": (r < r_mean) & (vr < vr_mean),
    }

    quadrant_results = []
    for name, mask in quadrants.items():
        if mask.sum() > 30:
            reg_q = linear_regression(res[mask], cosD[mask])
            quadrant_results.append(
                {
                    "quadrant": name,
                    "n": int(mask.sum()),
                    "eta": float(reg_q["eta"]),
                    "eta_error": float(reg_q["eta_error"]),
                    "snr": float(abs(reg_q["eta"]) / reg_q["eta_error"])
                    if reg_q["eta_error"] > 0
                    else 0.0,
                }
            )
            print_status(
                f"  {name}: η = {reg_q['eta']:.4e} ± {reg_q['eta_error']:.4e} "
                f"({abs(reg_q['eta'])/reg_q['eta_error']:.2f}σ, N={mask.sum()})",
                "CALC",
            )

    # ------------------------------------------------------------------
    # 8. TEP prediction test: η should increase with |v_r| if dynamical
    # ------------------------------------------------------------------
    print_status("═══ TEP dynamical prediction test ═══", "TITLE")

    # TEP predicts: experienced temporal shear ∝ |v_r| * |∇φ|
    # If |∇φ| ∝ 1/r, then shear ∝ |v_r|/r
    # But v_r and r are weakly correlated, so |v_r|/r is a distinct predictor.
    shear_proxy = np.abs(vr) / r
    shear_p85 = np.percentile(shear_proxy, 85)
    shear_p15 = np.percentile(shear_proxy, 15)
    mask_high_shear = shear_proxy >= shear_p85
    mask_low_shear = shear_proxy <= shear_p15

    reg_high = linear_regression(res[mask_high_shear], cosD[mask_high_shear])
    reg_low = linear_regression(res[mask_low_shear], cosD[mask_low_shear])

    diff_shear = reg_high["eta"] - reg_low["eta"]
    diff_shear_err = np.sqrt(reg_high["eta_error"] ** 2 + reg_low["eta_error"] ** 2)
    diff_shear_sig = abs(diff_shear) / diff_shear_err if diff_shear_err > 0 else 0.0

    print_status(
        f"High dynamical shear (|v_r|/r > {shear_p85:.3f}): "
        f"η = {reg_high['eta']:.4e} ± {reg_high['eta_error']:.4e}",
        "CALC",
    )
    print_status(
        f"Low dynamical shear (|v_r|/r < {shear_p15:.3f}):  "
        f"η = {reg_low['eta']:.4e} ± {reg_low['eta_error']:.4e}",
        "CALC",
    )
    print_status(
        f"Δη = {diff_shear:.4e} ± {diff_shear_err:.4e} ({diff_shear_sig:.2f}σ)",
        "CALC",
    )

    # ------------------------------------------------------------------
    # 9. Compile results
    # ------------------------------------------------------------------
    if cmb_controlled_result.get("available"):
        # Primary criterion: CMB-controlled model with velocity significance
        if cmb_controlled_result.get("eta_vr_p", 1.0) < 0.05:
            velocity_modulation_result = "SIGNIFICANT_CMB_CONTROLLED"
            print_status(
                "RESULT: Significant velocity-dependent modulation detected "
                "with CMB orientation controlled. TEP temporal topology is "
                "dynamical and cosmologically oriented.",
                "SUCCESS",
            )
        elif diff_v_sig > 2.0 or diff_shear_sig > 2.0:
            velocity_modulation_result = "MARGINAL"
            print_status(
                "RESULT: Marginal velocity modulation detected.",
                "SUCCESS",
            )
        else:
            velocity_modulation_result = "NOT_SIGNIFICANT_CMB_CONTROLLED"
            print_status(
                "RESULT: No significant velocity modulation with CMB control. "
                "Temporal topology may be static or CMB-only.",
                "INFO",
            )
        pipeline_status = "PASS"
    elif joint_result.get("available"):
        if joint_result.get("eta_vr_p", 1.0) < 0.05:
            velocity_modulation_result = "SIGNIFICANT_JOINT_MODEL"
            print_status(
                "RESULT: Significant velocity-dependent modulation detected. "
                "TEP temporal topology is dynamical.",
                "SUCCESS",
            )
        else:
            velocity_modulation_result = "NOT_SIGNIFICANT_JOINT_MODEL"
            print_status(
                "RESULT: No significant velocity modulation. "
                "Temporal topology appears static (distance-only).",
                "INFO",
            )
        pipeline_status = "PASS"
    else:
        velocity_modulation_result = "MODEL_FAILED"
        pipeline_status = "FAIL"
        print_status("RESULT: Joint model failed to converge.", "ERROR")

    results = {
        "step_id": "step_047",
        "status": pipeline_status,
        "velocity_modulation_result": velocity_modulation_result,
        "n_observations": n_clean,
        "n_outliers_removed": int(n - n_clean),
        "kinematic_ranges": {
            "speed_kms_min": float(vel_data["speed_kms"].min()),
            "speed_kms_max": float(vel_data["speed_kms"].max()),
            "radial_velocity_kms_min": float(vel_data["radial_velocity_kms"].min()),
            "radial_velocity_kms_max": float(vel_data["radial_velocity_kms"].max()),
            "distance_au_min": float(vel_data["distance_au"].min()),
            "distance_au_max": float(vel_data["distance_au"].max()),
        },
        "correlation_r_vr": float(corr_r_vr),
        "model_A_distance_only": {
            "perihelion_eta": float(eta_r_peri),
            "perihelion_eta_error": float(reg_peri["eta_error"]),
            "aphelion_eta": float(eta_r_aph),
            "aphelion_eta_error": float(reg_aph["eta_error"]),
            "delta_eta": float(diff_r),
            "delta_eta_error": float(diff_r_err),
            "significance_sigma": float(diff_r_sig),
        },
        "model_B_velocity_only": {
            "fast_recession_eta": float(eta_v_away),
            "fast_recession_eta_error": float(reg_away["eta_error"]),
            "fast_approach_eta": float(eta_v_in),
            "fast_approach_eta_error": float(reg_in["eta_error"]),
            "delta_eta": float(diff_v),
            "delta_eta_error": float(diff_v_err),
            "significance_sigma": float(diff_v_sig),
        },
        "model_C_speed_only": {
            "high_speed_eta": float(eta_speed_fast),
            "high_speed_eta_error": float(reg_fast["eta_error"]),
            "low_speed_eta": float(eta_speed_slow),
            "low_speed_eta_error": float(reg_slow["eta_error"]),
            "delta_eta": float(diff_speed),
            "delta_eta_error": float(diff_speed_err),
            "significance_sigma": float(diff_speed_sig),
        },
        "model_D_joint_fit": joint_result,
        "model_E_cmb_controlled": cmb_controlled_result,
        "binned_velocity_trend": vel_trend,
        "quadrant_analysis": quadrant_results,
        "dynamical_shear_test": {
            "high_shear_eta": float(reg_high["eta"]),
            "high_shear_eta_error": float(reg_high["eta_error"]),
            "low_shear_eta": float(reg_low["eta"]),
            "low_shear_eta_error": float(reg_low["eta_error"]),
            "delta_eta": float(diff_shear),
            "delta_eta_error": float(diff_shear_err),
            "significance_sigma": float(diff_shear_sig),
        },
    }

    return results


if __name__ == "__main__":
    log_dir = PROJECT_ROOT / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    logger = TEPLogger(
        "step_047", str(log_dir / "step_047_velocity_modulation.log")
    )
    set_step_logger(logger)

    data_path = (
        PROJECT_ROOT / "data" / "processed" / "INPOP19a_all_stations_residuals.csv"
    )
    if not data_path.exists():
        print_status(f"No processed INPOP19a residuals at {data_path}", "ERROR")
        sys.exit(1)

    df = pd.read_csv(data_path)

    results = velocity_modulation_analysis(df, verbose=True)

    if results:
        logger.save_step_results(
            results, PROJECT_ROOT, "step_047_velocity_modulation"
        )
        print_status("Orbital Velocity Modulation analysis complete.", "SUCCESS")
    else:
        print_status("Analysis failed.", "ERROR")
        sys.exit(1)
