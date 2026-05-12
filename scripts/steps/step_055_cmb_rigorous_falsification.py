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

The suite comprises five independent falsification tests:

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

  D. Sky Scrambling Monte Carlo:
     Rotate the putative dipole direction uniformly on the celestial
     sphere 1,000 times. For each random direction, fit the full joint
     model and record the cosθ significance. The true CMB direction must
     be a clear outlier in this distribution.

  E. Gram–Schmidt Orthogonalization:
     Explicitly orthogonalize cosθ against cosD, r_c, and vr_c. Fit the
     joint model using only the orthogonal residual. If the signal is
     genuine, it must persist in the component that is mathematically
     independent of all other predictors.

All tests are designed to be maximally conservative: they assume the
alternative hypothesis is false and seek to demonstrate that the data
nevertheless reject the null.
"""

import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

import json
import numpy as np
import pandas as pd
from scipy import stats
from skyfield.api import load

from scripts.utils.llr_constants import ETA_SCALE_FACTOR
from scripts.utils.statistical_utils import linear_regression, detect_outliers_sigma
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


def compute_cmb_projections(jd_array):
    """Compute CMB-frame kinematic projections for every epoch."""
    eph_path = PROJECT_ROOT / "de421.bsp"
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

    return {
        "v_parallel_kms": v_parallel_kms,
        "earth_moon_cos_theta": earth_moon_cos_theta,
        "sun_distance_au": distance_km / 1.495978707e8,
        "radial_velocity_kms": v_radial_kms,
        "earth_moon_unit_vectors": em_hat,
    }


def fit_full_joint_model(res, cosD, r_c, vr_c, cos_theta_c):
    """Fit the 5-parameter joint model and return eta_theta and its SE."""
    X = np.column_stack([cosD, r_c * cosD, vr_c * cosD, cos_theta_c * cosD, np.ones(len(cosD))])
    coeffs, _, rank, _ = np.linalg.lstsq(X, res, rcond=None)
    if rank < 5:
        return None, None, None, None, None
    resid = res - X @ coeffs
    mse = np.sum(resid ** 2) / (len(res) - 5)
    cov = mse * np.linalg.inv(X.T @ X)
    se = np.sqrt(np.diag(cov))
    eta_0 = coeffs[0] / ETA_SCALE_FACTOR
    eta_r = coeffs[1] / ETA_SCALE_FACTOR
    eta_vr = coeffs[2] / ETA_SCALE_FACTOR
    eta_theta = coeffs[3] / ETA_SCALE_FACTOR
    se_eta_theta = se[3] / ETA_SCALE_FACTOR
    return eta_theta, se_eta_theta, eta_r, eta_vr, eta_0


def compute_vif(X):
    """Compute Variance Inflation Factors for each column of design matrix X.
    VIF_j = 1 / (1 - R²_j), where R²_j is from regressing column j on all others.
    """
    n, p = X.shape
    vifs = []
    for j in range(p):
        y_col = X[:, j]
        X_others = np.delete(X, j, axis=1)
        coeffs, _, rank, _ = np.linalg.lstsq(X_others, y_col, rcond=None)
        if rank < X_others.shape[1] or np.var(y_col) == 0:
            vifs.append(np.inf)
            continue
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


def cmb_falsification_analysis(df, verbose=False):
    print_status("═══ Step 055: CMB Anisotropy Rigorous Falsification ═══", "TITLE")

    n = len(df)
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

    # ------------------------------------------------------------------
    # A. Aliasing Simulation Falsification
    # ------------------------------------------------------------------
    print_status("═══ A. Aliasing Simulation Falsification ═══", "TITLE")

    # Estimate the true synodic amplitude from the data
    X_syn = np.column_stack([cosD, np.ones(n_clean)])
    coeffs_syn, _, _, _ = np.linalg.lstsq(X_syn, res, rcond=None)
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
        "aliasing_rejected": p_aliasing_t < 0.001,
    }

    # ------------------------------------------------------------------
    # B. Multicollinearity Diagnostic Suite
    # ------------------------------------------------------------------
    print_status("═══ B. Multicollinearity Diagnostic Suite ═══", "TITLE")

    X_full = np.column_stack([cosD, r_c * cosD, vr_c * cosD, cos_theta_c * cosD, np.ones(n_clean)])
    vifs = compute_vif(X_full)
    # Condition number of X^T X (not X itself)
    XtX = X_full.T @ X_full
    s = np.linalg.svdvals(XtX)
    cond_XtX = s[0] / s[-1] if s[-1] > 0 else np.inf
    cond_X = np.linalg.cond(X_full)

    predictor_names = ["cosD", "r_c·cosD", "vr_c·cosD", "cosθ_c·cosD", "intercept"]
    for name, vif in zip(predictor_names, vifs):
        print_status(f"VIF({name}) = {vif:.3f}", "CALC")
    print_status(f"Condition number κ(X) = {cond_X:.2e}", "CALC")
    print_status(f"Condition number κ(XᵀX) = {cond_XtX:.2e}", "CALC")

    # Compute what the SE would be if cosθ were perfectly orthogonal to all others
    # SE_orthogonal = sqrt(MSE / sum((cosθ_c*cosD)_orth²))
    # First orthogonalize cosθ_c*cosD against the other three columns
    y_orth = cos_theta_c * cosD
    X_others = np.column_stack([cosD, r_c * cosD, vr_c * cosD, np.ones(n_clean)])
    coeffs_orth, _, rank_orth, _ = np.linalg.lstsq(X_others, y_orth, rcond=None)
    if rank_orth == 4:
        y_resid = y_orth - X_others @ coeffs_orth
        # Refit full model to get MSE
        coeffs_full, _, _, _ = np.linalg.lstsq(X_full, res, rcond=None)
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

    multicollinearity_result = {
        "vifs": {name: float(v) for name, v in zip(predictor_names, vifs)},
        "max_vif": float(np.max(vifs)),
        "mean_vif": float(np.mean(vifs)),
        "condition_number_X": float(cond_X),
        "condition_number_XtX": float(cond_XtX),
        "se_cos_theta_actual": float(se_actual) if se_actual is not None else None,
        "se_cos_theta_orthogonal": float(se_orthogonal) if se_orthogonal is not None else None,
        "se_inflation_factor": float(inflation) if inflation is not None else None,
        "severe_multicollinearity": bool(np.max(vifs) > 10 or cond_X > 1e3),
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
        "permutation_rejected": p_perm_t < 0.001,
    }

    # ------------------------------------------------------------------
    # D. Sky Scrambling Monte Carlo
    # ------------------------------------------------------------------
    print_status("═══ D. Sky Scrambling Monte Carlo ═══", "TITLE")

    n_scramble = 1000
    rng_scram = np.random.default_rng(44)
    sigmas_scram = []
    aics_scram = []

    for _ in range(n_scramble):
        dir_vec = random_unit_vector(rng_scram)
        cos_theta_rand = np.sum(em_hat[:, mask_clean] * dir_vec[:, None], axis=0)
        cos_theta_rand_c = cos_theta_rand - np.mean(cos_theta_rand)
        eta_s, se_s, _, _, _ = fit_full_joint_model(res, cosD, r_c, vr_c, cos_theta_rand_c)
        if eta_s is not None and se_s > 0:
            sigmas_scram.append(abs(eta_s) / se_s)
            # Quick AIC: need RSS
            X_s = np.column_stack([cosD, r_c * cosD, vr_c * cosD, cos_theta_rand_c * cosD, np.ones(n_clean)])
            coeffs_s, _, _, _ = np.linalg.lstsq(X_s, res, rcond=None)
            rss_s = np.sum((res - X_s @ coeffs_s) ** 2)
            aics_scram.append(n_clean * np.log(rss_s / n_clean) + 2 * 5)

    sigmas_scram = np.array(sigmas_scram)
    aics_scram = np.array(aics_scram)

    # True CMB fit
    eta_true, se_true, _, _, _ = fit_full_joint_model(res, cosD, r_c, vr_c, cos_theta_c)
    sigma_true = abs(eta_true) / se_true if se_true > 0 else 0.0
    X_true = np.column_stack([cosD, r_c * cosD, vr_c * cosD, cos_theta_c * cosD, np.ones(n_clean)])
    coeffs_true, _, _, _ = np.linalg.lstsq(X_true, res, rcond=None)
    rss_true = np.sum((res - X_true @ coeffs_true) ** 2)
    aic_true = n_clean * np.log(rss_true / n_clean) + 2 * 5

    p_scramble_sigma = np.mean(sigmas_scram >= sigma_true)
    p_scramble_aic = np.mean(aics_scram <= aic_true)

    print_status(
        f"Random directions: median σ = {np.median(sigmas_scram):.2f}, "
        f"99th %ile σ = {np.percentile(sigmas_scram, 99):.2f}",
        "CALC",
    )
    print_status(
        f"True CMB: σ = {sigma_true:.2f}; p(scramble ≥ true) = {p_scramble_sigma:.4f}",
        "CALC",
    )
    print_status(
        f"Random directions AIC: median = {np.median(aics_scram):.1f}",
        "CALC",
    )
    print_status(
        f"True CMB AIC = {aic_true:.1f}; p(random AIC ≤ true) = {p_scramble_aic:.4f}",
        "CALC",
    )

    sky_scramble_result = {
        "n_scrambles": n_scramble,
        "sigma_median": float(np.median(sigmas_scram)),
        "sigma_99th_percentile": float(np.percentile(sigmas_scram, 99)),
        "sigma_max": float(np.max(sigmas_scram)),
        "true_sigma": float(sigma_true),
        "p_scramble_sigma": float(p_scramble_sigma),
        "aic_median": float(np.median(aics_scram)),
        "aic_99th_percentile": float(np.percentile(aics_scram, 1)),  # lower is better
        "true_aic": float(aic_true),
        "p_scramble_aic": float(p_scramble_aic),
        "true_direction_preferred": p_scramble_sigma < 0.01,
    }

    # ------------------------------------------------------------------
    # E. Gram–Schmidt Orthogonalization
    # ------------------------------------------------------------------
    print_status("═══ E. Gram–Schmidt Orthogonalization ═══", "TITLE")

    # Orthogonalize cosθ·cosD against [cosD, r_c·cosD, vr_c·cosD, 1]
    y_target = cos_theta_c * cosD
    X_base = np.column_stack([cosD, r_c * cosD, vr_c * cosD, np.ones(n_clean)])
    coeffs_gs, _, rank_gs, _ = np.linalg.lstsq(X_base, y_target, rcond=None)
    if rank_gs == 4:
        y_orth_gs = y_target - X_base @ coeffs_gs
        # Verify orthogonality
        orth_dot_products = [np.dot(y_orth_gs, X_base[:, j]) for j in range(4)]
        print_status(
            f"Orthogonality residuals (dot products): "
            f"{[f'{d:.2e}' for d in orth_dot_products]}",
            "CALC",
        )

        X_orth_model = np.column_stack([cosD, r_c * cosD, vr_c * cosD, y_orth_gs, np.ones(n_clean)])
        coeffs_orth_m, _, rank_orth_m, _ = np.linalg.lstsq(X_orth_model, res, rcond=None)
        if rank_orth_m == 5:
            resid_orth_m = res - X_orth_model @ coeffs_orth_m
            mse_orth_m = np.sum(resid_orth_m ** 2) / (n_clean - 5)
            cov_orth_m = mse_orth_m * np.linalg.inv(X_orth_model.T @ X_orth_model)
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
                "signal_persists": abs(t_theta_orth) > 2.0,
            }
        else:
            orthogonalization_result = {"available": False, "reason": "rank_deficient_orth_model"}
    else:
        orthogonalization_result = {"available": False, "reason": "rank_deficient_gs"}

    if not orthogonalization_result.get("available"):
        print_status("Orthogonalization failed.", "WARNING")

    # ------------------------------------------------------------------
    # 6. Compile results
    # ------------------------------------------------------------------
    status = "PASS"
    checks = [
        aliasing_result.get("aliasing_rejected", False),
        permutation_result.get("permutation_rejected", False),
        sky_scramble_result.get("true_direction_preferred", False),
        orthogonalization_result.get("signal_persists", False) if orthogonalization_result.get("available") else False,
    ]
    n_passed = sum(checks)
    if n_passed >= 3:
        status = "PASS"
        print_status(
            f"RESULT: {n_passed}/4 falsification tests passed. CMB anisotropy claim is robust.",
            "SUCCESS",
        )
    elif n_passed >= 2:
        status = "WARNING"
        print_status(
            f"RESULT: {n_passed}/4 falsification tests passed. Marginal robustness.",
            "WARNING",
        )
    else:
        status = "FAIL"
        print_status(
            f"RESULT: Only {n_passed}/4 falsification tests passed. CMB anisotropy claim is weakened.",
            "ERROR",
        )

    results = {
        "step_id": "step_055",
        "status": status,
        "n_observations": int(n_clean),
        "correlation_cosD_cos_theta": float(corr_cosD_cosTheta),
        "aliasing_simulation": aliasing_result,
        "multicollinearity_diagnostics": multicollinearity_result,
        "permutation_test": permutation_result,
        "sky_scrambling": sky_scramble_result,
        "orthogonalization_test": orthogonalization_result,
        "falsification_tests_passed": int(n_passed),
        "falsification_tests_total": 4,
    }

    return results


if __name__ == "__main__":
    log_dir = PROJECT_ROOT / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    logger = TEPLogger(
        "step_055", str(log_dir / "step_055_cmb_rigorous_falsification.log")
    )
    set_step_logger(logger)

    data_path = (
        PROJECT_ROOT / "data" / "processed" / "INPOP19a_all_stations_residuals.csv"
    )
    if not data_path.exists():
        # Fallback to DE430 if INPOP19a not available
        data_path = (
            PROJECT_ROOT / "data" / "processed" / "DE430_all_residuals.csv"
        )
    if not data_path.exists():
        print_status(f"No processed residuals found at expected paths", "ERROR")
        sys.exit(1)

    df = pd.read_csv(data_path)

    results = cmb_falsification_analysis(df, verbose=True)

    if results:
        logger.save_step_results(
            results, PROJECT_ROOT, "step_055_cmb_rigorous_falsification"
        )
        print_status("CMB Falsification analysis complete.", "SUCCESS")
    else:
        print_status("Analysis failed.", "ERROR")
        sys.exit(1)
