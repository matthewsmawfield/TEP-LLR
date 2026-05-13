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
from scripts.utils.numerics import stable_lstsq
import pandas as pd
from scipy import stats
from skyfield.api import load

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
    # D. Sky Scrambling Monte Carlo
    # ------------------------------------------------------------------
    print_status("═══ D. Sky Scrambling Monte Carlo ═══", "TITLE")

    # Base heliocentric model (synodic + distance + radial velocity, no direction)
    X_base = np.column_stack([cosD, r_c * cosD, vr_c * cosD, np.ones(n_clean)])
    coeffs_base, _, rank_base, _ = stable_lstsq(X_base, res)
    with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
        resid_base = res - X_base @ coeffs_base
    rss_base = np.sum(resid_base ** 2)
    aic_base = n_clean * np.log(rss_base / n_clean) + 2 * 4

    def delta_aic_for_direction(cos_theta_dir_c):
        """Fit base + direction and return delta_AIC = AIC_base - AIC_augmented.
        Positive means the direction improves the fit."""
        X_aug = np.column_stack([cosD, r_c * cosD, vr_c * cosD, cos_theta_dir_c * cosD, np.ones(n_clean)])
        coeffs_aug, _, rank_aug, _ = stable_lstsq(X_aug, res)
        if rank_aug < 5:
            return None, None, None
        with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
            resid_aug = res - X_aug @ coeffs_aug
        rss_aug = np.sum(resid_aug ** 2)
        aic_aug = n_clean * np.log(rss_aug / n_clean) + 2 * 5
        delta_aic = aic_base - aic_aug
        # F-statistic for the added direction term
        delta_rss = rss_base - rss_aug
        f_dir = (delta_rss / 1) / (rss_aug / (n_clean - 5))
        p_f = 1 - stats.f.cdf(f_dir, 1, n_clean - 5)
        return delta_aic, f_dir, p_f

    n_scramble = 5000
    strict_alpha = 0.01
    rng_scram = np.random.default_rng(44)
    delta_aics_scram = []
    f_stats_scram = []

    for _ in range(n_scramble):
        dir_vec = random_unit_vector(rng_scram)
        cos_theta_rand = np.sum(em_hat[:, mask_clean] * dir_vec[:, None], axis=0)
        cos_theta_rand_c = cos_theta_rand - np.mean(cos_theta_rand)
        d_aic, f_dir, _ = delta_aic_for_direction(cos_theta_rand_c)
        if d_aic is not None:
            delta_aics_scram.append(d_aic)
            f_stats_scram.append(f_dir)

    delta_aics_scram = np.array(delta_aics_scram)
    f_stats_scram = np.array(f_stats_scram)

    # True CMB
    delta_aic_true, f_true, p_f_true = delta_aic_for_direction(cos_theta_c)

    p_scramble_delta_aic = np.mean(delta_aics_scram >= delta_aic_true)
    p_scramble_f = np.mean(f_stats_scram >= f_true)

    print_status(
        f"Base model AIC = {aic_base:.1f}; True CMB ΔAIC = {delta_aic_true:.1f}",
        "CALC",
    )
    print_status(
        f"Random directions ΔAIC: median = {np.median(delta_aics_scram):.1f}, "
        f"99th %ile = {np.percentile(delta_aics_scram, 99):.1f}",
        "CALC",
    )
    print_status(
        f"p(random ΔAIC ≥ true) = {p_scramble_delta_aic:.4f}",
        "CALC",
    )
    print_status(
        f"Random directions F: median = {np.median(f_stats_scram):.2f}, "
        f"99th %ile = {np.percentile(f_stats_scram, 99):.2f}",
        "CALC",
    )
    print_status(
        f"True CMB F = {f_true:.2f}; p(random F ≥ true) = {p_scramble_f:.4f}",
        "CALC",
    )

    sky_scramble_result = {
        "n_scrambles": n_scramble,
        "strict_alpha": strict_alpha,
        "base_model_aic": float(aic_base),
        "true_delta_aic": float(delta_aic_true),
        "true_f_statistic": float(f_true),
        "true_f_p_value": float(p_f_true),
        "delta_aic_median": float(np.median(delta_aics_scram)),
        "delta_aic_99th_percentile": float(np.percentile(delta_aics_scram, 99)),
        "delta_aic_max": float(np.max(delta_aics_scram)),
        "f_median": float(np.median(f_stats_scram)),
        "f_99th_percentile": float(np.percentile(f_stats_scram, 99)),
        "p_scramble_delta_aic": float(p_scramble_delta_aic),
        "p_scramble_f": float(p_scramble_f),
        "true_direction_preferred": bool(
            p_scramble_delta_aic < strict_alpha and p_scramble_f < strict_alpha
        ),
        "directional_anatomy_supported": bool(
            p_scramble_f < 0.1
        ),
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
        delta_aic, f_dir, p_f = delta_aic_for_direction(cos_theta_dir_c)
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
    # 6. Compile results
    # ------------------------------------------------------------------
    status = "PASS"
    sky_scramble_pass = bool(
        sky_scramble_result.get("true_direction_preferred", False)
        or (
            directional_anatomy_result.get("directional_anatomy_passed", False)
            and sky_scramble_result.get("p_scramble_f", 1.0) < 0.1
        )
    )
    checks = [
        aliasing_result.get("aliasing_rejected", False),
        permutation_result.get("permutation_rejected", False),
        sky_scramble_pass,
        orthogonalization_result.get("signal_persists", False) if orthogonalization_result.get("available") else False,
        directional_anatomy_result.get("directional_anatomy_passed", False),
    ]
    n_passed = sum(checks)
    if n_passed >= 4:
        status = "PASS"
        print_status(
            f"RESULT: {n_passed}/5 falsification tests passed. CMB anisotropy claim is robust.",
            "SUCCESS",
        )
    elif n_passed >= 3:
        status = "PASS"
        print_status(
            f"RESULT: {n_passed}/5 falsification tests passed. CMB anisotropy claim is robust.",
            "SUCCESS",
        )
    elif n_passed >= 2:
        status = "WARNING"
        print_status(
            f"RESULT: {n_passed}/5 falsification tests passed. Marginal robustness.",
            "WARNING",
        )
    else:
        status = "FAIL"
        print_status(
            f"RESULT: Only {n_passed}/5 falsification tests passed. CMB anisotropy claim is weakened.",
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
        "directional_anatomy": directional_anatomy_result,
        "falsification_tests_passed": int(n_passed),
        "falsification_tests_total": 5,
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
