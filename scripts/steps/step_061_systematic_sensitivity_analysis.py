#!/usr/bin/env python3
"""
Step 061: Systematic Amplitude Sensitivity Analysis

Quantifies whether each known systematic could produce the observed eta.
For each systematic source from Step 044, we:

  1. Compute the "required amplitude": the amplitude a systematic would need
     to have to fully explain the observed pooled eta.

  2. Compute the "exclusion ratio": required_amplitude / known_amplitude.
     If > 1, the known systematic is too small to explain the signal.

  3. Monte Carlo falsification: Generate N_MC synthetic datasets where the
     systematic is the ONLY signal (no true eta), at the known amplitude.
     Fit the full-systematic model (including cosD). Count how often
     |eta| >= observed_eta.

Extended systematics beyond the known list:

  A. Adversarial nuisance search — learn dominant nuisance directions from
     residuals via PCA on a rich ephemeris-like basis and via a non-parametric
     GP elongation surface; test whether they absorb cos(D).

  B. Era x station x lunation interaction grid — per-cell eta with multiplicity
     control registered for Step 042 (headline-claim coverage).

  C. Blind year hold-out — hold out 20% of calendar years at random; fit
     nuisances on train, estimate eta on test only.
"""

import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

import json
import argparse
import hashlib
import warnings
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.decomposition import PCA
from sklearn.exceptions import ConvergenceWarning
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, WhiteKernel, ConstantKernel as C

from scripts.steps.step_065_high_dimensional_absorption_test import (
    build_ephemeris_like_basis,
    SYNODIC_MONTH_DAYS,
)
from scripts.utils.statistical_utils import robust_regression, detect_outliers_sigma
from scripts.utils.llr_constants import ETA_SCALE_FACTOR, TEMPORAL_DRIFT_ERA_SPLIT_YEAR
from scripts.utils.config import get_config
from scripts.utils.logger import TEPLogger, set_step_logger, print_status

TEP_CONFIG = get_config()
BLIND_HOLDOUT_YEAR_FRACTION = 0.20
BLIND_HOLDOUT_N_SPLITS = 20
INTERACTION_GRID_MIN_N = 80
ERA_TIER_MID_YEAR = 1990.0
PCA_K_VALUES = (1, 3, 5, 10, 20)
GP_N_EL_BINS = 24
GP_N_TIME_BINS = 30
GP_MIN_PER_BIN = 8


def load_cleaned_inpop_df():
    """Load canonical INPOP19a residuals with 6σ cleaning."""
    input_path = PROJECT_ROOT / "data" / "processed" / "INPOP19a_all_stations_residuals.csv"
    if not input_path.exists():
        raise FileNotFoundError(f"Input file not found: {input_path}")
    df_raw = pd.read_csv(input_path)
    outlier_mask = detect_outliers_sigma(df_raw["residual_m"].values, sigma_threshold=6.0)
    df = df_raw.loc[~outlier_mask].copy()
    return df_raw, df, outlier_mask


def build_full_systematic_design(df):
    """Build full-systematic design matrix (includes cosD as first column)."""
    year = df["date_julian"].values / 365.25
    month = df["date_julian"].values / 27.32
    elongation = df["elongation_rad"].values
    return np.column_stack([
        np.cos(elongation),
        np.cos(2.0 * elongation),
        np.sin(2.0 * np.pi * month),
        np.cos(2.0 * np.pi * month),
        np.sin(2.0 * np.pi * year),
        np.cos(2.0 * np.pi * year),
        np.ones(len(df)),
    ])


def build_nuisance_only_design(df):
    """Full systematic nuisance terms without cos(D)."""
    return build_full_systematic_design(df)[:, 1:]


def fit_eta_from_design(y, design, scale_errors_by_birge=False):
    """Fit design matrix and return eta (coefficient 0 / ETA_SCALE_FACTOR)."""
    fit = robust_regression(y, design, scale_errors_by_birge=scale_errors_by_birge)
    if fit.get("status") == "SINGULAR" or not np.isfinite(fit["coefficients"][0]):
        return None
    eta = float(fit["coefficients"][0] / ETA_SCALE_FACTOR)
    eta_err = float(fit["errors"][0] / ETA_SCALE_FACTOR)
    snr = float(abs(eta) / max(eta_err, 1e-20))
    p_two_sided = float(2 * (1 - stats.norm.cdf(abs(eta) / max(eta_err, 1e-20))))
    return {
        "eta": eta,
        "eta_err": eta_err,
        "snr": snr,
        "p_two_sided": p_two_sided,
        "n_obs": int(fit["n_obs"]),
        "amplitude_m": float(fit["coefficients"][0]),
    }


def compute_station_noise_params(df):
    """Compute per-station RMS noise and sample sizes."""
    stations = df["station"].unique()
    params = {}
    for s in stations:
        mask = df["station"].values == s
        res = df["residual_m"].values[mask]
        if len(res) > 10:
            params[s] = {
                "rms": float(np.std(res, ddof=1)),
                "n": int(len(res)),
                "mean": float(np.mean(res)),
            }
    return params


def simulate_systematic_only(df, systematic_type, amplitude_m, n_mc=2000, seed=61):
    """
    Simulate data where the systematic is the ONLY signal.
    Fit full-systematic model, extract eta distribution.
    """
    rng = np.random.default_rng(seed)
    design = build_full_systematic_design(df)
    noise_params = compute_station_noise_params(df)

    elongation = df["elongation_rad"].values
    year = df["date_julian"].values / 365.25

    stations = df["station"].values
    if systematic_type == "ephemeris":
        signal = np.full(len(df), amplitude_m)
    elif systematic_type == "atmospheric":
        signal = amplitude_m * np.sin(2.0 * np.pi * year)
    elif systematic_type == "instrumental":
        signal = np.zeros(len(df))
        for s, params in noise_params.items():
            mask = stations == s
            signal[mask] = params["mean"] * (amplitude_m / max(abs(params["mean"]), 1e-10))
    elif systematic_type == "tidal":
        signal = amplitude_m * np.cos(2.0 * elongation)
    elif systematic_type == "thermal":
        day_frac = df["date_julian"].values % 1.0
        signal = amplitude_m * np.sin(2.0 * np.pi * day_frac)
    else:
        raise ValueError(f"Unknown systematic type: {systematic_type}")

    signal = np.nan_to_num(signal, nan=0.0)

    rms_vec = np.zeros(len(df), dtype=np.float64)
    for s, params in noise_params.items():
        mask = stations == s
        rms_vec[mask] = params["rms"]

    etas = []
    for _ in range(n_mc):
        noise = rng.normal(0.0, rms_vec)
        y = signal + noise
        outlier_mask = detect_outliers_sigma(y, sigma_threshold=6.0)
        kept = ~outlier_mask
        if kept.sum() < 100:
            etas.append(0.0)
            continue
        try:
            fit = robust_regression(y[kept], design[kept], scale_errors_by_birge=False)
            eta = fit["coefficients"][0] / ETA_SCALE_FACTOR
            etas.append(float(eta))
        except Exception:
            etas.append(0.0)

    return np.array(etas)


def _standardize_columns(X: np.ndarray) -> np.ndarray:
    """Column-wise z-score; leave near-constant columns at zero."""
    Xs = np.zeros_like(X, dtype=float)
    for j in range(X.shape[1]):
        col = X[:, j]
        std = float(np.std(col, ddof=1))
        if std < 1e-12:
            continue
        Xs[:, j] = (col - np.mean(col)) / std
    return Xs


def adversarial_pca_absorption(df):
    """
    Learn nuisance directions from residuals via PCA on a rich basis,
    then test whether they absorb cos(D).

    PCs are ranked by |correlation| with nuisance-stripped residuals so the
    search targets directions that actually explain residual structure.
    """
    y = df["residual_m"].values
    cos_d = np.cos(df["elongation_rad"].values)

    baseline = fit_eta_from_design(
        y, np.column_stack([cos_d, np.ones(len(df))]), scale_errors_by_birge=False
    )
    full_sys = fit_eta_from_design(y, build_full_systematic_design(df), scale_errors_by_birge=False)

    X_basis, basis_names, cond = build_ephemeris_like_basis(df)
    if not np.isfinite(cond):
        raise RuntimeError(f"Adversarial PCA basis has non-finite condition number: {cond}")

    nuis_fit = robust_regression(y, build_nuisance_only_design(df), scale_errors_by_birge=False)
    if nuis_fit.get("status") == "SINGULAR":
        raise RuntimeError("Nuisance-only fit singular in adversarial PCA block.")
    with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
        residual_after_nuisance = y - build_nuisance_only_design(df) @ nuis_fit["coefficients"]

    X_std = _standardize_columns(X_basis)
    n_components = min(20, X_std.shape[1], X_std.shape[0] - 1)
    pca = PCA(n_components=n_components, svd_solver="full")
    scores_full = pca.fit_transform(X_std)

    pc_corr = np.array([
        float(np.corrcoef(scores_full[:, j], residual_after_nuisance)[0, 1])
        if np.std(scores_full[:, j]) > 0
        else 0.0
        for j in range(scores_full.shape[1])
    ])
    pc_order = np.argsort(-np.abs(pc_corr))
    scores = scores_full[:, pc_order]

    design_all_pcs = np.column_stack([cos_d, scores, np.ones(len(df))])
    joint_all = fit_eta_from_design(y, design_all_pcs, scale_errors_by_birge=False)

    by_k = {}
    for k in PCA_K_VALUES:
        if k > n_components:
            continue
        design_k = np.column_stack([cos_d, scores[:, :k], np.ones(len(df))])
        fit_k = fit_eta_from_design(y, design_k, scale_errors_by_birge=False)
        if fit_k is None or baseline is None:
            continue
        absorption = 1.0 - abs(fit_k["eta"]) / max(abs(baseline["eta"]), 1e-20)
        by_k[str(k)] = {
            **fit_k,
            "absorption_fraction": float(absorption),
            "survival_fraction": float(1.0 - absorption),
            "variance_explained_cumulative": float(
                np.sum(pca.explained_variance_ratio_[pc_order[:k]])
            ),
        }

    return {
        "basis_n_terms": int(X_basis.shape[1]),
        "basis_condition_number": float(cond),
        "pca_n_components_used": int(n_components),
        "pc_ranked_by_residual_correlation": True,
        "top_pc_residual_correlations": [
            float(pc_corr[pc_order[i]]) for i in range(min(5, n_components))
        ],
        "baseline_cosd_only": baseline,
        "full_systematic": full_sys,
        "joint_cosd_plus_all_pcs": joint_all,
        "pca_by_n_components": by_k,
        "cosd_survives_rich_pc_basis": bool(
            joint_all is not None
            and baseline is not None
            and abs(joint_all["eta"]) >= 0.5 * abs(baseline["eta"])
        ),
    }


def compute_2d_bin_means(df, n_el=GP_N_EL_BINS, n_time=GP_N_TIME_BINS, min_per_bin=GP_MIN_PER_BIN):
    """Precision-weighted mean residuals on an elongation × time grid (Step 060-style)."""
    t_years = df["date_julian"].values / 365.25
    t_norm = t_years - float(np.min(t_years))
    elong = df["elongation_rad"].values
    residuals = df["residual_m"].values

    el_edges = np.linspace(0.0, np.pi, n_el + 1)
    t_edges = np.linspace(0.0, float(np.max(t_norm)), n_time + 1)

    elong_centers = []
    time_centers = []
    means = []
    errs = []
    ns = []

    for i in range(n_el):
        for j in range(n_time):
            mask = (
                (elong >= el_edges[i])
                & (elong < el_edges[i + 1])
                & (t_norm >= t_edges[j])
                & (t_norm < t_edges[j + 1])
            )
            n_bin = int(mask.sum())
            if n_bin < min_per_bin:
                continue
            res = residuals[mask]
            elong_centers.append(float(0.5 * (el_edges[i] + el_edges[i + 1])))
            time_centers.append(float(0.5 * (t_edges[j] + t_edges[j + 1])))
            means.append(float(np.mean(res)))
            errs.append(max(float(np.std(res, ddof=1) / np.sqrt(n_bin)), 1e-10))
            ns.append(n_bin)

    if len(means) < 20:
        raise RuntimeError(
            f"2D GP binning produced only {len(means)} valid cells; need >= 20."
        )

    return {
        "elongation_centers": elong_centers,
        "time_centers": time_centers,
        "means": means,
        "errs": errs,
        "ns": ns,
        "n_valid_bins": len(means),
    }


def adversarial_gp_absorption(df, seed=61):
    """
    Fit a 2D GP nuisance surface on elongation × time bin means (Step 060-style),
    subtract it, and test whether cos(D) is absorbed.
    """
    y = df["residual_m"].values
    cos_d = np.cos(df["elongation_rad"].values)
    t_years = df["date_julian"].values / 365.25
    t_norm = t_years - float(np.min(t_years))
    elong = df["elongation_rad"].values

    baseline = fit_eta_from_design(
        y, np.column_stack([cos_d, np.ones(len(df))]), scale_errors_by_birge=False
    )

    bins = compute_2d_bin_means(df)
    X_train = np.column_stack([bins["elongation_centers"], bins["time_centers"]])
    y_train = np.array(bins["means"])
    dy = np.array(bins["errs"])

    kernel = (
        C(1.0, (1e-2, 1e3))
        * RBF(
            length_scale=[0.6, 8.0],
            length_scale_bounds=(1e-2, 50.0),
        )
        + WhiteKernel(noise_level=1e-4, noise_level_bounds=(1e-8, 1.0))
    )
    gp = GaussianProcessRegressor(
        kernel=kernel,
        n_restarts_optimizer=5,
        random_state=seed,
        normalize_y=False,
        alpha=dy**2,
    )
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        warnings.simplefilter("ignore", category=ConvergenceWarning)
        gp.fit(X_train, y_train)
        X_pred = np.column_stack([elong, t_norm])
        gp_mean = gp.predict(X_pred)
    if not np.all(np.isfinite(gp_mean)):
        raise RuntimeError("2D GP prediction returned non-finite values.")

    y_adj = y - gp_mean
    post_gp = fit_eta_from_design(
        y_adj, np.column_stack([cos_d, np.ones(len(df))]), scale_errors_by_birge=False
    )

    joint = fit_eta_from_design(
        y, np.column_stack([cos_d, gp_mean, np.ones(len(df))]), scale_errors_by_birge=False
    )

    absorption = None
    if post_gp is not None and baseline is not None:
        absorption = 1.0 - abs(post_gp["eta"]) / max(abs(baseline["eta"]), 1e-20)

    return {
        "method": "2D binned elongation x time GP (Step 060-style bin means)",
        "n_valid_bins": int(bins["n_valid_bins"]),
        "gp_train_n": int(bins["n_valid_bins"]),
        "kernel": str(gp.kernel_),
        "log_marginal_likelihood": float(gp.log_marginal_likelihood_value_),
        "baseline_cosd_only": baseline,
        "cosd_after_gp_subtraction": post_gp,
        "joint_cosd_plus_gp_mean": joint,
        "absorption_fraction_gp_subtract": absorption,
        "survival_fraction_gp_subtract": float(1.0 - absorption) if absorption is not None else None,
        "cosd_survives_gp_nuisance": bool(
            post_gp is not None
            and baseline is not None
            and abs(post_gp["eta"]) >= 0.5 * abs(baseline["eta"])
        ),
    }


def _era_label(year_value):
    """Three-tier era partition for interaction grid coverage."""
    if year_value < ERA_TIER_MID_YEAR:
        return "pre_1990"
    if year_value < TEMPORAL_DRIFT_ERA_SPLIT_YEAR:
        return "era_1990s"
    return "post_2000"


def _lunation_quartile(date_julian):
    phase = (np.asarray(date_julian) % SYNODIC_MONTH_DAYS) / SYNODIC_MONTH_DAYS
    q = np.floor(phase * 4).astype(int)
    return np.clip(q, 0, 3)


def era_station_lunation_grid(df):
    """
    Fit full-systematic eta in each era x station x lunation-quartile cell.
    Returns per-cell results and pooled interaction model for headline claims.
    """
    df = df.copy()
    df["era"] = df["date_julian_year"].map(_era_label)
    df["lunation_q"] = _lunation_quartile(df["date_julian"].values)

    cells = []
    stations = sorted(df["station"].unique())
    eras = ["pre_1990", "era_1990s", "post_2000"]
    lunation_labels = [f"Q{i}" for i in range(4)]

    for era in eras:
        for station in stations:
            for lq in range(4):
                mask = (
                    (df["era"] == era)
                    & (df["station"] == station)
                    & (df["lunation_q"] == lq)
                )
                sub = df.loc[mask]
                n_cell = len(sub)
                entry = {
                    "era": era,
                    "station": station,
                    "lunation_quartile": lunation_labels[lq],
                    "n_obs": int(n_cell),
                    "sufficient_n": bool(n_cell >= INTERACTION_GRID_MIN_N),
                    "fit": None,
                }
                if n_cell >= INTERACTION_GRID_MIN_N:
                    fit = fit_eta_from_design(
                        sub["residual_m"].values,
                        build_full_systematic_design(sub),
                        scale_errors_by_birge=False,
                    )
                    entry["fit"] = fit
                cells.append(entry)

    sufficient = [c for c in cells if c["fit"] is not None]
    n_negative = sum(1 for c in sufficient if c["fit"]["eta"] < 0)
    n_significant = sum(1 for c in sufficient if c["fit"]["p_two_sided"] < 0.05)

    # Pooled interaction: cosD modulated by era, station, lunation phase
    elong = df["elongation_rad"].values
    cos_d = np.cos(elong)
    sin_lun = np.sin(2.0 * np.pi * (df["date_julian"].values % SYNODIC_MONTH_DAYS) / SYNODIC_MONTH_DAYS)
    cos_lun = np.cos(2.0 * np.pi * (df["date_julian"].values % SYNODIC_MONTH_DAYS) / SYNODIC_MONTH_DAYS)
    era_1990s = (df["era"] == "era_1990s").astype(float).values
    era_post = (df["era"] == "post_2000").astype(float).values
    station_refs = stations[1:]
    station_mats = [
        (df["station"].values == s).astype(float) for s in station_refs
    ]

    interact_cols = [
        cos_d,
        cos_d * era_1990s,
        cos_d * era_post,
        cos_d * sin_lun,
        cos_d * cos_lun,
    ]
    interact_names = [
        "cosD",
        "cosD_x_era_1990s",
        "cosD_x_era_post",
        "cosD_x_sin_lunation",
        "cosD_x_cos_lunation",
    ]
    for s, mat in zip(station_refs, station_mats):
        interact_cols.append(cos_d * mat)
        interact_names.append(f"cosD_x_station_{s}")

    nuisance = build_nuisance_only_design(df)
    pooled_design = np.column_stack(interact_cols + [nuisance])
    pooled_fit = robust_regression(df["residual_m"].values, pooled_design, scale_errors_by_birge=False)
    pooled_eta = float(pooled_fit["coefficients"][0] / ETA_SCALE_FACTOR)
    pooled_eta_err = float(pooled_fit["errors"][0] / ETA_SCALE_FACTOR)

    headline_claims = {
        "pooled_full_systematic": {
            "claim": "Primary pooled synodic coefficient (full systematic)",
            "source_cell": "all_data",
            "eta": pooled_eta,
            "eta_err": pooled_eta_err,
            "snr": float(abs(pooled_eta) / max(pooled_eta_err, 1e-20)),
            "p_two_sided": float(2 * (1 - stats.norm.cdf(abs(pooled_eta) / max(pooled_eta_err, 1e-20)))),
        },
        "per_era_pooled": {},
        "per_station_pooled": {},
    }
    for era in eras:
        sub = df[df["era"] == era]
        fit = fit_eta_from_design(
            sub["residual_m"].values, build_full_systematic_design(sub), scale_errors_by_birge=False
        )
        if fit:
            headline_claims["per_era_pooled"][era] = fit

    for station in stations:
        sub = df[df["station"] == station]
        fit = fit_eta_from_design(
            sub["residual_m"].values, build_full_systematic_design(sub), scale_errors_by_birge=False
        )
        if fit:
            headline_claims["per_station_pooled"][station] = fit

    # Multiplicity summary (Step 042 ingests cell tests)
    cell_tests = [
        {
            "name": (
                f"Grid {c['era']} / {c['station']} / lunation {c['lunation_quartile']}"
            ),
            "category": "interaction_grid",
            "sigma": c["fit"]["snr"],
            "p": c["fit"]["p_two_sided"],
            "n_obs": c["n_obs"],
            "eta": c["fit"]["eta"],
        }
        for c in sufficient
    ]

    return {
        "min_cell_n": INTERACTION_GRID_MIN_N,
        "n_cells_total": len(cells),
        "n_cells_fitted": len(sufficient),
        "n_cells_negative_eta": n_negative,
        "n_cells_significant_05": n_significant,
        "fraction_negative_eta": float(n_negative / max(len(sufficient), 1)),
        "cells": cells,
        "headline_claims": headline_claims,
        "pooled_interaction_model": {
            "eta_cosD": pooled_eta,
            "eta_cosD_err": pooled_eta_err,
            "interaction_term_names": interact_names,
        },
        "cell_tests_for_step_042": cell_tests,
    }


def blind_year_holdout(df, seed=61, n_splits=BLIND_HOLDOUT_N_SPLITS):
    """
    Hold out a random 20% of calendar years; fit nuisances on train,
    estimate eta from cos(D) on test only.
    """
    years = np.sort(df["date_julian_year"].astype(int).unique())
    n_holdout = max(1, int(round(len(years) * BLIND_HOLDOUT_YEAR_FRACTION)))
    rng = np.random.default_rng(seed)

    cos_d = np.cos(df["elongation_rad"].values)
    y = df["residual_m"].values
    year_vals = df["date_julian_year"].astype(int).values

    splits = []
    for rep in range(n_splits):
        perm = rng.permutation(years)
        test_years = set(perm[:n_holdout].tolist())
        train_mask = np.array([yv not in test_years for yv in year_vals])
        test_mask = ~train_mask

        if train_mask.sum() < 500 or test_mask.sum() < 100:
            continue

        X_nuisance_train = build_nuisance_only_design(df.loc[train_mask])
        y_train = y[train_mask]
        nuis_fit = robust_regression(y_train, X_nuisance_train, scale_errors_by_birge=False)
        if nuis_fit.get("status") == "SINGULAR":
            continue

        beta_nuis = nuis_fit["coefficients"]
        X_nuisance_test = build_nuisance_only_design(df.loc[test_mask])
        y_test = y[test_mask]
        cos_d_test = cos_d[test_mask]

        with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
            y_test_adj = y_test - X_nuisance_test @ beta_nuis
        if not np.all(np.isfinite(y_test_adj)):
            continue
        blind_test = fit_eta_from_design(
            y_test_adj,
            np.column_stack([cos_d_test, np.ones(test_mask.sum())]),
            scale_errors_by_birge=False,
        )

        # In-sample train eta (not blind — for contrast only)
        cos_d_train = cos_d[train_mask]
        train_eta = fit_eta_from_design(
            y_train,
            np.column_stack([cos_d_train, np.ones(train_mask.sum())]),
            scale_errors_by_birge=False,
        )

        splits.append({
            "replicate": rep,
            "n_train_years": int(len(years) - n_holdout),
            "n_test_years": int(n_holdout),
            "n_train_obs": int(train_mask.sum()),
            "n_test_obs": int(test_mask.sum()),
            "test_years": sorted(test_years),
            "train_eta_cosd_only": train_eta,
            "blind_test_eta": blind_test,
        })

    if not splits:
        raise RuntimeError("Blind year hold-out produced no valid splits.")

    blind_etas = np.array([s["blind_test_eta"]["eta"] for s in splits if s["blind_test_eta"]])
    blind_errs = np.array([s["blind_test_eta"]["eta_err"] for s in splits if s["blind_test_eta"]])
    blind_snrs = np.array([s["blind_test_eta"]["snr"] for s in splits if s["blind_test_eta"]])

    combined_eta = float(np.average(blind_etas, weights=1.0 / blind_errs**2))
    combined_err = float(1.0 / np.sqrt(np.sum(1.0 / blind_errs**2)))

    return {
        "holdout_year_fraction": BLIND_HOLDOUT_YEAR_FRACTION,
        "n_calendar_years": int(len(years)),
        "n_holdout_years_per_split": int(n_holdout),
        "n_replicates": len(splits),
        "blind_test_eta_mean": float(np.mean(blind_etas)),
        "blind_test_eta_std": float(np.std(blind_etas, ddof=1)),
        "blind_test_snr_mean": float(np.mean(blind_snrs)),
        "blind_test_snr_min": float(np.min(blind_snrs)),
        "inverse_variance_combined": {
            "eta": combined_eta,
            "eta_err": combined_err,
            "snr": float(abs(combined_eta) / max(combined_err, 1e-20)),
            "p_two_sided": float(2 * (1 - stats.norm.cdf(abs(combined_eta) / max(combined_err, 1e-20)))),
        },
        "n_splits_negative_eta": int(np.sum(blind_etas < 0)),
        "n_splits_positive_eta": int(np.sum(blind_etas > 0)),
        "all_splits_negative_eta": bool(np.all(blind_etas < 0)),
        "splits": splits,
        "cell_tests_for_step_042": [
            {
                "name": "Blind 20% year hold-out (inverse-variance combined)",
                "category": "validation",
                "sigma": float(abs(combined_eta) / max(combined_err, 1e-20)),
                "p": float(2 * (1 - stats.norm.cdf(abs(combined_eta) / max(combined_err, 1e-20)))),
                "n_obs": int(sum(s["n_test_obs"] for s in splits) / len(splits)),
                "eta": combined_eta,
            }
        ],
    }


def main():
    parser = argparse.ArgumentParser(description="Step 061: Systematic Sensitivity Analysis")
    parser.parse_args()

    log_dir = PROJECT_ROOT / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    logger = TEPLogger(
        "step_061", str(log_dir / "step_061_systematic_sensitivity_analysis.log")
    )
    set_step_logger(logger)

    print_status("Starting Step 061: Systematic Amplitude Sensitivity Analysis", "TITLE")

    df_raw, df, outlier_mask = load_cleaned_inpop_df()
    print_status(f"Loaded data: {len(df_raw):,} observations", "INFO")
    print_status(
        f"Canonical 6σ-cleaned sample: {len(df):,} observations "
        f"({int(outlier_mask.sum()):,} outliers removed)",
        "INFO",
    )

    step_044_path = PROJECT_ROOT / "results/outputs/step_044_systematic_projection_analysis.json"
    if not step_044_path.exists():
        print_status("Step 044 results not found. Run Step 044 first.", "ERROR")
        sys.exit(1)

    with open(step_044_path, "r") as f:
        step_044 = json.load(f)

    step_040_path = PROJECT_ROOT / "results/outputs/step_040_unified_results_table.json"
    observed_eta = None
    if step_040_path.exists():
        with open(step_040_path, "r") as f:
            step_040 = json.load(f)
        observed_eta = (
            step_040.get("primary_estimands", {})
            .get("full_systematic_ols", {})
            .get("eta")
        )

    if observed_eta is None:
        step_003_path = PROJECT_ROOT / "results/outputs/step_003_statistical_analysis.json"
        if not step_003_path.exists():
            print_status("Neither Step 040 nor Step 003 results are available.", "ERROR")
            sys.exit(1)
        with open(step_003_path, "r") as f:
            step_003 = json.load(f)
        observed_eta = step_003.get("eta_full_systematic", step_003.get("eta_ols"))

    if observed_eta is None:
        print_status("Unable to resolve observed eta from upstream outputs.", "ERROR")
        sys.exit(1)

    observed_abs_eta = abs(observed_eta)
    print_status(f"Observed |eta| = {observed_abs_eta:.3e}", "INFO")

    systematics = step_044.get("systematic_projections", {})
    results = {}
    required_amp_m = observed_abs_eta * ETA_SCALE_FACTOR

    print_status("", "INFO")
    print_status(">>> Computing known-systematic sensitivity...", "PROCESS")

    for sys_name, sys_data in systematics.items():
        known_amp_m = sys_data.get("bias_A_m", 0)
        known_eta = sys_data.get("bias_eta", 0)
        ratio = required_amp_m / max(abs(known_amp_m), 1e-20)

        print_status(
            f"  {sys_name}: known_amp={known_amp_m*100:.3f} cm, "
            f"required={required_amp_m*100:.3f} cm, ratio={ratio:.1f}x",
            "CALC",
        )

        stable_seed = 61 + int(hashlib.sha256(sys_name.encode("utf-8")).hexdigest()[:8], 16) % 1000
        etas_mc = simulate_systematic_only(df, sys_name, known_amp_m, n_mc=2000, seed=stable_seed)
        p_exceed = float(np.mean(np.abs(etas_mc) >= observed_abs_eta))
        mean_eta_mc = float(np.mean(np.abs(etas_mc)))
        std_eta_mc = float(np.std(np.abs(etas_mc)))
        max_eta_mc = float(np.max(np.abs(etas_mc)))
        percentile_observed = float(np.mean(np.abs(etas_mc) <= observed_abs_eta) * 100)

        results[sys_name] = {
            "known_amplitude_m": float(known_amp_m),
            "known_eta": float(known_eta),
            "required_amplitude_m": float(required_amp_m),
            "ratio_required_to_known": float(ratio),
            "monte_carlo": {
                "n_mc": 2000,
                "mean_abs_eta": mean_eta_mc,
                "std_abs_eta": std_eta_mc,
                "max_abs_eta": max_eta_mc,
                "p_exceed_observed": p_exceed,
                "percentile_of_observed": percentile_observed,
            },
        }

    print_status("", "INFO")
    print_status(">>> Adversarial PCA nuisance search...", "PROCESS")
    pca_adv = adversarial_pca_absorption(df)
    print_status(
        f"  cosD survives all PCs: {pca_adv['cosd_survives_rich_pc_basis']}, "
        f"joint eta={pca_adv['joint_cosd_plus_all_pcs']['eta']:.3e}",
        "CALC",
    )

    print_status(">>> Adversarial 2D GP (elongation x time bins)...", "PROCESS")
    gp_adv = adversarial_gp_absorption(df, seed=TEP_CONFIG.get("RANDOM_SEED", 42))
    print_status(
        f"  cosD survives GP subtraction: {gp_adv['cosd_survives_gp_nuisance']}, "
        f"eta after GP={gp_adv['cosd_after_gp_subtraction']['eta']:.3e}",
        "CALC",
    )

    print_status(">>> Era x station x lunation interaction grid...", "PROCESS")
    grid = era_station_lunation_grid(df)
    print_status(
        f"  Fitted {grid['n_cells_fitted']}/{grid['n_cells_total']} cells; "
        f"{grid['n_cells_negative_eta']} negative, {grid['n_cells_significant_05']} significant",
        "CALC",
    )

    print_status(">>> Blind 20% year hold-out...", "PROCESS")
    blind = blind_year_holdout(df, seed=61)
    print_status(
        f"  Combined blind test eta={blind['inverse_variance_combined']['eta']:.3e} "
        f"({blind['inverse_variance_combined']['snr']:.2f}σ), "
        f"all splits negative={blind['all_splits_negative_eta']}",
        "CALC",
    )

    all_exceed = all(res["monte_carlo"]["p_exceed_observed"] < 0.05 for res in results.values())
    min_ratio = min(res["ratio_required_to_known"] for res in results.values())

    interpretation = (
        f"Known systematics: to explain |eta|={observed_abs_eta:.2e} requires "
        f"{required_amp_m*100:.2f} cm; smallest ratio to a known systematic is {min_ratio:.1f}x. "
        f"Adversarial PCA (20 residual-ranked PCs): joint eta="
        f"{pca_adv['joint_cosd_plus_all_pcs']['eta']:.2e} "
        f"(survival={pca_adv['cosd_survives_rich_pc_basis']}). "
        f"GP elongation nuisance absorbs "
        f"{gp_adv['absorption_fraction_gp_subtract']*100:.0f}% of cos(D) amplitude "
        f"(eta after subtraction={gp_adv['cosd_after_gp_subtraction']['eta']:.2e}). "
        f"Interaction grid: {grid['n_cells_negative_eta']}/{grid['n_cells_fitted']} fitted cells "
        f"have negative eta; blind 20% year hold-out combined eta="
        f"{blind['inverse_variance_combined']['eta']:.2e} "
        f"({blind['inverse_variance_combined']['snr']:.1f}σ). "
    )
    if all_exceed:
        interpretation += (
            "Modeled known-systematic-only MC falsification remains below observed |eta|. "
            "Unmodeled or data-driven nuisances are additionally stress-tested above."
        )
    else:
        interpretation += "At least one known systematic MC run approaches observed |eta|; see per-source table."

    output = {
        "step_id": "step_061",
        "status": "PASS",
        "method": (
            "Systematic amplitude sensitivity with Monte Carlo falsification; "
            "adversarial PCA/GP nuisance search; era×station×lunation grid; "
            "blind year hold-out"
        ),
        "n_raw": int(len(df_raw)),
        "n_outliers_removed": int(outlier_mask.sum()),
        "n_observations": int(len(df)),
        "observed_abs_eta": float(observed_abs_eta),
        "required_amplitude_m": float(required_amp_m),
        "required_amplitude_cm": float(required_amp_m * 100),
        "systematics": results,
        "adversarial_pca": pca_adv,
        "adversarial_gp": gp_adv,
        "interaction_grid": grid,
        "blind_year_holdout": blind,
        "tests_for_step_042": (
            grid["cell_tests_for_step_042"] + blind["cell_tests_for_step_042"]
        ),
        "interpretation": interpretation,
    }

    logger.save_step_results(output, PROJECT_ROOT, "step_061_systematic_sensitivity_analysis")
    print_status("Systematic Sensitivity Analysis Complete.", "SUCCESS")


if __name__ == "__main__":
    main()
