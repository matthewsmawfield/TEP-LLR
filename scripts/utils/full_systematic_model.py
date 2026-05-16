#!/usr/bin/env python3
"""Shared full-systematic regression utilities for station-leverage analyses."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from scripts.utils.llr_constants import ETA_SCALE_FACTOR
from scripts.utils.statistical_utils import (
    cluster_robust_variance,
    detect_outliers_sigma,
    robust_regression,
)

FULL_SYSTEMATIC_NAMES = [
    "cosD",
    "cos2D",
    "sin_m",
    "cos_m",
    "sin_y",
    "cos_y",
    "const",
]
STATIONS = ["Grasse", "APO", "McDonald2", "Matera", "Haleakala"]
POWERED_SNR_THRESHOLD = 3.0


def load_canonical_clean_df(project_root: Path) -> tuple[pd.DataFrame, pd.DataFrame, np.ndarray]:
    """Load INPOP19a residuals with canonical 6σ cleaning and time ordering."""
    input_path = project_root / "data" / "processed" / "INPOP19a_all_stations_residuals.csv"
    if not input_path.exists():
        raise FileNotFoundError(f"Input file not found: {input_path}")
    df_raw = pd.read_csv(input_path)
    if "date_julian" not in df_raw.columns or "station" not in df_raw.columns:
        raise ValueError("INPOP residuals CSV must include date_julian and station.")
    outlier_mask = detect_outliers_sigma(df_raw["residual_m"].values, sigma_threshold=6.0)
    df_clean = df_raw.loc[~outlier_mask].copy()
    df_clean = df_clean.sort_values(["date_julian", "station"], kind="mergesort").reset_index(
        drop=True
    )
    return df_raw, df_clean, outlier_mask


def build_full_systematic_matrix(df: pd.DataFrame) -> np.ndarray:
    """Build full-systematic design matrix (cosD as column 0)."""
    year = df["date_julian"].values / 365.25
    month = df["date_julian"].values / 27.32
    elongation = df["elongation_rad"].values
    return np.column_stack(
        [
            np.cos(elongation),
            np.cos(2.0 * elongation),
            np.sin(2.0 * np.pi * month),
            np.cos(2.0 * np.pi * month),
            np.sin(2.0 * np.pi * year),
            np.cos(2.0 * np.pi * year),
            np.ones(len(df)),
        ]
    )


def fit_eta_from_design(
    y: np.ndarray,
    X: np.ndarray,
    cluster_ids: np.ndarray | None = None,
    scale_errors_by_birge: bool = False,
) -> dict[str, Any]:
    """Fit design matrix; return η with OLS and optional cluster-robust errors."""
    reg = robust_regression(y, X, scale_errors_by_birge=scale_errors_by_birge)
    if reg.get("status") == "SINGULAR" or not np.isfinite(reg["coefficients"][0]):
        raise RuntimeError("Singular full-systematic regression.")

    coeffs = reg["coefficients"]
    errs_ols = reg["errors"]
    with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
        resid = y - X @ coeffs

    eta = float(coeffs[0] / ETA_SCALE_FACTOR)
    eta_err = float(errs_ols[0] / ETA_SCALE_FACTOR)
    snr = float(abs(eta) / max(eta_err, 1e-20))

    eta_err_cluster = None
    snr_cluster = None
    n_clusters = None
    if cluster_ids is not None and len(cluster_ids) == len(y):
        cr = cluster_robust_variance(
            X, resid, cluster_ids, small_sample_correction=True
        )
        eta_err_cluster = float(cr["se_cluster"][0] / ETA_SCALE_FACTOR)
        snr_cluster = float(abs(eta) / max(eta_err_cluster, 1e-20))
        n_clusters = int(cr["n_clusters"])

    return {
        "eta": eta,
        "eta_err": eta_err,
        "eta_err_cluster": eta_err_cluster,
        "snr": snr,
        "snr_cluster": snr_cluster,
        "n": int(len(y)),
        "n_clusters": n_clusters,
        "amplitude_m": float(coeffs[0]),
        "amplitude_err_m": float(errs_ols[0]),
    }


def fit_full_systematic_on_df(
    df: pd.DataFrame,
    cluster_by_station: bool = True,
) -> dict[str, Any]:
    """Fit full-systematic model on a dataframe subset."""
    if len(df) < 100:
        raise ValueError(f"Insufficient observations for fit: N={len(df)}")
    y = df["residual_m"].values
    X = build_full_systematic_matrix(df)
    cluster_ids = df["station"].values if cluster_by_station else None
    return fit_eta_from_design(y, X, cluster_ids=cluster_ids)


def fit_common_eta_station_systematics(df: pd.DataFrame) -> dict[str, Any]:
    """
    Mixed model: common η on cos(D) + station-specific nuisance systematics.

    H0: residual = η_common·cosD + Σ_s [α_s·cos2D + …] · I(station=s)
    """
    res = df["residual_m"].values
    st = df["station"].values
    el = df["elongation_rad"].values
    jd = df["date_julian"].values
    cos_c = np.cos(el)
    cos2d = np.cos(2.0 * el)
    year = jd / 365.25
    sin_y = np.sin(2.0 * np.pi * year)
    cos_y = np.cos(2.0 * np.pi * year)
    month = jd / 27.32
    sin_m = np.sin(2.0 * np.pi * month)
    cos_m = np.cos(2.0 * np.pi * month)

    n = len(res)
    stations = [s for s in STATIONS if s in set(st)]
    ns = len(stations)

    X0 = np.zeros((n, 1 + 6 * ns))
    X0[:, 0] = cos_c
    for i, station in enumerate(stations):
        mask = st == station
        X0[mask, 1 + i] = cos2d[mask]
        X0[mask, 1 + ns + i] = sin_m[mask]
        X0[mask, 1 + 2 * ns + i] = cos_m[mask]
        X0[mask, 1 + 3 * ns + i] = sin_y[mask]
        X0[mask, 1 + 4 * ns + i] = cos_y[mask]
        X0[mask, 1 + 5 * ns + i] = 1.0

    fit = fit_eta_from_design(res, X0, cluster_ids=st)
    fit["method"] = "common_eta_station_specific_systematics"
    fit["stations_modeled"] = stations
    return fit


def fit_non_grasse_with_grasse_nuisance(df: pd.DataFrame) -> dict[str, Any]:
    """
    Grasse-conditioned estimand: common η with explicit Grasse nuisance block.

    All stations contribute to cos(D); Grasse receives an additional cos(D)×I(Grasse)
    interaction so differential Grasse leverage on the synodic term is absorbed before
    reading the pooled η coefficient.
    """
    res = df["residual_m"].values
    st = df["station"].values
    el = df["elongation_rad"].values
    jd = df["date_julian"].values
    cos_c = np.cos(el)
    is_grasse = (st == "Grasse").astype(float)
    cos2d = np.cos(2.0 * el)
    year = jd / 365.25
    month = jd / 27.32
    sin_m = np.sin(2.0 * np.pi * month)
    cos_m = np.cos(2.0 * np.pi * month)
    sin_y = np.sin(2.0 * np.pi * year)
    cos_y = np.cos(2.0 * np.pi * year)

    # cosD, cosD×Grasse, cos2D, monthly, annual, const
    X = np.column_stack(
        [
            cos_c,
            cos_c * is_grasse,
            cos2d,
            sin_m,
            cos_m,
            sin_y,
            cos_y,
            np.ones(len(res)),
        ]
    )
    fit = fit_eta_from_design(res, X, cluster_ids=st)
    fit["method"] = "pooled_cosD_with_grasse_interaction"
    return fit


def summarize_powered(eta: float, eta_err: float, threshold: float = POWERED_SNR_THRESHOLD) -> dict[str, Any]:
    snr = abs(eta) / max(eta_err, 1e-20)
    return {
        "snr": float(snr),
        "powered": bool(snr > threshold),
        "power_label": "powered" if snr > threshold else "underpowered",
    }
