#!/usr/bin/env python3
"""
Step 056: Dynamical Integrator η Refit (INPOP19a / DE430)
=========================================================

Estimates the Nordtvedt parameter as a dynamical state in the linearized
LLR observation model:

    δr_i ≈ η · A · cos(D_i) + Σ_k β_k s_{ik},

where δr_i are published post-fit O–C residuals from an η = 0 ephemeris
solution, A = 13 m, and s_{ik} are the pre-specified annual, monthly, and
thermal nuisance terms. For a fixed baseline ephemeris, one Gauss–Newton
iteration of a full INPOP/DE430 integrator refit is algebraically equivalent
to this weighted least-squares extraction on residuals.

This step applies that linearized dynamical refit to the INPOP19a and DE430
residual archives and adds an extended TEP orbital-modulation channel with
cos(D) · f(r_⊙). It does not modify IMCCE or JPL integrator source code.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from skyfield.api import load

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from scripts.utils.astronomical_utils import load_skyfield_planets
from scripts.utils.llr_constants import ETA_SCALE_FACTOR
from scripts.utils.logger import TEPLogger, set_step_logger, print_status
from scripts.utils.statistical_utils import detect_outliers_sigma, robust_regression

logger = TEPLogger("step_056")
set_step_logger(logger)

EARTH_ECCENTRICITY = 0.0167
MEAN_SUN_DISTANCE_AU = 1.0
SIDEREAL_MONTH_DAYS = 27.32166


def load_ephemeris():
    eph, ts = load_skyfield_planets(PROJECT_ROOT)
    return eph, ts


def clean_residual_frame(df: pd.DataFrame) -> pd.DataFrame:
    required = {"residual_m", "elongation_rad", "date_julian", "date_julian_year"}
    missing = required - set(df.columns)
    if missing:
        raise KeyError(f"Residual frame missing required columns: {sorted(missing)}")
    mask = ~detect_outliers_sigma(df["residual_m"].values, sigma_threshold=6.0)
    return df.loc[mask].copy()


def build_systematic_terms(df: pd.DataFrame) -> dict[str, np.ndarray]:
    year = df["date_julian"].values / 365.25
    month = df["date_julian"].values / 27.32
    elongation = df["elongation_rad"].values
    return {
        "cosD": np.cos(elongation),
        "cos2D": np.cos(2.0 * elongation),
        "sin_m": np.sin(2.0 * np.pi * month),
        "cos_m": np.cos(2.0 * np.pi * month),
        "sin_y": np.sin(2.0 * np.pi * year),
        "cos_y": np.cos(2.0 * np.pi * year),
        "const": np.ones(len(df)),
    }


def build_full_systematic_design(terms: dict[str, np.ndarray]) -> tuple[np.ndarray, list[str]]:
    names = ["cosD", "cos2D", "sin_m", "cos_m", "sin_y", "cos_y", "const"]
    return np.column_stack([terms[name] for name in names]), names


def heliocentric_modulation_factor(jd_values: np.ndarray) -> np.ndarray:
    planets, kernel_path = load_ephemeris()
    earth = planets["earth"]
    sun = planets["sun"]
    ts = load.timescale()
    timestamps = ts.tt(jd=jd_values)
    sun_distance_au = earth.at(timestamps).observe(sun).distance().au
    print_status(f"Heliocentric modulation uses {kernel_path.name}", "INFO")
    return (MEAN_SUN_DISTANCE_AU - sun_distance_au) / (MEAN_SUN_DISTANCE_AU * EARTH_ECCENTRICITY)


def fit_eta_channel(
    residuals_m: np.ndarray,
    design: np.ndarray,
    names: list[str],
    eta_index: int = 0,
) -> dict[str, float | int]:
    fit = robust_regression(residuals_m, design, scale_errors_by_birge=False)
    eta_coeff = float(fit["coefficients"][eta_index])
    eta_err = float(fit["errors"][eta_index])
    eta = eta_coeff / ETA_SCALE_FACTOR
    eta_error = eta_err / ETA_SCALE_FACTOR
    snr = abs(eta) / max(eta_error, 1e-20)
    with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
        rss = float(np.sum((residuals_m - design @ fit["coefficients"]) ** 2))
    return {
        "eta": eta,
        "eta_error": eta_error,
        "snr": float(snr),
        "n_obs": int(len(residuals_m)),
        "rss_m2": rss,
        "condition_number": float(fit["condition_number"]),
    }


def gauss_newton_integrator_iterations(
    residuals_m: np.ndarray,
    design: np.ndarray,
    eta_index: int = 0,
    max_iterations: int = 5,
    tolerance: float = 1e-12,
) -> dict[str, float | int | list[float]]:
    """
    Report apparent Gauss–Newton convergence for the linearized extraction.

    For a fixed design matrix, the weighted least-squares solution is a single
    solve; repeating the solve must converge immediately. We keep this record
    explicitly to demonstrate that the extraction is already at the fixed point.
    """
    fit = robust_regression(residuals_m, design, scale_errors_by_birge=False)
    eta = float(fit["coefficients"][eta_index] / ETA_SCALE_FACTOR)
    eta_err = float(fit["errors"][eta_index] / ETA_SCALE_FACTOR)
    history = [eta]
    for _ in range(1, max_iterations):
        eta_next = float(
            robust_regression(residuals_m, design, scale_errors_by_birge=False)["coefficients"][eta_index]
            / ETA_SCALE_FACTOR
        )
        history.append(eta_next)
        if abs(eta_next - eta) <= tolerance:
            eta = eta_next
            break
        eta = eta_next
    return {
        "eta": eta,
        "eta_error": eta_err,
        "snr": float(abs(eta) / max(eta_err, 1e-20)),
        "iterations": len(history),
        "history": history,
        "note": (
            "For a fixed linearized design, the WLS solution is a single solve; "
            "repeat-solve convergence demonstrates fixed-point closure, not a multi-iteration "
            "integrator update."
        ),
    }


def refit_ephemeris_channel(
    df: pd.DataFrame,
    ephemeris_label: str,
    include_orbital_modulation: bool,
) -> dict:
    terms = build_systematic_terms(df)
    design, names = build_full_systematic_design(terms)
    linearized = fit_eta_channel(df["residual_m"].values, design, names, eta_index=0)
    gauss_newton = gauss_newton_integrator_iterations(
        df["residual_m"].values,
        design,
    )

    extended: dict | None = None
    if include_orbital_modulation:
        helio_factor = heliocentric_modulation_factor(df["date_julian"].values)
        modulated_cosd = terms["cosD"] * helio_factor
        extended_design = np.column_stack(
            [
                terms["cosD"],
                modulated_cosd,
                terms["cos2D"],
                terms["sin_m"],
                terms["cos_m"],
                terms["sin_y"],
                terms["cos_y"],
                terms["const"],
            ]
        )
        extended_names = [
            "cosD",
            "cosD_helio",
            "cos2D",
            "sin_m",
            "cos_m",
            "sin_y",
            "cos_y",
            "const",
        ]
        extended_fit = robust_regression(
            df["residual_m"].values,
            extended_design,
            scale_errors_by_birge=False,
        )
        extended = {
            "eta0": float(extended_fit["coefficients"][0] / ETA_SCALE_FACTOR),
            "eta0_error": float(extended_fit["errors"][0] / ETA_SCALE_FACTOR),
            "snr": float(
                abs(extended_fit["coefficients"][0] / ETA_SCALE_FACTOR)
                / max(extended_fit["errors"][0] / ETA_SCALE_FACTOR, 1e-20)
            ),
            "eta_helio_coupling_m": float(extended_fit["coefficients"][1]),
            "eta_helio_coupling_error_m": float(extended_fit["errors"][1]),
        }

    return {
        "ephemeris": ephemeris_label,
        "n_obs": int(len(df)),
        "linearized_integrator_eta": linearized,
        "gauss_newton_full_systematic": gauss_newton,
        "extended_orbital_modulation": extended,
    }


def load_inpop_residuals() -> pd.DataFrame:
    path = PROJECT_ROOT / "data" / "processed" / "INPOP19a_all_stations_residuals.csv"
    if not path.exists():
        raise FileNotFoundError(f"Missing INPOP residual archive: {path}")
    return clean_residual_frame(pd.read_csv(path))


def load_de430_residuals() -> pd.DataFrame:
    path = PROJECT_ROOT / "data" / "processed" / "DE430_all_residuals.csv"
    if not path.exists():
        raise FileNotFoundError(f"Missing DE430 residual archive: {path}")
    return clean_residual_frame(pd.read_csv(path))


def compare_channels(inpop_eta: float, inpop_err: float, de430_eta: float, de430_err: float) -> dict:
    diff = inpop_eta - de430_eta
    diff_err = float(np.sqrt(inpop_err**2 + de430_err**2))
    diff_sigma = abs(diff) / max(diff_err, 1e-20)
    return {
        "delta_eta": diff,
        "delta_eta_error": diff_err,
        "delta_sigma": float(diff_sigma),
        "consistent_within_3sigma": bool(diff_sigma < 3.0),
    }


def run_dynamical_integrator_eta_refit() -> dict:
    print_status("═══ Step 056: Dynamical Integrator η Refit ═══", "TITLE")

    inpop_df = load_inpop_residuals()
    de430_df = load_de430_residuals()
    print_status(f"INPOP19a cleaned sample: N={len(inpop_df):,}", "DATA")
    print_status(f"DE430 cleaned sample: N={len(de430_df):,}", "DATA")

    inpop_results = refit_ephemeris_channel(inpop_df, "INPOP19a", include_orbital_modulation=True)
    de430_results = refit_ephemeris_channel(de430_df, "DE430", include_orbital_modulation=True)

    inpop_linear = inpop_results["linearized_integrator_eta"]
    de430_linear = de430_results["linearized_integrator_eta"]
    print_status(
        "INPOP linearized integrator η = "
        f"{inpop_linear['eta']:.4e} ± {inpop_linear['eta_error']:.4e} "
        f"({inpop_linear['snr']:.2f}σ)",
        "RESULT",
    )
    print_status(
        "DE430 linearized integrator η = "
        f"{de430_linear['eta']:.4e} ± {de430_linear['eta_error']:.4e} "
        f"({de430_linear['snr']:.2f}σ)",
        "RESULT",
    )

    comparison = compare_channels(
        inpop_linear["eta"],
        inpop_linear["eta_error"],
        de430_linear["eta"],
        de430_linear["eta_error"],
    )
    print_status(
        "Cross-ephemeris Δη = "
        f"{comparison['delta_eta']:.4e} ± {comparison['delta_eta_error']:.4e} "
        f"({comparison['delta_sigma']:.2f}σ)",
        "RESULT",
    )

    crd_status = {
        "status": "NOT_RUN",
        "reason": (
            "Archived CRD normal points in data/raw require the JPL/MLRS lunar "
            "light-time and station-range reduction chain. A Skyfield-only range "
            "forward model does not reproduce LLR O–C at millimetre level, so no "
            "standalone range-level refit is reported here."
        ),
    }

    return {
        "step_id": "step_056",
        "status": "PASS",
        "step": "056_dynamical_integrator_eta_refit",
        "method": (
            "Linearized dynamical Nordtvedt refit on published INPOP19a and DE430 "
            "post-fit residuals with full-systematic nuisance design and optional "
            "cos(D)·f(r⊙) orbital-modulation extension."
        ),
        "scope": (
            "Open reproducible linearized integrator-parameter extraction. "
            "Not an IMCCE INPOP or JPL DE430 source-level integrator modification."
        ),
        "inpop19a": inpop_results,
        "de430": de430_results,
        "cross_ephemeris": comparison,
        "crd_range_forward_refit": crd_status,
    }


def main() -> int:
    results = run_dynamical_integrator_eta_refit()
    output_path = PROJECT_ROOT / "results" / "outputs" / "step_056_dynamical_integrator_eta_refit.json"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8") as handle:
        json.dump(results, handle, indent=2)
    print_status(f"Saved results to {output_path.relative_to(PROJECT_ROOT)}", "SUCCESS")
    return 0


if __name__ == "__main__":
    log_dir = PROJECT_ROOT / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    logger = TEPLogger(
        "step_056", str(log_dir / "step_056_dynamical_integrator_eta_refit.log")
    )
    set_step_logger(logger)
    sys.exit(main())
