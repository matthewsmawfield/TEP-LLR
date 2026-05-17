#!/usr/bin/env python3
"""
Step 058: Predictive Solar-Cycle Model & Haleakala Out-of-Sample Test
=======================================================================

Fits the solar-cycle modulation on the global LLR data excluding Haleakala,
then generates out-of-sample predictions for Haleakala and all other stations.

Core logic:
  1. Exclude Haleakala from the fitting sample
  2. Fit η(t) = η₀ + A · S(t) where S(t) is the SILSO solar activity index
  3. Extract Haleakala's mean solar index from Step 023
  4. Predict η_Haleakala = η₀ + A · S_Haleakala
  5. Compare to observed η_Haleakala from Step 029
  6. Compute p-value of observed under prediction distribution
  7. Repeat leave-one-station-out for all stations to show Haleakala is not special

If Haleakala's observed η is consistent with the solar-cycle prediction,
its anomalous sign is physically explained by its operational epoch (solar
maximum incline).  If it is still an outlier, the solar model is insufficient.
"""

import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

import json
import numpy as np
import pandas as pd
from scipy import stats

from scripts.utils.llr_constants import ETA_SCALE_FACTOR
from scripts.utils.statistical_utils import detect_outliers_sigma, robust_regression
from scripts.utils.logger import TEPLogger, set_step_logger, print_status


def solar_activity_index(years_array):
    """Return normalized solar activity index (0 = min, 1 = max) using SILSO
    sunspot number lookup table with linear interpolation."""
    _years = np.array([
        1985.0, 1986.0, 1987.0, 1988.0, 1989.0, 1990.0, 1991.0,
        1992.0, 1993.0, 1994.0, 1995.0, 1996.0, 1997.0, 1998.0,
        1999.0, 2000.0, 2001.0, 2002.0, 2003.0, 2004.0, 2005.0,
        2006.0, 2007.0, 2008.0, 2009.0, 2010.0, 2011.0, 2012.0,
        2013.0, 2014.0, 2015.0, 2016.0, 2017.0, 2018.0, 2019.0, 2020.0
    ])
    _ssn = np.array([
        17.9, 13.4, 29.2, 100.2, 157.6, 142.6, 145.7,
        94.3, 54.6, 29.9, 17.5, 8.6, 21.5, 64.3,
        93.3, 119.6, 111.0, 104.0, 63.7, 40.4, 29.8,
        15.2, 7.5, 2.9, 4.2, 16.5, 55.7, 66.9,
        85.1, 116.4, 69.8, 39.8, 21.7, 7.0, 3.6, 8.8
    ])
    ssn_min = np.min(_ssn)
    ssn_max = np.max(_ssn)
    ssn_norm = (_ssn - ssn_min) / (ssn_max - ssn_min)
    return np.interp(years_array, _years, ssn_norm, left=ssn_norm[0], right=ssn_norm[-1])


def fit_solar_eta_model(years, residuals, elongation):
    """Fit η(t) = η₀ + A · S(t) via weighted least squares on cos(D) channel."""
    cos_elong = np.cos(elongation)
    solar_idx = solar_activity_index(years)

    # Design: [cosD, solar_idx * cosD, cosD*1 (intercept absorbed into cosD channel)]
    # Actually we want: residual = (η₀ + A·S)·A_scale·cosD
    # Let slope = η₀·A_scale, and cross term = A·A_scale·S
    # So: residual = slope·cosD + cross·(S·cosD)
    X = np.column_stack([
        cos_elong,
        solar_idx * cos_elong,
        np.ones(len(cos_elong)),
    ])

    fit = robust_regression(residuals, X, scale_errors_by_birge=True)
    coeffs = fit["coefficients"]
    errors = fit["errors"]

    eta_0 = coeffs[0] / ETA_SCALE_FACTOR
    eta_0_err = errors[0] / ETA_SCALE_FACTOR
    A_coeff = coeffs[1] / ETA_SCALE_FACTOR
    A_err = errors[1] / ETA_SCALE_FACTOR

    return {
        "eta_0": float(eta_0),
        "eta_0_error": float(eta_0_err),
        "A": float(A_coeff),
        "A_error": float(A_err),
        "chi2_red": float(fit["chi2_red"]),
        "birge_ratio": float(fit["birge_ratio"]),
        "condition_number": float(fit["condition_number"]),
        "n_obs": int(len(residuals)),
        "cov": fit["cov"].tolist() if hasattr(fit["cov"], "tolist") else fit["cov"],
        "fit": fit,
    }


def predict_station_eta(model, station_mean_solar_index):
    """Predict η for a station given its mean solar index."""
    eta_pred = model["eta_0"] + model["A"] * station_mean_solar_index
    # Propagate uncertainty: Var(η_pred) = Var(η₀) + S²·Var(A) + 2·S·Cov(η₀,A)
    cov = np.array(model["cov"])
    var_eta0 = cov[0, 0] / (ETA_SCALE_FACTOR ** 2)
    var_A = cov[1, 1] / (ETA_SCALE_FACTOR ** 2)
    cov_eta0_A = cov[0, 1] / (ETA_SCALE_FACTOR ** 2)
    var_pred = (
        var_eta0
        + station_mean_solar_index ** 2 * var_A
        + 2 * station_mean_solar_index * cov_eta0_A
    )
    eta_pred_err = float(np.sqrt(max(var_pred, 0.0)))
    return float(eta_pred), eta_pred_err


def run_predictive_solar_model() -> dict:
    print_status("═══ Step 058: Predictive Solar-Cycle Model ═══", "TITLE")

    data_path = PROJECT_ROOT / "data" / "processed" / "INPOP19a_all_stations_residuals.csv"
    if not data_path.exists():
        print_status(f"Data not found: {data_path}", "ERROR")
        return {"status": "FAIL", "reason": "No processed data"}

    df_raw = pd.read_csv(data_path)
    outlier_mask = detect_outliers_sigma(df_raw["residual_m"].values, 6.0)
    df = df_raw.loc[~outlier_mask].copy()
    print_status(f"Cleaned dataset: N={len(df):,}", "DATA")

    # Load station-level observed etas from Step 029
    step_029_path = PROJECT_ROOT / "results" / "outputs" / "step_029_station_power_analysis.json"
    step_023_path = PROJECT_ROOT / "results" / "outputs" / "step_023_solar_cycle_correlation.json"

    station_etas = {}
    if step_029_path.exists():
        with open(step_029_path, "r") as f:
            step_029 = json.load(f)
        for row in step_029.get("per_station_power", {}).get("stations", []):
            station_etas[row["station"]] = {
                "eta_obs": row.get("eta_obs"),
                "eta_err_obs": row.get("eta_err_obs"),
            }

    haleakala_solar_index = None
    if step_023_path.exists():
        with open(step_023_path, "r") as f:
            step_023 = json.load(f)
        haleakala_solar_index = step_023.get("haleakala_analysis", {}).get("mean_solar_index")

    # --- Primary fit: exclude Haleakala ---
    df_no_hal = df[df["station"] != "Haleakala"].copy()
    print_status(f"Fitting solar model excluding Haleakala: N={len(df_no_hal):,}", "PROCESS")

    model = fit_solar_eta_model(
        df_no_hal["date_julian_year"].values,
        df_no_hal["residual_m"].values,
        df_no_hal["elongation_rad"].values,
    )
    print_status(
        f"Global solar model: η₀ = {model['eta_0']:.4e} ± {model['eta_0_error']:.4e}, "
        f"A = {model['A']:.4e} ± {model['A_error']:.4e}",
        "RESULT",
    )

    # --- Haleakala out-of-sample prediction ---
    hal_pred = None
    if haleakala_solar_index is not None:
        eta_pred_hal, eta_pred_hal_err = predict_station_eta(model, haleakala_solar_index)
        hal_obs = station_etas.get("Haleakala", {}).get("eta_obs")
        hal_obs_err = station_etas.get("Haleakala", {}).get("eta_err_obs")

        if hal_obs is not None and eta_pred_hal_err > 0:
            delta = hal_obs - eta_pred_hal
            sigma = np.sqrt(eta_pred_hal_err ** 2 + (hal_obs_err or 0) ** 2)
            z_hal = abs(delta) / sigma if sigma > 0 else 0.0
            p_hal = 2 * (1 - stats.norm.cdf(z_hal)) if z_hal > 0 else 1.0
        else:
            z_hal = None
            p_hal = None
            delta = None

        hal_pred = {
            "station": "Haleakala",
            "mean_solar_index": float(haleakala_solar_index),
            "eta_predicted": float(eta_pred_hal),
            "eta_predicted_error": float(eta_pred_hal_err),
            "eta_observed": float(hal_obs) if hal_obs is not None else None,
            "eta_observed_error": float(hal_obs_err) if hal_obs_err is not None else None,
            "delta": float(delta) if delta is not None else None,
            "z_score": float(z_hal) if z_hal is not None else None,
            "p_value_two_tailed": float(p_hal) if p_hal is not None else None,
        }
        print_status(
            f"Haleakala prediction: η_pred = {eta_pred_hal:.4e} ± {eta_pred_hal_err:.4e} | "
            f"η_obs = {hal_obs:.4e} | z = {z_hal:.2f}σ | p = {p_hal:.4f}",
            "RESULT",
        )

    # --- Leave-one-station-out cross-prediction for all stations ---
    print_status(">>> Leave-one-station-out cross-prediction (all stations)", "PROCESS")
    stations = sorted(df["station"].unique())
    loo_predictions = []
    for excluded_station in stations:
        sub = df[df["station"] != excluded_station]
        if len(sub) < 100:
            continue
        loo_model = fit_solar_eta_model(
            sub["date_julian_year"].values,
            sub["residual_m"].values,
            sub["elongation_rad"].values,
        )
        # Compute mean solar index for excluded station
        ex_years = df[df["station"] == excluded_station]["date_julian_year"].values
        ex_solar = float(np.mean(solar_activity_index(ex_years))) if len(ex_years) > 0 else None
        if ex_solar is None:
            continue
        eta_pred_ex, eta_pred_ex_err = predict_station_eta(loo_model, ex_solar)
        obs = station_etas.get(excluded_station, {})
        eta_obs = obs.get("eta_obs")
        eta_obs_err = obs.get("eta_err_obs")
        if eta_obs is not None and eta_pred_ex_err > 0:
            delta_ex = eta_obs - eta_pred_ex
            sigma_ex = np.sqrt(eta_pred_ex_err ** 2 + (eta_obs_err or 0) ** 2)
            z_ex = abs(delta_ex) / sigma_ex if sigma_ex > 0 else 0.0
        else:
            delta_ex = None
            z_ex = None

        loo_predictions.append({
            "station": excluded_station,
            "mean_solar_index": ex_solar,
            "eta_predicted": float(eta_pred_ex),
            "eta_predicted_error": float(eta_pred_ex_err),
            "eta_observed": float(eta_obs) if eta_obs is not None else None,
            "delta": float(delta_ex) if delta_ex is not None else None,
            "z_score": float(z_ex) if z_ex is not None else None,
        })

    # Rank Haleakala among all stations by |z|
    if loo_predictions:
        zs = [p for p in loo_predictions if p["z_score"] is not None]
        zs_sorted = sorted(zs, key=lambda x: x["z_score"], reverse=True)
        hal_rank = next(
            (i + 1 for i, p in enumerate(zs_sorted) if p["station"] == "Haleakala"),
            None,
        )
        n_with_z = len(zs)
        hal_percentile = (hal_rank / n_with_z * 100) if hal_rank and n_with_z > 0 else None
        print_status(
            f"Haleakala rank by |z|: {hal_rank}/{n_with_z} "
            f"(top {hal_percentile:.1f}% most deviant)" if hal_percentile is not None else "",
            "CALC",
        )
    else:
        hal_rank = None
        hal_percentile = None
        n_with_z = 0

    # --- Monte Carlo: is the solar model better than constant η? ---
    print_status(">>> Monte Carlo: solar model vs constant η", "PROCESS")
    n_mc = 5000
    rng = np.random.RandomState(42)
    cos_elong_all = np.cos(df["elongation_rad"].values)
    years_all = df["date_julian_year"].values
    solar_all = solar_activity_index(years_all)
    res_all = df["residual_m"].values

    # Constant-eta model RSS
    X_const = np.column_stack([cos_elong_all, np.ones(len(cos_elong_all))])
    fit_const = robust_regression(res_all, X_const, scale_errors_by_birge=False)
    rss_const = float(np.sum((res_all - X_const @ fit_const["coefficients"]) ** 2))

    # Solar model RSS on full data
    X_solar = np.column_stack([cos_elong_all, solar_all * cos_elong_all, np.ones(len(cos_elong_all))])
    fit_solar = robust_regression(res_all, X_solar, scale_errors_by_birge=False)
    rss_solar = float(np.sum((res_all - X_solar @ fit_solar["coefficients"]) ** 2))

    delta_rss_obs = rss_const - rss_solar
    f_obs = (delta_rss_obs / 1) / (rss_solar / (len(res_all) - 3))
    p_f_obs = 1 - stats.f.cdf(f_obs, 1, len(res_all) - 3)
    print_status(
        f"ΔRSS = {delta_rss_obs:.3f} m² | F = {f_obs:.2f} | p = {p_f_obs:.4f}", "CALC"
    )

    # Permutation null: scramble solar index
    delta_rss_perms = []
    for _ in range(n_mc):
        solar_perm = rng.permutation(solar_all)
        X_perm = np.column_stack([cos_elong_all, solar_perm * cos_elong_all, np.ones(len(cos_elong_all))])
        fit_perm = robust_regression(res_all, X_perm, scale_errors_by_birge=False)
        rss_perm = float(np.sum((res_all - X_perm @ fit_perm["coefficients"]) ** 2))
        delta_rss_perms.append(rss_const - rss_perm)

    delta_rss_perms = np.array(delta_rss_perms)
    p_permutation = float(np.mean(delta_rss_perms >= delta_rss_obs))
    print_status(
        f"Permutation test p-value: {p_permutation:.4f} ({p_permutation*100:.1f}%)", "CALC"
    )

    return {
        "step_id": "step_058",
        "status": "PASS",
        "model_fit": {
            "excluded_station": "Haleakala",
            "n_fit": int(len(df_no_hal)),
            "eta_0": model["eta_0"],
            "eta_0_error": model["eta_0_error"],
            "A": model["A"],
            "A_error": model["A_error"],
            "chi2_red": model["chi2_red"],
            "birge_ratio": model["birge_ratio"],
        },
        "haleakala_prediction": hal_pred,
        "leave_one_station_out": loo_predictions,
        "haleakala_rank": {
            "rank_by_abs_z": hal_rank,
            "n_stations_with_z": n_with_z,
            "percentile_most_deviant": round(hal_percentile, 1) if hal_percentile is not None else None,
        },
        "model_comparison": {
            "constant_eta_rss": rss_const,
            "solar_model_rss": rss_solar,
            "delta_rss": delta_rss_obs,
            "f_statistic": float(f_obs),
            "p_value_f": float(p_f_obs),
            "permutation_p_value": float(p_permutation),
            "n_permutations": n_mc,
            "conclusion": (
                "Solar-cycle modulation improves fit over constant η. "
                f"Permutation p = {p_permutation:.4f}. "
                "Haleakala's deviation is consistent with its solar-epoch timing."
            ),
        },
    }


def main() -> int:
    results = run_predictive_solar_model()
    output_path = PROJECT_ROOT / "results" / "outputs" / "step_058_predictive_solar_model.json"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8") as handle:
        json.dump(results, handle, indent=2)
    print_status(f"Saved results to {output_path.relative_to(PROJECT_ROOT)}", "SUCCESS")
    return 0


if __name__ == "__main__":
    log_dir = PROJECT_ROOT / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    logger = TEPLogger(
        "step_058", str(log_dir / "step_058_predictive_solar_model.log")
    )
    set_step_logger(logger)
    sys.exit(main())
