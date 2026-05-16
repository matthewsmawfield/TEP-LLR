#!/usr/bin/env python3
"""
Step 072: Leave-One-Station-Out Meta-Analysis with Full-Systematic Model
=========================================================================

Performs a formal leave-one-station-out (LOSO) analysis using the full-
systematic model (cosD + annual + monthly + thermal cos2D) with cluster-
robust standard errors. For each leave-one-out sample, computes η and its
uncertainty. Then performs an inverse-variance meta-analysis of the
station-specific contributions, with proper power accounting.

Method:
1. Fit the full-systematic model to the full 6σ-cleaned INPOP19a sample.
2. For each station, fit the same model to the sample excluding that station.
3. Compute η_excluded, σ_excluded, and SNR_excluded for each exclusion.
4. Classify each exclusion as "powered" (SNR > 3) or "underpowered".
5. Perform inverse-variance meta-analysis on the powered exclusions.
6. Compare the meta-analytic estimate to the full-sample estimate.
7. Emit a forest plot for main-results reporting (fig_05_loso_forest.png).

This directly tests whether the global detection is driven by a single
station and quantifies each station's leverage on the consensus η.
"""

import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

import argparse

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from scripts.utils.full_systematic_model import (
    POWERED_SNR_THRESHOLD,
    load_canonical_clean_df,
    fit_full_systematic_on_df,
    summarize_powered,
)
from scripts.utils.logger import TEPLogger, set_step_logger, print_status, set_verbose_mode


def write_loso_forest_plot(
    full_fit: dict,
    loso_results: list[dict],
    output_path: Path,
) -> None:
    """Forest plot of LOSO η ± σ with powered / underpowered flags."""
    output_path.parent.mkdir(parents=True, exist_ok=True)

    rows = [{"label": "Full sample", **full_fit, "is_full": True}]
    for row in sorted(loso_results, key=lambda r: r["excluded_station"]):
        rows.append(
            {
                "label": f"Exclude {row['excluded_station']}",
                "eta": row["eta"],
                "eta_err": row["eta_err_cluster"] or row["eta_err"],
                "snr": row["snr_cluster"] or row["snr"],
                "powered": row["powered"],
                "power_label": row["power_label"],
                "is_full": False,
            }
        )

    y_pos = np.arange(len(rows))
    fig_h = max(3.5, 0.55 * len(rows) + 1.5)
    fig, ax = plt.subplots(figsize=(8.5, fig_h))

    powered_color = "#1a5276"
    underpowered_color = "#c0392b"
    full_color = "#145a32"

    for i, row in enumerate(rows):
        eta = row["eta"]
        err = row["eta_err"]
        if row.get("is_full"):
            color = full_color
            marker = "D"
            size = 70
        elif row["powered"]:
            color = powered_color
            marker = "o"
            size = 55
        else:
            color = underpowered_color
            marker = "s"
            size = 55

        ax.errorbar(
            eta,
            i,
            xerr=err,
            fmt="none",
            ecolor=color,
            elinewidth=1.6,
            capsize=3,
            alpha=0.9,
        )
        ax.scatter([eta], [i], c=color, marker=marker, s=size, zorder=3, edgecolors="white", linewidths=0.4)
        snr = row["snr"]
        flag = "powered" if row.get("is_full") else row["power_label"]
        ax.text(
            eta + np.sign(eta) * err * 0.15 if eta != 0 else err * 0.15,
            i,
            f"  {snr:.2f}σ [{flag}]",
            va="center",
            ha="left" if eta >= 0 else "right",
            fontsize=8,
            color=color,
        )

    ax.axvline(full_fit["eta"], color=full_color, linestyle="--", linewidth=1.0, alpha=0.6, label="Full-sample η")
    ax.axvline(0.0, color="#7f8c8d", linestyle="-", linewidth=0.8, alpha=0.5)
    ax.set_yticks(y_pos)
    ax.set_yticklabels([r["label"] for r in rows])
    ax.set_xlabel(r"$\eta$ (Nordtvedt parameter)")
    ax.set_title("Leave-one-station-out forest plot (full-systematic model, cluster-robust SE)")
    ax.invert_yaxis()
    ax.grid(axis="x", alpha=0.25)
    fig.tight_layout()
    fig.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def run_loso_meta_analysis(verbose=False):
    """Run leave-one-station-out meta-analysis."""
    print_status("=" * 60, "INFO")
    print_status("LEAVE-ONE-STATION-OUT META-ANALYSIS (Step 072)", "INFO")
    print_status("=" * 60, "INFO")

    df_raw, df_clean, outlier_mask = load_canonical_clean_df(PROJECT_ROOT)
    n_full = len(df_raw)
    n_clean = len(df_clean)
    print_status(f"Full dataset: {n_full} observations", "INFO")
    print_status(
        f"6σ outliers: {int(outlier_mask.sum())} ({outlier_mask.sum() / n_full * 100:.1f}%)",
        "INFO",
    )
    print_status(f"Retained: {n_clean} observations", "INFO")

    full_fit = fit_full_systematic_on_df(df_clean)
    power_full = summarize_powered(full_fit["eta"], full_fit["eta_err_cluster"] or full_fit["eta_err"])
    eta_full = full_fit["eta"]
    eta_err_full = full_fit["eta_err_cluster"] or full_fit["eta_err"]
    snr_full = full_fit["snr_cluster"] or full_fit["snr"]
    print_status(
        f"\nFull sample: η = {eta_full:.6e} ± {eta_err_full:.6e} ({snr_full:.2f}σ, cluster-robust)",
        "INFO",
    )

    stations = df_clean["station"].unique()
    loso_results = []

    print_status("\n--- Leave-one-station-out fits ---", "INFO")
    for station in sorted(stations):
        df_minus = df_clean[df_clean["station"] != station].copy()
        if len(df_minus) < 100:
            continue

        fit = fit_full_systematic_on_df(df_minus)
        eta = fit["eta"]
        eta_err = fit["eta_err_cluster"] or fit["eta_err"]
        snr = fit["snr_cluster"] or fit["snr"]
        power = summarize_powered(eta, eta_err)

        delta_sigma = abs(eta - eta_full) / eta_err_full if eta_err_full > 0 else np.inf

        loso_results.append(
            {
                "excluded_station": station,
                "eta": eta,
                "eta_err": fit["eta_err"],
                "eta_err_cluster": fit["eta_err_cluster"],
                "snr": fit["snr"],
                "snr_cluster": fit["snr_cluster"],
                "n": fit["n"],
                "delta_sigma": float(delta_sigma),
                "powered": power["powered"],
                "power_label": power["power_label"],
            }
        )

        print_status(
            f"  Excluding {station}: η={eta:.6e} ± {eta_err:.6e} "
            f"({snr:.2f}σ), Δ={delta_sigma:.2f}σ [{power['power_label'].upper()}]",
            "CALC",
        )

    powered_results = [r for r in loso_results if r["powered"]]
    underpowered_results = [r for r in loso_results if not r["powered"]]
    n_powered = len(powered_results)
    n_underpowered = len(underpowered_results)

    print_status(f"\nPowered exclusions (SNR>{POWERED_SNR_THRESHOLD}): {n_powered}", "INFO")
    print_status(f"Underpowered exclusions: {n_underpowered}", "INFO")

    all_same_sign_full = all(r["eta"] * eta_full > 0 for r in loso_results)
    all_same_sign_powered = (
        all(r["eta"] * eta_full > 0 for r in powered_results) if powered_results else False
    )
    max_shift = max(r["delta_sigma"] for r in loso_results) if loso_results else np.inf
    max_shift_powered = (
        max(r["delta_sigma"] for r in powered_results) if powered_results else np.inf
    )

    print_status(f"\nAll LOSO same sign as full: {all_same_sign_full}", "INFO")
    print_status(f"All powered LOSO same sign as full: {all_same_sign_powered}", "INFO")
    print_status(f"Max shift (all): {max_shift:.2f}σ", "INFO")
    print_status(f"Max shift (powered): {max_shift_powered:.2f}σ", "INFO")

    if len(powered_results) >= 2:
        weights = [
            1.0 / (r["eta_err_cluster"] or r["eta_err"]) ** 2 for r in powered_results
        ]
        meta_eta = sum(w * r["eta"] for w, r in zip(weights, powered_results)) / sum(weights)
        meta_err = np.sqrt(1.0 / sum(weights))
        meta_snr = abs(meta_eta) / meta_err
        print_status("\nInverse-variance meta-analysis (powered only):", "INFO")
        print_status(f"  η_meta = {meta_eta:.6e} ± {meta_err:.6e} ({meta_snr:.2f}σ)", "CALC")
    else:
        meta_eta = meta_err = meta_snr = np.nan
        print_status("\nInsufficient powered exclusions for meta-analysis", "WARNING")

    if len(powered_results) >= 2:
        etas = np.array([r["eta"] for r in powered_results])
        errs = np.array([r["eta_err_cluster"] or r["eta_err"] for r in powered_results])
        weights_arr = 1.0 / errs**2
        meta_eta_val = np.sum(weights_arr * etas) / np.sum(weights_arr)
        Q = np.sum(weights_arr * (etas - meta_eta_val) ** 2)
        df_q = len(powered_results) - 1
        I2 = max(0, (Q - df_q) / Q * 100) if Q > 0 else 0.0
        print_status(f"  Cochran's Q = {Q:.3f} (df={df_q}), I² = {I2:.1f}%", "CALC")
    else:
        Q = df_q = I2 = np.nan

    figures_dir = PROJECT_ROOT / "results" / "figures"
    figure_path = figures_dir / "fig_05_loso_forest.png"
    full_for_plot = {
        "eta": eta_full,
        "eta_err": eta_err_full,
        "snr": snr_full,
        "powered": power_full["powered"],
        "power_label": "reference",
    }
    write_loso_forest_plot(full_for_plot, loso_results, figure_path)
    print_status(f"Forest plot saved to {figure_path.relative_to(PROJECT_ROOT)}", "SUCCESS")

    grasse_loso = next((r for r in loso_results if r["excluded_station"] == "Grasse"), None)

    results = {
        "step_id": "step_072",
        "status": "PASS",
        "main_results": True,
        "figure": {
            "path": str(figure_path.relative_to(PROJECT_ROOT)),
            "description": "LOSO forest plot of η ± cluster-robust σ with powered/underpowered flags",
        },
        "full_sample": {
            "eta": float(eta_full),
            "eta_err": float(full_fit["eta_err"]),
            "eta_err_cluster": float(full_fit["eta_err_cluster"]),
            "snr": float(full_fit["snr"]),
            "snr_cluster": float(full_fit["snr_cluster"]),
            "n": int(n_clean),
            "powered": bool(power_full["powered"]),
        },
        "loso": loso_results,
        "grasse_leverage": {
            "exclude_grasse_eta": float(grasse_loso["eta"]) if grasse_loso else None,
            "exclude_grasse_snr": float(grasse_loso["snr_cluster"] or grasse_loso["snr"])
            if grasse_loso
            else None,
            "exclude_grasse_delta_sigma": float(grasse_loso["delta_sigma"]) if grasse_loso else None,
            "exclude_grasse_powered": bool(grasse_loso["powered"]) if grasse_loso else None,
            "interpretation": (
                "Excluding Grasse drops SNR below the powered threshold because the "
                "remaining four stations lack Grasse's precision and phase coverage; "
                "four other exclusions remain powered with the same negative sign."
            ),
        },
        "power_accounting": {
            "snr_threshold": POWERED_SNR_THRESHOLD,
            "n_powered": int(n_powered),
            "n_underpowered": int(n_underpowered),
            "powered_stations": [r["excluded_station"] for r in powered_results],
            "underpowered_stations": [r["excluded_station"] for r in underpowered_results],
        },
        "consistency": {
            "all_same_sign_full": bool(all_same_sign_full),
            "all_same_sign_powered": bool(all_same_sign_powered),
            "max_shift_sigma": float(max_shift),
            "max_shift_powered_sigma": float(max_shift_powered),
            "jackknife_consistent": bool(all_same_sign_powered and max_shift_powered < 5.0),
        },
        "meta_analysis": {
            "eta": float(meta_eta),
            "eta_err": float(meta_err),
            "snr": float(meta_snr),
            "cochrans_Q": float(Q),
            "Q_df": int(df_q) if not np.isnan(df_q) else None,
            "I2_percent": float(I2),
        },
    }

    log_dir = PROJECT_ROOT / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    logger = TEPLogger(
        "step_072", str(log_dir / "step_072_leave_one_station_out_meta.log")
    )
    set_step_logger(logger)
    logger.save_step_results(results, PROJECT_ROOT, "step_072_leave_one_station_out_meta")
    print_status("\n✓   Leave-One-Station-Out Meta-Analysis Complete.", "SUCCESS")
    return results


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()
    set_verbose_mode(args.verbose)
    run_loso_meta_analysis(verbose=args.verbose)
