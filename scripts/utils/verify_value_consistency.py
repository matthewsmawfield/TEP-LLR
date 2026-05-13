#!/usr/bin/env python3
"""Verify manuscript source values against canonical pipeline outputs."""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]
OUTPUTS_DIR = PROJECT_ROOT / "results" / "outputs"
COMPONENTS_DIR = PROJECT_ROOT / "site" / "components"


def load_json(path: Path) -> dict:
    if not path.exists():
        raise FileNotFoundError(f"Missing required pipeline output: {path}")
    with path.open(encoding="utf-8") as handle:
        return json.load(handle)


def format_sci_latex(value: float, decimals: int = 2) -> str:
    if value == 0:
        return "0"
    sign = "-" if value < 0 else ""
    magnitude = abs(value)
    exponent = int(math.floor(math.log10(magnitude)))
    mantissa = magnitude / (10**exponent)
    mantissa_str = f"{mantissa:.{decimals}f}"
    if float(mantissa_str) >= 10:
        mantissa_str = f"{float(mantissa_str) / 10:.{decimals}f}"
        exponent += 1
    return f"{sign}{mantissa_str} \\times 10^{{{exponent}}}"


def format_sigma(value: float, decimals: int = 2) -> str:
    return f"{value:.{decimals}f}\\sigma"


def read_manuscript_sources() -> str:
    if not COMPONENTS_DIR.is_dir():
        raise FileNotFoundError(f"Missing manuscript components directory: {COMPONENTS_DIR}")
    return "\n".join(
        path.read_text(encoding="utf-8")
        for path in sorted(COMPONENTS_DIR.glob("*.html"))
    )


def require_substring(text: str, needle: str, label: str, errors: list[str]) -> None:
    if needle not in text:
        errors.append(f"{label}: expected manuscript source to contain {needle!r}")


def forbid_substring(text: str, needle: str, label: str, errors: list[str]) -> None:
    if needle in text:
        errors.append(f"{label}: found stale manuscript value {needle!r}")


def main() -> int:
    errors: list[str] = []

    step_001 = load_json(OUTPUTS_DIR / "step_001_data_preprocessing.json")
    step_003 = load_json(OUTPUTS_DIR / "step_003_statistical_analysis.json")
    step_040 = load_json(OUTPUTS_DIR / "step_040_unified_results_table.json")
    step_050 = load_json(OUTPUTS_DIR / "step_050_corrected_tep_analysis.json")
    step_006b = load_json(OUTPUTS_DIR / "step_006b_de430_outlier_robustness.json")
    step_042 = load_json(OUTPUTS_DIR / "step_042_multiple_testing_correction.json")
    step_046 = load_json(OUTPUTS_DIR / "step_046_station_balanced_tep.json")
    step_055 = load_json(OUTPUTS_DIR / "step_055_cmb_rigorous_falsification.json")
    step_056 = load_json(OUTPUTS_DIR / "step_056_dynamical_integrator_eta_refit.json")
    step_029 = load_json(OUTPUTS_DIR / "step_029_station_power_analysis.json")

    primary = step_040["primary_estimands"]["full_systematic_ols"]
    cluster = step_050["models"]["m5_full_corrected"]["cluster_robust"]
    cosd_only = step_040["primary_estimands"]["full_sample_ols"]
    de430 = step_006b["threshold_sweep"]["6sigma"]
    phase_chi2 = step_006b["phase_chi_square"]
    inpop_integrator = step_056["inpop19a"]["linearized_integrator_eta"]
    de430_integrator = step_056["de430"]["linearized_integrator_eta"]
    delta_integrator = step_056["cross_ephemeris"]

    primary_multiplicity = next(
        test
        for test in step_042["tests"]
        if test["name"] == "Full-Systematic OLS (primary estimand)"
    )
    step_046_full = next(test for test in step_046["tests"] if test["name"] == "full_sample")
    step_046_equal_n = next(
        test for test in step_046["tests"] if test["name"] == "equal_n_subsample"
    )
    step_046_grasse_capped = next(
        test for test in step_046["tests"] if test["name"] == "grasse_capped"
    )
    sky_scramble = step_055["sky_scrambling"]
    wls_eta = step_029["summary"]["precision_weighted_eta"]
    wls_snr = step_029["summary"]["precision_weighted_snr"]

    total_obs = step_001["combined"]["total_obs"]
    n_clean = primary["n_obs"]

    eta = primary["eta"]
    eta_error = primary["eta_error"]
    snr = primary["snr"]
    cluster_snr = cluster["snr_cluster"]
    cosd_eta = cosd_only["eta"]
    cosd_snr = cosd_only["snr"]
    keplerian = step_050["keplerian_inclusion_proxy"]
    kepler_eta = keplerian["eta_after_kepler_partialing"]
    joint_eta = keplerian["eta_joint_with_kepler_terms"]

    manuscript = read_manuscript_sources()

    require_substring(
        manuscript,
        format_sci_latex(eta),
        "primary eta",
        errors,
    )
    require_substring(
        manuscript,
        format_sci_latex(eta_error),
        "primary eta uncertainty",
        errors,
    )
    require_substring(
        manuscript,
        format_sigma(snr),
        "primary SNR",
        errors,
    )
    require_substring(
        manuscript,
        format_sigma(cluster_snr),
        "cluster-robust SNR",
        errors,
    )
    require_substring(
        manuscript,
        format_sci_latex(kepler_eta["eta"]),
        "Keplerian partialing eta",
        errors,
    )
    require_substring(
        manuscript,
        format_sigma(kepler_eta["snr"]),
        "Keplerian partialing SNR",
        errors,
    )
    require_substring(
        manuscript,
        format_sci_latex(joint_eta["eta"]),
        "Keplerian joint-model eta",
        errors,
    )
    require_substring(
        manuscript,
        format_sigma(joint_eta["snr"]),
        "Keplerian joint-model SNR",
        errors,
    )
    require_substring(
        manuscript,
        "R^2 \\approx 0.001",
        "Keplerian-only R-squared",
        errors,
    )
    require_substring(
        manuscript,
        format_sci_latex(cosd_eta),
        "cosD-only eta",
        errors,
    )
    require_substring(
        manuscript,
        format_sigma(cosd_snr),
        "cosD-only SNR",
        errors,
    )
    require_substring(
        manuscript,
        f"{total_obs:,}",
        "raw observation count",
        errors,
    )
    require_substring(
        manuscript,
        f"{n_clean:,}",
        "cleaned observation count",
        errors,
    )
    require_substring(
        manuscript,
        format_sci_latex(de430["eta"]),
        "DE430 cleaned eta",
        errors,
    )
    require_substring(
        manuscript,
        format_sigma(de430["snr"]),
        "DE430 cleaned SNR",
        errors,
    )
    require_substring(
        manuscript,
        f"\\chi^2 = {phase_chi2['chi2']:.1f}",
        "DE430 phase chi-square",
        errors,
    )
    require_substring(
        manuscript,
        format_sci_latex(inpop_integrator["eta"]),
        "Step 056 INPOP linearized integrator eta",
        errors,
    )
    require_substring(
        manuscript,
        format_sigma(inpop_integrator["snr"]),
        "Step 056 INPOP linearized integrator SNR",
        errors,
    )
    require_substring(
        manuscript,
        format_sci_latex(de430_integrator["eta"]),
        "Step 056 DE430 linearized integrator eta",
        errors,
    )
    require_substring(
        manuscript,
        format_sigma(de430_integrator["snr"]),
        "Step 056 DE430 linearized integrator SNR",
        errors,
    )
    require_substring(
        manuscript,
        format_sci_latex(delta_integrator["delta_eta"]),
        "Step 056 cross-ephemeris delta eta",
        errors,
    )
    require_substring(
        manuscript,
        format_sigma(primary_multiplicity["bonferroni_sigma"]),
        "Step 042 Bonferroni-adjusted primary SNR",
        errors,
    )
    require_substring(
        manuscript,
        format_sigma(step_046_full["full_systematic"]["snr"]),
        "Step 046 full-archive full-systematic SNR",
        errors,
    )
    require_substring(
        manuscript,
        format_sigma(step_046_equal_n["full_systematic"]["snr"]),
        "Step 046 equal-N full-systematic SNR",
        errors,
    )
    require_substring(
        manuscript,
        format_sigma(step_046_grasse_capped["full_systematic"]["snr"]),
        "Step 046 Grasse-capped full-systematic SNR",
        errors,
    )
    require_substring(
        manuscript,
        format_sci_latex(wls_eta),
        "Step 029 WLS eta",
        errors,
    )
    require_substring(
        manuscript,
        format_sigma(wls_snr),
        "Step 029 WLS SNR",
        errors,
    )
    require_substring(
        manuscript,
        f"{sky_scramble['n_scrambles']:,}",
        "Step 055 sky-scramble draw count",
        errors,
    )
    p_marginal = float(sky_scramble["p_scramble_delta_aic"])
    p_marginal_str = f"{p_marginal:.3f}"

    require_substring(
        manuscript,
        f"p = {p_marginal_str}",
        "Step 055 marginal sky-scramble p-value",
        errors,
    )
    marginal_needle = f"marginal ($p = {p_marginal_str}$)"
    marginal_alt = f"marginally specific (uniform $p = {p_marginal_str}$)"
    if marginal_needle not in manuscript and marginal_alt not in manuscript:
        errors.append(
            "Step 055 marginal sky-scramble interpretation: expected manuscript "
            f"source to contain {marginal_needle!r} or {marginal_alt!r}"
        )

    forbid_substring(manuscript, "figure-placeholder", "figure placeholders", errors)
    forbid_substring(manuscript, r"\chi^2 = 50.7", "DE430 phase chi-square stale value", errors)

    if errors:
        print("Manuscript value consistency check failed:")
        for error in errors:
            print(f"  - {error}")
        return 1

    print("Manuscript value consistency check passed.")
    print(
        "Verified primary estimand "
        f"eta={format_sci_latex(eta)} ({format_sigma(snr)}), "
        f"Keplerian partialing eta={format_sci_latex(kepler_eta['eta'])} "
        f"({format_sigma(kepler_eta['snr'])}), "
        f"N={n_clean:,}, DE430 chi^2={phase_chi2['chi2']:.1f}."
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
