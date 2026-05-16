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


def compact_whitespace(text: str) -> str:
    return " ".join(text.split())


def require_substring(text: str, needle: str, label: str, errors: list[str]) -> None:
    if needle not in text:
        errors.append(f"{label}: expected manuscript source to contain {needle!r}")


def require_phrase(text: str, phrase: str, label: str, errors: list[str]) -> None:
    if compact_whitespace(phrase) not in compact_whitespace(text):
        errors.append(f"{label}: expected manuscript source to contain phrase {phrase!r}")


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
    step_064_pi = load_json(OUTPUTS_DIR / "step_064_prediction_coverage.json")
    step_029 = load_json(OUTPUTS_DIR / "step_029_station_power_analysis.json")
    step_059 = load_json(OUTPUTS_DIR / "step_059_grasse_systematic_sufficiency.json")
    step_061 = load_json(OUTPUTS_DIR / "step_061_systematic_sensitivity_analysis.json")
    step_065 = load_json(OUTPUTS_DIR / "step_065_high_dimensional_absorption_test.json")
    step_073_path = OUTPUTS_DIR / "step_073_laplace_bayes_factor.json"
    step_073 = load_json(step_073_path) if step_073_path.is_file() else None
    step_072 = load_json(OUTPUTS_DIR / "step_072_leave_one_station_out_meta.json")
    step_044 = load_json(OUTPUTS_DIR / "step_044_systematic_projection_analysis.json")

    primary = step_040["primary_estimands"]["precision_weighted_full_systematic"]
    cooks_excised = step_040["primary_estimands"]["cooks_excised_full_systematic"]
    leverage_cosd = step_040["primary_estimands"]["leverage_excised_ols"]
    mcmc_est = step_040["primary_estimands"]["bayesian_mcmc"]
    sensitivity_full = step_040["primary_estimands"]["full_systematic_ols"]
    cosd_only = step_040["primary_estimands"]["full_sample_ols"]
    de430 = step_006b["threshold_sweep"]["6sigma"]
    phase_chi2 = step_006b["phase_chi_square"]
    inpop_integrator = step_056["inpop19a"]["linearized_integrator_eta"]
    de430_integrator = step_056["de430"]["linearized_integrator_eta"]
    delta_integrator = step_056["cross_ephemeris"]

    primary_multiplicity = next(
        (
            test
            for test in step_042["tests"]
            if test.get("category") == "primary"
            and test.get("analysis_type") == "independent_hypothesis"
            and "primary headline" in test.get("name", "").lower()
        ),
        None,
    )
    if primary_multiplicity is None:
        primary_multiplicity = next(
            test
            for test in step_042["tests"]
            if test.get("category") == "primary"
            and test.get("analysis_type") == "independent_hypothesis"
        )
    step_046_full = next(test for test in step_046["tests"] if test["name"] == "full_sample")
    step_046_equal_n = next(
        test for test in step_046["tests"] if test["name"] == "equal_n_subsample"
    )
    step_046_grasse_capped = next(
        test for test in step_046["tests"] if test["name"] == "grasse_capped"
    )
    sky_scramble = step_055["sky_scrambling"]
    sky_matched = sky_scramble["correlation_matched"]
    wls_eta = step_029["summary"]["precision_weighted_eta"]
    wls_snr = step_029["summary"]["precision_weighted_snr"]

    total_obs = step_001["combined"]["total_obs"]
    n_clean = primary["n_obs"]
    n_cooks_excised = cooks_excised["n_obs"]

    eta = primary["eta"]
    eta_error = primary["eta_error"]
    snr = primary["snr"]
    cluster_snr = primary["snr_cluster"]
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
        format_sci_latex(sensitivity_full["eta"]),
        "non-excised full-systematic sensitivity eta",
        errors,
    )
    require_substring(
        manuscript,
        format_sigma(sensitivity_full["snr"]),
        "non-excised full-systematic sensitivity SNR",
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
        "cleaned observation count (headline full archive)",
        errors,
    )
    require_substring(
        manuscript,
        "23{,}837",
        "Cook-excised observation count",
        errors,
    )
    require_substring(
        manuscript,
        format_sci_latex(cooks_excised["eta"]),
        "Cook-excised full-systematic eta",
        errors,
    )
    require_substring(
        manuscript,
        format_sci_latex(cooks_excised["eta_error"]),
        "Cook-excised full-systematic eta uncertainty",
        errors,
    )
    require_substring(
        manuscript,
        format_sigma(cooks_excised["snr"]),
        "Cook-excised OLS SNR",
        errors,
    )
    require_substring(
        manuscript,
        format_sigma(cooks_excised["snr_cluster"]),
        "Cook-excised cluster-robust SNR",
        errors,
    )
    require_substring(
        manuscript,
        format_sci_latex(leverage_cosd["eta"]),
        "Cook-excised cosD-only eta (Step 017 diagnostic)",
        errors,
    )
    require_substring(
        manuscript,
        format_sigma(leverage_cosd["snr"]),
        "Cook-excised cosD-only SNR (Step 017 diagnostic)",
        errors,
    )
    require_substring(
        manuscript,
        format_sci_latex(mcmc_est["eta"]),
        "Bayesian MCMC eta",
        errors,
    )
    require_substring(
        manuscript,
        format_sci_latex(mcmc_est["eta_error"]),
        "Bayesian MCMC eta uncertainty",
        errors,
    )
    require_substring(
        manuscript,
        format_sigma(mcmc_est["snr"]),
        "Bayesian MCMC SNR",
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
    p_matched = float(sky_matched["p_scramble_f_eff_matched"])
    p_matched_str = f"{p_matched:.3f}"

    # Manuscript reports sky-scramble as p_F ≈ 0.095 (not always "p = 0.095").
    marginal_needles = [
        f"p = {p_marginal_str}",
        f"p_F = {p_marginal_str}",
        f"p_F \\approx {p_marginal_str}",
        f"($p_F \\approx {p_marginal_str}$",
        f"marginally specific ($p_F \\approx {p_marginal_str}$",
    ]
    if not any(needle in manuscript for needle in marginal_needles):
        errors.append(
            "Step 055 marginal sky-scramble p-value: expected manuscript source to "
            f"contain one of {marginal_needles!r}"
        )
    matched_needles = [
        f"correlation-matched $p = {p_matched_str}$",
        f"correlation-matched $p_F \\approx {p_matched_str}$",
        f"correlation-matched $p_F \\approx {p_matched_str}",
    ]
    if not any(needle in manuscript for needle in matched_needles):
        errors.append(
            "Step 055 correlation-matched p-value: expected manuscript source to "
            f"contain one of {matched_needles!r}"
        )
    require_substring(
        manuscript,
        f"{step_055['falsification_tests_passed']}/5 diagnostics",
        "Step 055 diagnostic count",
        errors,
    )
    refined = sky_scramble.get("refined_directional_nulls", {})
    if refined:
        gs_p = float(refined["orthogonal_scramble_null"]["p_f_eff"])
        require_substring(
            manuscript,
            "6 \\times 10^{-4}",
            "Step 055 phase-null p-value",
            errors,
        )
        require_substring(
            manuscript,
            f"{gs_p:.3f}",
            "Step 055 GS orthogonal scramble p-value",
            errors,
        )
    require_phrase(
        manuscript,
        "strong residual-channel evidence consistent with TEP",
        "positive TEP framing (claim hierarchy)",
        errors,
    )
    require_phrase(
        manuscript,
        "modeled known-systematic-only explanations fail",
        "known-systematic framing",
        errors,
    )
    require_phrase(
        manuscript,
        "bounded station-leverage risk",
        "bounded station leverage framing",
        errors,
    )
    require_phrase(
        manuscript,
        "no catalogued systematic at its observed amplitude reproduces the headline",
        "abstract caveat containment",
        errors,
    )
    require_phrase(
        manuscript,
        "integrator refit remains the open closure",
        "abstract integrator closure",
        errors,
    )

    for step_id, status in (
        ("step_029", step_029["status"]),
        ("step_059", step_059["status"]),
        ("step_065", step_065["status"]),
    ):
        if status == "WARNING":
            errors.append(
                f"{step_id}: legacy top-level WARNING status; use explicit risk fields with PASS"
            )
    if step_061["status"] != "PASS":
        errors.append(f"step_061: expected PASS status, got {step_061['status']!r}")

    blind = step_061.get("blind_year_holdout", {}).get("inverse_variance_combined", {})
    pca_joint = step_061.get("adversarial_pca", {}).get("joint_cosd_plus_all_pcs", {})
    require_substring(
        manuscript,
        format_sci_latex(float(blind["eta"]), 2),
        "Step 061 blind hold-out eta",
        errors,
    )
    require_substring(
        manuscript,
        format_sigma(float(blind["snr"]), 1),
        "Step 061 blind hold-out SNR",
        errors,
    )
    require_substring(
        manuscript,
        format_sci_latex(float(pca_joint["eta"]), 2),
        "Step 061 adversarial PCA joint eta",
        errors,
    )
    require_substring(
        manuscript,
        str(step_042["n_tests"]),
        "Step 042 n_tests in manuscript",
        errors,
    )
    forbid_substring(
        manuscript,
        "20 significance measures, of which four",
        "stale Step 042 test count (20)",
        errors,
    )
    forbid_substring(
        manuscript,
        "sixteen as sensitivity analyses",
        "stale Step 042 sensitivity count (16)",
        errors,
    )

    ledger_path = OUTPUTS_DIR / "tep_evidence_ledger.json"
    if ledger_path.is_file():
        ledger = load_json(ledger_path)
        pillars = ledger.get("positive_evidence", [])
        if not pillars:
            errors.append("tep_evidence_ledger.json: missing positive_evidence pillars")
        else:
            first = pillars[0]
            if "precision-weighted" not in first.get("pillar", "").lower():
                errors.append(
                    "tep_evidence_ledger.json: first pillar must be precision-weighted WLS headline"
                )
            ledger_eta = first.get("result", {}).get("eta")
            if ledger_eta is None or abs(float(ledger_eta) - float(primary["eta"])) > 1e-9:
                errors.append(
                    "tep_evidence_ledger.json: primary pillar eta does not match step_040 PW headline"
                )
        directional = next(
            (p for p in pillars if "directional" in p.get("pillar", "").lower()),
            None,
        )
        if directional is None:
            errors.append("tep_evidence_ledger.json: missing directional residual structure pillar")
        elif "dual_axis_r_cmb_gal" not in directional.get("result", {}):
            errors.append("tep_evidence_ledger.json: directional pillar missing dual_axis_r_cmb_gal")
        loso_pillar = next(
            (p for p in pillars if "leave-one-station-out" in p.get("pillar", "").lower()),
            None,
        )
        if loso_pillar is None:
            errors.append("tep_evidence_ledger.json: missing LOSO meta pillar")
    else:
        errors.append("tep_evidence_ledger.json: missing (run generate_evidence_ledger.py)")

    pi_cov = step_064_pi["prediction_interval_coverage"]["wls_published_sigma"]["coverage"]
    pi_chi2 = step_064_pi["prediction_interval_coverage"]["wls_published_sigma"]["chi2_reduced"]
    require_substring(manuscript, "Step 064-PI", "Step 064-PI uncertainty calibration", errors)
    require_substring(manuscript, "Step 064-SRP", "Step 064-SRP SRP check", errors)
    require_substring(
        manuscript,
        f"{pi_cov['68pct']['observed_coverage'] * 100:.1f}%",
        "Step 064-PI 68% observed coverage",
        errors,
    )
    require_substring(
        manuscript,
        f"{pi_chi2:.2f}".rstrip("0").rstrip("."),
        "Step 064-PI chi2_reduced",
        errors,
    )
    if step_064_pi.get("step_id") not in ("step_064_pi", "step_064"):
        errors.append(
            f"step_064_prediction_coverage.json: unexpected step_id {step_064_pi.get('step_id')!r}"
        )

    forbid_substring(
        manuscript,
        "Orbital SRP Scaling (Step 064)",
        "ambiguous Step 064 SRP heading (use 064-SRP)",
        errors,
    )
    forbid_substring(manuscript, "pipeline of 73 steps", "stale 73-step pipeline count", errors)
    forbid_substring(manuscript, "figure-placeholder", "figure placeholders", errors)
    forbid_substring(manuscript, r"\chi^2 = 50.7", "DE430 phase chi-square stale value", errors)
    forbid_substring(manuscript, "composite sky-scramble scoring rule", "stale CMB scoring language", errors)
    forbid_substring(manuscript, "formally excluded", "overclaim exclusion language", errors)
    forbid_substring(manuscript, "quantitatively falsified", "overclaim falsification language", errors)
    forbid_substring(manuscript, "exhaust the nominated conventional explanations", "overclaim exhaustive language", errors)
    forbid_substring(
        manuscript,
        "BF}_{10} = 1.2 \\times 10^{14}",
        "stale Laplace headline Bayes factor in abstract",
        errors,
    )
    forbid_substring(manuscript, "4.8.6 Leave-One-Station-Out", "duplicate LOSO section 4.8.6", errors)
    forbid_substring(manuscript, "P(\\eta<0)=0.943", "stale bootstrap P(eta<0) (use 0.945)", errors)
    forbid_substring(manuscript, "Sections 4.0, 4.26, and 5.5", "wrong Grasse cross-ref to 4.26", errors)
    forbid_substring(
        manuscript,
        "not alternative explanations that reproduce the observed signal",
        "stale abstract overclaim on alternatives",
        errors,
    )

    meta = step_072["meta_analysis"]
    meta_eta_str = format_sci_latex(float(meta["eta"]), 2)
    require_substring(manuscript, meta_eta_str, "Step 072 meta eta", errors)
    require_substring(manuscript, "12.8\\sigma", "Step 072 meta SNR (rounded)", errors)
    require_substring(manuscript, "fig_05_loso_forest.png", "LOSO forest figure path", errors)
    require_substring(manuscript, "Step 072", "Step 072 LOSO reference", errors)
    require_substring(manuscript, "0.984", "dual-axis collinearity r", errors)
    forbid_substring(
        manuscript,
        "\\mathcal{B}_{\\rm TEP,GR} = 74.3",
        "stale Savage-Dickey headline value",
        errors,
    )

    phase_diff = step_044.get("phase_locked_differential", {})
    if phase_diff:
        require_substring(
            manuscript,
            format_sci_latex(phase_diff["eta_differential"]),
            "Step 044 phase-locked differential eta",
            errors,
        )
        require_substring(
            manuscript,
            format_sigma(phase_diff["snr_differential"]),
            "Step 044 phase-locked differential SNR",
            errors,
        )

    if step_073 is not None:
        cross = step_073.get("evidence_cross_checks", {})
        grid_bf = cross.get("grid_quadrature", {}).get("bf10")
        bridge_bf = cross.get("bridge_sampling", {}).get("bf10")
        if grid_bf is not None:
            grid_rounded = round(float(grid_bf))
            require_substring(
                manuscript,
                f"BF}}_{{10}} \\approx {grid_rounded}",
                "Step 073 grid quadrature BF (rounded)",
                errors,
            )
        if bridge_bf is not None:
            bridge_rounded = round(float(bridge_bf))
            require_substring(
                manuscript,
                f"\\approx {bridge_rounded}",
                "Step 073 bridge sampling BF (rounded)",
                errors,
            )

    if errors:
        print("Manuscript value consistency check failed:")
        for error in errors:
            print(f"  - {error}")
        return 1

    print("Manuscript value consistency check passed.")
    print(
        "Verified headline precision-weighted estimand "
        f"eta={format_sci_latex(eta)} ({format_sigma(snr)}; "
        f"{format_sigma(cluster_snr)} cluster-robust), "
        f"Cook-excised diagnostic eta={format_sci_latex(cooks_excised['eta'])} "
        f"({format_sigma(cooks_excised['snr'])}), "
        f"N={n_clean:,} (Cook-excised N={n_cooks_excised:,}), "
        f"Keplerian partialing eta={format_sci_latex(kepler_eta['eta'])} "
        f"({format_sigma(kepler_eta['snr'])}), "
        f"DE430 chi^2={phase_chi2['chi2']:.1f}."
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
