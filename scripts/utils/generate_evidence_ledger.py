#!/usr/bin/env python3
"""Generate a compact evidence ledger from canonical pipeline outputs."""

from __future__ import annotations

import json
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]
OUTPUTS_DIR = PROJECT_ROOT / "results" / "outputs"


def load(name: str) -> dict:
    with (OUTPUTS_DIR / name).open(encoding="utf-8") as handle:
        return json.load(handle)


def sci(value: float) -> str:
    return f"{value:.3e}"


def sigma(value: float) -> str:
    return f"{value:.2f}sigma"


def main() -> int:
    step_040 = load("step_040_unified_results_table.json")
    step_050 = load("step_050_corrected_tep_analysis.json")
    step_029 = load("step_029_station_power_analysis.json")
    step_055 = load("step_055_cmb_rigorous_falsification.json")
    step_056 = load("step_056_dynamical_integrator_eta_refit.json")
    step_059 = load("step_059_grasse_systematic_sufficiency.json")
    step_071 = load("step_071_stratified_equal_n.json")
    step_072 = load("step_072_leave_one_station_out_meta.json")
    step_061 = load("step_061_systematic_sensitivity_analysis.json")
    step_064_pi = load("step_064_prediction_coverage.json")
    step_065 = load("step_065_high_dimensional_absorption_test.json")

    headline = step_040["primary_estimands"]["precision_weighted_full_systematic"]
    cooks_excised = step_040["primary_estimands"]["cooks_excised_full_systematic"]
    common_eta = step_050["station_univ"]["common_eta_station_systematics"]
    f_test = step_050["station_univ"]["f_test"]
    systematics = step_061["systematics"]
    min_ratio_name, min_ratio = min(
        (
            (name, entry["ratio_required_to_known"])
            for name, entry in systematics.items()
        ),
        key=lambda item: item[1],
    )
    max_p_exceed = max(
        entry["monte_carlo"]["p_exceed_observed"]
        for entry in systematics.values()
    )

    ledger = {
        "summary": "Strong residual-channel evidence for TEP with bounded station and source-level-refit risks.",
        "positive_evidence": [
            {
                "pillar": "Primary synodic estimand (precision-weighted WLS)",
                "source": "step_040_unified_results_table.json / step_050_corrected_tep_analysis.json",
                "result": {
                    "eta": headline["eta"],
                    "eta_error": headline["eta_error"],
                    "snr": headline["snr"],
                    "snr_cluster": headline["snr_cluster"],
                    "n_obs": headline["n_obs"],
                },
                "interpretation": "Headline Nordtvedt estimate from precision-weighted full-systematic WLS on all cleaned shots (no row deletion).",
            },
            {
                "pillar": "Uncertainty calibration (prediction intervals and η brackets)",
                "source": "step_064_prediction_coverage.json",
                "result": {
                    "wls_68pct_observed": step_064_pi["abstract_linkage"]["wls_68pct_observed"],
                    "wls_95pct_observed": step_064_pi["abstract_linkage"]["wls_95pct_observed"],
                    "chi2_reduced_wls": step_064_pi["abstract_linkage"]["chi2_reduced_wls"],
                    "headline_snr_wls": step_064_pi["abstract_linkage"]["headline_snr_wls"],
                    "sigma_calibration_scale": step_064_pi["headline_eta_intervals"][
                        "sigma_calibration_scale_68pct"
                    ],
                    "loso_conformal_95_excludes_zero": step_064_pi["headline_eta_intervals"][
                        "loso_conformal"
                    ]["intervals"]["95pct"]["ci_upper"]
                    < 0.0,
                },
                "interpretation": step_064_pi["interpretation"],
            },
            {
                "pillar": "Leverage diagnostic and cross-station support",
                "source": "step_040_unified_results_table.json / step_050_corrected_tep_analysis.json",
                "result": {
                    "cooks_excised_eta": cooks_excised["eta"],
                    "cooks_excised_snr": cooks_excised["snr"],
                    "cooks_excised_n_obs": cooks_excised["n_obs"],
                    "common_eta": common_eta["eta"],
                    "common_eta_snr": common_eta["snr"],
                    "station_deviation_p": f_test["p"],
                },
                "interpretation": "Cook's-Distance excised OLS and common-eta mixed model remain negative and significant; sign stable under leverage removal.",
            },
            {
                "pillar": "Known-systematic-only simulations",
                "source": "step_061_systematic_sensitivity_analysis.json",
                "result": {
                    "required_amplitude_cm": step_061["required_amplitude_cm"],
                    "minimum_required_to_known_ratio": min_ratio,
                    "minimum_ratio_systematic": min_ratio_name,
                    "max_p_exceed_observed": max_p_exceed,
                    "n_mc_per_systematic": next(iter(systematics.values()))["monte_carlo"]["n_mc"],
                },
                "interpretation": "No modeled known systematic at observed amplitude reproduces the full-systematic eta.",
            },
            {
                "pillar": "Adversarial nuisance + blind hold-out (Step 061)",
                "source": "step_061_systematic_sensitivity_analysis.json",
                "result": {
                    "pca_joint_eta": step_061["adversarial_pca"]["joint_cosd_plus_all_pcs"]["eta"],
                    "pca_joint_snr": step_061["adversarial_pca"]["joint_cosd_plus_all_pcs"]["snr"],
                    "gp_absorption_fraction": step_061["adversarial_gp"]["absorption_fraction_gp_subtract"],
                    "blind_holdout_eta": step_061["blind_year_holdout"]["inverse_variance_combined"]["eta"],
                    "blind_holdout_snr": step_061["blind_year_holdout"]["inverse_variance_combined"]["snr"],
                    "interaction_cells_fitted": step_061["interaction_grid"]["n_cells_fitted"],
                    "interaction_cells_negative_eta": step_061["interaction_grid"]["n_cells_negative_eta"],
                },
                "interpretation": (
                    "Data-driven PCA does not absorb cos(D); blind year hold-out "
                    "recovers significant negative eta with nuisances trained only on "
                    "non-held-out years."
                ),
            },
            {
                "pillar": "Linearized post-fit extraction",
                "source": "step_056_dynamical_integrator_eta_refit.json",
                "result": {
                    "inpop_eta": step_056["inpop19a"]["linearized_integrator_eta"]["eta"],
                    "inpop_snr": step_056["inpop19a"]["linearized_integrator_eta"]["snr"],
                    "de430_eta": step_056["de430"]["linearized_integrator_eta"]["eta"],
                    "de430_snr": step_056["de430"]["linearized_integrator_eta"]["snr"],
                },
                "interpretation": "Published residual archives recover a consistent negative eta under the same full-systematic nuisance design.",
            },
            {
                "pillar": "Directional residual structure",
                "source": "step_055_cmb_rigorous_falsification.json",
                "result": {
                    "diagnostics_passed": step_055["falsification_tests_passed"],
                    "diagnostics_total": step_055["falsification_tests_total"],
                    "sky_scramble_p": step_055["sky_scrambling"]["p_scramble_f"],
                    "correlation_matched_p": step_055["sky_scrambling"]["correlation_matched"]["p_scramble_f_matched"],
                    "phase_null_p_eff": step_055["sky_scrambling"]["refined_directional_nulls"]["phase_null"]["p_f_eff"],
                    "orthogonal_scramble_p_eff": step_055["sky_scrambling"]["refined_directional_nulls"]["orthogonal_scramble_null"]["p_f_eff"],
                    "dual_axis_r_cmb_gal": step_055["dual_axis_identifiability"]["cmb_galactic_cos_correlation"],
                    "dual_eta_cmb_t": step_055["dual_axis_identifiability"]["dual_axis"]["eta_cmb_t"],
                    "dual_eta_gal_t": step_055["dual_axis_identifiability"]["dual_axis"]["eta_gal_t"],
                },
                "interpretation": (
                    "Corroborative fixed-sky directional anatomy on the residual channel, not a replacement eta estimator. "
                    "Synodic phase coupling is decisively rejected; uniform-axis uniqueness remains marginal. "
                    "Dual-axis fit: r(cos theta_CMB, cos theta_gal)=0.984; Planck term absorbed when both axes are included."
                ),
            },
            {
                "pillar": "Leave-one-station-out meta (Step 072)",
                "source": "step_072_leave_one_station_out_meta.json",
                "result": {
                    "eta_meta": step_072["meta_analysis"]["eta"],
                    "eta_meta_err": step_072["meta_analysis"]["eta_err"],
                    "snr_meta": step_072["meta_analysis"]["snr"],
                    "cochrans_Q": step_072["meta_analysis"]["cochrans_Q"],
                    "I2_percent": step_072["meta_analysis"]["I2_percent"],
                    "n_powered_loso": step_072["power_accounting"]["n_powered"],
                },
                "interpretation": (
                    "Four of five LOSO exclusions remain powered with the same negative sign; "
                    "excluding Grasse is underpowered. Inverse-variance meta on powered exclusions: "
                    f"{sigma(step_072['meta_analysis']['snr'])}, I2=0%."
                ),
            },
        ],
        "bounded_risks": [
            {
                "risk": "Station leverage (LOSO + Grasse-conditioned)",
                "source": (
                    "step_029_station_power_analysis.json / "
                    "step_059_grasse_systematic_sufficiency.json / "
                    "step_072_leave_one_station_out_meta.json / "
                    "step_071_stratified_equal_n.json"
                ),
                "result": {
                    "loso_powered_exclusions": step_072["power_accounting"]["n_powered"],
                    "loso_meta_snr": step_072["meta_analysis"]["snr"],
                    "exclude_grasse_snr": step_072["grasse_leverage"]["exclude_grasse_snr"],
                    "non_grasse_eta": step_059["grasse_conditioned_estimand"]["non_grasse_direct"]["eta"],
                    "non_grasse_snr": step_059["grasse_conditioned_estimand"]["non_grasse_direct"]["snr_cluster"]
                    or step_059["grasse_conditioned_estimand"]["non_grasse_direct"]["snr"],
                    "equal_n_sign_retained": step_071["stability_report"]["eta_snr_stable_under_balance"][
                        "sign_retained_all_subsamples"
                    ],
                },
                "status": [step_029["status"], step_059["status"], step_072["status"], step_071["status"]],
                "interpretation": (
                    "Grasse leverage is material: excluding Grasse is underpowered, but four other "
                    "LOSO exclusions remain powered with the same negative sign (meta-analysis "
                    f"{sigma(step_072['meta_analysis']['snr'])}). Non-Grasse and Grasse-conditioned "
                    "estimands retain negative eta; equal-N balance reduces SNR without sign reversal."
                ),
            },
            {
                "risk": "Source-level absorption/refit",
                "source": "step_065_high_dimensional_absorption_test.json",
                "status": step_065["status"],
                "interpretation": "High-dimensional residual-basis stress test supports sideband survival but preserves the need for source-level INPOP/DE430 eta-free refits.",
            },
        ],
    }

    json_path = OUTPUTS_DIR / "tep_evidence_ledger.json"
    md_path = OUTPUTS_DIR / "tep_evidence_ledger.md"
    json_path.write_text(json.dumps(ledger, indent=2) + "\n", encoding="utf-8")

    lines = [
        "# TEP Evidence Ledger",
        "",
        ledger["summary"],
        "",
        "## Positive Evidence",
        "",
    ]
    for item in ledger["positive_evidence"]:
        lines.extend([
            f"### {item['pillar']}",
            f"- Source: `{item['source']}`",
            f"- Interpretation: {item['interpretation']}",
        ])
        for key, value in item["result"].items():
            if isinstance(value, float):
                value_text = sigma(value) if "snr" in key else sci(value)
            else:
                value_text = str(value)
            lines.append(f"- {key}: {value_text}")
        lines.append("")

    lines.extend(["## Bounded Risks", ""])
    for item in ledger["bounded_risks"]:
        lines.extend([
            f"### {item['risk']}",
            f"- Source: `{item['source']}`",
            f"- Status: {item['status']}",
            f"- Interpretation: {item['interpretation']}",
            "",
        ])

    md_path.write_text("\n".join(lines), encoding="utf-8")
    print(f"Wrote {json_path.relative_to(PROJECT_ROOT)}")
    print(f"Wrote {md_path.relative_to(PROJECT_ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
