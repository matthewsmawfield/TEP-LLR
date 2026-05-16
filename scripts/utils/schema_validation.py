#!/usr/bin/env python3
"""Validate canonical pipeline JSON outputs and required result fields."""

from __future__ import annotations

import json
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]
OUTPUTS_DIR = PROJECT_ROOT / "results" / "outputs"

REQUIRED_OUTPUTS = {
    "step_001_data_preprocessing.json": ["combined"],
    "step_040_unified_results_table.json": ["primary_estimands"],
    "step_046_station_balanced_tep.json": [
        "status",
        "stress_test_result",
        "primary_full_systematic_significant",
        "balanced_full_systematic_all_significant",
    ],
    "step_050_corrected_tep_analysis.json": [
        "models",
        "station_univ",
        "precision_weighted_full_systematic",
        "cooks_excised_full_systematic",
    ],
    "step_055_cmb_rigorous_falsification.json": ["status", "sky_scrambling"],
    "step_056_dynamical_integrator_eta_refit.json": ["inpop19a", "de430", "cross_ephemeris"],
    "step_057_haleakala_null_fluctuation.json": [
        "status",
        "haleakala_simulation_tep",
        "haleakala_simulation_gr",
        "interpretation",
    ],
}

# Pipeline steps that must record a top-level PASS for CI / manuscript integrity.
REQUIRE_TOP_LEVEL_PASS = frozenset(
    {
        "step_055_cmb_rigorous_falsification.json",
        "step_056_dynamical_integrator_eta_refit.json",
        "step_057_haleakala_null_fluctuation.json",
        "step_040_unified_results_table.json",
        "step_046_station_balanced_tep.json",
        "step_050_corrected_tep_analysis.json",
        "step_042_multiple_testing_correction.json",
    }
)

OPTIONAL_OUTPUTS = {
    "step_012_station_dominance.json": ["full_sample_eta", "dominant_station_jackknife"],
}

VALID_TOP_LEVEL_STATUSES = {"PASS", "WARNING", "FAIL", "ERROR"}


def load_json(path: Path) -> dict:
    with path.open(encoding="utf-8") as handle:
        return json.load(handle)


def validate_required_keys(payload: dict, required_keys: list[str], label: str, errors: list[str]) -> None:
    for key in required_keys:
        if key not in payload:
            errors.append(f"{label}: missing required key {key!r}")


def main() -> int:
    errors: list[str] = []

    for filename, required_keys in REQUIRED_OUTPUTS.items():
        path = OUTPUTS_DIR / filename
        if not path.exists():
            errors.append(f"Missing required pipeline output: {path}")
            continue
        payload = load_json(path)
        validate_required_keys(payload, required_keys, filename, errors)
        if filename == "step_040_unified_results_table.json":
            pw = payload.get("primary_estimands", {}).get(
                "precision_weighted_full_systematic"
            )
            if pw is None:
                errors.append(
                    f"{filename}: primary_estimands missing "
                    "precision_weighted_full_systematic (headline estimand)"
                )
            elif "PRIMARY HEADLINE" not in str(pw.get("status", "")):
                errors.append(
                    f"{filename}: precision_weighted_full_systematic.status "
                    "should mark PRIMARY HEADLINE estimand"
                )
        status_raw = payload.get("status")
        if status_raw is not None:
            status = str(status_raw).strip().upper()
            if status not in VALID_TOP_LEVEL_STATUSES:
                errors.append(
                    f"{filename}: non-standard top-level status {status_raw!r}; "
                    f"expected one of {sorted(VALID_TOP_LEVEL_STATUSES)}"
                )
        if filename in REQUIRE_TOP_LEVEL_PASS:
            status = str(payload.get("status", "")).strip().upper()
            if status != "PASS":
                errors.append(
                    f"{filename}: expected top-level status 'PASS', found {payload.get('status', '')!r}"
                )

    for filename, required_keys in OPTIONAL_OUTPUTS.items():
        path = OUTPUTS_DIR / filename
        if not path.exists():
            continue
        payload = load_json(path)
        validate_required_keys(payload, required_keys, filename, errors)

    for path in sorted(OUTPUTS_DIR.glob("step_*.json")):
        payload = load_json(path)
        if not isinstance(payload, dict) or "status" not in payload:
            continue
        status = str(payload["status"]).strip().upper()
        if status not in VALID_TOP_LEVEL_STATUSES:
            errors.append(
                f"{path.name}: non-standard top-level status {payload['status']!r}; "
                f"expected one of {sorted(VALID_TOP_LEVEL_STATUSES)}"
            )

    if errors:
        print("Pipeline schema validation failed:")
        for error in errors:
            print(f"  - {error}")
        return 1

    print("Pipeline schema validation passed.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
