#!/usr/bin/env python3
"""Load required canonical quantities from pipeline JSON outputs (no silent fallbacks)."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

PROJECT_ROOT = Path(__file__).resolve().parents[2]
OUTPUTS_DIR = PROJECT_ROOT / "results" / "outputs"


def load_json_required(filename: str) -> dict[str, Any]:
    path = OUTPUTS_DIR / filename
    if not path.is_file():
        raise FileNotFoundError(f"Required pipeline output missing: {path}")
    with path.open(encoding="utf-8") as handle:
        return json.load(handle)


def load_headline_eta() -> float:
    """Precision-weighted full-systematic η from Step 050 (primary manuscript estimand)."""
    step_050 = load_json_required("step_050_corrected_tep_analysis.json")
    try:
        return float(step_050["precision_weighted_full_systematic"]["eta"])
    except (KeyError, TypeError, ValueError) as exc:
        raise KeyError(
            "step_050_corrected_tep_analysis.json missing "
            "precision_weighted_full_systematic.eta"
        ) from exc


def load_full_systematic_ols_eta() -> float:
    """Unweighted full-systematic OLS η (sensitivity bound; Step 050 m5 or Step 040)."""
    step_050 = load_json_required("step_050_corrected_tep_analysis.json")
    try:
        return float(step_050["models"]["m5_full_corrected"]["eta"])
    except (KeyError, TypeError, ValueError) as exc:
        step_040 = load_json_required("step_040_unified_results_table.json")
        try:
            return float(step_040["primary_estimands"]["full_systematic_ols"]["eta"])
        except (KeyError, TypeError, ValueError) as exc2:
            raise KeyError(
                "Cannot resolve full-systematic OLS eta from step_050 or step_040"
            ) from exc2
