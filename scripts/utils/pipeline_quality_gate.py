#!/usr/bin/env python3
"""Run the core manuscript and pipeline quality gates.

This utility is intentionally small and reviewer-facing: it collects the
checks that defend the headline TEP evidence spine from stale values, schema
drift, accidental overclaim wording, and unexpected legacy WARNING statuses in step JSON.
"""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]
OUTPUTS_DIR = PROJECT_ROOT / "results" / "outputs"


def run_check(label: str, command: list[str]) -> bool:
    print(f"\n=== {label} ===")
    result = subprocess.run(command, cwd=PROJECT_ROOT, text=True)
    if result.returncode == 0:
        print(f"PASS: {label}")
        return True
    print(f"FAIL: {label} (exit {result.returncode})")
    return False


def load_status(path: Path) -> str | None:
    with path.open(encoding="utf-8") as handle:
        return json.load(handle).get("status")


def check_warning_ledger() -> bool:
    """Fail if any step output still uses legacy top-level status=WARNING."""
    print("\n=== Warning Ledger (legacy status scan) ===")
    ok = True
    warned: list[str] = []

    for path in sorted(OUTPUTS_DIR.glob("step_*.json")):
        try:
            status = load_status(path)
        except Exception as exc:
            print(f"FAIL: could not read {path.name}: {exc}")
            ok = False
            continue
        if status == "WARNING":
            warned.append(path.name)

    if warned:
        print(
            "FAIL: top-level status=WARNING is retired for bounded-risk reporting; "
            f"use explicit *_result / risk_flags fields with PASS. Found: {', '.join(warned)}"
        )
        ok = False
    else:
        print("PASS: no step_*.json uses legacy top-level WARNING status")

    return ok


def main() -> int:
    checks = [
        run_check("Schema validation", [sys.executable, "scripts/utils/schema_validation.py"]),
        run_check("Evidence ledger generation", [sys.executable, "scripts/utils/generate_evidence_ledger.py"]),
        run_check("Manuscript value/framing consistency", [sys.executable, "scripts/utils/verify_value_consistency.py"]),
        run_check("Python compile check", [sys.executable, "-m", "compileall", "scripts"]),
        check_warning_ledger(),
    ]

    if all(checks):
        print("\nQUALITY GATE PASS: pipeline and manuscript evidence spine are internally consistent.")
        return 0

    print("\nQUALITY GATE FAIL: resolve the items above before relying on the manuscript.")
    return 1


if __name__ == "__main__":
    sys.exit(main())
