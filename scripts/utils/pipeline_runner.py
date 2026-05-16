#!/usr/bin/env python3
"""
Shared pipeline runner utility for TEP-LLR pipeline.

Provides run_step() and run_pipeline() for executing the canonical
60-step Lunar Laser Ranging analysis pipeline.
"""

import datetime
import os
import subprocess
import sys
import time
import json
import platform
import hashlib
import psutil
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]
STEPS_DIR = PROJECT_ROOT / "scripts" / "steps"
OUTPUTS_DIR = PROJECT_ROOT / "results" / "outputs"

OUTPUT_NAME_OVERRIDES = {
    # Historical output stem retained for manuscript compatibility.
    "step_043_temporal_bin_variation_analysis.py": "step_043_temporal_bin_variation.json",
}


def _fmt_elapsed(seconds: float) -> str:
    if seconds < 60:
        return f"{seconds:.1f}s"
    m, s = divmod(int(seconds), 60)
    return f"{m}m{s:02d}s"


def _get_file_hash(path: Path) -> str:
    """Compute SHA-256 hash of a file for reproducibility auditing."""
    if not path.exists():
        return "MISSING"
    sha256 = hashlib.sha256()
    with open(path, "rb") as f:
        while chunk := f.read(8192):
            sha256.update(chunk)
    return sha256.hexdigest()[:12]


def _expected_output_path(script_name: str) -> Path:
    """Return the canonical JSON output path for a pipeline step."""
    output_name = OUTPUT_NAME_OVERRIDES.get(script_name, script_name.replace(".py", ".json"))
    return OUTPUTS_DIR / output_name


def _validate_step_output(output_path: Path, pre_run_mtime: float | None) -> tuple[bool, str]:
    """
    Verify that a step wrote a fresh, non-failing JSON output.

    This prevents false PASS audit entries when a subprocess exits 0 after
    returning None or writing no output, which would otherwise leave a stale
    artifact for downstream steps.
    """
    if not output_path.exists():
        return False, f"missing expected output {output_path.relative_to(PROJECT_ROOT)}"

    post_run_mtime = output_path.stat().st_mtime
    if pre_run_mtime is not None and post_run_mtime <= pre_run_mtime:
        return False, f"stale output {output_path.relative_to(PROJECT_ROOT)} was not updated"

    try:
        with output_path.open(encoding="utf-8") as handle:
            payload = json.load(handle)
    except json.JSONDecodeError as exc:
        return False, f"invalid JSON in {output_path.relative_to(PROJECT_ROOT)}: {exc}"

    if payload is None:
        return False, f"{output_path.relative_to(PROJECT_ROOT)} contains JSON null"

    if isinstance(payload, dict):
        status = str(payload.get("status", "")).upper()
        if status in {"FAIL", "FAILED", "ERROR"}:
            return False, f"{output_path.relative_to(PROJECT_ROOT)} reports status={status}"
        if status == "WARNING":
            return True, "fresh output verified with status=WARNING"

    return True, "fresh output verified"


def generate_audit_report(pipeline_name: str, results: list, elapsed_total: float) -> Path:
    """
    Generate a formal Research-Grade Audit Report (JSON).
    Captures system state, dependency hashes, and pipeline telemetry.
    """
    if any(r["status"] == "FAIL" for r in results):
        audit_status = "FAIL"
    elif any(r["status"] == "WARNING" for r in results):
        audit_status = "WARNING"
    else:
        audit_status = "PASS"

    audit_data = {
        "report_metadata": {
            "pipeline_name": pipeline_name,
            "timestamp": datetime.datetime.now().isoformat(),
            "elapsed_total_s": round(elapsed_total, 2),
            "status": audit_status
        },
        "system_telemetry": {
            "os": platform.system(),
            "os_release": platform.release(),
            "machine": platform.machine(),
            "processor": platform.processor(),
            "cpu_count_logical": psutil.cpu_count(logical=True),
            "cpu_count_physical": psutil.cpu_count(logical=False),
            "memory_total_gb": round(psutil.virtual_memory().total / (1024**3), 2),
            "python_version": sys.version.split()[0],
            "cpu_architecture": platform.machine(),
        },
        "dependency_audit": {
            "scripts": {r["name"]: _get_file_hash(STEPS_DIR / r["name"]) for r in results},
            "utils": {
                path.name: _get_file_hash(path)
                for path in sorted((PROJECT_ROOT / "scripts" / "utils").glob("*.py"))
            },
            "config": {
                "config.json": _get_file_hash(PROJECT_ROOT / "config.json"),
            },
        },
        "step_telemetry": results
    }
    
    audit_path = PROJECT_ROOT / "results" / "audits" / f"RESEARCH_AUDIT_{datetime.datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
    audit_path.parent.mkdir(parents=True, exist_ok=True)
    
    with open(audit_path, "w") as f:
        json.dump(audit_data, f, indent=4)
        
    return audit_path


def run_step(script_name: str, step_idx: int, total: int, label: str = "STEP") -> dict:
    """
    Run a single pipeline step as a subprocess.

    Returns a dict with keys: name, status ("PASS" | "WARNING" | "FAIL"), elapsed_s, returncode.
    Never raises — failures are captured in the return dict.
    """
    script_path = STEPS_DIR / script_name
    output_path = _expected_output_path(script_name)
    pre_run_mtime = output_path.stat().st_mtime if output_path.exists() else None
    phase = f"[{step_idx:>3}/{total}]"

    print(f"\n{'─'*70}")
    print(f" {phase}  {label}: {script_name}")
    print(f"{'─'*70}\n")

    t0 = time.perf_counter()
    env = os.environ.copy()
    env["PYTHONPATH"] = str(PROJECT_ROOT) + os.pathsep + env.get("PYTHONPATH", "")
    result = subprocess.run(
        [sys.executable, str(script_path)],
        cwd=str(PROJECT_ROOT),
        capture_output=False,
        env=env,
    )
    elapsed = time.perf_counter() - t0

    output_ok = False
    output_check = "not checked"
    if result.returncode == 0:
        output_ok, output_check = _validate_step_output(output_path, pre_run_mtime)

    output_warning = output_check.endswith("status=WARNING")
    if result.returncode == 0 and output_ok:
        status = "WARNING" if output_warning else "PASS"
    else:
        status = "FAIL"
    if status == "PASS":
        print(f"\n✓  PASS  {script_name}  ({_fmt_elapsed(elapsed)})")
    elif status == "WARNING":
        print(f"\n⚠  WARNING  {script_name}  ({_fmt_elapsed(elapsed)})")
        print(f"   output check: {output_check}")
    else:
        print(f"\n✗  FAIL  {script_name}  rc={result.returncode}  ({_fmt_elapsed(elapsed)})")
        if result.returncode == 0:
            print(f"   output check: {output_check}")

    return {
        "name":        script_name,
        "status":      status,
        "elapsed_s":   round(elapsed, 2),
        "returncode":  result.returncode,
        "output":      str(output_path.relative_to(PROJECT_ROOT)),
        "output_check": output_check,
    }


def run_pipeline(
    pipeline_name: str,
    steps: list,
    description: str = "",
    stop_on_failure: bool = False,
) -> list:
    """
    Run a list of step scripts and print a summary table.

    Parameters
    ----------
    pipeline_name  : Human-readable name for the header banner.
    steps          : List of script filenames (e.g. "step_001_uncover_load.py").
    description    : Optional one-line description printed in the banner.
    stop_on_failure: If True, stop at first failure (default: continue all).

    Returns
    -------
    List of result dicts (one per step).
    """
    wall_start = time.perf_counter()
    n = len(steps)

    print("╔" + "═" * 68 + "╗")
    print(f"║  TEP-LLR: {pipeline_name:<58}  ║")
    if description:
        print(f"║  {description:<66}  ║")
    print(f"║  Steps: {n:<60}  ║")
    print(f"║  Started: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'):<58}  ║")
    print("╚" + "═" * 68 + "╝\n")

    results = []
    for i, step in enumerate(steps, start=1):
        r = run_step(step, i, n, label=pipeline_name)
        results.append(r)
        if stop_on_failure and r["status"] == "FAIL":
            print(f"\n  ⚠  Pipeline stopped at {step} (stop_on_failure=True)")
            break

    total_elapsed = time.perf_counter() - wall_start
    n_pass = sum(1 for r in results if r["status"] == "PASS")
    n_warn = sum(1 for r in results if r["status"] == "WARNING")
    n_fail = sum(1 for r in results if r["status"] == "FAIL")

    # Generate Audit Report
    audit_path = generate_audit_report(pipeline_name, results, total_elapsed)

    print()
    print("╔" + "═" * 68 + "╗")
    print(f"║  {pipeline_name} — COMPLETE" + " " * (42 - len(pipeline_name)) + "  ║")
    print("╠" + "═" * 68 + "╣")
    for r in results:
        icon = "✓" if r["status"] == "PASS" else "⚠" if r["status"] == "WARNING" else "✗"
        t = _fmt_elapsed(r["elapsed_s"])
        name = r["name"][:52]
        stat = r["status"]
        print(f"║  {icon} {name:<52}  {stat:<8}  {t:>5}  ║")
    print("╠" + "═" * 68 + "╣")
    pct = 100 * n_pass // len(results) if results else 0
    print(f"║  PASS {n_pass}/{len(results)} ({pct}%)   WARNING: {n_warn}   FAIL: {n_fail}   "
          f"Total: {_fmt_elapsed(total_elapsed):<37}  ║")
    print(f"║  AUDIT: {str(audit_path.relative_to(PROJECT_ROOT)):<58}  ║")
    print("╚" + "═" * 68 + "╝")

    return results
