#!/usr/bin/env python3
"""
Shared pipeline runner utility for TEP-LLR pipeline.

Provides run_step() and run_pipeline() for executing the 33-step
Lunar Laser Ranging analysis pipeline.
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


def generate_audit_report(pipeline_name: str, results: list, elapsed_total: float) -> Path:
    """
    Generate a formal Research-Grade Audit Report (JSON).
    Captures system state, dependency hashes, and pipeline telemetry.
    """
    audit_data = {
        "report_metadata": {
            "pipeline_name": pipeline_name,
            "timestamp": datetime.datetime.now().isoformat(),
            "elapsed_total_s": round(elapsed_total, 2),
            "status": "PASS" if all(r["status"] == "PASS" for r in results) else "FAIL"
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
            "is_m4_pro": "Apple M4 Pro" in platform.processor() or os.uname().machine == 'arm64' # Heuristic for Apple Silicon
        },
        "dependency_audit": {
            "scripts": {r["name"]: _get_file_hash(STEPS_DIR / r["name"]) for r in results},
            "utils": {
                "statistical_utils.py": _get_file_hash(PROJECT_ROOT / "scripts" / "utils" / "statistical_utils.py"),
                "llr_constants.py": _get_file_hash(PROJECT_ROOT / "scripts" / "utils" / "llr_constants.py"),
            }
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

    Returns a dict with keys: name, status ("PASS" | "FAIL"), elapsed_s, returncode.
    Never raises — failures are captured in the return dict.
    """
    script_path = STEPS_DIR / script_name
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

    status = "PASS" if result.returncode == 0 else "FAIL"
    if status == "PASS":
        print(f"\n✓  PASS  {script_name}  ({_fmt_elapsed(elapsed)})")
    else:
        print(f"\n✗  FAIL  {script_name}  rc={result.returncode}  ({_fmt_elapsed(elapsed)})")

    return {
        "name":        script_name,
        "status":      status,
        "elapsed_s":   round(elapsed, 2),
        "returncode":  result.returncode,
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
        if stop_on_failure and r["status"] != "PASS":
            print(f"\n  ⚠  Pipeline stopped at {step} (stop_on_failure=True)")
            break

    total_elapsed = time.perf_counter() - wall_start
    n_pass = sum(1 for r in results if r["status"] == "PASS")
    n_fail = sum(1 for r in results if r["status"] == "FAIL")

    # Generate Audit Report
    audit_path = generate_audit_report(pipeline_name, results, total_elapsed)

    print()
    print("╔" + "═" * 68 + "╗")
    print(f"║  {pipeline_name} — COMPLETE" + " " * (42 - len(pipeline_name)) + "  ║")
    print("╠" + "═" * 68 + "╣")
    for r in results:
        icon = "✓" if r["status"] == "PASS" else "✗"
        t = _fmt_elapsed(r["elapsed_s"])
        name = r["name"][:52]
        stat = r["status"]
        print(f"║  {icon} {name:<52}  {stat:<8}  {t:>5}  ║")
    print("╠" + "═" * 68 + "╣")
    pct = 100 * n_pass // len(results) if results else 0
    print(f"║  PASS {n_pass}/{len(results)} ({pct}%)   FAIL: {n_fail}   "
          f"Total: {_fmt_elapsed(total_elapsed):<37}  ║")
    print(f"║  AUDIT: {str(audit_path.relative_to(PROJECT_ROOT)):<58}  ║")
    print("╚" + "═" * 68 + "╝")

    return results
