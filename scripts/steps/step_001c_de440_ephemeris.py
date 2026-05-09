#!/usr/bin/env python3
"""
Step 001c: DE440 Ephemeris Validation and Preparation
Validates the presence and integrity of the DE440 BSP file.
"""

import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from scripts.utils.logger import TEPLogger, set_step_logger, print_status, set_verbose_mode

import argparse

# Add project root to path

def validate_de440(verbose=False):
    raw_dir = PROJECT_ROOT / "data" / "raw"
    de440_file = raw_dir / "de440.bsp"

    status = "MISSING"
    file_size = 0

    if de440_file.exists():
        file_size = de440_file.stat().st_size
        if file_size > 100 * 1024 * 1024:  # > 100MB
            status = "PASS"
            if verbose:
                print_status(
                    f"DE440 found: {de440_file.name} ({file_size/1e6:.1f} MB)", "INFO")
                print_status("DE440 file check:", "CALC")
                print_status(f"  Size: {file_size} bytes", "CALC")
                print_status(f"  Status: {status}", "CALC")
        else:
            status = "INVALID"
            print_status(
                f"DE440 file too small: {file_size/1e6:.1f} MB", "WARNING")
    else:
        print_status(f"DE440 file not found at {de440_file}", "WARNING")

    return {
        "status": status,
        "file_path": str(de440_file.relative_to(PROJECT_ROOT)) if de440_file.exists() else None,
        "file_size_mb": round(file_size / 1e6, 2),
        "coverage": "1550-2650" if status == "PASS" else "N/A"
    }

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Step 001c: DE440 Ephemeris Processing")
    args = parser.parse_args()

    log_dir = PROJECT_ROOT / "logs"
    logger = TEPLogger("step_001c", str(
        log_dir / "step_001c_de440_ephemeris.log"))
    set_step_logger(logger)
    set_verbose_mode(True)

    print_status("Checking DE440 Ephemeris...", "TITLE")
    summary = validate_de440(verbose=True)

    results = {
        "step_id": "step_001c",
        "validation": summary,
        "status": "PASS" if summary["status"] == "PASS" else "FAIL"
    }

    logger.save_step_results(results, PROJECT_ROOT,
                             "step_001c_de440_ephemeris")
    if summary["status"] == "PASS":
        print_status("DE440 Validation Complete.", "SUCCESS")