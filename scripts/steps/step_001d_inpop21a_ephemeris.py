#!/usr/bin/env python3
"""
Step 001d: INPOP21a Ephemeris Validation and Preparation
Validates the presence and integrity of the INPOP21a binary file.
"""

import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from scripts.utils.logger import TEPLogger, set_step_logger, print_status, set_verbose_mode

import argparse

# Add project root to path

def validate_inpop21a(verbose=False):
    raw_dir = PROJECT_ROOT / "data" / "raw"
    inpop_file = raw_dir / "inpop21a_TDB_m1000_p1000_littleendian.dat"

    status = "MISSING"
    file_size = 0

    if inpop_file.exists():
        file_size = inpop_file.stat().st_size
        if file_size > 150 * 1024 * 1024:  # > 150MB
            status = "PASS"
            if verbose:
                print_status(
                    f"INPOP21a found: {inpop_file.name} ({file_size/1e6:.1f} MB)", "INFO")
                print_status("INPOP21a file check:", "CALC")
                print_status(f"  Size: {file_size} bytes", "CALC")
                print_status(f"  Status: {status}", "CALC")
        else:
            status = "INVALID"
            print_status(
                f"INPOP21a file too small: {file_size/1e6:.1f} MB", "WARNING")
    else:
        print_status(f"INPOP21a file not found at {inpop_file}", "WARNING")

    return {
        "status": status,
        "file_path": str(inpop_file.relative_to(PROJECT_ROOT)) if inpop_file.exists() else None,
        "file_size_mb": round(file_size / 1e6, 2),
        "coverage": "1000-2000" if status == "PASS" else "N/A"
    }

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Step 001d: INPOP21a Ephemeris Processing")
    args = parser.parse_args()

    log_dir = PROJECT_ROOT / "logs"
    logger = TEPLogger("step_001d", str(
        log_dir / "step_001d_inpop21a_ephemeris.log"))
    set_step_logger(logger)
    set_verbose_mode(True)

    print_status("Checking INPOP21a Ephemeris...", "TITLE")
    summary = validate_inpop21a(verbose=True)

    results = {
        "step_id": "step_001d",
        "validation": summary,
        "status": "PASS" if summary["status"] == "PASS" else "FAIL"
    }

    logger.save_step_results(results, PROJECT_ROOT,
                             "step_001d_inpop21a_ephemeris")
    if summary["status"] == "PASS":
        print_status("INPOP21a Validation Complete.", "SUCCESS")