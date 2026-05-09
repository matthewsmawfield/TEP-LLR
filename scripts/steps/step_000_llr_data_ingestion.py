#!/usr/bin/env python3
"""
Step 000: LLR Data Ingestion (Audit Mode)
Checks for presence of required raw data files.
"""

import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from scripts.utils.logger import TEPLogger, set_step_logger, print_status

import argparse
from datetime import datetime

# Add project root to path

def check_data(verbose=False):
    raw_dir = PROJECT_ROOT / "data" / "raw"
    required = [
        "INPOP19a_APO_residuals.txt",
        "INPOP19a_Grasse_residuals.txt",
        "INPOP19a_Matera_residuals.txt",
        "INPOP19a_McDonald2_residuals.txt",
        "INPOP19a_Haleakala_residuals.txt"
    ]

    found = []
    missing = []

    for f in required:
        file_path = raw_dir / f
        if file_path.exists():
            size_kb = file_path.stat().st_size / 1024
            found.append(f)
            if verbose:
                print_status(f"Found: {f:<35} ({size_kb:>8.1f} KB)", "INFO")
        else:
            missing.append(f)

    if verbose:
        print_status(
            f"Data audit complete. Station coverage: {len(found)}/{len(required)}", "INFO")

    return found, missing

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Step 000: LLR Data Ingestion Audit")
    args = parser.parse_args()

    log_dir = PROJECT_ROOT / "logs"
    logger = TEPLogger("step_000", str(
        log_dir / "step_000_llr_data_ingestion.log"))
    set_step_logger(logger)

    print_status("Starting data integrity check...", "TITLE")

    found, missing = check_data(verbose=True)

    results = {
        "step_id": "step_000",
        "timestamp": datetime.now().isoformat(),
        "files_found": found,
        "files_missing": missing,
        "status": "PASS" if not missing else "FAIL"
    }

    if missing:
        print_status(f"Missing required data files: {missing}", "WARNING")
        # Fail if core stations are missing to ensure high academic standard
        core_stations = ["INPOP19a_APO_residuals.txt", "INPOP19a_Grasse_residuals.txt"]
        for core in core_stations:
            if core in missing:
                print_status(f"CRITICAL: Core station data missing ({core}). Cannot proceed with robust analysis.", "ERROR")
                sys.exit(1)
        
        # If other stations are missing, log as a major weakness
        print_status(f"MAJOR WEAKNESS: Analysis proceeding without {len(missing)} stations.", "WARNING")
    else:
        print_status(
            f"All {len(found)} core LLR residual files found.", "SUCCESS")

    logger.save_step_results(results, PROJECT_ROOT,
                             "step_000_llr_data_ingestion")