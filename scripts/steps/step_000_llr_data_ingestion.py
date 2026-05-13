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
import hashlib
import json
from datetime import datetime


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        while chunk := handle.read(8192):
            digest.update(chunk)
    return digest.hexdigest()


def verify_manifest(raw_dir: Path, filenames: list[str], verbose: bool = False) -> tuple[list[dict], list[str]]:
    manifest_path = raw_dir / "data_manifest.json"
    if not manifest_path.exists():
        raise FileNotFoundError(
            f"Missing data integrity manifest: {manifest_path}. "
            "Cannot verify raw INPOP19a files against recorded checksums."
        )

    with open(manifest_path, "r", encoding="utf-8") as handle:
        manifest = json.load(handle)

    manifest_files = manifest.get("files", {})
    verified = []
    failures = []

    for filename in filenames:
        entry = manifest_files.get(filename)
        if entry is None:
            failures.append(f"{filename}: missing manifest entry")
            continue

        expected_sha = entry.get("sha256")
        if not expected_sha:
            failures.append(f"{filename}: manifest entry has no sha256")
            continue

        file_path = raw_dir / filename
        if not file_path.exists():
            failures.append(f"{filename}: file missing")
            continue

        actual_sha = _sha256_file(file_path)
        if actual_sha != expected_sha:
            failures.append(
                f"{filename}: sha256 mismatch (expected {expected_sha}, got {actual_sha})"
            )
            continue

        verified.append({
            "filename": filename,
            "sha256": actual_sha,
            "size_bytes": file_path.stat().st_size,
            "source": entry.get("source"),
            "source_url": entry.get("source_url"),
        })
        if verbose:
            print_status(
                f"Verified: {filename:<35} sha256={actual_sha[:12]}...",
                "SUCCESS",
            )

    return verified, failures


def check_data(verbose=False):
    raw_dir = PROJECT_ROOT / "data" / "raw"
    required = [
        "INPOP19a_APO_residuals.txt",
        "INPOP19a_Grasse_residuals.txt",
        "INPOP19a_Matera_residuals.txt",
        "INPOP19a_McDonald2_residuals.txt",
        "INPOP19a_Haleakala_residuals.txt",
        "DE430_2014-2018_residuals.dat",
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

    if missing:
        print_status(
            f"CRITICAL: Missing required raw data files ({len(missing)}): {missing}",
            "ERROR",
        )
        print_status(
            "Pipeline requires the full five-station INPOP19a archive and DE430 residuals.",
            "ERROR",
        )
        results = {
            "step_id": "step_000",
            "timestamp": datetime.now().isoformat(),
            "files_found": found,
            "files_missing": missing,
            "files_verified": [],
            "verification_failures": [],
            "status": "FAIL",
        }
        logger.save_step_results(results, PROJECT_ROOT,
                                 "step_000_llr_data_ingestion")
        sys.exit(1)

    verified_files, verification_failures = verify_manifest(
        PROJECT_ROOT / "data" / "raw",
        found,
        verbose=True,
    )

    results = {
        "step_id": "step_000",
        "timestamp": datetime.now().isoformat(),
        "files_found": found,
        "files_missing": missing,
        "files_verified": verified_files,
        "verification_failures": verification_failures,
        "status": "PASS" if not verification_failures else "FAIL"
    }

    print_status(
        f"All {len(found)} required LLR residual files present.", "SUCCESS")

    if verification_failures:
        print_status(
            f"Data integrity verification failed for {len(verification_failures)} file(s).",
            "ERROR",
        )
        for failure in verification_failures:
            print_status(failure, "ERROR")
        logger.save_step_results(results, PROJECT_ROOT,
                                 "step_000_llr_data_ingestion")
        sys.exit(1)

    logger.save_step_results(results, PROJECT_ROOT,
                             "step_000_llr_data_ingestion")
