#!/usr/bin/env python3
"""
TEP-LLR Step 001: Data Preprocessing
"""

import sys
from pathlib import Path

# Add project root to path
PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from scripts.utils.parse_inpop_mini import process_all_stations
from scripts.utils.logger import TEPLogger, set_step_logger, print_status

if __name__ == "__main__":
    log_dir = PROJECT_ROOT / "logs"
    log_dir.mkdir(exist_ok=True)
    logger = TEPLogger("step_001", str(log_dir / "step_001_data_preprocessing.log"))
    set_step_logger(logger)

    print_status("Starting INPOP MINI format preprocessing...", "TITLE")
    summary = process_all_stations()

    if summary:
        summary["step_id"] = "step_001"
        summary["status"] = "PASS"
        logger.save_step_results(summary, PROJECT_ROOT, "step_001_data_preprocessing")
        print_status("Preprocessing complete.", "SUCCESS")
