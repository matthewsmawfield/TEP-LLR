import logging
import sys
from pathlib import Path

# Module-level logger used by print_status when set via set_step_logger
_active_logger = None


class TEPLogger:
    def __init__(self, name, log_file_path=None):
        self.logger = logging.getLogger(name)
        self.logger.setLevel(logging.INFO)

        # Remove existing handlers to avoid duplicates
        if self.logger.hasHandlers():
            self.logger.handlers.clear()

        formatter = logging.Formatter('[%(asctime)s] %(message)s', datefmt='%H:%M:%S')

        # Console handler — plain, no timestamp (timestamps go to file only)
        ch = logging.StreamHandler(sys.stdout)
        ch.setFormatter(logging.Formatter('%(message)s'))
        self.logger.addHandler(ch)

        # File handler — full timestamp format
        if log_file_path:
            log_path = Path(log_file_path)
            log_path.parent.mkdir(parents=True, exist_ok=True)
            fh = logging.FileHandler(str(log_path), mode='w')
            fh.setFormatter(formatter)
            self.logger.addHandler(fh)

    def info(self, msg):
        self.logger.info(msg)

    def warning(self, msg):
        self.logger.warning(msg)

    def error(self, msg):
        self.logger.error(msg)

    def debug(self, msg):
        self.logger.debug(msg)

    def exception(self, msg):
        self.logger.exception(msg)

    def save_step_results(self, results: dict, project_root, step_name: str):
        import json
        import os
        from pathlib import Path
        root = Path(project_root)
        out_path = root / "results" / "outputs" / f"{step_name}.json"
        out_path.parent.mkdir(parents=True, exist_ok=True)
        with open(out_path, 'w') as f:
            json.dump(results, f, indent=4)
        # Show relative path for cleaner logs
        rel_path = out_path.relative_to(root) if out_path.is_relative_to(root) else out_path
        self.info(f"Results saved to {rel_path}")


def set_step_logger(logger: TEPLogger):
    """Register a TEPLogger so that print_status routes through it."""
    global _active_logger
    _active_logger = logger

_verbose_mode = False
def set_verbose_mode(verbose: bool):
    global _verbose_mode
    _verbose_mode = verbose

def get_verbose_mode() -> bool:
    global _verbose_mode
    return _verbose_mode


_LEVEL_PREFIXES = {
    "TITLE":   "═══",
    "PROCESS": ">>>",
    "SUCCESS": "✓  ",
    "ERROR":   "✗  ",
    "WARNING": "⚠  ",
    "INFO":    "   ",
    None:      "   ",
}


def print_status(msg: str, level=None):
    """
    Print a formatted status message.

    If a step logger has been registered via set_step_logger(), the message
    is routed through it (→ both console and log file). Otherwise falls back
    to a plain print.
    """
    prefix = _LEVEL_PREFIXES.get(level, "   ")
    formatted = f"{prefix} {msg}" if msg else ""

    if _active_logger is not None:
        if level == "ERROR":
            _active_logger.error(formatted)
        elif level == "WARNING":
            _active_logger.warning(formatted)
        else:
            _active_logger.info(formatted)
    else:
        print(formatted)
