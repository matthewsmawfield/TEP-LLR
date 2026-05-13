import json
from pathlib import Path
import multiprocessing as mp

# Default configuration values
DEFAULT_CONFIG = {
    "N_WORKERS": min(mp.cpu_count(), 12),
    "RANDOM_SEED": 42,
    "MCMC_LIGHT_WALKERS": 24,
    "MCMC_LIGHT_STEPS": 400,
    "MCMC_HEAVY_WALKERS": 32,
    "MCMC_HEAVY_STEPS": 1000,
    "MCMC_STANDARD_WALKERS": 32,
    "MCMC_STANDARD_STEPS": 3000,
    "MCMC_BURN_IN": 1000,
    "BOOTSTRAP_ITERATIONS": 10000,
    "PERMUTATION_ITERATIONS": 10000,
    "THEIL_SEN_SAMPLES": 10000,
    "THEIL_SEN_BOOTSTRAPS": 1000,
    "LEVERAGE_BOOTSTRAPS": 1000
}

# Cache the loaded config
_CONFIG = None

def get_config():
    global _CONFIG
    if _CONFIG is not None:
        return _CONFIG
        
    config_path = Path(__file__).resolve().parent.parent.parent / "config.json"

    if not config_path.exists():
        raise FileNotFoundError(
            f"Missing required pipeline configuration: {config_path}. "
            "Restore config.json rather than falling back to implicit defaults."
        )

    try:
        with open(config_path, "r") as f:
            config = json.load(f)
    except Exception as e:
        raise RuntimeError(f"Failed to load required config.json: {e}. Fix the file before running the pipeline.")

    missing_keys = sorted(set(DEFAULT_CONFIG) - set(config))
    if missing_keys:
        raise RuntimeError(f"config.json missing required keys: {', '.join(missing_keys)}")
            
    _CONFIG = config.copy()
    return _CONFIG.copy()
