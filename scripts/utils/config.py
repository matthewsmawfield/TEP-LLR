import json
from pathlib import Path
import multiprocessing as mp

# Default configuration values
DEFAULT_CONFIG = {
    "N_WORKERS": min(mp.cpu_count(), 12),
    "MCMC_LIGHT_WALKERS": 24,
    "MCMC_LIGHT_STEPS": 400,
    "MCMC_HEAVY_WALKERS": 32,
    "MCMC_HEAVY_STEPS": 1000,
    "MCMC_STANDARD_WALKERS": 32,
    "MCMC_STANDARD_STEPS": 2000,
    "MCMC_BURN_IN": 500,
    "BOOTSTRAP_ITERATIONS": 10000,
    "PERMUTATION_ITERATIONS": 10000,
    "THEIL_SEN_SAMPLES": 10000,
    "THEIL_SEN_BOOTSTRAPS": 1000,
    "LEVERAGE_BOOTSTRAPS": 200
}

# Cache the loaded config
_CONFIG = None

def get_config():
    global _CONFIG
    if _CONFIG is not None:
        return _CONFIG
        
    config = DEFAULT_CONFIG.copy()
    config_path = Path(__file__).resolve().parent.parent.parent / "config.json"
    
    if config_path.exists():
        try:
            with open(config_path, "r") as f:
                user_config = json.load(f)
                config.update(user_config)
        except Exception as e:
            raise RuntimeError(f"Failed to load config.json: {e}. Fix or remove the file to use defaults.")
            
    _CONFIG = config
    return _CONFIG
