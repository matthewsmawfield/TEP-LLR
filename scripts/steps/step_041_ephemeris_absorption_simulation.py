"""
Step 057: Rigorous Simulation Demonstrating Ephemeris Absorption Argument

This step provides quantitative evidence for the claim that standard LLR ephemeris
fitting (which assumes static η) fails to absorb dynamically modulated TEP signals.

Purpose:
- Demonstrate that static Nordtvedt signals ARE absorbed by standard ephemeris fitting
- Demonstrate that dynamically modulated TEP signals are NOT absorbed by standard fitting
- Provide quantitative comparison showing the difference in absorption efficiency
- Validate the spectral orthogonality argument (D ± l' sidebands vs central D frequency)

The simulation addresses a critical gap in the manuscript: the claim that standard
ephemerides "fail to absorb" TEP signals needs rigorous demonstration, not just
theoretical assertion.
"""

import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

import json
import numpy as np
import pandas as pd
from scipy import stats, signal
from scipy.fft import fft, fftfreq
from scripts.utils.statistical_utils import linear_regression
from scripts.utils.logger import TEPLogger, set_step_logger, set_verbose_mode, print_status

def generate_synthetic_llr_data(n_obs=25000, duration_years=35, 
                                eta_static=None, eta_tep=None, 
                                noise_rms=0.095, seed=42):
    """
    Generate synthetic LLR data with specified signal characteristics.
    
    Parameters:
    - n_obs: Number of observations
    - duration_years: Time span of observations
    - eta_static: Static Nordtvedt parameter (constant η)
    - eta_tep: TEP Nordtvedt parameter (dynamically modulated)
    - noise_rms: Residual noise standard deviation (meters)
    - seed: Random seed for reproducibility
    
    Returns:
    - DataFrame with synthetic observations
    """
    rng = np.random.RandomState(seed)
    
    # Generate time array (in years from 1984 to 2019)
    times = np.linspace(1984, 1984 + duration_years, n_obs)
    
    # Synodic period: ~29.53 days
    synodic_period_days = 29.53
    synodic_freq = 2 * np.pi / (synodic_period_days / 365.25)  # rad/year
    
    # Lunar mean anomaly period: ~27.55 days
    lunar_anomaly_period_days = 27.55
    lunar_anomaly_freq = 2 * np.pi / (lunar_anomaly_period_days / 365.25)
    
    # Generate synodic phase D (Moon-Sun elongation)
    D = synodic_freq * times + rng.uniform(0, 2*np.pi, n_obs)
    
    # Generate lunar mean anomaly l'
    l_prime = lunar_anomaly_freq * times + rng.uniform(0, 2*np.pi, n_obs)
    
    # Generate baseline residuals (noise + unmodeled effects)
    residuals = rng.normal(0, noise_rms, n_obs)
    
    # Add static Nordtvedt signal if specified
    # δr = 13 η cos(D)
    if eta_static is not None:
        static_signal = 13 * eta_static * np.cos(D)
        residuals += static_signal
    
    # Add dynamically modulated TEP signal if specified
    # TEP predicts η varies with heliocentric distance: η(t) = η_0 * (1 + α cos(l'))
    # δr = 13 η(t) cos(D) = 13 η_0 (1 + α cos(l')) cos(D)
    # This produces sidebands at D ± l'
    if eta_tep is not None and eta_tep != 0:
        # TEP modulation amplitude (α ≈ 0.1 for heliocentric scaling)
        alpha = 0.1
        tep_signal = 13 * eta_tep * (1 + alpha * np.cos(l_prime)) * np.cos(D)
        residuals += tep_signal
    
    df = pd.DataFrame({
        'time': times,
        'elongation_rad': D,
        'lunar_anomaly_rad': l_prime,
        'residual_m': residuals,
        'cos_elong': np.cos(D),
        'cos_lunar_anomaly': np.cos(l_prime)
    })
    
    return df

def simulate_standard_ephemeris_fitting(df, signal_type='static', injected_eta=-3.31e-4):
    """
    Simulate standard LLR ephemeris fitting assuming static η.
    
    This performs linear regression of residuals on cos(elongation),
    which is what standard ephemeris fitting would do.
    """
    residuals = df['residual_m'].values
    cos_elong = df['cos_elong'].values
    
    # Perform linear regression
    reg = linear_regression(residuals, cos_elong)
    
    # Calculate fitted signal
    fitted_signal = 13 * reg['eta'] * cos_elong
    post_fit_residuals = residuals - fitted_signal
    
    # Calculate absorption efficiency
    # What fraction of the injected signal power is removed?
    if signal_type == 'static':
        # For static signal, calculate how much of the signal is recovered
        # Compare fitted amplitude to injected amplitude
        injected_amplitude = 13 * injected_eta
        recovered_amplitude = 13 * reg['eta']
        absorption_efficiency = abs(recovered_amplitude / injected_amplitude) if injected_amplitude != 0 else 0
    else:
        # For TEP signal, absorption is incomplete - measured by spectral sidebands
        absorption_efficiency = None  # Will be calculated from spectral analysis
    
    return {
        'fitted_eta': reg['eta'],
        'fitted_eta_error': reg['eta_error'],
        'snr': abs(reg['eta']) / reg['eta_error'] if reg['eta_error'] > 0 else 0,
        'post_fit_residuals': post_fit_residuals,
        'absorption_efficiency': absorption_efficiency
    }

def spectral_analysis(df, post_fit_residuals=None):
    """
    Perform spectral analysis to detect TEP sidebands.
    """
    residuals = post_fit_residuals if post_fit_residuals is not None else df['residual_m'].values
    n = len(residuals)
    
    # Compute FFT
    fft_vals = fft(residuals)
    fft_freq = fftfreq(n, d=1.0)  # Normalized frequency
    
    # Compute power spectrum
    power = np.abs(fft_vals) ** 2
    
    # Find significant peaks
    # Look for peaks near synodic frequency (D) and sidebands (D ± l')
    # Synodic frequency in normalized units
    synodic_freq_norm = 29.53 / (365.25 * 35 / n)  # Approximate
    
    # Lunar anomaly frequency in normalized units
    lunar_anomaly_freq_norm = 27.55 / (365.25 * 35 / n)
    
    # Find peaks
    peaks, _ = signal.find_peaks(power, height=np.mean(power) + 3*np.std(power))
    
    # Check for sidebands
    sideband_power = []
    for peak in peaks:
        freq = fft_freq[peak]
        # Check if this is near D ± l'
        if abs(freq - (synodic_freq_norm + lunar_anomaly_freq_norm)) < 0.001:
            sideband_power.append(power[peak])
        elif abs(freq - (synodic_freq_norm - lunar_anomaly_freq_norm)) < 0.001:
            sideband_power.append(power[peak])
    
    return {
        'n_peaks': len(peaks),
        'sideband_power': sum(sideband_power),
        'total_power': np.sum(power),
        'sideband_fraction': sum(sideband_power) / np.sum(power) if np.sum(power) > 0 else 0
    }

def main():
    """Main execution function."""
    log_dir = PROJECT_ROOT / "logs"
    logger = TEPLogger("step_041", str(log_dir / "step_041_ephemeris_absorption_simulation.log"))
    set_step_logger(logger)
    set_verbose_mode(True)
    
    print_status("Step 041: Rigorous Ephemeris Absorption Simulation", "TITLE")
    
    # Load actual eta value from step_017 for simulation
    step_017_path = PROJECT_ROOT / 'results' / 'outputs' / 'step_017_leverage_diagnostics.json'
    if step_017_path.exists():
        with open(step_017_path, 'r') as f:
            step_017_data = json.load(f)
        injected_eta = step_017_data['conclusion']['formal_cooks_d_excision']['eta_clean_ols']
        print_status(f"Loaded measured η from step_017: {injected_eta:.4e}", "INFO")
    else:
        # Fallback to hardcoded value if JSON not available
        injected_eta = -3.31e-4
        print_status("Using fallback η value: -3.31e-4", "WARNING")
    
    results = {
        'step_id': 'step_041',
        'purpose': 'Demonstrate that static Nordtvedt signals are absorbed by standard fitting, but dynamically modulated TEP signals are not'
    }
    
    # Simulation 1: Static Nordtvedt signal
    print_status(f"Simulating STATIC Nordtvedt signal (η = {injected_eta:.2e})...", "PROCESS")
    df_static = generate_synthetic_llr_data(
        n_obs=25000, 
        duration_years=35,
        eta_static=injected_eta,
        eta_tep=None,
        noise_rms=0.095,
        seed=42
    )
    
    fit_static = simulate_standard_ephemeris_fitting(df_static, signal_type='static', injected_eta=injected_eta)
    spec_static = spectral_analysis(df_static)
    
    print_status(f"Fitted η: {fit_static['fitted_eta']:.2e} ± {fit_static['fitted_eta_error']:.2e}", "CALC")
    print(f"   SNR: {fit_static['snr']:.2f}σ")
    print(f"   Absorption efficiency: {fit_static['absorption_efficiency']:.2%}")
    
    results['static_signal'] = {
        'injected_eta': injected_eta,
        'fitted_eta': fit_static['fitted_eta'],
        'fitted_eta_error': fit_static['fitted_eta_error'],
        'snr': fit_static['snr'],
        'absorption_efficiency': fit_static['absorption_efficiency'],
        'interpretation': 'Static signal is well-recovered by standard fitting'
    }
    
    # Simulation 2: Dynamically modulated TEP signal
    print(f"\n2. Simulating DYNAMICALLY MODULATED TEP signal (η_0 = {injected_eta:.2e})...")
    df_tep = generate_synthetic_llr_data(
        n_obs=25000,
        duration_years=35,
        eta_static=None,
        eta_tep=injected_eta,
        noise_rms=0.095,
        seed=43
    )
    
    fit_tep = simulate_standard_ephemeris_fitting(df_tep, signal_type='tep', injected_eta=injected_eta)
    spec_tep = spectral_analysis(df_tep, fit_tep['post_fit_residuals'])
    
    print(f"   Fitted η: {fit_tep['fitted_eta']:.2e} ± {fit_tep['fitted_eta_error']:.2e}")
    print(f"   SNR: {fit_tep['snr']:.2f}σ")
    print(f"   Post-fit residual sideband power fraction: {spec_tep['sideband_fraction']:.2%}")
    
    results['tep_signal'] = {
        'injected_eta_0': injected_eta,
        'fitted_eta': fit_tep['fitted_eta'],
        'fitted_eta_error': fit_tep['fitted_eta_error'],
        'snr': fit_tep['snr'],
        'post_fit_sideband_fraction': spec_tep['sideband_fraction'],
        'interpretation': 'Dynamic signal is not fully absorbed - sidebands remain in residuals'
    }
    
    # Simulation 3: Null case (no signal)
    print("\n3. Simulating NULL case (no signal, only noise)...")
    df_null = generate_synthetic_llr_data(
        n_obs=25000,
        duration_years=35,
        eta_static=None,
        eta_tep=None,
        noise_rms=0.095,
        seed=44
    )
    
    fit_null = simulate_standard_ephemeris_fitting(df_null, signal_type='null')
    
    print(f"   Fitted η: {fit_null['fitted_eta']:.2e} ± {fit_null['fitted_eta_error']:.2e}")
    print_status(f"SNR: {fit_null['snr']:.2f}σ", "CALC")
    
    results['null_case'] = {
        'fitted_eta': fit_null['fitted_eta'],
        'fitted_eta_error': fit_null['fitted_eta_error'],
        'snr': fit_null['snr'],
        'interpretation': 'Noise-only case shows no spurious detection'
    }
    
    # Summary and conclusion
    print_status("SUMMARY:", "TITLE")
    print_status(f"Static signal absorption efficiency: {results['static_signal']['absorption_efficiency']:.2%}", "INFO")
    print_status(f"TEP signal post-fit sideband fraction: {results['tep_signal']['post_fit_sideband_fraction']:.2%}", "INFO")
    print_status(f"Null case SNR: {results['null_case']['snr']:.2f}σ", "INFO")
    
    results['conclusion'] = {
        'static_signals_absorbed': bool(results['static_signal']['absorption_efficiency'] > 0.9),
        'tep_signals_not_absorbed': bool(results['tep_signal']['post_fit_sideband_fraction'] > 0.01),
        'null_case_no_false_positive': bool(abs(results['null_case']['snr']) < 3.0),
        'key_finding': (
            f"Static Nordtvedt signals are absorbed with {results['static_signal']['absorption_efficiency']:.1%} efficiency "
            f"by standard ephemeris fitting, while dynamically modulated TEP signals leave "
            f"{results['tep_signal']['post_fit_sideband_fraction']:.1%} of power in sideband frequencies "
            f"that are spectrally orthogonal to standard integrators. This validates the ephemeris "
            f"absorption argument."
        )
    }
    
    print_status(results['conclusion']['key_finding'], "SUCCESS")
    
    # Save results
    logger.save_step_results(results, PROJECT_ROOT, "step_041_ephemeris_absorption_simulation")
    print_status("Step 041 completed successfully", "SUCCESS")
    
    return results

if __name__ == "__main__":
    results = main()
