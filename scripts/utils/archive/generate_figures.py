#!/usr/bin/env python3
"""
Generate figures for TEP-LLR manuscript.

This script creates publication-quality figures for the TEP-LLR analysis,
including residual plots, station comparisons, and phase-binned analyses.
All figures are saved to results/figures/ directory.

Author: TEP-LLR Analysis Pipeline
Date: 2024
"""

import sys
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib

# Add the project root to the Python path
PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from scripts.utils.llr_constants import ETA_SCALE_FACTOR

matplotlib.use('Agg')


def load_data():
    """
    Load processed LLR data from CSV file.

    Returns:
        pd.DataFrame: DataFrame with LLR residual data including
                     date_julian, date_julian_year, residual_m,
                     elongation_rad, station columns
    """
    data_file = PROJECT_ROOT / "data" / "processed" / \
        "INPOP19a_all_stations_residuals.csv"
    df = pd.read_csv(data_file)
    return df


def plot_residuals_vs_cos_elongation(df, output_dir):
    """
    Create scatter plot of residuals vs cos(elongation).

    Args:
        df: DataFrame with LLR residual data
        output_dir: Directory to save the figure

    The plot includes:
    - Scatter plot of residuals vs cos(elongation) with transparency
    - Linear fit line showing the TEP signal
    - Fit equation in the legend
    """
    fig, ax = plt.subplots(figsize=(10, 6))

    cos_elong = np.cos(df['elongation_rad'].values)
    residuals = df['residual_m'].values

    # Scatter plot with transparency
    ax.scatter(cos_elong, residuals, alpha=0.3, s=1, color='blue')

    # Add fit line
    A = np.sum(residuals * cos_elong) / np.sum(cos_elong**2)
    x_fit = np.linspace(-1, 1, 100)
    y_fit = A * x_fit
    ax.plot(x_fit, y_fit, 'r-', linewidth=2, label=f'Fit: δr = {A:.4f} cos(D)')

    ax.set_xlabel('cos(D)', fontsize=12)
    ax.set_ylabel('Residual (m)', fontsize=12)
    ax.set_title('LLR Residuals vs cos(Moon-Sun Elongation)', fontsize=14)
    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    output_file = output_dir / "residuals_vs_cos_elongation.png"
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved: {output_file}")


def plot_residuals_vs_elongation(df, output_dir):
    """
    Create scatter plot of residuals vs elongation angle.

    Args:
        df: DataFrame with LLR residual data
        output_dir: Directory to save the figure

    The plot includes:
    - Scatter plot of residuals vs elongation (0 to 2π radians)
    - Measured TEP signal curve (computed from data, not hardcoded)
    - Legend with the measured eta value
    """
    fig, ax = plt.subplots(figsize=(10, 6))

    elongation = df['elongation_rad'].values
    residuals = df['residual_m'].values

    # Scatter plot with transparency
    ax.scatter(elongation, residuals, alpha=0.3, s=1, color='green')

    # Add measured TEP signal (computed from data)
    A = np.sum(residuals * np.cos(elongation)) / np.sum(np.cos(elongation)**2)
    eta = A / ETA_SCALE_FACTOR
    elongation_fit = np.linspace(0, 2*np.pi, 100)
    predicted = 13 * eta * np.cos(elongation_fit)
    ax.plot(elongation_fit, predicted, 'r-', linewidth=2,
            label=f'Measured TEP: η = {eta:.2e}')

    ax.set_xlabel('Elongation (rad)', fontsize=12)
    ax.set_ylabel('Residual (m)', fontsize=12)
    ax.set_title('LLR Residuals vs Moon-Sun Elongation', fontsize=14)
    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    output_file = output_dir / "residuals_vs_elongation.png"
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved: {output_file}")


def plot_station_comparison(df, output_dir):
    """
    Create bar plot comparing eta values across stations.

    Args:
        df: DataFrame with LLR residual data
        output_dir: Directory to save the figure

    The plot includes:
    - Bar chart of eta values for each station with error bars
    - Horizontal line at eta=0 (GR prediction)
    - Only includes stations with >=50 observations
    """
    # Station-by-station analysis
    stations = df['station'].unique()
    station_results = {}

    for station in stations:
        station_data = df[df['station'] == station]
        residuals = station_data['residual_m'].values
        elongation = station_data['elongation_rad'].values
        cos_elong = np.cos(elongation)

        if len(residuals) < 50:
            continue

        A = np.sum(residuals * cos_elong) / np.sum(cos_elong**2)
        eta = A / ETA_SCALE_FACTOR

        # Calculate uncertainty
        A_err = np.sqrt(np.sum(residuals**2) /
                        (len(residuals) * np.sum(cos_elong**2)))
        eta_err = A_err / ETA_SCALE_FACTOR

        station_results[station] = {'eta': eta,
                                    'eta_err': eta_err, 'n': len(residuals)}

    # Plot
    fig, ax = plt.subplots(figsize=(10, 6))

    station_names = list(station_results.keys())
    etas = [station_results[s]['eta'] for s in station_names]
    eta_errs = [station_results[s]['eta_err'] for s in station_names]

    x_pos = np.arange(len(station_names))
    ax.bar(x_pos, etas, yerr=eta_errs, capsize=5, alpha=0.7, color='steelblue')
    ax.axhline(y=0, color='r', linestyle='--',
               linewidth=2, label='GR prediction (η=0)')

    ax.set_xlabel('Station', fontsize=12)
    ax.set_ylabel('Nordtvedt Parameter η', fontsize=12)
    ax.set_title(
        'Station-by-Station Nordtvedt Parameter Estimates', fontsize=14)
    ax.set_xticks(x_pos)
    ax.set_xticklabels(station_names, rotation=45)
    ax.legend()
    ax.grid(True, alpha=0.3, axis='y')

    plt.tight_layout()
    output_file = output_dir / "station_comparison.png"
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved: {output_file}")


def plot_phase_binned_residuals(df, output_dir):
    """
    Create phase-binned mean residuals plot.

    Args:
        df: DataFrame with LLR residual data
        output_dir: Directory to save the figure

    The plot includes:
    - Error bar plot of mean residuals in 12 phase bins
    - Measured TEP signal curve (computed from data, not hardcoded)
    - Bins with <10 observations are shown as NaN
    """
    n_bins = 12
    phase_bins = np.linspace(0, 2*np.pi, n_bins + 1)
    bin_centers = (phase_bins[:-1] + phase_bins[1:]) / 2

    bin_means = []
    bin_errors = []

    for i in range(n_bins):
        mask = (df['elongation_rad'] >= phase_bins[i]) & (
            df['elongation_rad'] < phase_bins[i+1])
        if np.sum(mask) < 10:
            bin_means.append(np.nan)
            bin_errors.append(np.nan)
        else:
            bin_residuals = df.loc[mask, 'residual_m']
            bin_means.append(np.mean(bin_residuals))
            bin_errors.append(
                np.std(bin_residuals, ddof=1) / np.sqrt(np.sum(mask)))

    fig, ax = plt.subplots(figsize=(10, 6))

    ax.errorbar(bin_centers, bin_means, yerr=bin_errors,
                fmt='o', capsize=5, markersize=8, color='blue')

    # Add measured TEP signal (computed from data)
    A = np.sum(df['residual_m'].values * np.cos(df['elongation_rad'].values)
               ) / np.sum(np.cos(df['elongation_rad'].values)**2)
    eta = A / ETA_SCALE_FACTOR
    predicted = 13 * eta * np.cos(bin_centers)
    ax.plot(bin_centers, predicted, 'r-', linewidth=2,
            label=f'Measured TEP: η = {eta:.2e}')

    ax.set_xlabel('Elongation (rad)', fontsize=12)
    ax.set_ylabel('Mean Residual (m)', fontsize=12)
    ax.set_title('Phase-Binned Mean LLR Residuals', fontsize=14)
    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    output_file = output_dir / "phase_binned_residuals.png"
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved: {output_file}")


def main():
    """
    Main function to generate all figures for TEP-LLR manuscript.

    This function:
    1. Creates the output directory (results/figures/)
    2. Loads processed LLR data
    3. Generates four publication-quality figures:
       - residuals_vs_cos_elongation.png
       - residuals_vs_elongation.png
       - station_comparison.png
       - phase_binned_residuals.png
    4. Prints status messages for each figure
    """

    # Create output directory
    output_dir = PROJECT_ROOT / "results" / "figures"
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load data
    df = load_data()
    print(f"Loaded {len(df)} observations")

    # Generate figures
    plot_residuals_vs_cos_elongation(df, output_dir)
    plot_residuals_vs_elongation(df, output_dir)
    plot_station_comparison(df, output_dir)
    plot_phase_binned_residuals(df, output_dir)

    print("\nAll figures generated successfully!")


if __name__ == "__main__":
    main()
