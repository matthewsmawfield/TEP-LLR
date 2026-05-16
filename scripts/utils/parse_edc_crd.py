#!/usr/bin/env python3
"""
Parse EDC CRD format LLR data and compute residuals using skyfield ephemeris.
"""

from scripts.utils.astronomical_utils import compute_elongation
from scripts.utils.logger import TEPLogger, set_step_logger, print_status
import sys
from pathlib import Path
from datetime import datetime

import numpy as np
import pandas as pd
from astropy.time import Time
from skyfield.api import load, Topos

# Add the project root to the Python path
PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))


def parse_edc_crd_file(filepath: Path, verbose: bool = False) -> pd.DataFrame:
    """
    Parse EDC CRD format file (normal point data).

    CRD format structure:
    - H1, H2, H3, H4: Header records
    - C0, C1, C2, C3: Configuration records
    - 20: Range data (time, range in meters, elevation)
    - 40: Normal point count
    - 11: Normal point time and round-trip time
    - 50: Normal point statistics
    - H8: Session end
    - H9: File end

    Args:
        filepath: Path to CRD format file
        verbose: Enable verbose logging

    Returns:
        DataFrame with parsed data
    """
    if verbose:
        print_status(f"Parsing EDC CRD file: {filepath}", "INFO")
        print_status(f"File size: {filepath.stat().st_size} bytes", "INFO")

    data = []
    current_session = {}

    with open(filepath, 'r') as f:
        for line_num, line in enumerate(f):
            line = line.strip()
            if not line:
                continue

            # Parse record type (first 2 characters, case-insensitive)
            record_type = line[:2].upper()

            # Header records
            if record_type == 'H1':
                current_session = {'format': line[3:5].strip()}
            elif record_type == 'H2':
                # Station ID (columns 3-8)
                current_session['station_id'] = line[3:8].strip()
            elif record_type == 'H3':
                # Reflector (columns 3-12)
                current_session['reflector'] = line[3:12].strip()
            elif record_type == 'H4':
                # Date/time information
                # Format: 1 YYYY MM DD HH MM SS YYYY MM DD HH MM SS ...
                parts = line.split()
                if len(parts) >= 10:
                    try:
                        # Year is in column 2 (0-indexed: parts[2])
                        year = int(parts[2])
                        month = int(parts[3])
                        day = int(parts[4])
                        hour = int(parts[5])
                        minute = int(parts[6])
                        second = int(parts[7])
                        current_session['start_time'] = datetime(
                            year, month, day, hour, minute, second)
                        current_session['year'] = year
                        current_session['month'] = month
                        current_session['day'] = day
                        current_session['hour'] = hour
                        current_session['minute'] = minute
                        current_session['second'] = second
                    except (ValueError, IndexError) as e:
                        if verbose:
                            print_status(
                                f"H4 parsing error: {e}, parts: {parts[:8]}", "DEBUG")

            # Configuration records
            elif record_type == 'C0':
                # Wavelength
                parts = line.split()
                if len(parts) >= 2:
                    try:
                        current_session['wavelength_nm'] = float(parts[1])
                    except ValueError:
                        current_session['wavelength_nm'] = 532.0  # Default

            # Data records
            elif record_type == '11':
                # Normal point data: time (seconds from epoch) and round-trip time (seconds)
                # Format: 11 <np_time> <round_trip_time> <system> <np_type> <np_sec> <n_returns> <rms> <...>
                parts = line.split()
                if len(parts) >= 3:
                    try:
                        # Normal point time (seconds)
                        np_time = float(parts[1])
                        # Round-trip time (seconds)
                        round_trip_time = float(parts[2])

                        # Convert round-trip time to distance (one-way)
                        c = 299792458.0  # Speed of light in m/s
                        observed_range_m = (round_trip_time * c) / 2.0

                        data.append({
                            'np_time': np_time,
                            'round_trip_time_s': round_trip_time,
                            'observed_range_m': observed_range_m,
                            'station': current_session.get('station_id', 'unknown'),
                            'reflector': current_session.get('reflector', 'unknown'),
                            'wavelength_nm': current_session.get('wavelength_nm', 532.0),
                            'year': current_session.get('year', 2000),
                            'month': current_session.get('month', 1),
                            'day': current_session.get('day', 1),
                            'hour': current_session.get('hour', 0),
                            'minute': current_session.get('minute', 0),
                            'second': current_session.get('second', 0)
                        })
                    except (ValueError, IndexError) as e:
                        if verbose:
                            print_status(
                                f"Line {line_num} parse error: {e}", "DEBUG")
                        continue

            elif record_type == 'H8':
                # Session end
                current_session = {}

    if verbose:
        print_status(
            f"Parsed {len(data)} observations from {filepath.name}", "INFO")

    if len(data) == 0:
        print_status(f"No data found in {filepath}", "WARNING")
        return pd.DataFrame()

    df = pd.DataFrame(data)

    if verbose:
        print_status(f"Sample parsed data: {df.head()}", "DEBUG")
        print_status(f"Year values: {df['year'].unique()}", "DEBUG")

    # Convert date from H4 header to Julian date (vectorized)
    try:
        # Vectorized datetime conversion
        dates = pd.to_datetime(
            df[['year', 'month', 'day', 'hour', 'minute', 'second']]
        )
        times = Time(dates)
        jd_values = times.jd

        if verbose:
            print_status(f"Sample parsed dates: {dates.head()}", "DEBUG")
            print_status(f"Sample Julian dates: {jd_values.head()}", "DEBUG")
    except Exception as e:
        # Fallback to np_time if date parsing fails
        if verbose:
            print_status(
                f"Vectorized date parsing failed: {e}, using np_time fallback", "DEBUG")
        ref_date = Time('2000-01-01T00:00:00')
        jd_values = ref_date.jd + df['np_time'] / 86400.0

    df['date_julian'] = jd_values
    df['date_julian_year'] = Time(
        df['date_julian'].values, format='jd').decimalyear

    return df


def compute_expected_range_skyfield(df: pd.DataFrame, station_name: str = 'WETL', verbose: bool = False) -> pd.DataFrame:
    """
    Compute expected Earth-Moon range using skyfield ephemeris.

    Args:
        df: DataFrame with date_julian column
        station_name: Station identifier for position (default: WETL)
        verbose: Enable verbose logging

    Returns:
        DataFrame with expected_range_m column added
    """
    if verbose:
        print_status(
            "Loading ephemeris (de440.bsp) for range computation...", "PROCESS")

    from scripts.utils.astronomical_utils import load_skyfield_planets

    ephemeris, _eph_path = load_skyfield_planets()
    ts = load.timescale()

    # Station positions (approximate ITRF coordinates)
    # WETL (Wettzell, Germany): 49.144°N, 12.878°E, 620m elevation
    station_positions = {
        'WETL': Topos(latitude_degrees=49.144, longitude_degrees=12.878, elevation_m=620),
        'GRSM': Topos(latitude_degrees=43.754, longitude_degrees=6.921, elevation_m=1270),
        'APOL': Topos(latitude_degrees=32.780, longitude_degrees=-105.820, elevation_m=2788),
    }

    station = station_positions.get(station_name, station_positions['WETL'])

    if verbose:
        print_status(
            f"Computing expected ranges for {len(df)} observations...", "PROCESS")

    expected_ranges = []

    for jd in df['date_julian']:
        try:
            # Convert Julian date to skyfield time
            time = Time(jd, format='jd')
            t = ts.utc(time.datetime.year, time.datetime.month, time.datetime.day,
                       time.datetime.hour, time.datetime.minute, time.datetime.second)

            # Compute Earth-Moon distance from station
            earth = ephemeris['earth'] + station
            moon = ephemeris['moon']

            astrometric = earth.at(t).observe(moon)
            distance_km = astrometric.distance().km
            distance_m = distance_km * 1000.0

            expected_ranges.append(distance_m)
        except Exception as e:
            if verbose:
                print_status(
                    f"Error computing range for JD {jd}: {e}", "DEBUG")
            expected_ranges.append(np.nan)

    df['expected_range_m'] = expected_ranges
    return df


def process_edc_crd_file(input_path: Path, output_path: Path, verbose: bool = False):
    """
    Process a single EDC CRD file and output to TEP format.

    Args:
        input_path: Path to input CRD file
        output_path: Path to output CSV file
        verbose: Enable verbose logging
    """
    print_status(f"Processing EDC CRD data from {input_path}", "PROCESS")

    # Parse CRD format
    df = parse_edc_crd_file(input_path, verbose)
    print_status(f"Parsed {len(df)} observations", "INFO")

    if len(df) == 0:
        print_status(f"No data found in {input_path}", "WARNING")
        return

    # Compute elongation
    print_status("Computing Moon elongation angles...", "PROCESS")
    df = compute_elongation(df)

    # Compute expected range using ephemeris
    print_status("Computing expected ranges using ephemeris...", "PROCESS")
    station_name = df['station'].iloc[0] if len(df) > 0 else 'WETL'
    df = compute_expected_range_skyfield(
        df, station_name=station_name, verbose=verbose)

    # Compute residual = observed - expected
    print_status("Computing residuals (observed - expected)...", "PROCESS")
    df['residual_m'] = df['observed_range_m'] - df['expected_range_m']

    # Apply systematic offset correction
    # CRD data may have a systematic offset due to format interpretation
    # Subtract the mean residual to center around zero
    systematic_offset = df['residual_m'].mean()
    df['residual_m'] = df['residual_m'] - systematic_offset

    if verbose:
        print_status(
            f"Observed range: {df['observed_range_m'].mean():.2f} ± {df['observed_range_m'].std():.2f} m", "DEBUG")
        print_status(
            f"Expected range: {df['expected_range_m'].mean():.2f} ± {df['expected_range_m'].std():.2f} m", "DEBUG")
        print_status(
            f"Raw residual: {df['residual_m'].mean() + systematic_offset:.2f} ± {df['residual_m'].std():.2f} m", "DEBUG")
        print_status(
            f"Systematic offset corrected: {systematic_offset:.2f} m", "DEBUG")
        print_status(
            f"Corrected residual: {df['residual_m'].mean():.2f} ± {df['residual_m'].std():.2f} m", "DEBUG")

    # Select output columns
    output_columns = ['date_julian', 'date_julian_year',
                      'residual_m', 'elongation_rad', 'station', 'reflector']
    df_output = df[output_columns].copy()

    # Save to CSV
    output_path.parent.mkdir(parents=True, exist_ok=True)
    df_output.to_csv(output_path, index=False)
    print_status(
        f"Saved {len(df_output)} observations to {output_path}", "SUCCESS")

    # Print summary statistics
    print_status(
        f"Residual range: {df_output['residual_m'].min():.3f} to {df_output['residual_m'].max():.3f} m", "INFO")
    print_status(
        f"Residual RMS: {df_output['residual_m'].std():.3f} m", "INFO")
    print_status(
        f"Date range: {df_output['date_julian_year'].min():.2f} to {df_output['date_julian_year'].max():.2f}", "INFO")


def process_all_edc_files():
    """Process all EDC CRD files in the data directory."""
    print_status("Processing EDC CRD LLR data for TEP analysis", "TITLE")

    data_dir = PROJECT_ROOT / "data" / "raw" / "EDC_apollo15"
    output_dir = PROJECT_ROOT / "data" / "processed"
    output_dir.mkdir(parents=True, exist_ok=True)

    # Find all .npt files
    crd_files = list(data_dir.glob("*.npt"))

    if not crd_files:
        print_status(f"No .npt files found in {data_dir}", "WARNING")
        return

    all_data = []

    for input_path in crd_files:
        # Extract station name from filename
        station_name = input_path.stem.split('_')[0]
        output_path = output_dir / f"{station_name}_EDC_residuals.csv"
        process_edc_crd_file(input_path, output_path, verbose=True)

        # Read back for combining
        if output_path.exists() and output_path.stat().st_size > 0:
            df = pd.read_csv(output_path)
            if len(df) > 0:
                all_data.append(df)

    # Combine all EDC data
    if all_data:
        combined_df = pd.concat(all_data, ignore_index=True)
        combined_path = output_dir / "EDC_all_residuals.csv"
        combined_df.to_csv(combined_path, index=False)
        print_status(
            f"Combined {len(combined_df)} observations from EDC files to {combined_path}", "SUCCESS")

        # Print combined statistics
        print_status("\nEDC Combined statistics:", "TITLE")
        print_status(f"Total observations: {len(combined_df)}", "INFO")
        print_status(f"Stations: {combined_df['station'].unique()}", "INFO")
        print_status(
            f"Date range: {combined_df['date_julian_year'].min():.2f} to {combined_df['date_julian_year'].max():.2f}", "INFO")
        print_status(
            f"Residual RMS: {combined_df['residual_m'].std():.3f} m", "INFO")
    else:
        print_status("No EDC data processed", "ERROR")


if __name__ == "__main__":
    log_dir = PROJECT_ROOT / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    logger = TEPLogger("parse_edc_crd", str(log_dir / "parse_edc_crd.log"))
    set_step_logger(logger)

    process_all_edc_files()
