#!/usr/bin/env python3
"""
Parse DE430 format LLR residuals and convert to TEP analysis format.

DE430 format from Paris Observatory (Geoazur):
Column 1: Modified Julian date (MJD) or similar (e.g., 23825.699490)
Column 2: Residual in meters (O-C, e.g., 0.04321 m)
Column 3: Timestamp/ID (long string)
Column 4: Additional metadata

Output format for TEP analysis:
- residual_m: Residual in meters
- elongation_rad: Moon elongation angle in radians
- date_julian: Julian date
"""

from scripts.utils.astronomical_utils import compute_elongation
from scripts.utils.logger import TEPLogger, set_step_logger, print_status
import sys
from pathlib import Path

import pandas as pd
from astropy.time import Time

# Add the project root to the Python path
PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))


def parse_de430_file(filepath: Path) -> pd.DataFrame:
    """
    Parse DE430 format file.

    Args:
        filepath: Path to DE430 format file

    Returns:
        DataFrame with parsed data
    """
    data = []

    with open(filepath, 'r') as f:
        for line_num, line in enumerate(f):
            parts = line.strip().split()
            if len(parts) < 3:
                continue

            try:
                # Column 1: Date (appears to be some other format, not used)
                float(parts[0])

                # Column 2: Residual in meters
                residual_m = float(parts[1])

                # Column 3: Timestamp/ID - this appears to contain the actual date
                # Format: 5120150326164716328590825462182877449301910037002793631
                # Decoding: 5 1 2015 03 26 → year=2015, month=03, day=26
                timestamp = parts[2]

                # Decode date from timestamp
                # The timestamp appears to be: [station_code?] [year] [month] [day] [time]...
                # Format: 5 1 2015 03 26 16 47 16...
                #         ^ ^ ^    ^  ^  ^  ^  ^
                #         | | |    |  |  |  |  |
                #         | | |    |  |  |  |  +-- hour?
                #         | | |    |  |  |  +----- minute?
                #         | | |    |  |  +-------- second?
                #         | | |    |  +----------- day
                #         | | |    +-------------- month
                #         | | +------------------- year
                #         | +--------------------- ??? (maybe 1 for first digit of year?)
                #         +----------------------- ??? (maybe station code?)

                # Try to extract year, month, day from timestamp
                # Looking at pattern: 5 1 2015 03 26
                # Position 1-2: "51" or "5" and "1"
                # Position 3-6: "2015"
                # Position 7-8: "03"
                # Position 9-10: "26"

                if len(timestamp) >= 10:
                    try:
                        year = int(timestamp[2:6])  # 2015
                        month = int(timestamp[6:8])  # 03
                        day = int(timestamp[8:10])   # 26

                        # Create astropy Time object
                        time = Time(
                            f'{year}-{month:02d}-{day:02d}', format='iso')

                        data.append({
                            'date_julian': time.jd,
                            'date_julian_year': time.decimalyear,
                            'residual_m': residual_m,
                            'timestamp': timestamp,
                            'station': 'DE430'
                        })
                    except (ValueError, IndexError) as e:
                        # If date decoding fails, skip this line
                        print_status(
                            f"Line {line_num} date decode error: {e}", "DEBUG")
                        continue
                else:
                    print_status(
                        f"Line {line_num}: timestamp too short: {timestamp[:20]}...", "DEBUG")
                    continue

            except (ValueError, IndexError) as e:
                print_status(f"Line {line_num} parse error: {e}", "DEBUG")
                continue

    df = pd.DataFrame(data)
    print_status(f"Parsed {len(df)} observations from {filepath.name}", "INFO")
    return df


def process_de430_file(input_path: Path, output_path: Path):
    """
    Process DE430 file and output to TEP format.

    Args:
        input_path: Path to input DE430 file
        output_path: Path to output CSV file
    """
    print_status(f"Processing DE430 data from {input_path}", "PROCESS")

    # Parse DE430 format
    df = parse_de430_file(input_path)

    if len(df) == 0:
        print_status(f"No data found in {input_path}", "WARNING")
        return

    # Compute elongation
    print_status("Computing Moon elongation angles...", "PROCESS")
    df = compute_elongation(df)

    # Select output columns
    output_columns = ['date_julian', 'date_julian_year',
                      'residual_m', 'elongation_rad', 'station']
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


if __name__ == "__main__":
    log_dir = PROJECT_ROOT / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    logger = TEPLogger("parse_de430", str(log_dir / "parse_de430.log"))
    set_step_logger(logger)

    data_dir = PROJECT_ROOT / "data" / "raw"
    output_dir = PROJECT_ROOT / "data" / "processed"
    output_dir.mkdir(parents=True, exist_ok=True)

    input_path = data_dir / "DE430_2014-2018_residuals.dat"
    output_path = output_dir / "DE430_residuals.csv"

    process_de430_file(input_path, output_path)
