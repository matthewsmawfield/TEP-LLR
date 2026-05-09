#!/usr/bin/env python3
"""
CDDIS LLR Data Download Script

Downloads LLR residuals from NASA Crustal Dynamics Data Information System (CDDIS).

REQUIRES NASA Earthdata Credentials:
1. Register at https://urs.earthdata.nasa.gov/
2. Create .netrc file with credentials:
   machine urs.earthdata.nasa.gov login <username> password <password>
3. Set .netrc permissions: chmod 600 ~/.netrc

Alternative: Use environment variables:
export NASA_USERNAME=<username>
export NASA_PASSWORD=<password>

Author: TEP-LLR Analysis Pipeline
Date: 2026-05-09
"""

import os
import sys
from pathlib import Path
import requests
from datetime import datetime, timedelta
import logging

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from scripts.utils.logger import TEPLogger, set_step_logger, print_status

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class CDDISDownloader:
    """Download LLR data from CDDIS archive."""

    def __init__(self, username=None, password=None):
        """
        Initialize CDDIS downloader.

        Args:
            username: NASA Earthdata username (or None to use env var)
            password: NASA Earthdata password (or None to use env var)
        """
        self.username = username or os.environ.get('NASA_USERNAME')
        self.password = password or os.environ.get('NASA_PASSWORD')
        self.session = None

        if not self.username or not self.password:
            print_status(
                "NASA Earthdata credentials required. "
                "Set NASA_USERNAME and NASA_PASSWORD environment variables, "
                "or provide them as arguments.", "WARNING")
            print_status(
                "Register at: https://urs.earthdata.nasa.gov/", "INFO")

    def authenticate(self):
        """Authenticate with NASA Earthdata."""
        if not self.username or not self.password:
            raise ValueError("NASA Earthdata credentials not provided")

        # NASA Earthdata authentication URL
        auth_url = "https://urs.earthdata.nasa.gov/oauth/tokens"

        print_status(f"Authenticating as {self.username}...", "PROCESS")
        self.session = requests.Session()
        self.session.auth = (self.username, self.password)

        # Test authentication
        try:
            response = self.session.get("https://urs.earthdata.nasa.gov/profile")
            if response.status_code == 200:
                print_status("Authentication successful", "SUCCESS")
            else:
                raise Exception(f"Authentication failed: {response.status_code}")
        except Exception as e:
            raise Exception(f"Authentication error: {e}")

    def download_llr_residuals(self, start_year=2019, end_year=2024, output_dir=None):
        """
        Download LLR residuals from CDDIS for specified year range.

        Args:
            start_year: Start year (default: 2019)
            end_year: End year (default: 2024)
            output_dir: Output directory (default: data/raw/cddis)

        Returns:
            List of downloaded file paths
        """
        if not self.session:
            self.authenticate()

        if output_dir is None:
            output_dir = PROJECT_ROOT / "data" / "raw" / "cddis"

        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        downloaded_files = []

        # CDDIS LLR data directory structure
        # Base URL for LLR data
        base_url = "https://cddis.nasa.gov/archive/slr/data"

        print_status(f"Downloading LLR residuals for {start_year}-{end_year}...", "PROCESS")

        for year in range(start_year, end_year + 1):
            print_status(f"Processing year {year}...", "INFO")

            # CDDIS LLR data is organized by year and station
            # This is a placeholder - actual URL structure needs to be verified
            # Typical LLR data format: https://cddis.nasa.gov/archive/slr/data/llr/YYYY/station/

            # Example URL (needs verification):
            year_url = f"{base_url}/llr/{year}/"

            try:
                response = self.session.get(year_url)
                if response.status_code == 200:
                    # Parse response to find available files
                    # This is a placeholder - actual parsing depends on CDDIS structure
                    print_status(f"Found data for year {year}", "INFO")

                    # Download files (placeholder)
                    # Actual implementation depends on CDDIS file structure
                    # year_files = self._parse_cddis_directory(response.text)
                    # for file_info in year_files:
                    #     file_path = self._download_file(file_info, output_dir)
                    #     if file_path:
                    #         downloaded_files.append(file_path)

                else:
                    print_status(f"No data found for year {year} (HTTP {response.status_code})", "WARNING")

            except Exception as e:
                print_status(f"Error processing year {year}: {e}", "ERROR")
                continue

        print_status(f"Downloaded {len(downloaded_files)} files", "SUCCESS")
        return downloaded_files

    def download_de430_residuals(self, output_dir=None):
        """
        Attempt to download extended DE430 residuals if available.

        Args:
            output_dir: Output directory (default: data/raw/cddis)

        Returns:
            Path to downloaded file or None if not available
        """
        if not self.session:
            self.authenticate()

        if output_dir is None:
            output_dir = PROJECT_ROOT / "data" / "raw" / "cddis"

        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        print_status("Searching for extended DE430 residuals...", "PROCESS")

        # Search for DE430 residuals in CDDIS archive
        # This is exploratory - actual availability needs verification
        search_urls = [
            "https://cddis.nasa.gov/archive/slr/products/ephemerides/",
            "https://cddis.nasa.gov/archive/slr/data/llr/ephemerides/",
        ]

        for url in search_urls:
            try:
                response = self.session.get(url)
                if response.status_code == 200:
                    print_status(f"Found ephemeris directory: {url}", "INFO")
                    # Parse and download if DE430 residuals found
                    # This requires knowledge of CDDIS file structure
                else:
                    print_status(f"No data at {url} (HTTP {response.status_code})", "DEBUG")
            except Exception as e:
                print_status(f"Error accessing {url}: {e}", "DEBUG")

        print_status("Extended DE430 residuals search complete", "INFO")
        return None

    def _download_file(self, file_url, output_dir):
        """Download a single file from CDDIS."""
        try:
            filename = file_url.split('/')[-1]
            output_path = output_dir / filename

            print_status(f"Downloading {filename}...", "PROCESS")
            response = self.session.get(file_url, stream=True)

            if response.status_code == 200:
                with open(output_path, 'wb') as f:
                    for chunk in response.iter_content(chunk_size=8192):
                        f.write(chunk)

                print_status(f"Downloaded {filename}", "SUCCESS")
                return output_path
            else:
                print_status(f"Failed to download {filename} (HTTP {response.status_code})", "WARNING")
                return None

        except Exception as e:
            print_status(f"Error downloading {file_url}: {e}", "ERROR")
            return None


def main():
    """Main execution function."""
    log_dir = PROJECT_ROOT / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)

    logger_instance = TEPLogger("cddis_download", str(log_dir / "cddis_download.log"))
    set_step_logger(logger_instance)

    print_status("CDDIS LLR Data Download", "TITLE")
    print_status("="*80, "TITLE")

    # Initialize downloader
    downloader = CDDISDownloader()

    # Check if credentials are available
    if not downloader.username or not downloader.password:
        print_status("NASA Earthdata credentials not found", "WARNING")
        print_status("To use this script:", "INFO")
        print_status("1. Register at https://urs.earthdata.nasa.gov/", "INFO")
        print_status("2. Set environment variables:", "INFO")
        print_status("   export NASA_USERNAME=<your_username>", "INFO")
        print_status("   export NASA_PASSWORD=<your_password>", "INFO")
        print_status("3. Or create .netrc file:", "INFO")
        print_status("   machine urs.earthdata.nasa.gov login <username> password <password>", "INFO")
        print_status("   chmod 600 ~/.netrc", "INFO")
        print_status("="*80, "TITLE")
        return

    try:
        # Attempt to download extended DE430 residuals
        print_status("\nAttempting to download extended DE430 residuals...", "TITLE")
        de430_file = downloader.download_de430_residuals()

        if de430_file:
            print_status(f"Extended DE430 residuals downloaded: {de430_file}", "SUCCESS")
        else:
            print_status("Extended DE430 residuals not found on CDDIS", "INFO")
            print_status("Attempting to download general LLR residuals 2019-2024...", "INFO")

            # Download LLR residuals for extended baseline
            downloaded_files = downloader.download_llr_residuals(
                start_year=2019,
                end_year=2024
            )

            if downloaded_files:
                print_status(f"Downloaded {len(downloaded_files)} LLR data files", "SUCCESS")
            else:
                print_status("No LLR data files downloaded", "WARNING")

    except Exception as e:
        print_status(f"Error during download: {e}", "ERROR")
        print_status("Check credentials and network connection", "INFO")


if __name__ == "__main__":
    main()
