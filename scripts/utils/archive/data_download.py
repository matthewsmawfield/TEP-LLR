#!/usr/bin/env python3
"""
LLR Data Download Utilities

This module provides utilities for downloading LLR data and ephemerides from
various public sources. Uses existing TEP-SLR authentication infrastructure.

Data Sources:
1. NASA CDDIS/Earthdata - LLR station data (uses ~/.netrc or CDDIS_USER/CDDIS_PASS)
2. JPL SSD - DE440/DE441 ephemerides (publicly accessible)
3. IMCCE - INPOP21 ephemerides (publicly accessible)
4. ILRS - Station-specific data (requires ILRS access)

Author: TEP-LLR Analysis Pipeline
Date: 2024
"""

import os
import netrc
import requests
import logging
from pathlib import Path
from urllib.parse import urljoin
import json

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def get_auth():
    """
    Get NASA Earthdata/CDDIS authentication credentials.

    Uses the same authentication method as TEP-SLR data acquisition.
    Checks ~/.netrc file first, then environment variables.

    Returns:
        Tuple of (username, password) or None if not found
    """
    try:
        auth = netrc.netrc().authenticators("urs.earthdata.nasa.gov")
        if auth:
            logger.info(
                "Using credentials from ~/.netrc for urs.earthdata.nasa.gov")
            return (auth[0], auth[2])
    except Exception as e:
        logger.debug(f"No .netrc file found: {e}")

    user = os.getenv("CDDIS_USER")
    passwd = os.getenv("CDDIS_PASS")
    if user and passwd:
        logger.info(
            "Using credentials from CDDIS_USER/CDDIS_PASS environment variables")
        return (user, passwd)

    logger.warning("No NASA Earthdata credentials found")
    return None


class LLRDataDownloader:
    """Download LLR data from various sources."""

    def __init__(self, download_dir: str = 'data/raw'):
        self.download_dir = Path(download_dir)
        self.download_dir.mkdir(parents=True, exist_ok=True)

    def download_jpl_ephemeris(self, ephemeris: str = 'DE440') -> dict:
        """
        Download JPL ephemeris data.

        JPL ephemerides are publicly accessible from:
        https://ssd.jpl.nasa.gov/ftp/eph/planets/Linux/

        Args:
            ephemeris: Ephemeris name (e.g., 'DE440', 'DE441')

        Returns:
            Dictionary with download results
        """
        logger.info(
            f"Attempting to download {ephemeris} ephemeris from JPL...")

        # JPL NAIF ephemeris base URL
        base_url = "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/"

        # Map ephemeris names to file patterns
        ephemeris_files = {
            'DE440': 'de440.bsp',
            'DE441': 'de441.bsp',
            'DE430': 'de430.bsp',
            'DE421': 'de421.bsp'
        }

        if ephemeris not in ephemeris_files:
            return {
                'success': False,
                'error': f'Unknown ephemeris: {ephemeris}',
                'available': list(ephemeris_files.keys())
            }

        file_name = ephemeris_files[ephemeris]
        url = urljoin(base_url, file_name)
        output_path = self.download_dir / file_name

        try:
            response = requests.get(url, stream=True, timeout=30)
            response.raise_for_status()

            total_size = int(response.headers.get('content-length', 0))
            logger.info(f"File size: {total_size / (1024*1024):.2f} MB")

            with open(output_path, 'wb') as f:
                for chunk in response.iter_content(chunk_size=8192):
                    if chunk:
                        f.write(chunk)

            logger.info(
                f"Successfully downloaded {file_name} to {output_path}")

            return {
                'success': True,
                'file': str(output_path),
                'size_bytes': output_path.stat().st_size,
                'url': url
            }

        except requests.exceptions.RequestException as e:
            logger.error(f"Failed to download {ephemeris}: {e}")
            return {
                'success': False,
                'error': str(e),
                'manual_instructions': f"Manually download from {url}"
            }

    def download_inpop_ephemeris(self, ephemeris: str = 'INPOP21a') -> dict:
        """
        Download INPOP ephemeris from IMCCE.

        INPOP ephemerides are publicly accessible from:
        https://ftp.imcce.fr/pub/ephem/planets/inpop21a/

        Args:
            ephemeris: Ephemeris name (e.g., 'INPOP21a', 'INPOP19a')

        Returns:
            Dictionary with download results
        """
        logger.info(
            f"Attempting to download {ephemeris} ephemeris from IMCCE...")

        # IMCCE ephemeris base URL
        base_url = "https://ftp.imcce.fr/pub/ephem/planets/"

        # Map ephemeris names to file patterns (use .dat files for JPL-compatible format)
        ephemeris_files = {
            'INPOP21a': 'inpop21a/inpop21a_TDB_m1000_p1000_littleendian.dat',
            'INPOP19a': 'inpop19a/inpop19a_TDB_m1000_p1000_littleendian.dat',
            'INPOP17a': 'inpop17a/inpop17a_TDB_m1000_p1000_littleendian.dat'
        }

        if ephemeris not in ephemeris_files:
            return {
                'success': False,
                'error': f'Unknown ephemeris: {ephemeris}',
                'available': list(ephemeris_files.keys())
            }

        file_name = ephemeris_files[ephemeris]
        url = urljoin(base_url, file_name)

        try:
            # Directly download the ephemeris file
            logger.info(f"Downloading {file_name}...")
            response = requests.get(url, timeout=120)
            response.raise_for_status()

            output_path = self.download_dir / file_name.split('/')[-1]
            with open(output_path, 'wb') as f:
                f.write(response.content)

            logger.info(f"Downloaded {file_name} to {output_path}")
            return {
                'success': True,
                'file': str(output_path),
                'size_bytes': output_path.stat().st_size,
                'url': url
            }

        except requests.exceptions.RequestException as e:
            logger.error(f"Failed to access {ephemeris}: {e}")
            return {
                'success': False,
                'error': str(e),
                'manual_instructions': f"Manually browse {url}"
            }

    def download_llr_data_cddis(self, station: str = None, years: list = None) -> dict:
        """
        Download LLR data from NASA CDDIS using existing authentication.

        Args:
            station: Station name (e.g., 'GRASSE', 'APO')
            years: List of years to download

        Returns:
            Dictionary with download results
        """
        logger.info("Attempting to download LLR data from NASA CDDIS...")

        auth = get_auth()
        if not auth:
            logger.error("NASA Earthdata credentials not found")
            logger.info(
                "Configure ~/.netrc for urs.earthdata.nasa.gov or set CDDIS_USER/CDDIS_PASS")
            return {
                'success': False,
                'error': 'Authentication required',
                'instructions': 'Configure ~/.netrc or CDDIS_USER/CDDIS_PASS environment variables'
            }

        session = requests.Session()
        session.auth = auth
        session.headers.update({"User-Agent": "TEP-LLR/data_download"})

        # CDDIS LLR data directory
#         base_url = "https://cddis.nasa.gov/archive/slr/products/resource/" # autofixed F841 base_url

        try:
            # Download ILRS Data Handling File
            logger.info("Downloading ILRS Data Handling File...")
            ilrs_url = "https://cddis.nasa.gov/archive/slr/products/resource/ILRS_Data_Handling_File.snx"
            try:
                ilrs_resp = session.get(ilrs_url, timeout=60)
                if ilrs_resp.status_code == 200:
                    output_path = self.download_dir / 'ILRS_Data_Handling_File.snx'
                    with open(output_path, 'wb') as f:
                        f.write(ilrs_resp.content)
                    logger.info(
                        f"Downloaded ILRS Data Handling File to {output_path}")
            except Exception as e:
                logger.warning(
                    f"Could not download ILRS Data Handling File: {e}")

            # Download SLRF2020 station coordinates
            logger.info("Downloading SLRF2020 station coordinates...")
            slrf_url = "https://cddis.nasa.gov/archive/slr/products/resource/SLRF2020_POS+VEL.snx"
            try:
                slrf_resp = session.get(slrf_url, timeout=60)
                if slrf_resp.status_code == 200:
                    output_path = self.download_dir / 'SLRF2020_POS+VEL.snx'
                    with open(output_path, 'wb') as f:
                        f.write(slrf_resp.content)
                    logger.info(f"Downloaded SLRF2020 to {output_path}")
            except Exception as e:
                logger.warning(f"Could not download SLRF2020: {e}")

            # Download LLR normal point data for key stations
            # First, list available SLR data directories
            logger.info("Listing CDDIS SLR archive structure...")
            slr_base_url = "https://cddis.nasa.gov/archive/slr/"
            try:
                slr_resp = session.get(slr_base_url, timeout=30)
                if slr_resp.status_code == 200:
                    logger.info("Successfully accessed SLR archive base")
                    # List available directories
                    import re
                    dirs = re.findall(r'href="([^"]*/)"', slr_resp.text)
                    logger.info(f"Available directories: {dirs[:10]}")
            except Exception as e:
                logger.warning(f"Could not list SLR directories: {e}")

            # Download LLR normal point data from CDDIS CRD archive
            logger.info(
                "Downloading LLR normal point data from CDDIS CRD archive...")

            # CRD normal point data directories
            crd_urls = [
                "https://cddis.nasa.gov/archive/slr/data/npt_crd/",
                "https://cddis.nasa.gov/archive/slr/data/npt_crd_v2/",
            ]

            results = {'success': True, 'downloaded': [], 'crd_downloaded': {}}

            for crd_url in crd_urls:
                try:
                    logger.info(f"Accessing CRD directory: {crd_url}")
                    crd_resp = session.get(crd_url, timeout=30)
                    if crd_resp.status_code == 200:
                        logger.info("Successfully accessed CRD directory")

                        # Look for Moon/Luna directories
                        import re
                        directories = re.findall(
                            r'href="([^"]*/)"', crd_resp.text)
                        moon_dirs = [
                            d for d in directories if 'luna' in d.lower() or 'moon' in d.lower()]

                        if moon_dirs:
                            logger.info(
                                f"Found Moon directories: {moon_dirs[:5]}")

                            # Download recent CRD files from Moon directories
                            # Limit to first 3 directories
                            for moon_dir in moon_dirs[:3]:
                                dir_url = f"{crd_url.rstrip('/')}/{moon_dir}"
                                try:
                                    dir_resp = session.get(dir_url, timeout=30)
                                    if dir_resp.status_code == 200:
                                        # Look for CRD files
                                        crd_files = re.findall(
                                            r'href="([^"]*\.crd)"', dir_resp.text)
                                        if crd_files:
                                            logger.info(
                                                f"  Found {len(crd_files)} CRD files in {moon_dir}")
                                            downloaded_count = 0
                                            # Last 10 files
                                            for crd_file in crd_files[-10:]:
                                                file_url = f"{dir_url}/{crd_file}"
                                                try:
                                                    file_resp = session.get(
                                                        file_url, timeout=60)
                                                    if file_resp.status_code == 200:
                                                        output_path = self.download_dir / \
                                                            f"llr_{moon_dir}_{crd_file}"
                                                        with open(output_path, 'wb') as f:
                                                            f.write(
                                                                file_resp.content)
                                                        results['downloaded'].append(
                                                            str(output_path))
                                                        downloaded_count += 1
                                                        logger.info(
                                                            f"  Downloaded {crd_file} to {output_path}")
                                                except Exception as e:
                                                    logger.warning(
                                                        f"  Could not download {crd_file}: {e}")
                                            results['crd_downloaded'][moon_dir] = downloaded_count
                                except Exception as e:
                                    logger.warning(
                                        f"Could not access directory {moon_dir}: {e}")
                        else:
                            logger.info(
                                f"No Moon directories found in {crd_url}")
                        break
                except Exception as e:
                    logger.warning(f"Could not access {crd_url}: {e}")

            return results

        except Exception as e:
            logger.error(f"Error downloading from CDDIS: {e}")
            return {'success': False, 'error': str(e)}

    def download_wettzell_data(self) -> dict:
        """
        Download Wettzell station LLR data.

        Wettzell is an active LLR station in Germany. Data may be available
        through ILRS or directly from the station.

        Returns:
            Dictionary with download instructions
        """
        logger.info("Wettzell LLR data access information")

        instructions = {
            'station': 'Wettzell (WETL)',
            'location': 'Germany',
            'status': 'Active',
            'sources': [
                {
                    'name': 'ILRS Data Center',
                    'url': 'https://ilrs.gsfc.nasa.gov/data_and_products/data/index.html',
                    'authentication': 'Required'
                },
                {
                    'name': 'BKG (Bundesamt für Kartographie und Geodäsie)',
                    'url': 'https://www.bkg.bund.de/',
                    'authentication': 'Contact required'
                }
            ],
            'instructions': [
                '1. Contact ILRS or BKG for data access',
                '2. Request Wettzell LLR normal point data',
                '3. Place downloaded files in data/raw/WETZELL_residuals.txt',
                '4. Run step_007 to ingest the data'
            ]
        }

        return instructions


def main():
    """
    Main function to demonstrate data download capabilities.
    """
    logger.info("=" * 80)
    logger.info("LLR DATA DOWNLOAD UTILITIES")
    logger.info("=" * 80)
    logger.info("")

    downloader = LLRDataDownloader()

    # Attempt to download DE440 ephemeris
    logger.info("Attempting to download DE440 ephemeris...")
    de440_result = downloader.download_jpl_ephemeris('DE440')
    logger.info(f"DE440 download result: {de440_result['success']}")
    if not de440_result['success']:
        logger.info(
            f"Manual download required: {de440_result.get('manual_instructions', 'N/A')}")
    logger.info("")

    # Attempt to access INPOP21a ephemeris
    logger.info("Attempting to access INPOP21a ephemeris...")
    inpop_result = downloader.download_inpop_ephemeris('INPOP21a')
    logger.info(f"INPOP21a access result: {inpop_result['success']}")
    if inpop_result.get('manual_download_required'):
        logger.info(
            f"Manual download required: {inpop_result['instructions']}")
    logger.info("")

    # LLR data from CDDIS
    logger.info("LLR data from NASA CDDIS...")
    cddis_result = downloader.download_llr_data_cddis()
    if cddis_result.get('success'):
        logger.info(
            f"Downloaded {len(cddis_result['downloaded'])} files from CDDIS")
        for file_path in cddis_result['downloaded']:
            logger.info(f"  - {file_path}")
    else:
        logger.info(
            f"Download failed: {cddis_result.get('error', 'Unknown error')}")
    logger.info("")

    # Wettzell data
    logger.info("Wettzell station data...")
    wettzell_result = downloader.download_wettzell_data()
    logger.info(f"Station: {wettzell_result['station']}")
    logger.info(f"Status: {wettzell_result['status']}")
    logger.info("")

    # Save summary
    summary = {
        'download_attempts': {
            'DE440': de440_result,
            'INPOP21a': inpop_result,
            'CDDIS_LL': cddis_result,
            'Wettzell': wettzell_result
        },
        'recommendations': [
            'Register for NASA Earthdata to access CDDIS LLR data',
            'Download DE440 ephemeris manually if automated download fails',
            'Contact ILRS for Wettzell station data',
            'Place all downloaded data in data/raw/ directory'
        ]
    }

    output_path = Path('results/outputs/data_download_summary.json')
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with open(output_path, 'w') as f:
        json.dump(summary, f, indent=2)

    logger.info(f"Download summary saved to {output_path}")
    logger.info("")
    logger.info("DATA DOWNLOAD UTILITIES DEMONSTRATION COMPLETE")


if __name__ == '__main__':
    main()
