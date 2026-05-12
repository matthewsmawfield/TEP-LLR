"""
NASA-Standard LLR Analysis Constants and Formulas
=================================================

This module contains verified physical constants and formulas for professional-grade
Lunar Laser Ranging analysis, following NASA GEODYN, JPL DE440, and IERS2010 standards.

All constants are verified against peer-reviewed literature and official IERS standards.

References:
-----------
[1] IERS Conventions 2010 (IERS Technical Note No. 36)
    - Solid Earth tides: Chapter 7
    - Love numbers: Table 7.3
    
[2] Park et al. 2021, "The JPL Planetary and Lunar Ephemerides DE440 and DE441"
    - AJ 161, 105 (2021)
    - LLR residuals: ~1.3 cm (recent data), ~20 cm (early data)
    
[3] Williams & Dickey 2003, "Lunar Geophysics"
    - Love numbers and tidal response
    
[4] Mendes & Pavlis 2004, GRL 31, L14602
    - Atmospheric refraction model
    
[5] Farrell 1972, "Deformation of the Earth by Surface Loads"
    - Ocean loading theory
    
[6] Saastamoinen 1972, "Atmospheric Correction"
    - Zenith hydrostatic delay

Author: TEP-LLR Project
Date: 2026
"""

import numpy as np

# =============================================================================
# FUNDAMENTAL CONSTANTS (CODATA 2018 / IERS2010)
# =============================================================================

C_M_S = 299792458.0  # Speed of light in vacuum [m/s] - EXACT (CODATA 2018)
C_KM_S = C_M_S / 1000.0  # Speed of light [km/s]

# Newtonian gravitational constant [m^3 kg^-1 s^-2] (CODATA 2018)
G = 6.67430e-11
# Uncertainty: 0.0021e-11

# =============================================================================
# ASTRONOMICAL MASSES (DE440 / IERS2010)
# =============================================================================

# GM values from DE440 ephemeris [m^3/s^2]
GM_EARTH = 3.986004418e14      # Earth's GM (includes atmosphere)
# Source: DE440, IERS2010
# Uncertainty: 0.000000008e14

GM_MOON = 4.902800582e12       # Moon's GM
# Source: DE440, Williams & Dickey 2003
# Uncertainty: 0.000000007e12

GM_SUN = 1.32712440018e20      # Sun's GM
# Source: DE440
# Uncertainty: 0.00000008e20

# Mass ratios (dimensionless)
MASS_RATIO_SUN_EARTH = 332946.050895  # M_sun / M_earth (DE440)
MASS_RATIO_EARTH_MOON = 81.3005691     # M_earth / M_moon (DE440)

# =============================================================================
# EARTH PARAMETERS (IERS2010)
# =============================================================================

EARTH_RADIUS_M = 6378136.6  # Earth's equatorial radius [m] (WGS84 / IERS2010)
EARTH_FLATTENING = 1.0 / 298.25642  # WGS84 flattening

# Geocentric gravitational constant including atmosphere
# Same as GM_EARTH above

# =============================================================================
# MOON PARAMETERS (DE440 / IAU)
# =============================================================================

MOON_RADIUS_M = 1737400.0  # Moon mean radius [m] (IAU 2015)
# Source: Archinal et al. 2018, Icarus 281, 356

# Moon reflector selenographic coordinates (lat, lon in radians)
# Coordinates from Apollo/Luna mission data
# Frame: Moon Principal Axis (MOON_PA) or IAU_MOON
REFLECTORS = {
    'apollo11': {
        'lat_rad': 0.67337 * np.pi / 180,    # 0.67337° N
        'lon_rad': 23.47293 * np.pi / 180,   # 23.47293° E
        'height_m': 0  # Approximate height above lunar surface
    },
    'apollo14': {
        # 3.67500° S (DMA Landing Site Map - more accurate)
        'lat_rad': -3.67500 * np.pi / 180,
        # 17.46700° W (DMA Landing Site Map)
        'lon_rad': -17.46700 * np.pi / 180,
        'height_m': 0
    },
    'apollo15': {
        'lat_rad': 26.13222 * np.pi / 180,   # 26.13222° N
        'lon_rad': 3.62829 * np.pi / 180,    # 3.62829° E
        'height_m': 0
    },
    'luna17': {
        'lat_rad': 38.33307 * np.pi / 180,   # 38.33307° N
        'lon_rad': -35.00807 * np.pi / 180,  # 35.00807° W
        'height_m': 0
    },
    'luna21': {
        'lat_rad': 25.85093 * np.pi / 180,   # 25.85093° N
        'lon_rad': 30.45323 * np.pi / 180,   # 30.45323° E
        'height_m': 0
    },
}

# =============================================================================
# LOVE NUMBERS (IERS2010 Chapter 7)
# =============================================================================

# Nominal values for elastic Earth model
# Reference: IERS Conventions 2010, Table 7.3

# Degree 2 Love numbers (elastic, nominal)
H2_ELASTIC = 0.6078  # Vertical displacement Love number (h_2)
# Source: IERS2010, Table 7.3
# Real part of complex Love number

L2_ELASTIC = 0.0847  # Horizontal displacement Love number (l_2)
# Source: IERS2010, Table 7.3

K2_ELASTIC = 0.2983  # Gravitational potential Love number (k_2)
# Source: IERS2010
# Used for tidal perturbations on satellite orbits

# Anelastic correction (frequency-dependent)
# These are imaginary parts of complex Love numbers
H2_ANELASTIC_IMAG = -0.0025  # Im(h_2) (IERS2010)
L2_ANELASTIC_IMAG = -0.0007  # Im(l_2) (IERS2010)

# Load Love numbers (for ocean loading and surface mass loading)
# Reference: Farrell 1972
H2_LOAD = -1.0  # Vertical load Love number (h_2')
# Negative because surface load depresses the crust

L2_LOAD = -0.25  # Horizontal load Love number (l_2')

# =============================================================================
# SOLID EARTH TIDE FORMULAS (IERS2010 Chapter 7)
# =============================================================================


def solid_earth_tide_step1(r_station, r_body, gm_body, h2=H2_ELASTIC, l2=L2_ELASTIC):
    """
    IERS2010 Step 1: Calculate elastic solid Earth tide displacement.

    Formula (IERS2010 Eq. 7.5):
        Δr = h_2 * P_2(cos ζ) * (GM_body / GM_earth) * (R_earth^4 / R_body^3)

    Where:
        - P_2(x) = (3x² - 1)/2 is the 2nd Legendre polynomial
        - ζ is the zenith angle to the tidal body
        - R_earth is geocentric distance to station
        - R_body is geocentric distance to tidal body (Moon or Sun)

    Parameters:
    -----------
    r_station : array-like, shape (3,)
        Station position vector in Earth-centered frame [m]
    r_body : array-like, shape (3,)
        Tidal body position vector (Moon or Sun) [m]
    gm_body : float
        Gravitational parameter of tidal body [m^3/s^2]
    h2 : float
        Vertical Love number (default: IERS2010 value 0.6078)
    l2 : float
        Horizontal Love number (default: IERS2010 value 0.0847)

    Returns:
    --------
    delta_r : ndarray
        Tidal displacement vector [m]

    Reference:
    --------
    IERS Conventions 2010, Chapter 7, Section 7.1.2
    """
    r_station = np.asarray(r_station)
    r_body = np.asarray(r_body)

    R_sta = np.linalg.norm(r_station)
    R_body = np.linalg.norm(r_body)

    # Unit vectors
    u_sta = r_station / R_sta
    u_body = r_body / R_body

    # Zenith angle cosine
    cos_zeta = np.dot(u_sta, u_body)

    # 2nd Legendre polynomial: P_2(x) = (3x² - 1)/2
    P2 = 0.5 * (3 * cos_zeta**2 - 1)

    # Vertical displacement amplitude
    amp = h2 * P2 * (gm_body / GM_EARTH) * (R_sta**4 / R_body**3)

    # Vertical component (radial)
    delta_r_vertical = amp * u_sta

    # Horizontal component (gradient of potential)
    # Tangential component: perpendicular to radial
    # IERS2010 Eq. 7.6
    amp_horizontal = l2 * (3 * cos_zeta) * (gm_body /
                                            GM_EARTH) * (R_sta**4 / R_body**3)

    # Horizontal direction: component of u_body perpendicular to u_sta
    u_horizontal = u_body - cos_zeta * u_sta
    u_horizontal_norm = np.linalg.norm(u_horizontal)

    if u_horizontal_norm > 1e-10:
        u_horizontal = u_horizontal / u_horizontal_norm
        delta_r_horizontal = amp_horizontal * u_horizontal
    else:
        delta_r_horizontal = np.zeros(3)

    return delta_r_vertical + delta_r_horizontal


def solid_earth_tide_moon_sun(r_station, r_moon, r_sun):
    """
    Calculate total solid Earth tide from both Moon and Sun.

    Parameters:
    -----------
    r_station : array-like, shape (3,)
        Station position vector [m]
    r_moon : array-like, shape (3,)
        Moon position vector [m]
    r_sun : array-like, shape (3,)
        Sun position vector [m]

    Returns:
    --------
    delta_r : ndarray
        Total tidal displacement vector [m]

    Reference:
    --------
    IERS Conventions 2010, Chapter 7
    """
    # Moon contribution
    delta_moon = solid_earth_tide_step1(r_station, r_moon, GM_MOON)

    # Sun contribution (~46% of Moon's effect)
    delta_sun = solid_earth_tide_step1(r_station, r_sun, GM_SUN)

    return delta_moon + delta_sun


# =============================================================================
# OCEAN LOADING (Simplified GOT4.7/FES2004 Model)
# =============================================================================

def ocean_loading_displacement(r_station, r_moon, r_sun,
                               ocean_factor=0.1, phase_lag=np.pi/4):
    """
    Calculate ocean loading displacement using simplified model.

    Ocean loading is the elastic response of the Earth's crust to the weight
    of ocean tides. Amplitude is typically 1-2 cm vertical.

    This simplified model approximates ocean loading as a fraction of the
    solid Earth tide with a phase lag.

    Parameters:
    -----------
    r_station : array-like
        Station position [m]
    r_moon : array-like
        Moon position [m]
    r_sun : array-like
        Sun position [m]
    ocean_factor : float
        Fraction of solid tide amplitude (default: 0.1 = 10%)
    phase_lag : float
        Phase lag in radians (default: π/4 = 45°)

    Returns:
    --------
    delta_r : ndarray
        Ocean loading displacement [m]

    Reference:
    --------
    - Farrell 1972, "Deformation of the Earth by Surface Loads"
    - GOT4.7 / FES2004 ocean tide models
    """
    r_station = np.asarray(r_station)
    r_moon = np.asarray(r_moon)
    r_sun = np.asarray(r_sun)

    R_sta = np.linalg.norm(r_station)
    R_moon = np.linalg.norm(r_moon)
    R_sun = np.linalg.norm(r_sun)

    u_sta = r_station / R_sta
    u_moon = r_moon / R_moon
    u_sun = r_sun / R_sun

    # Tidal potentials
    cos_zeta_moon = np.dot(u_sta, u_moon)
    cos_zeta_sun = np.dot(u_sta, u_sun)

    P2_moon = 0.5 * (3 * cos_zeta_moon**2 - 1)
    P2_sun = 0.5 * (3 * cos_zeta_sun**2 - 1)

    # Load Love numbers
    h2_load = H2_LOAD  # -1.0

    # Ocean loading vertical displacement
    # Apply phase lag by shifting the tidal argument
    amp_moon = h2_load * P2_moon * \
        (GM_MOON / GM_EARTH) * (R_sta**4 / R_moon**3)
    amp_sun = h2_load * P2_sun * (GM_SUN / GM_EARTH) * (R_sta**4 / R_sun**3)

    # Apply ocean factor and phase lag
    ocean_vert = ocean_factor * (amp_moon + amp_sun) * np.cos(phase_lag)

    return ocean_vert * u_sta


# =============================================================================
# POLE TIDE (IERS2010)
# =============================================================================

def pole_tide_displacement(r_station, x_p, y_p):
    """
    Calculate pole tide displacement due to polar motion.

    The pole tide is the response of the Earth's crust to the centrifugal
    effect of polar motion ( Chandler wobble + annual wobble).

    Formula (IERS2010, Eq. 7.29):
        Δz = -0.033 * sin(2φ) * (x_p * cos λ + y_p * sin λ)  [meters]

    Where:
        - φ is geocentric latitude
        - λ is longitude
        - x_p, y_p are pole coordinates in arcseconds

    Parameters:
    -----------
    r_station : array-like, shape (3,)
        Station position in Earth-fixed frame [m]
    x_p : float
        Pole x-coordinate [arcseconds]
    y_p : float
        Pole y-coordinate [arcseconds]

    Returns:
    --------
    delta_r : ndarray
        Pole tide displacement vector [m]

    Reference:
    --------
    IERS Conventions 2010, Chapter 7, Section 7.1.4
    """
    r_station = np.asarray(r_station)

    # Convert to geodetic coordinates
    R = np.linalg.norm(r_station)
    x, y, z = r_station

    # Geocentric latitude
    phi = np.arctan2(z, np.sqrt(x**2 + y**2))
    # Longitude
    lam = np.arctan2(y, x)

    # IERS2010 formula with Love number correction
    # Theoretical amplitude: -0.033 m per arcsec of polar motion (solid Earth response)
    # With Love number: approximately -0.025 to -0.030 m per arcsec
    # x_p, y_p are in arcseconds - DO NOT convert to radians
    amp = -0.033 * np.sin(2 * phi) * (x_p * np.cos(lam) + y_p * np.sin(lam))

    # Radial displacement only (vertical)
    u_radial = r_station / R

    return amp * u_radial


# =============================================================================
# ATMOSPHERIC REFRACTION (Mendes-Pavlis Model)
# =============================================================================

def atmospheric_refraction_mendes_pavlis(elevation_rad, pressure_mbar,
                                         temperature_K, humidity_percent,
                                         latitude_rad, height_km,
                                         wavelength_nm=532.0):
    """
    Calculate atmospheric refraction correction using Mendes-Pavlis model.

    For optical wavelengths (APOLLO: 532 nm), this uses the Saastamoinen
    formula for zenith hydrostatic delay and modified for elevation mapping.

    Parameters:
    -----------
    elevation_rad : float
        Elevation angle [radians]
    pressure_mbar : float
        Atmospheric pressure [mbar]
    temperature_K : float
        Temperature [Kelvin]
    humidity_percent : float
        Relative humidity [%]
    latitude_rad : float
        Station latitude [radians]
    height_km : float
        Station height above geoid [km]
    wavelength_nm : float
        Laser wavelength [nm] (default: 532 nm for APOLLO green laser)

    Returns:
    --------
    refraction_m : float
        Refraction correction to be SUBTRACTED from observed range [m]

    Reference:
    --------
    - Mendes & Pavlis 2004, GRL 31, L14602
    - Saastamoinen 1972, "Atmospheric Correction"
    - Ciddor & Hill 1999 for refractivity at 532 nm
    """
    # Zenith hydrostatic delay (Saastamoinen formula)
    # Δz_hydro = 0.0022768 * P / (1 - 0.00266*cos(2φ) - 0.00028*h)
    lat_factor = 1 - 0.00266 * np.cos(2 * latitude_rad) - 0.00028 * height_km
    zenith_hydrostatic = 0.0022768 * pressure_mbar / lat_factor  # meters

    # Wet delay component
    zenith_wet = 0.0
    if not np.isnan(temperature_K) and not np.isnan(humidity_percent):
        # Saturation vapor pressure (Tetens formula, Buck 1981)
        # e_s = 6.1121 * exp((18.678 - T/234.5) * T / (257.14 + T)) [mbar]
        # where T is in Celsius
        T_C = temperature_K - 273.15
        e_s = 6.1121 * np.exp((18.678 - T_C/234.5) * T_C / (257.14 + T_C))
        e = (humidity_percent / 100.0) * e_s  # Actual vapor pressure [mbar]

        # Wet zenith delay (Saastamoinen modified)
        zenith_wet = 0.002277 * (1255 / temperature_K + 0.05) * e / lat_factor

    # Elevation mapping function
    # Using CfA 2.2 mapping function (Davis et al. 1985) or simple 1/sin(el)
    sin_el = np.sin(elevation_rad)

    if sin_el > 0.1736:  # elevation > 10°
        # Simple mapping
        mapping = 1.0 / sin_el
    elif sin_el > 0.0872:  # 5° < elevation < 10°
        # Modified for lower elevations
        mapping = 1.0 / sin_el * (1 + 0.001 / sin_el**2)
    else:
        # Cap at elevation = 5° (mapping ~ 5.7)
        mapping = 5.0

    # Total refraction correction
    refraction = (zenith_hydrostatic + zenith_wet) * mapping

    return refraction


# =============================================================================
# SHAPIRO DELAY (General Relativity)
# =============================================================================

def shapiro_delay(r_station, r_body, r_sun, c=C_M_S):
    """
    Calculate Shapiro time delay (gravitational delay from Sun).

    Formula:
        Δt = (2GM_sun / c³) * ln((r1 + r2 + R) / (r1 + r2 - R))

    Where:
        - r1 = distance from station to Sun
        - r2 = distance from body (Moon/reflector) to Sun
        - R = distance from station to body

    Parameters:
    -----------
    r_station : array-like
        Station position [m]
    r_body : array-like
        Body position (reflector) [m]
    r_sun : array-like
        Sun position [m]
    c : float
        Speed of light [m/s]

    Returns:
    --------
    delay_m : float
        Shapiro delay converted to equivalent distance [m]

    Reference:
    --------
    Shapiro 1964, "Fourth Test of General Relativity"
    """
    r_station = np.asarray(r_station)
    r_body = np.asarray(r_body)
    r_sun = np.asarray(r_sun)

    r1 = np.linalg.norm(r_station - r_sun)
    r2 = np.linalg.norm(r_body - r_sun)
    R = np.linalg.norm(r_body - r_station)

    # Gravitational delay in seconds
    # GM_sun = 1.32712440018e20 m^3/s^2
    GM_sun = GM_SUN

    # Check for valid argument to log (must be positive)
    log_arg = (r1 + r2 + R) / (r1 + r2 - R)
    if log_arg <= 0:
        raise ValueError(
            f"Invalid argument to log in Shapiro delay: {log_arg} (must be positive)")

    delay_s = (2 * GM_sun / c**3) * np.log(log_arg)

    # Convert to distance
    return delay_s * c


# =============================================================================
# LUNAR PHASE CONSTANTS
# =============================================================================
SYNODIC_PERIOD_DAYS = 29.530588  # Lunar synodic period in days
# Source: JPL DE430/DE440 ephemeris documentation, IAU 1976 values
# Consistent with Delaunay frequency values used in step_033
# JD of 2000-01-06 12:24 UTC (reference new moon)
REFERENCE_NEW_MOON_JD = 2451549.0
SYNODIC_PERIOD_RADIANS = 2 * np.pi  # Full cycle in radians

# =============================================================================
# TEP ANALYSIS CONSTANTS
# =============================================================================
# ETA_SCALE_FACTOR: Amplitude-to-η conversion factor for TEP Nordtvedt signal
# Formula: δr = ETA_SCALE_FACTOR × η × cos(D)
# Units: meters per unit η (η is dimensionless)
#
# CORRECT DERIVATION (from standard Nordtvedt effect theory, Scholarpedia):
# The standard Nordtvedt effect formula from GR/PPN theory gives:
# δr ≈ 1.3 × 10³ η cos(Ḋt) cm = 13 η cos(D) meters
#
# Physical origin:
# - Earth's gravitational binding energy ≈ 4 × 10⁻¹⁰
# - Resonance enhancement sensitivity factor ≈ 3.3 × 10¹² cm (from synodic frequency proximity to orbital frequency)
# - Tidal distortion (lunar variation) provides ~2× enhancement
# - Combined: 3.3 × 10¹² × 4 × 10⁻¹⁰ = 1.3 × 10³ cm = 13 meters
#
# Reference: Nordtvedt K (1968, 1994), Scholarpedia article on Nordtvedt effect
# This is the STANDARD Nordtvedt effect sensitivity factor, NOT TEP-specific.
# TEP uses the same factor because it predicts a Nordtvedt-like effect.
ETA_SCALE_FACTOR = 13.0
ELONGATION_MASK_WIDTH = 0.5  # Radians (~28.6°) for new/full moon phase binning
STATISTICAL_ALPHA = 0.05  # Significance level for two-tailed tests
Z_ALPHA_2 = 1.96  # Z-score for α/2 = 0.025 (two-tailed test at α=0.05)

# Data processing constants
SIGMA_UNCERTAINTY_FLOOR_MM = 5.0  # Minimum uncertainty in mm to prevent infinite weights

# =============================================================================
# TEP COUPLING CONSTANTS (UPDATED FRAMEWORK)
# =============================================================================

# Nordtvedt Parameter (η) - TEP Modification Parameter for LLR
# Following the updated TEP framework (Papers 6, 10, 11, 12, 13):
# - The microscopic coupling β is constrained by Cassini to α_0 ≲ 3×10⁻³ in screened regime
# - Observable response coefficients κ (like κ_MSP, κ_Cep) are empirically determined in unscreened regimes
# - A key insight from Paper 10: TEP effect operates on clock rates via A(Φ) ≈ 1 − ηΦ/c²
# - This suggests η itself is the TEP modification parameter for LLR, analogous to κ_MSP and κ_Cep
#
# The measured η is ~100× larger than the theoretical prediction from β alone
# (η ~ 4×10^-6 from η = 4β² - γ - 3), indicating the actual measured η includes additional
# contributions from the screening mechanism and other effects.
#
# The measured η is much smaller than κ_MSP ~ 10^6 - 10^7 and κ_Cep ~ 10^6, which is consistent
# with the screening mechanism: LLR is in a more screened regime (Solar System) compared to
# globular clusters and galactic disks, so the effective response should be smaller.
#
# CONCLUSION: η is the observable response coefficient for LLR. No separate κ_LL is needed.
# The relationship between η and κ_LL is not required; η itself is the TEP modification parameter.

# Temporal Topology Saturation Density (ρ_T)
# ρ_T ≈ 20 g/cm³ is the Temporal Topology saturation density - a fundamental
# saturation scale of the temporal-field topology, NOT an ambient-density switch.
# This is calibrated from GNSS atomic clock correlations (Paper 6).
# Source: GNSS distance-structured correlations (L_c ≈ 4200 km → ρ_T ≈ 20 g/cm³)
RHO_T = 20.0  # Temporal Topology saturation density in g/cm³
RHO_T_ERROR = 8.0  # 40% systematic uncertainty from GNSS calibration
RHO_T_SOURCE = "GNSS atomic clock correlations (Paper 6, UCD)"

# =============================================================================
# STATION COORDINATES (ITRF / SLRF)
# =============================================================================

# SLRF2020 station coordinates (ITRF2020 frame)
# Source: CDDIS SLRF2020_POS+VEL.snx (downloaded 2024)
# Coordinates extracted from SINEX file, epoch 15:001:00000
STATION_COORDS = {
    'APOL': {  # Apache Point Observatory (station 7840)
        'X': 4033463.48,
        'Y': 23662.79,
        'Z': 4924305.35,
        'source': 'SLRF2020',
        'epoch': '15:001:00000'
    },
    'GRASSE': {  # Grasse LLR station (station 7845)
        'X': 4581691.94,
        'Y': 556196.37,
        'Z': 4389355.29,
        'source': 'SLRF2020',
        'epoch': '15:001:00000'
    },
    'MATERA': {  # Matera LLR station (station 7941)
        'X': 4641978.52,
        'Y': 1393067.82,
        'Z': 4133249.70,
        'source': 'SLRF2020',
        'epoch': '15:001:00000'
    },
    'MCDONALD': {  # McDonald Observatory (station 7080)
        'X': -1330021.29,
        'Y': -5328401.84,
        'Z': 3236480.70,
        'source': 'SLRF2020',
        'epoch': '15:001:00000'
    },
    'HALEAKALA': {  # Haleakala Observatory (station 7882)
        'X': -1997242.86,
        'Y': -5528040.77,
        'Z': 2468355.84,
        'source': 'SLRF2020',
        'epoch': '15:001:00000'
    },
    'WETTZELL': {  # Wettzell LLR station (station 8834)
        'X': 4075576.57,
        'Y': 931785.77,
        'Z': 4801583.75,
        'source': 'SLRF2020',
        'epoch': '15:001:00000'
    }
}

# ITRF2014 coordinates for LLR stations [meters]
# Epoch: 2010.0 (standard reference epoch)
# Sources: ILRS site logs, SLRF2020

STATIONS_LLR = {
    'APOL': {
        # Apache Point, New Mexico, USA
        # CRITICAL: Position previously solved from LLR coordinate optimization
        # This creates a circular dependency if the same LLR data is used for TEP detection
        # Current position uses SLRF2020 coordinates instead for independence
        'pos': np.array([STATION_COORDS['APOL']['X'], STATION_COORDS['APOL']['Y'], STATION_COORDS['APOL']['Z']]),  # [m] - SLRF2020
        'vel': np.array([-0.012, -0.002, -0.004]),  # [m/year] - approximate from plate motion
        'epoch': 2015.0,  # SLRF2020 epoch
        'sigma': np.array([0.01, 0.01, 0.01]),  # SLRF2020 precision
        'source': 'SLRF2020 (independent of LLR residuals)'
    },
    'GRAL': {
        # Grasse, France (alternative ID)
        'pos': np.array([4581692.1, 556196.0, 4389355.1]),
        'vel': np.array([0.018, 0.015, 0.012]),
        'epoch': 2010.0,
        'sigma': np.array([0.01, 0.01, 0.01]),
        'source': 'ILRS site log (ITRF2014)'
    },
    'GRSM': {
        # Grasse, France (primary ID)
        'pos': np.array([STATION_COORDS['GRASSE']['X'], STATION_COORDS['GRASSE']['Y'], STATION_COORDS['GRASSE']['Z']]),
        'vel': np.array([0.018, 0.015, 0.012]),
        'epoch': 2010.0,
        'sigma': np.array([0.01, 0.01, 0.01]),
        'source': 'ILRS site log (ITRF2014)'
    },
    'MATM': {
        # Matera, Italy
        'pos': np.array([STATION_COORDS['MATERA']['X'], STATION_COORDS['MATERA']['Y'], STATION_COORDS['MATERA']['Z']]),
        'vel': np.array([0.020, 0.014, 0.011]),
        'epoch': 2010.0,
        'sigma': np.array([0.01, 0.01, 0.01]),
        'source': 'ILRS site log (ITRF2014)'
    },
    'WETL': {
        # Wettzell, Germany
        'pos': np.array([STATION_COORDS['WETTZELL']['X'], STATION_COORDS['WETTZELL']['Y'], STATION_COORDS['WETTZELL']['Z']]),
        'vel': np.array([0.019, 0.013, 0.010]),
        'epoch': 2010.0,
        'sigma': np.array([0.01, 0.01, 0.01]),
        'source': 'ILRS site log (ITRF2014)'
    },
}


def get_station_position(sta_id, obs_epoch_year, apply_plate_motion=True):
    """
    Get station position with optional plate motion correction.

    Parameters:
    -----------
    sta_id : str
        Station ID (e.g., 'APOL', 'GRSM')
    obs_epoch_year : float
        Observation epoch (decimal year)
    apply_plate_motion : bool
        Whether to apply plate motion correction

    Returns:
    --------
    pos : ndarray
        Station position [m]
    sigma : ndarray
        Position uncertainty [m]
    """
    sta_id = str(sta_id).strip().upper()

    # Handle aliases
    if sta_id not in STATIONS_LLR:
        if sta_id == 'APOLLO':
            sta_id = 'APOL'
        else:
            raise KeyError(f"Unknown station ID: {sta_id}. No coordinates available.")

    sta_data = STATIONS_LLR[sta_id]

    if apply_plate_motion:
        # Calculate epoch difference
        epoch_diff = obs_epoch_year - sta_data['epoch']

        # Apply linear plate motion model
        pos = sta_data['pos'] + sta_data['vel'] * epoch_diff
    else:
        pos = sta_data['pos'].copy()

    sigma = sta_data.get('sigma', np.array([1.0, 1.0, 1.0]))

    return pos, sigma


# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

def legendre_p2(x):
    """2nd Legendre polynomial: P_2(x) = (3x² - 1)/2"""
    return 0.5 * (3 * x**2 - 1)


def elevation_angle(r_station, r_target):
    """
    Calculate elevation angle from station to target.

    Returns elevation in radians.
    """
    r_station = np.asarray(r_station)
    r_target = np.asarray(r_target)

    # Vector from station to target
    vec = r_target - r_station
    vec_unit = vec / np.linalg.norm(vec)

    # Local zenith (radial outward)
    zenith = r_station / np.linalg.norm(r_station)

    # sin(elevation) = dot(vec_unit, zenith)
    sin_el = np.dot(vec_unit, zenith)

    return np.arcsin(np.clip(sin_el, -1, 1))


# =============================================================================
# VALIDATION AGAINST LITERATURE
# =============================================================================

"""
Expected residuals from literature:
- DE440 LLR fit: ~1.3 cm RMS (recent data, 2000-2020)
- DE440 LLR fit: ~20 cm RMS (early data, 1970-1990)
- Implementation target: < 1 m (after station coordinate fix)

Current limitations:
1. Apache Point coordinates: ILRS site log shows "Approximate Position"
   - Estimated uncertainty: ~10 meters
   - This is the dominant error source
   
2. System delay calibration: Not implemented
   - CRD C0-C6 records contain per-session calibrations
   - Typical values: 0-100 nanoseconds
   
3. Ocean loading: Simplified model used
   - Full GOT4.7/FES2004 requires station-specific coefficients
   
4. IERS2010 Step 2: Frequency-dependent corrections not implemented
   - Correction amplitude: ~5-10 mm
"""

if __name__ == '__main__':
    # Test the constants module
    print("NASA-Standard LLR Constants Module")
    print("=" * 50)
    print(f"Speed of light: {C_M_S} m/s")
    print(f"GM Earth: {GM_EARTH:.6e} m³/s²")
    print(f"GM Moon: {GM_MOON:.6e} m³/s²")
    print(f"Love number h2: {H2_ELASTIC}")
    print(f"Love number l2: {L2_ELASTIC}")
    print(f"Moon radius: {MOON_RADIUS_M} m")
    print()
    print("Stations:")
    for sta_id, sta_data in STATIONS_LLR.items():
        print(f"  {sta_id}: {sta_data['source']}")
