# Temporal Equivalence Principle: Lunar Laser Ranging and the Nordtvedt Effect

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19446028.svg)](https://doi.org/10.5281/zenodo.19446028)
[![License: CC BY 4.0](https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)

![TEP-LLR: Lunar Laser Ranging](site/public/image.webp)

**Author:** Matthew Lukin Smawfield  
**Version:** v0.3 (Lucknow)  
**First published:** 17 May 2026 · **Last updated:** 16 September 2026
**Status:** Preprint (Open for Collaboration)  
**DOI:** [10.5281/zenodo.19446028](https://doi.org/10.5281/zenodo.19446028)  
**Website:** [https://mlsmawfield.com/tep/llr/](https://mlsmawfield.com/tep/llr/)  
**Paper Series:** TEP Series: Paper 17 (Lunar Laser Ranging)    

## Abstract

The Temporal Equivalence Principle (TEP) proposes that proper time is a dynamical scalar field $\phi$ rather than a geometric parameter. This field couples universally to all matter via a conformal metric $\tilde{g}_{\mu\nu} = A^2(\phi)\,g_{\mu\nu}$; in deep gravitational potential wells, high ambient density suppresses the local gradient of $\phi$ — a mechanism termed Temporal Shear screening — with the degree of suppression scaling with a body's gravitational compactness $\Phi/c^2$. Because Earth and Moon have different compactness and interior shielding profiles, TEP motivates a compactness-dependent Strong Equivalence Principle response that would appear in Lunar Laser Ranging as a synodic Earth-Moon range modulation $\delta r = 13\eta_{\rm resid}\cos D$, where $D$ is the Moon-Sun elongation angle and $\eta_{\rm resid}$ is the residual-channel amplitude (an Observable Response Coefficient, distinct from the standard metric Nordtvedt parameter $\eta$). The expected residual-channel amplitude is at the millimetre level for $\eta_{\rm resid} \sim 10^{-4}$.

This work analyses 26,207 raw Lunar Laser Ranging O–C residuals from the public INPOP19a ephemeris archives (Paris Observatory, Geoazur), comprising $N = 25{,}445$ measurements from five international stations (1984–2019) after standard $6\sigma$ MAD outlier cleaning. The primary estimand is a precision-weighted synodic regression on the post-fit residual channel, with inverse station-variance weighting and a full systematic model with annual, monthly, and thermal $\cos(2D)$ regressors and a common intercept.

The primary precision-weighted residual-channel estimate is $\eta_{\rm resid} = -3.91 \times 10^{-4} \pm 5.63 \times 10^{-5}$ ($6.94\sigma$; $6.78\sigma$ cluster-robust). The amplitude stabilises in the band $-3.2$ to $-4.1 \times 10^{-4}$ across the estimator hierarchy, from cosD-only ($5.25\sigma$) through full-systematic OLS ($6.17\sigma$) to precision-weighted ($6.94\sigma$). Robustness checks — common-$\eta_{\rm resid}$ mixed model, leave-one-station-out meta-analysis, wild cluster bootstrap, phase-locked differential, cross-ephemeris validation on cleaned DE430 ($\eta_{\rm resid} = -7.04 \times 10^{-4}$, $5.96\sigma$), parametric GR-null bootstrap, and a frequency null scan — all support residual-channel survival of the synodic component (Section 4). If interpreted within the static Parametrized Post-Newtonian framework — an interpretation the residual-channel estimand does not yet license (Section 4.14) — the measured amplitude would correspond to $\beta_{\rm PPN} = 0.999902 \pm 1.07 \times 10^{-5}$ from the joint $(\beta_{\rm PPN}, \gamma)$ contour; this mapping is reported for scale only and is not advanced as a PPN measurement pending the integrator-level refit. A full-sky directional scan on the residual channel (2,664 right-ascension/declination grid points at 5° intervals) places the Planck CMB dipole axis at rank 226/2664 (top 8.5%).

The result is therefore framed as a high-significance residual-channel candidate with a TEP interpretation, not as a completed replacement for direct-fit LLR bounds. Source-level numerical refits of the INPOP or DE430 integrators with $\eta_{\rm resid}$ left free remain the critical open closure test. Ephemeris-absorption stress tests bound the residual-channel survival amplitude but do not replace a full dynamical integrator-level confirmation.

Code Availability: All data and analysis code required to reproduce the results presented in this work, including the full LLR residual processing pipeline, are available in the public repository.

## Key Results

**TEP Signal in INPOP19a LLR Residuals:**
- **26,207 observations** from 5 stations (APO, Grasse, Matera, McDonald2, Haleakala); **25,445** after $6\sigma$ MAD cleaning
- **Date range:** 1984–2019
- **Residual precision:** 9.5 cm RMS
- **Cook's-Distance-excised full-systematic OLS (leverage diagnostic):** $\eta_{\rm resid} = -3.87 \times 10^{-4} \pm 4.95 \times 10^{-5}$ ($N = 23{,}837$ after excision), 7.82σ; cluster-robust 8.65σ
- **Precision-weighted full-systematic (consensus):** $\eta_{\rm resid} = -3.91 \times 10^{-4} \pm 5.63 \times 10^{-5}$, 6.94σ; cluster-robust 6.78σ
- **Full-systematic OLS without excision (sensitivity upper bound):** $\eta_{\rm resid} = -4.06 \times 10^{-4} \pm 6.58 \times 10^{-5}$ ($N = 25{,}445$), 6.17σ; cluster-robust 6.52σ
- **Common-$\eta_{\rm resid}$ mixed model with station systematics (pooling):** $\eta_{\rm resid} = -4.32 \times 10^{-4} \pm 6.76 \times 10^{-5}$, 6.40σ; $F(4, 25{,}410) = 1.19$, $p = 0.31$
- **Phase-locked new/full-moon differential (robustness):** $\eta_{\rm resid} = -5.95 \times 10^{-4} \pm 1.01 \times 10^{-4}$, 5.91σ
- **cosD-only OLS (baseline):** $\eta_{\rm resid} = -3.18 \times 10^{-4} \pm 6.05 \times 10^{-5}$, 5.25σ

---

## The TEP Research Program

| Paper | Repository | Title | DOI |
|-------|-----------|-------|-----|
| **Paper 0** | [TEP](https://github.com/matthewsmawfield/TEP) | Temporal Equivalence Principle: Dynamic Time & Emergent Light Speed | [10.5281/zenodo.16921911](https://doi.org/10.5281/zenodo.16921911) |
| **Paper 1** | [TEP-GNSS](https://github.com/matthewsmawfield/TEP-GNSS) | Global Time Echoes: Distance-Structured Correlations in GNSS Clocks | [10.5281/zenodo.17127229](https://doi.org/10.5281/zenodo.17127229) |
| **Paper 2** | [TEP-GNSS-II](https://github.com/matthewsmawfield/TEP-GNSS-II) | Global Time Echoes: 25-Year Temporal Evolution | [10.5281/zenodo.17517141](https://doi.org/10.5281/zenodo.17517141) |
| **Paper 3** | [TEP-GNSS-RINEX](https://github.com/matthewsmawfield/TEP-GNSS-RINEX) | Global Time Echoes: Raw RINEX Validation of Distance-Structured Correlations in GNSS Clocks | [10.5281/zenodo.17860166](https://doi.org/10.5281/zenodo.17860166) |
| **Paper 4** | [TEP-GL](https://github.com/matthewsmawfield/TEP-GL) | Temporal-Spatial Coupling in Gravitational Lensing: A Reinterpretation of Dark Matter Observations | [10.5281/zenodo.17982540](https://doi.org/10.5281/zenodo.17982540) |
| **Paper 5** | [TEP-GTE](https://github.com/matthewsmawfield/TEP-GTE) | Global Time Echoes: Empirical Validation of the Temporal Equivalence Principle | [10.5281/zenodo.18004832](https://doi.org/10.5281/zenodo.18004832) |
| **Paper 6** | [TEP-UCD](https://github.com/matthewsmawfield/TEP-UCD) | Universal Critical Density: Unifying Atomic, Galactic, and Compact Object Scales | [10.5281/zenodo.18064365](https://doi.org/10.5281/zenodo.18064365) |
| **Paper 7** | [TEP-RBH](https://github.com/matthewsmawfield/TEP-RBH) | The Soliton Wake: A Runaway Black Hole as a Gravitational Soliton | [10.5281/zenodo.18059250](https://doi.org/10.5281/zenodo.18059250) |
| **Paper 8** | [TEP-SLR](https://github.com/matthewsmawfield/TEP-SLR) | Global Time Echoes: Optical-Domain Consistency Test via Satellite Laser Ranging | [10.5281/zenodo.18064581](https://doi.org/10.5281/zenodo.18064581) |
| **Paper 9** | [TEP-EXP](https://github.com/matthewsmawfield/TEP-EXP) | What Do Precision Tests of General Relativity Actually Measure? | [10.5281/zenodo.18109760](https://doi.org/10.5281/zenodo.18109760) |
| **Paper 10** | [TEP-COS](https://github.com/matthewsmawfield/TEP-COS) | The Temporal Equivalence Principle: Suppressed Density Scaling in Globular Cluster Pulsars | [10.5281/zenodo.18165798](https://doi.org/10.5281/zenodo.18165798) |
| **Paper 11** | [TEP-H0](https://github.com/matthewsmawfield/TEP-H0) | The Cepheid Bias: Resolving the Hubble Tension | [10.5281/zenodo.18209702](https://doi.org/10.5281/zenodo.18209702) |
| **Paper 12** | [TEP-JWST](https://github.com/matthewsmawfield/TEP-JWST) | The Temporal Equivalence Principle: A Unified Resolution to the JWST High-Redshift Anomalies | [10.5281/zenodo.19000827](https://doi.org/10.5281/zenodo.19000827) |
| **Paper 13** | [TEP-WB](https://github.com/matthewsmawfield/TEP-WB) | The Temporal Equivalence Principle: Temporal Shear Recovery in Gaia DR3 Wide Binaries | [10.5281/zenodo.19102061](https://doi.org/10.5281/zenodo.19102061) |
| **Paper 15** | [TEP-EFA](https://github.com/matthewsmawfield/TEP-EFA) | Temporal Equivalence Principle: Temporal Shear in the Earth Flyby Anomaly | [10.5281/zenodo.19454862](https://doi.org/10.5281/zenodo.19454862) |
| **Paper 16** | [TEP-J0437](https://github.com/matthewsmawfield/TEP-J0437) | Synchronization Holonomy in Pulsar Scintillation | [10.5281/zenodo.19454620](https://doi.org/10.5281/zenodo.19454620) |
| **Paper 17** | **TEP-LLR** (This repo) | Lunar Laser Ranging and the Nordtvedt Effect | [10.5281/zenodo.19446028](https://doi.org/10.5281/zenodo.19446028) |

## Directory Structure

```text
TEP-LLR/
├── archive/                # Archived old/unused scripts and files
├── data/
│   ├── raw/                 # INPOP19a residual files (MINI format)
│   └── processed/           # Parsed residuals with elongation angles
├── logs/                    # Execution logs
├── manuscripts/             # Generated PDF/Markdown outputs
├── results/                 # Analytical outputs and figures
├── scripts/
│   ├── steps/               # Sequential analysis pipeline
│   │   ├── step_000_llr_data_ingestion.py      # Verify raw data availability and hashes
│   │   ├── step_001_data_preprocessing.py      # Parse MINI format
│   │   ├── step_002_de430_preprocessing.py     # DE430 ephemeris processing
│   │   ├── step_003_statistical_analysis.py    # Basic TEP detection analysis
│   │   ├── step_004_detection_analysis_advanced.py  # Advanced analysis (M4 Pro optimized)
│   │   ├── ... (87 canonical scripts: step_000 through step_082, including 006b and 046b)
│   │   └── run_all_steps.py                     # Run complete pipeline
│   └── utils/               # Shared utilities
│       ├── llr_constants.py                       # Physical constants
│       ├── parse_de430.py                         # DE430 parser
│       ├── parse_inpop_mini.py                    # INPOP parser
│       ├── logger.py                              # Logging utilities
│       ├── pipeline_runner.py                     # Pipeline execution
│       ├── pipeline_quality_gate.py               # Reviewer-facing audit gate
│       ├── generate_evidence_ledger.py            # Evidence summary artifact
│       ├── schema_validation.py                   # Output schema checks
│       ├── statistical_utils.py                   # Statistical utilities
│       └── verify_value_consistency.py            # Manuscript value audit
├── site/
│   └── components/          # HTML source of truth for manuscript
├── README.md
└── requirements.txt         # Python dependencies
```

## Installation

```bash
# Clone repository
git clone https://github.com/matthewsmawfield/TEP-LLR.git
cd TEP-LLR

# Install Python dependencies
pip install -r requirements.txt

```

## Essential Data Files

- `data/processed/INPOP19a_all_stations_residuals.csv` - Main LLR residuals dataset
- `data/raw/INPOP19a_*_residuals.txt` - INPOP19a MINI format residual files verified by `data/raw/data_manifest.json`
- `data/raw/DE430_2014-2018_residuals.dat` - DE430 residual archive verified by `data/raw/data_manifest.json`

## Data Sources

- **Paris Observatory (Geoazur):** INPOP19a lunar ephemerides with O-C residuals (primary data source)
- **Stations:** Apache Point (APO), Grasse, Matera, McDonald2, Haleakala
- **Ephemeris:** INPOP19a (2019) - most recent INPOP release with LLR residuals

## Reproduction Pipeline

```bash
# Run complete pipeline (recommended)
python scripts/steps/run_all_steps.py

# Validate structured outputs and manuscript consistency
python scripts/utils/schema_validation.py
python scripts/utils/verify_value_consistency.py

# Run the reviewer-facing quality gate
python scripts/utils/pipeline_quality_gate.py

# Generate the evidence ledger directly
python scripts/utils/generate_evidence_ledger.py

# Or run steps individually:
# Step 0: Verify required raw residual files and checksums
python scripts/steps/step_000_llr_data_ingestion.py --verbose

# Step 1: Parse MINI format and compute elongation angles
python scripts/steps/step_001_data_preprocessing.py --verbose

# Step 2: DE430 ephemeris processing
python scripts/steps/step_002_de430_preprocessing.py --verbose

# Step 3: Run basic TEP detection analysis
python scripts/steps/step_003_statistical_analysis.py --verbose

# Step 4: Run advanced TEP detection analysis (M4 Pro optimized)
python scripts/steps/step_004_detection_analysis_advanced.py --verbose
```

## Citation

```bibtex
@article{smawfield2026llr,
  title={Temporal Equivalence Principle: Lunar Laser Ranging and the Nordtvedt Effect},
  author={Smawfield, Matthew Lukin},
  journal={Zenodo},
  year={2026},
  doi={10.5281/zenodo.19446028},
  note={Preprint v0.3 (Lucknow)}
}
```

---

## Open Science Statement

These are working preprints shared in the spirit of open science—all manuscripts, analysis code, and data products are openly available under Creative Commons and MIT licenses to encourage and facilitate replication. Feedback and collaboration are warmly invited and welcome.

---

**Contact:** matthew@mlsmawfield.com  
**ORCID:** [0009-0003-8219-3159](https://orcid.org/0009-0003-8219-3159)
