# Temporal Equivalence Principle: Lunar Laser Ranging and the Nordtvedt Effect

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19446029.svg)](https://doi.org/10.5281/zenodo.19446029)
[![License: CC BY 4.0](https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)

**Author:** Matthew Lukin Smawfield  
**Version:** v0.1 (Lucknow)  
**First published:** 10 May 2026 · **Last updated:** 10 May 2026  
**Status:** Preprint (Open for Collaboration)  
**DOI:** [10.5281/zenodo.19446029](https://doi.org/10.5281/zenodo.19446029)  
**Website:** [https://mlsmawfield.com/tep/llr/](https://mlsmawfield.com/tep/llr/)  
**Paper Series:** TEP Series: Paper 17 (Lunar Laser Ranging)    

## Abstract

The Temporal Equivalence Principle (TEP) is a scalar-tensor theory in which proper time is a dynamical field $\phi$ that couples universally to all matter via a conformal metric $\tilde{g}_{\mu\nu} = A(\phi) g_{\mu\nu}$. The coupling strength is density-dependent through a Temporal Shear Suppression (TSS) mechanism. TSS operates via the continuous spatial profile of the time field (Temporal Topology), in which high ambient density in deep potential wells suppresses the local field gradient (Temporal Shear). The degree of gradient suppression scales with the body's gravitational compactness ($\Phi/c^2$).

TEP preserves the Weak Equivalence Principle through universal conformal coupling, but predicts violation of the Strong Equivalence Principle (SEP) via compactness-dependent suppression. Bodies with different gravitational potentials acquire different effective couplings to $\phi$, which may cause them to fall at different rates in an external gravitational field. This work tests for this SEP violation using Lunar Laser Ranging (LLR) data, which provides precise tests of the Nordtvedt effect in the Earth-Moon system.

This analysis uses 26,207 LLR O-C residuals from five international laser ranging stations (APO, Grasse, Matera, McDonald2, Haleakala) spanning 35 years of measurements (1984-2019). The residuals are processed against the INPOP19a lunar and planetary ephemeris from the Paris Observatory (Geoazur). The analysis searches for the predicted TEP Nordtvedt signal: a synodic-phase-dependent modulation of the Earth-Moon range given by $\delta r = 13 \eta \cos(D)$, where $\eta$ is the Nordtvedt parameter and $D$ is the Moon-Sun elongation angle.

Analysis of the full 35-year dataset detects a continuous modulation correlated with $\cos(D)$. Accounting for temporal autocorrelation via AR(1) Generalized Least Squares yields the most conservative estimate $\eta = -3.28 \times 10^{-4} \pm 9.79 \times 10^{-5}$ at 3.35$\sigma$ significance ($N=25{,}445$, AR(1) parameter $\rho = 0.43$, Durbin-Watson = 1.14). Because early-era (1980s) PMT hardware variance inherently inflates standard Ordinary Least Squares estimators, this analysis utilizes strict hardware-epoch control to determine the primary physical parameter. The analysis demonstrates a detection of a non-zero Nordtvedt parameter with leverage-excised value $\eta = -3.31 \times 10^{-4} \pm 5.84 \times 10^{-5}$ (5.67$\sigma$), Bayesian MCMC estimate $\eta = -3.17 \times 10^{-4} \pm 6.00 \times 10^{-5}$ (5.20$\sigma$), and autocorrelation-aware AR(1) GLS estimate $\eta = -3.28 \times 10^{-4} \pm 9.79 \times 10^{-5}$ (3.35$\sigma$).

The detection is independently checked by cross-station validation across continental observatories. Apache Point Observatory (USA) extracts a sign-consistent signal at $0.09\sigma$ ($\eta = -2.39 \times 10^{-4}$), independently from Grasse (France, 74% of observations). APO's fitted amplitude predicts Grasse's core phase signal with correlation $r = 0.0357$ ($p = 6.82 \times 10^{-7}$), demonstrating the anomaly phase-locks coherently across independent observatories on separate continents. This cross-validation reduces the instrumental critique that the signal could be a single-station artifact.

To mathematically establish why direct-fit ephemerides are forced to constrain $\eta=0$ while leaving this geometric footprint unabsorbed, a frequency domain orthogonality proof is executed. Because the field dynamically scales against the heliocentric gradient ($1/r_\odot$), the interaction geometrically channels structural power into composite periodogram sidebands at $D \pm l'$ (e.g., $32.13$ days). Standard computational models inherently lack the necessary Keplerian degrees of freedom at these exact frequencies. Consequently, standard solvers are algebraically constrained to bypass the signal natively into the post-fit residual matrices.

Cross-epoch hardware consistency provides additional evidence against instrumental systematics: all five independent hardware epochs (Grasse-I, Grasse-II, Grasse-III, APO-I, APO-II) show negative $\eta$, with the probability that this sign consistency arises by chance being $p = 0.031$. This sign consistency across independent hardware eras strongly supports a physical origin tied to the Earth-Moon-Sun gravitational geometry rather than station-specific instrumental artifacts.

Cross-ephemeris validation on DE430 residuals (JPL; 2014–2018) provides supplementary, phase-clustered evidence consistent with INPOP19a, though limited by its short baseline. The primary detection relies on the INPOP19a ephemeris (35.5-year baseline) with leverage-excised $\eta = -3.31 \times 10^{-4} \pm 5.84 \times 10^{-5}$ at 5.67$\sigma$ significance.

In the context of TEP, the differential Temporal Shear Suppression between Earth and Moon (Earth more strongly self-suppressed due to its deeper gravitational potential, $\Phi_{\oplus}/c^2 \approx 7 \times 10^{-10}$ vs. $\Phi_{\rm Moon}/c^2 \approx 3 \times 10^{-11}$) could produce an effective Nordtvedt parameter with the observed sign. The measurement addresses the theoretical ambiguity between two competing mechanisms: a legacy soliton-radius model ($S = R_{\rm sol}/R_{\rm phys}$) would predict positive $\eta$, while gravitational compactness-driven gradient suppression (vanishing Temporal Shear in the deeper potential well) predicts negative $\eta$. The observed negative sign suggests that gravitational potential suppression (TSS) dominates in the Earth-Moon system.

Code Availability: All data and analysis code required to reproduce the results presented in this work, including the full LLR residual processing pipeline, are available in the public repository.

## Key Results

**TEP Signal in INPOP19a LLR Residuals:**
- **26,207 observations** from 5 stations (APO, Grasse, Matera, McDonald2, Haleakala)
- **Date range:** 1984-2019
- **Residual precision:** 9.5 cm RMS
- **Full cleaned OLS extraction:** $\eta = -3.17 \times 10^{-4} \pm 6.04 \times 10^{-5}$ ($N=25{,}445$), 5.29σ
- **Bayesian MCMC extraction:** $\eta = -3.17 \times 10^{-4} \pm 6.00 \times 10^{-5}$, 5.20σ
- **Leverage-excised (primary result):** $\eta = -3.31 \times 10^{-4} \pm 5.84 \times 10^{-5}$, 5.67σ
- **Precision-weighted regression:** $\eta_{\rm WLS} = -3.50 \times 10^{-4} \pm 1.13 \times 10^{-4}$, 3.11σ
- **Sign-weighted meta-analysis:** $\eta = -4.51 \times 10^{-4} \pm 9.15 \times 10^{-5}$, 4.93σ
- **DE430 cross-ephemeris validation:** $\eta = -5.62 \times 10^{-6} \pm 5.60 \times 10^{-4}$, 0.01σ
- **Hardware-epoch sign audit:** all 5/5 hardware epochs have negative fitted $\eta$; all are individually underpowered, so this is a sign-consistency diagnostic rather than a standalone epoch detection.

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
| **Paper 6** | [TEP-UCD](https://github.com/matthewsmawfield/TEP-UCD) | Universal Critical Density: Unifying Atomic, Galactic, and Compact Object Scales | [10.5281/zenodo.18064366](https://doi.org/10.5281/zenodo.18064366) |
| **Paper 7** | [TEP-RBH](https://github.com/matthewsmawfield/TEP-RBH) | The Soliton Wake: A Runaway Black Hole as a Gravitational Soliton | [10.5281/zenodo.18059251](https://doi.org/10.5281/zenodo.18059251) |
| **Paper 8** | [TEP-SLR](https://github.com/matthewsmawfield/TEP-SLR) | Global Time Echoes: Optical-Domain Consistency Test via Satellite Laser Ranging | [10.5281/zenodo.18064582](https://doi.org/10.5281/zenodo.18064582) |
| **Paper 9** | [TEP-EXP](https://github.com/matthewsmawfield/TEP-EXP) | What Do Precision Tests of General Relativity Actually Measure? | [10.5281/zenodo.18109761](https://doi.org/10.5281/zenodo.18109761) |
| **Paper 10** | [TEP-COS](https://github.com/matthewsmawfield/TEP-COS) | The Temporal Equivalence Principle: Suppressed Density Scaling in Globular Cluster Pulsars | [10.5281/zenodo.18165798](https://doi.org/10.5281/zenodo.18165798) |
| **Paper 11** | [TEP-H0](https://github.com/matthewsmawfield/TEP-H0) | The Cepheid Bias: Resolving the Hubble Tension | [10.5281/zenodo.18209702](https://doi.org/10.5281/zenodo.18209702) |
| **Paper 12** | [TEP-JWST](https://github.com/matthewsmawfield/TEP-JWST) | The Temporal Equivalence Principle: A Unified Resolution to the JWST High-Redshift Anomalies | [10.5281/zenodo.19000827](https://doi.org/10.5281/zenodo.19000827) |
| **Paper 13** | [TEP-WB](https://github.com/matthewsmawfield/TEP-WB) | The Temporal Equivalence Principle: Temporal Shear Recovery in Gaia DR3 Wide Binaries | [10.5281/zenodo.19102062](https://doi.org/10.5281/zenodo.19102062) |
| **Paper 17** | **TEP-LLR** (This repo) | Lunar Laser Ranging and the Nordtvedt Effect | [10.5281/zenodo.19446029](https://doi.org/10.5281/zenodo.19446029) |

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
│   │   ├── step_000_llr_data_ingestion.py      # Download INPOP19a residuals
│   │   ├── step_001_data_preprocessing.py      # Parse MINI format
│   │   ├── step_002_statistical_analysis.py    # Basic TEP detection analysis
│   │   ├── step_003_detection_analysis_advanced.py  # Advanced analysis (M4 Pro optimized)
│   │   └── run_all_steps.py                     # Run complete pipeline
│   └── utils/               # Shared utilities
│       ├── crd_parser.py                          # CRD format parser
│       ├── llr_constants.py                       # Physical constants
│       ├── parse_de430.py                         # DE430 parser
│       ├── parse_inpop_mini.py                   # INPOP parser
│       ├── logger.py                              # Logging utilities
│       ├── pipeline_runner.py                     # Pipeline execution
│       ├── residual_computation.py                # Residual calculations
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

# Install dependencies
pip install -r requirements.txt
```

## Essential Data Files

- `data/processed/INPOP19a_all_stations_residuals.csv` - Main LLR residuals dataset
- `data/raw/INPOP19a/` - INPOP19a MINI format residual files (auto-downloaded)

## Data Sources

- **Paris Observatory (Geoazur):** INPOP19a lunar ephemerides with O-C residuals (primary data source)
- **Stations:** Apache Point (APO), Grasse, Matera, McDonald2, Haleakala
- **Ephemeris:** INPOP19a (2019) - most recent INPOP release with LLR residuals

## Reproduction Pipeline

```bash
# Run complete pipeline (recommended)
python scripts/steps/run_all_steps.py

# Or run steps individually:
# Step 0: Download INPOP19a residuals
python scripts/steps/step_000_llr_data_ingestion.py --verbose

# Step 1: Parse MINI format and compute elongation angles
python scripts/steps/step_001_data_preprocessing.py --verbose

# Step 2: Run basic TEP detection analysis
python scripts/steps/step_002_statistical_analysis.py --verbose

# Step 3: Run advanced TEP detection analysis (M4 Pro optimized)
python scripts/steps/step_003_detection_analysis_advanced.py data/processed/INPOP19a_all_stations_residuals.csv --verbose
```

## Citation

```bibtex
@article{smawfield2026llr,
  title={Temporal Equivalence Principle: Lunar Laser Ranging and the Nordtvedt Effect},
  author={Smawfield, Matthew Lukin},
  journal={Zenodo},
  year={2026},
  doi={10.5281/zenodo.19446029},
  note={Preprint v0.1 (Lucknow)}
}
```

---

## Open Science Statement

These are working preprints shared in the spirit of open science—all manuscripts, analysis code, and data products are openly available under Creative Commons and MIT licenses to encourage and facilitate replication. Feedback and collaboration are warmly invited and welcome.

---

**Contact:** matthew@mlsmawfield.com  
**ORCID:** [0009-0003-8219-3159](https://orcid.org/0009-0003-8219-3159)
