# Temporal Equivalence Principle: Lunar Laser Ranging and the Nordtvedt Effect

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19446029.svg)](https://doi.org/10.5281/zenodo.19446029)
[![License: CC BY 4.0](https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)

**Author:** Matthew Lukin Smawfield  
**Version:** v0.1 (Lucknow)  
**First published:** 10 May 2026 · **Last updated:** 12 May 2026  
**Status:** Preprint (Open for Collaboration)  
**DOI:** [10.5281/zenodo.19446029](https://doi.org/10.5281/zenodo.19446029)  
**Website:** [https://mlsmawfield.com/tep/llr/](https://mlsmawfield.com/tep/llr/)  
**Paper Series:** TEP Series: Paper 17 (Lunar Laser Ranging)    

## Abstract

The Temporal Equivalence Principle (TEP) is a scalar-tensor theory in which proper time is a dynamical field $\phi$ that couples universally to all matter via a conformal metric $\tilde{g}_{\mu\nu} = A^2(\phi) g_{\mu\nu}$. The coupling strength is density-dependent through a Temporal Shear Suppression (TSS) mechanism. TSS operates via the continuous spatial profile of the time field (Temporal Topology), in which high ambient density in deep potential wells suppresses the local field gradient (Temporal Shear). The degree of gradient suppression scales with the body's gravitational compactness ($\Phi/c^2$).

TEP preserves the Weak Equivalence Principle through universal conformal coupling, but predicts violation of the Strong Equivalence Principle (SEP) via compactness-dependent suppression. Bodies with different gravitational potentials acquire different effective couplings to $\phi$, which may cause them to fall at different rates in an external gravitational field. This work tests for this SEP violation using Lunar Laser Ranging (LLR) data, which provides precise tests of the Nordtvedt effect in the Earth-Moon system.

This analysis uses 26,207 raw LLR O-C residuals from five international laser ranging stations (APO, Grasse, Matera, McDonald2, Haleakala) spanning 35 years of measurements (1984–2019), with 25,445 retained after standard $6\sigma$ MAD outlier cleaning. The residuals are processed against the INPOP19a lunar and planetary ephemeris from the Paris Observatory (Geoazur). The analysis searches for the predicted TEP Nordtvedt signal: a synodic-phase-dependent modulation of the Earth-Moon range given by $\delta r = 13 \eta \cos(D)$, where $\eta$ is the Nordtvedt parameter and $D$ is the Moon-Sun elongation angle.

Analysis of the full 35-year dataset detects a continuous modulation correlated with $\cos(D)$. The primary physical parameter is extracted using a full systematic model that controls for annual, monthly, and thermal $\cos(2D)$ aliases: $\eta = -4.05 \times 10^{-4} \pm 6.57 \times 10^{-5}$ ($6.17\sigma$). The signal strengthens as more systematics are controlled, from $5.25\sigma$ ($\cos D$-only) to $6.17\sigma$ (full model). Cluster-robust standard errors with Cameron-Miller finite-cluster correction across five stations yield $\eta = -4.05 \times 10^{-4}$ at $6.52\sigma$. Step 040 fixes a single canonical estimator hierarchy; nonparametric and precision-weighted checks bound the amplitude without treating every diagnostic as a co-equal discovery.

Cross-station validation demonstrates signal universality when systematic controls are properly applied. A common Nordtvedt parameter with station-specific annual, monthly, and thermal terms yields $\eta = -4.31 \times 10^{-4} \pm 6.74 \times 10^{-5}$ ($6.40\sigma$), with $F(4, 25,410) = 1.19$ ($p = 0.31$) showing no evidence for station-specific deviations. Independent synodic tests reinforce the pooled estimate: a phase-locked new/full-moon differential gives $\eta = -5.74 \times 10^{-4} \pm 9.59 \times 10^{-5}$ ($5.99\sigma$), and a multi-frequency null scan finds no significant power at 55 tested non-synodic factors after correction. A high-precision clean subset (Grasse C-SPAD era plus APO; Step 053) reaches $7.2\sigma$ cluster-robust significance on $N = 12{,}576$.

Joint regressions that add heliocentric radial velocity (Step 047) and CMB dipole orientation (Step 048) expose a hierarchy in which CMB orientation remains highly significant ($\eta_\theta = -9.76 \times 10^{-4}$, $t = -11.03$) while heliocentric distance becomes non-significant ($p = 0.171$). Step 055 falsification tests exclude aliasing, multicollinearity, and permutation artifacts; anti-dipole and joint-model direction rotations anchor the regression structure to the CMB dipole.

To mathematically establish why direct-fit ephemerides are forced to constrain $\eta=0$ while leaving this geometric footprint unabsorbed, a frequency domain orthogonality proof is executed. Because the field dynamically scales against the heliocentric gradient ($1/r_\odot$), the interaction geometrically channels structural power into composite periodogram sidebands at $D \pm l'$ (e.g., $32.13$ days). Standard computational models inherently lack the necessary Keplerian degrees of freedom at these exact frequencies. Consequently, standard solvers are algebraically constrained to bypass the signal natively into the post-fit residual matrices.

Cross-ephemeris validation on DE430 residuals (JPL; 2014–2018) provides supplementary evidence consistent with INPOP19a. On the matched 2014–2018 window, INPOP19a and DE430 both yield negative $\eta$ at $10.0\sigma$ and $5.96\sigma$ (cosD-only), with $\Delta\eta = +3.42 \times 10^{-4}$ ($2.77\sigma$), and the canonical full-systematic model gives $7.72\sigma$ and $5.04\sigma$ on the same span. The primary detection relies on the INPOP19a ephemeris (35.5-year baseline) with full-systematic OLS $\eta = -4.05 \times 10^{-4} \pm 6.57 \times 10^{-5}$ at $6.17\sigma$ significance, and cluster-robust $\eta = -4.05 \times 10^{-4}$ at $6.52\sigma$ with Cameron-Miller finite-cluster correction across 5 stations.

In the context of TEP, the differential Temporal Shear Suppression between Earth and Moon could produce an effective Nordtvedt parameter with the observed negative sign, consistent with gravitational compactness-driven gradient suppression (TSS) dominating in the Earth-Moon system. Orthogonality proofs (Steps 032, 041, 054) and the ephemeris-absorption simulation establish that a synodic $\cos(D)$ coupling of this magnitude cannot be absorbed by static Keplerian fitting, so survival in post-fit residuals is the structurally required channel for this component of TEP.

Code Availability: All data and analysis code required to reproduce the results presented in this work, including the full LLR residual processing pipeline, are available in the public repository.

## Key Results

**TEP Signal in INPOP19a LLR Residuals:**
- **26,207 observations** from 5 stations (APO, Grasse, Matera, McDonald2, Haleakala); **25,445** after $6\sigma$ MAD cleaning
- **Date range:** 1984-2019
- **Residual precision:** 9.5 cm RMS
- **Full-systematic OLS (primary):** $\eta = -4.05 \times 10^{-4} \pm 6.57 \times 10^{-5}$ ($N=25{,}445$), 6.17σ; cluster-robust 6.52σ
- **Common-$\eta$ mixed model (station-specific systematics):** $\eta = -4.31 \times 10^{-4} \pm 6.74 \times 10^{-5}$, 6.40σ
- **Phase-locked new/full-moon differential:** $\eta = -5.74 \times 10^{-4} \pm 9.59 \times 10^{-5}$, 5.99σ
- **Clean subset (Grasse C-SPAD + APO):** $\eta = -3.35 \times 10^{-4}$, 7.25σ cluster-robust on $N=12{,}576$
- **Full-model AR(1) GLS (robustness check):** $\eta = -4.46 \times 10^{-4} \pm 9.57 \times 10^{-5}$, 4.66σ
- **Bayesian MCMC extraction:** $\eta = -2.87 \times 10^{-4} \pm 6.61 \times 10^{-5}$, 4.35σ
- **Cook's Distance leverage-excised (diagnostic):** $\eta = -3.31 \times 10^{-4} \pm 5.84 \times 10^{-5}$, 5.67σ
- **Precision-weighted regression:** $\eta_{\rm WLS} = -3.50 \times 10^{-4} \pm 1.13 \times 10^{-4}$, 3.11σ
- **Matched-window ephemeris (2014–2018, cosD-only):** INPOP19a $\eta = -3.61 \times 10^{-4} \pm 3.59 \times 10^{-5}$, 10.0σ; DE430 $\eta = -7.03 \times 10^{-4} \pm 1.18 \times 10^{-4}$, 5.96σ; $\Delta\eta = +3.42 \times 10^{-4}$ (2.77σ)
- **Matched-window full-systematic:** INPOP19a 7.72σ; DE430 5.04σ; $\Delta\eta = +3.09 \times 10^{-4}$ (2.49σ)

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
| **Paper 9** | [TEP-EXP](https://github.com/matthewsmawfield/TEP-EXP) | What Do Precision Tests of General Relativity Actually Measure? | [10.5281/zenodo.18109760](https://doi.org/10.5281/zenodo.18109760) |
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
│   │   ├── step_000_llr_data_ingestion.py      # Verify raw data availability and hashes
│   │   ├── step_001_data_preprocessing.py      # Parse MINI format
│   │   ├── step_002_de430_preprocessing.py     # DE430 ephemeris processing
│   │   ├── step_003_statistical_analysis.py    # Basic TEP detection analysis
│   │   ├── step_004_detection_analysis_advanced.py  # Advanced analysis (M4 Pro optimized)
│   │   ├── ... (58 canonical steps: step_000 through step_056, including 006b)
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

# Install Python dependencies
pip install -r requirements.txt

# Install static-site build dependencies
cd site
npm ci
cd ..
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

# Rebuild the manuscript from site/components/ (runs node directly; avoids zsh
# `compdef` noise when npm is wired into an interactive shell)
./scripts/build_manuscript.sh

# Equivalent via npm from site/ (may print a harmless zsh completion warning locally)
npm --prefix site run build:markdown

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
