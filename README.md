# Temporal Equivalence Principle: Lunar Laser Ranging and the Nordtvedt Effect

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19446029.svg)](https://doi.org/10.5281/zenodo.19446029)
[![License: CC BY 4.0](https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)

**Author:** Matthew Lukin Smawfield  
**Version:** v0.1 (Lucknow)  
**First published:** 14 May 2026 · **Last updated:** 14 May 2026  
**Status:** Preprint (Open for Collaboration)  
**DOI:** [10.5281/zenodo.19446029](https://doi.org/10.5281/zenodo.19446029)  
**Website:** [https://mlsmawfield.com/tep/llr/](https://mlsmawfield.com/tep/llr/)  
**Paper Series:** TEP Series: Paper 17 (Lunar Laser Ranging)    

## Abstract

The Temporal Equivalence Principle (TEP) is a scalar-tensor theory in which proper time is a dynamical field $\phi$ that couples universally to all matter via a conformal metric $\tilde{g}_{\mu\nu} = A^2(\phi) g_{\mu\nu}$, with $A(\phi)=\exp(\beta\phi/M_{\rm Pl})$. The coupling strength is density-dependent through a screening of Temporal Shear governed by Temporal Topology, in which high ambient density in deep potential wells suppresses the local field gradient. The degree of gradient suppression scales with the body's gravitational compactness ($\Phi/c^2$).

TEP preserves the Weak Equivalence Principle through universal conformal coupling, but predicts violation of the Strong Equivalence Principle (SEP) via compactness-dependent suppression. Bodies with different gravitational potentials acquire different effective couplings to $\phi$, which may cause them to fall at different rates in an external gravitational field. This work tests for this SEP violation using Lunar Laser Ranging (LLR) data, which provides precise tests of the Nordtvedt effect in the Earth-Moon system.

This analysis uses 26,207 raw LLR O-C residuals from five international laser ranging stations (APO, Grasse, Matera, McDonald2, Haleakala) spanning 35 years of measurements (1984–2019), with 25,445 retained after standard $6\sigma$ MAD outlier cleaning. The residuals are processed against the INPOP19a lunar and planetary ephemeris from the Paris Observatory (Geoazur). The analysis searches for the predicted TEP Nordtvedt signal: a synodic-phase-dependent modulation of the Earth-Moon range given by $\delta r = 13 \eta \cos(D)$, where $\eta$ is the Nordtvedt parameter and $D$ is the Moon-Sun elongation angle.

Analysis of the full 35-year dataset detects a synodic modulation correlated with $\cos(D)$. The headline estimand is a precision-weighted full-systematic regression controlling for annual, monthly, and thermal $\cos(2D)$ aliases: $\eta = -3.91 \times 10^{-4} \pm 5.63 \times 10^{-5}$ at $6.94\sigma$ ($6.78\sigma$ cluster-robust). The unweighted full-systematic OLS sensitivity bound is $\eta = -4.06 \times 10^{-4} \pm 6.58 \times 10^{-5}$ ($6.17\sigma$ / $6.52\sigma$ cluster-robust). A Cook's-Distance leverage diagnostic returns a consistent $\eta = -3.87 \times 10^{-4}$ ($7.82\sigma$ / $8.65\sigma$ cluster-robust), confirming the detection is not driven by high-leverage outliers. The signal strengthens as more systematics are controlled and is stable across mixed-model and station-specific specifications.

Independent confirmation comes from a phase-locked new/full-moon differential ($\eta = -5.95 \times 10^{-4}$, $6.16\sigma$) and a frequency null scan with no significant non-synodic power. Cross-ephemeris validation on DE430 and orthogonality tests support residual-channel survival of the synodic component.

In the TEP framework, differential screening of Temporal Shear between Earth and Moon could produce an effective Nordtvedt parameter with the observed negative sign, consistent with gravitational compactness-driven Temporal Shear suppression dominating in the Earth-Moon system. Source-level INPOP or DE430 integrator refits with $\eta$ free remain the definitive external confirmation test.

Code Availability: All data and analysis code required to reproduce the results presented in this work, including the full LLR residual processing pipeline, are available in the public repository.

## Key Results

**TEP Signal in INPOP19a LLR Residuals:**
- **26,207 observations** from 5 stations (APO, Grasse, Matera, McDonald2, Haleakala); **25,445** after $6\sigma$ MAD cleaning
- **Date range:** 1984–2019
- **Residual precision:** 9.5 cm RMS
- **Cook's-Distance-excised full-systematic OLS (leverage diagnostic):** $\eta = -3.87 \times 10^{-4} \pm 4.95 \times 10^{-5}$ ($N = 23{,}837$ after excision), 7.82σ; cluster-robust 8.65σ
- **Precision-weighted full-systematic (consensus):** $\eta = -3.91 \times 10^{-4} \pm 5.63 \times 10^{-5}$, 6.94σ; cluster-robust 6.78σ
- **Full-systematic OLS without excision (sensitivity upper bound):** $\eta = -4.06 \times 10^{-4} \pm 6.58 \times 10^{-5}$ ($N = 25{,}445$), 6.17σ; cluster-robust 6.52σ
- **Common-$\eta$ mixed model with station systematics (pooling):** $\eta = -4.31 \times 10^{-4} \pm 6.74 \times 10^{-5}$, 6.40σ; $F(4, 25{,}410) = 1.19$, $p = 0.31$
- **Phase-locked new/full-moon differential (robustness):** $\eta = -5.95 \times 10^{-4} \pm 9.66 \times 10^{-5}$, 6.16σ
- **cosD-only OLS (baseline):** $\eta = -3.18 \times 10^{-4} \pm 6.05 \times 10^{-5}$, 5.25σ

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
│   │   ├── ... (74 canonical steps: step_000 through step_073, including 006b and 046b)
│   │   └── run_all_steps.py                     # Run complete pipeline
│   └── utils/               # Shared utilities
│       ├── crd_parser.py                          # CRD format parser
│       ├── llr_constants.py                       # Physical constants
│       ├── parse_de430.py                         # DE430 parser
│       ├── parse_inpop_mini.py                   # INPOP parser
│       ├── logger.py                              # Logging utilities
│       ├── pipeline_runner.py                     # Pipeline execution
│       ├── pipeline_quality_gate.py               # Reviewer-facing audit gate
│       ├── generate_evidence_ledger.py            # Evidence summary artifact
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

# Run the reviewer-facing quality gate
python scripts/utils/pipeline_quality_gate.py

# Generate the evidence ledger directly
python scripts/utils/generate_evidence_ledger.py

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
