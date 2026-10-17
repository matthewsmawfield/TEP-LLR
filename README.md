# Temporal Equivalence Principle: Lunar Laser Ranging and the Nordtvedt Effect

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19446028.svg)](https://doi.org/10.5281/zenodo.19446028)
[![License: CC BY 4.0](https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)
![TEP-LLR: Lunar Laser Ranging](site/public/image.webp)

**Author:** Matthew Lukin Smawfield  
**Version:** v0.5 (Lucknow)  
**First published:** 17 May 2026 · **Last updated:** 17 October 2026  
**Status:** Preprint  
**DOI:** [10.5281/zenodo.19446028](https://doi.org/10.5281/zenodo.19446028)  
**Website:** [https://mlsmawfield.com/tep/llr/](https://mlsmawfield.com/tep/llr/)  
**Paper Series:** TEP Series · Group C (III. Laser Ranging, Spacecraft Tracking & Solar-System Tests) · Paper 17 (`C17`)

## Abstract

<!-- Generated from site/components/0_abstract.html; update the manuscript source, then run tools/sync_current_corpus_readmes.py. -->

The Temporal Equivalence Principle (TEP) proposes that the local rate of proper time is governed by a dynamical scalar field $\phi=\bar\phi+M_{\rm Pl}u$ rather than treated as a fixed geometric parameter. This field couples universally to all matter via the canonical unified metric $\tilde{g}_{\mu\nu} = -A_{\rm loc}²(u_{\rm grav}(\mathbf x))U_\mu U_\nu+A_{\rm sp}²(u_{\rm drift}(t,\mathbf x))\Pi_{\mu\nu}+B(\phi)\Pi_\mu^\alpha\Pi_\nu^\beta\nabla_\alpha\phi\nabla_\beta\phi$ ($U_\mu\equiv-\nabla_\mu\bar\phi/\sqrt{2\bar X}$, $A_{\rm loc}(u_{\rm grav})=e^{\beta_A u_{\rm grav}}$, $A_{\rm sp}(u_{\rm drift})=e^{\beta_A u_{\rm drift}}=e^{-u_{\rm drift}}$; Paper 0, §2.2, Step 80). Temporal Shear is pinned where the ambient time-gradient is steep — weak-field language: $g \gg g_t \approx cH_0/2$ — by the environmental-state operator $\mathcal S_\Sigma(\mathcal E)$. The canonical pairwise benchmark gives an Earth–Moon response of order $10^{-14}$ and does not account for the reported range residual as a differential-fall signal; the analysis instead targets the synodic term in the post-fit residual channel, $\delta r = 13\eta_{\rm resid}\cos D$, where $D$ is the Moon-Sun elongation angle and $\eta_{\rm resid}$ is the residual-channel amplitude (an Observable Response Coefficient, distinct from the standard metric Nordtvedt parameter $\eta$). Evaluating both the Cassini Shapiro delay and the Earth–Moon laser path on the single solved two-branch solar profile $(k/\sqrt{2})q y²+q² y³=1$ ($k=16.03$, $g_t=cH_0/2$; Step 083) gives $y_{\rm tb}(1.6\,R_\odot)=2.16\times10^{-8}$ ($\gamma_{\rm PPN}-1=-4.32\times10^{-8}$, passing Cassini) and $y_{\rm tb}(1\text{ AU})=1.47\times10^{-5}$, which suppresses the direct scalar optical excursion across the Earth–Moon baseline to $\delta u_\odot^{\rm tb}=7.43\times10^{-16}\cos D$ ($\delta r_{\rm path}^{\rm tb}=-1.43\times10^{-4}$ mm, $\eta_{\rm tb}=-1.10\times10^{-8}$), $3.4\times10⁴$ times below the measured post-fit INPOP19a amplitude $\eta_{\rm resid} = -(3.91\pm0.56)\times10^{-4}$ ($-5.08\pm0.73$ mm) and consistent with direct-fit LLR bounds. Because at a fixed cosmic epoch inside the screened Solar System ($A_{\rm sp}(t_0)=1$) TEP's spatial 3-metric $\tilde\gamma_{ij}=\gamma_{ij}=(1-2\Phi_N/c²)\delta_{ij}$ is identical to GR's, the unscreened solar Newtonian potential differential $\Delta\Phi_{N,\odot}/c²=2.54\times10^{-11}\cos D$ (which defines an unscreened reference ceiling $-\tfrac12 L\,\delta\phi_{\rm eff}\cos D=-4.85$ mm $\cos D$, $\eta_{\rm ceiling}=-3.73\times10^{-4}$ at $T_{\rm Moon}=0.9969$; Paper 0 Step 76, Step 083) cancels in the BCRS$\to$GCRS reduction identically to GR and would give $\gamma_{\rm PPN}-1=-2/3$ at Cassini if attributed to an unscreened solar scalar field; the measured $\eta_{\rm resid}$ is therefore reported as an empirical post-fit residual-channel detection whose dynamical or ephemeris-reduction transfer map awaits a source-level integrator refit.

This work analyses 26,207 raw Lunar Laser Ranging O–C residuals from the public INPOP19a ephemeris archives (Paris Observatory, Geoazur), comprising $N = 25{,}445$ measurements from five international stations (1984–2019) after standard $6\sigma$ MAD outlier cleaning. The primary estimand is a precision-weighted synodic regression on the post-fit residual channel, with inverse station-variance weighting and a full systematic model with annual, monthly, and thermal $\cos(2D)$ regressors and a common intercept.

The primary precision-weighted residual-channel estimate is $\eta_{\rm resid} = -3.91 \times 10^{-4} \pm 5.63 \times 10^{-5}$ (nominal WLS $6.94\sigma$; $6.78\sigma$ cluster-robust; wild cluster bootstrap 95% CI $[-5.37, -2.73]\times 10^{-4}$ across the $G=5$ station clusters). The amplitude stabilises in the band $-3.2$ to $-4.1 \times 10^{-4}$ across the estimator hierarchy, from cosD-only ($5.25\sigma$) through full-systematic OLS ($6.17\sigma$) to precision-weighted ($6.94\sigma$). Robustness checks — common-$\eta_{\rm resid}$ mixed model, leave-one-station-out meta-analysis, wild cluster bootstrap, phase-locked differential, cross-ephemeris validation on cleaned DE430 ($\eta_{\rm resid} = -7.04 \times 10^{-4}$, $5.96\sigma$; confirming the sign but differing in amplitude by $2.77\sigma$ from INPOP19a in the matched window), parametric GR-null bootstrap, and a frequency null scan — all support residual-channel survival of the synodic component (Section 4). If interpreted within the static Parametrized Post-Newtonian framework — an interpretation the residual-channel estimand does not yet license (Section 4.14) — the measured amplitude would correspond to $\beta_{\rm PPN} = 0.999902 \pm 1.07 \times 10^{-5}$ from the joint $(\beta_{\rm PPN}, \gamma)$ contour; this mapping is reported for scale only and is not advanced as a PPN measurement pending the integrator-level refit.

The result is therefore framed as a high-significance residual-channel candidate with a TEP interpretation, not as a completed replacement for direct-fit LLR bounds. The $6.94\sigma$ figure is a detection of a synodic pattern in post-fit residuals; it is not a measurement of the Nordtvedt parameter, and the quantitative relationship between $\eta_{\rm resid}$ and the direct-fit $\eta$ is not yet established. Source-level numerical refits of the INPOP or DE430 integrators with $\eta_{\rm resid}$ left free remain the critical open closure test. Ephemeris-absorption stress tests bound the residual-channel survival amplitude but do not replace a full dynamical integrator-level confirmation.

Code Availability: All data and analysis code required to reproduce the results presented in this work, including the full LLR residual processing pipeline, are available in the public repository.

Keywords: temporal equivalence principle, lunar laser ranging, LLR, equivalence principle, Nordtvedt effect, post-Newtonian, scalar-tensor gravity, strong equivalence principle

Current manuscript source: [0_abstract.html](site/components/0_abstract.html).

## Key Results

The current results and their comparison assumptions are reported in [4. Results](site/components/4_results.html), [6. Conclusion](site/components/6_conclusion.html). The abstract above is synchronized from the current manuscript.

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
| **Paper 16** | [TEP-J0437](https://github.com/matthewsmawfield/TEP-J0437) | Synchronization Holonomy in Pulsar Scintillation | [10.5281/zenodo.19454619](https://doi.org/10.5281/zenodo.19454619) |
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
│   │   ├── ... (91 canonical scripts: step_000 through step_086, including 006b and 046b)
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
  note={Preprint v0.5 (Lucknow)}
}
```

---

## Open Science Statement

These are working preprints shared in the spirit of open science—all manuscripts, analysis code, and data products are openly available under Creative Commons and MIT licenses to encourage and facilitate replication. Feedback and collaboration are warmly invited and welcome.

---

**Contact:** matthew@mlsmawfield.com  
**ORCID:** [0009-0003-8219-3159](https://orcid.org/0009-0003-8219-3159)
