#!/usr/bin/env python3
"""
Renumber Section 4 subsections sequentially and update cross-references
across all site/components/*.html files.
"""

import re
from pathlib import Path

COMPONENTS_DIR = Path("/Users/matthewsmawfield/www/Temporal Equivalence Principle/TEP-LLR/site/components")

# =============================================================================
# 1. Exact header replacements in 4_results.html
#    (old substring -> new substring)
# =============================================================================
HEADER_REPLACEMENTS = [
    # h3 headers
    ("<h3>4.0 Independent Evidence Spine (Steps 040, 045, 048, 055)</h3>", "<h3>4.0 Independent Evidence Spine</h3>"),
    ("<h3>4.27 Multiple-Testing Consolidation (Step 042)</h3>", "<h3>4.1 Multiple-Testing Consolidation</h3>"),
    ("<h3>4.1 Cross-Station Validation: Defense Against Single-Station Dominance</h3>", "<h3>4.2 Cross-Station Validation: Defense Against Single-Station Dominance</h3>"),
    ("<h3>4.26 Station-Balanced TEP Stress Test (Step 046)</h3>", "<h3>4.3 Station-Balanced TEP Stress Test</h3>"),
    ("<h3>4.27 Grasse-Specific Systematic Sufficiency (Step 059)</h3>", "<h3>4.4 Grasse-Specific Systematic Sufficiency</h3>"),
    ("<h3>4.28 Gaussian Process Non-Parametric Extraction (Step 060)</h3>", "<h3>4.5 Gaussian Process Non-Parametric Extraction</h3>"),
    ("<h3>4.29 Systematic Amplitude Sensitivity (Step 061)</h3>", "<h3>4.6 Systematic Amplitude Sensitivity</h3>"),
    ("<h3>4.2 Correlation Analysis</h3>", "<h3>4.7 Correlation Analysis</h3>"),
    ("<h3>4.3 Kinematic Signal Extraction and Robust Regression</h3>", "<h3>4.8 Kinematic Signal Extraction and Robust Regression</h3>"),
    ("<h3>4.4 Differential Analysis</h3>", "<h3>4.9 Differential Analysis</h3>"),
    ("<h3>4.5 Station-by-Station Consistency</h3>", "<h3>4.10 Station-by-Station Consistency</h3>"),
    ("<h3>4.6 Hardware Epoch Analysis (Step 030)</h3>", "<h3>4.11 Hardware Epoch Analysis</h3>"),
    ("<h3>4.7 Cross-Ephemeris Validation: INPOP19a Robustness and DE430 Fragility</h3>", "<h3>4.12 Cross-Ephemeris Validation: INPOP19a Robustness and DE430 Fragility</h3>"),
    ("<h3>4.28 Linearized Post-Fit η Extraction (Step 056)</h3>", "<h3>4.13 Linearized Post-Fit η Extraction</h3>"),
    ("<h3>4.8 Physical Interpretation</h3>", "<h3>4.14 Physical Interpretation</h3>"),
    ("<h3>4.9 Robustness and Systematic Error Analysis</h3>", "<h3>4.15 Robustness and Systematic Error Analysis</h3>"),
    ("<h3>4.10 Power and Sensitivity Analysis</h3>", "<h3>4.16 Power and Sensitivity Analysis</h3>"),
    ("<h3>4.11 Temporal Evolution and Synodic Phase Coherence</h3>", "<h3>4.17 Temporal Evolution and Synodic Phase Coherence</h3>"),
    ("<h3>4.12 Limitations of the Analysis</h3>", "<h3>4.18 Limitations of the Analysis</h3>"),
    ("<h3>4.13 Extended Systematic Analysis (Steps 010–013)</h3>", "<h3>4.19 Extended Systematic Analysis</h3>"),
    ("<h3>4.14 Bayesian Inference and Evidence</h3>", "<h3>4.20 Bayesian Inference and Evidence</h3>"),
    ("<h3>4.15 Spectral Analysis and Frequency Specificity</h3>", "<h3>4.21 Spectral Analysis and Frequency Specificity</h3>"),
    ("<h3>4.16 Full Hard Audit Verification</h3>", "<h3>4.22 Full Hard Audit Verification</h3>"),
    ("<h3>4.17 Demonstration of Ephemeris Absorption Masking (Steps 023 & 036)</h3>", "<h3>4.23 Demonstration of Ephemeris Absorption Masking</h3>"),
    ("<h3>4.18 Differential Suppression (Environmental Amplitude Scaling) (Step 022)</h3>", "<h3>4.24 Differential Suppression (Environmental Amplitude Scaling)</h3>"),
    ("<h3>4.19 Quantitative η Prediction (Step 033)</h3>", "<h3>4.25 Quantitative η Prediction</h3>"),
    ("<h3>4.20 Decoupling Thermal Array Deformation (Step 024)</h3>", "<h3>4.26 Decoupling Thermal Array Deformation</h3>"),
    ("<h3>4.21 Leverage Temporal Clustering (Step 025)</h3>", "<h3>4.27 Leverage Temporal Clustering</h3>"),
    ("<h3>4.22 False-Positive Diagnostic Results (Steps 029–030)</h3>", "<h3>4.28 False-Positive Diagnostic Results</h3>"),
    ("<h3>4.23 Frequency-Specific Null Testing (Step 015)</h3>", "<h3>4.29 Frequency-Specific Null Testing</h3>"),
    ("<h3>4.24 Cross-Validation, Station Distribution, and Covariate Shift (Steps 051–052)</h3>", "<h3>4.30 Cross-Validation, Station Distribution, and Covariate Shift</h3>"),
    ("<h3>4.25 Clean-Subset High-SNR Analysis and Orbital Orthogonality (Steps 053–054)</h3>", "<h3>4.31 Clean-Subset High-SNR Analysis and Orbital Orthogonality</h3>"),

    # h4 headers
    ("<h4>Primary Robust Estimand: Cook's-Distance-Excised Full Systematic Model</h4>", "<h4>4.8.1 Primary Robust Estimand: Cook's-Distance-Excised Full Systematic Model</h4>"),
    ("<h4>Robustness Checks and Diagnostic Tests</h4>", "<h4>4.8.2 Robustness Checks and Diagnostic Tests</h4>"),
    ("<h4>Haleakala Audit</h4>", "<h4>4.10.1 Haleakala Audit</h4>"),
    ("<h4>Full-Systematic Robust Estimator Consensus: Resolving the Factor-of-Two Concern</h4>", "<h4>4.15.1 Full-Systematic Robust Estimator Consensus: Resolving the Factor-of-Two Concern</h4>"),
    ("<h4>Systematic Projection Analysis (Step 044)</h4>", "<h4>4.15.2 Systematic Projection Analysis</h4>"),
    ("<h4>Phase-Locked Differential Analysis</h4>", "<h4>4.15.3 Phase-Locked Differential Analysis</h4>"),
    ("<h4>4.13.1 Systematic Control Analysis (Step 010)</h4>", "<h4>4.19.1 Systematic Control Analysis</h4>"),
    ("<h4>4.13.2 Noise Injection and Signal Recovery (Step 011)</h4>", "<h4>4.19.2 Noise Injection and Signal Recovery</h4>"),
    ("<h4>4.13.3 Subsample Robustness (Step 012)</h4>", "<h4>4.19.3 Subsample Robustness</h4>"),
    ("<h4>4.13.4 Station Decomposition (Step 013)</h4>", "<h4>4.19.4 Station Decomposition</h4>"),
    ("<h4>4.18.1 Orbital Velocity Modulation of Temporal Shear (Step 047)</h4>", "<h4>4.24.1 Orbital Velocity Modulation of Temporal Shear</h4>"),
    ("<h4>4.18.2 CMB Dipole Anisotropy (Step 048)</h4>", "<h4>4.24.2 CMB Dipole Anisotropy</h4>"),
    ("<h4>4.18.3 Canonical Full-Systematic Extraction (Step 050)</h4>", "<h4>4.24.3 Canonical Full-Systematic Extraction</h4>"),
    ("<h4>4.18.4 Per-Station Power and Observed Extraction (Step 029)</h4>", "<h4>4.24.4 Per-Station Power and Observed Extraction</h4>"),
    ("<h4>4.18.5 Precision-Weighted Station Regression (Step 029)</h4>", "<h4>4.24.5 Precision-Weighted Station Regression</h4>"),
    ("<h4>4.18.6 Frequency Domain Orthogonality and Sideband Harmonics (Step 032)</h4>", "<h4>4.24.6 Frequency Domain Orthogonality and Sideband Harmonics</h4>"),
    ("<h4>4.22.1 Day/Night Thermal Bias Null Test (Step 027)</h4>", "<h4>4.28.1 Day/Night Thermal Bias Null Test</h4>"),
    ("<h4>4.22.2 True Geometric Elongation Null Test (Step 030)</h4>", "<h4>4.28.2 True Geometric Elongation Null Test</h4>"),
]

# =============================================================================
# 2. General text replacements (section references, table labels, etc.)
#    Using regex with negative lookbehind/lookahead to avoid partial matches.
# =============================================================================
# Each tuple: (old_number, new_number)
# Order matters: longer old strings first to avoid shadowing.
SECTION_NUMBER_MAP = [
    # sub-sub-sections
    ("4.18.6", "4.24.6"),
    ("4.18.5", "4.24.5"),
    ("4.18.4", "4.24.4"),
    ("4.18.3", "4.24.3"),
    ("4.18.2", "4.24.2"),
    ("4.18.1", "4.24.1"),
    ("4.13.4", "4.19.4"),
    ("4.13.3", "4.19.3"),
    ("4.13.2", "4.19.2"),
    ("4.13.1", "4.19.1"),
    ("4.22.2", "4.28.2"),
    ("4.22.1", "4.28.1"),
    # sub-sections
    ("4.29", "4.6"),
    ("4.28", "4.13"),   # Linearized Post-Fit η Extraction (Step 056) — the *second* 4.28 in old doc
    ("4.27", "4.1"),    # Multiple-Testing Consolidation (Step 042) — the *first* 4.27 in old doc
    ("4.26", "4.3"),
    ("4.25", "4.31"),
    ("4.24", "4.30"),
    ("4.23", "4.29"),
    ("4.22", "4.28"),
    ("4.21", "4.27"),
    ("4.20", "4.26"),
    ("4.19", "4.25"),
    ("4.18", "4.24"),
    ("4.17", "4.23"),
    ("4.16", "4.22"),
    ("4.15", "4.21"),
    ("4.14", "4.20"),
    ("4.13", "4.19"),
    ("4.12", "4.18"),
    ("4.11", "4.17"),
    ("4.10", "4.16"),
    ("4.9", "4.15"),
    ("4.8", "4.14"),
    ("4.7", "4.12"),
    ("4.6", "4.11"),
    ("4.5", "4.10"),
    ("4.4", "4.9"),
    ("4.3", "4.8"),
    ("4.2", "4.7"),
    ("4.1", "4.2"),
]

# Additional exact replacements for table cells / ranges / captions
EXACT_TEXT_REPLACEMENTS = [
    # Table caption
    ("Table 4.27a: Required Grasse-only systematic vs. known",
     "Table 4.4a: Required Grasse-only systematic vs. known"),
    ("projection from Step 044 (Table 4.27a).",
     "projection from Step 044 (Table 4.4a)."),

    # Read-order table cells & caption
    ("Recommended read order for §4 (section numbers follow pipeline\n            grouping, not print order).",
     "Recommended read order for §4."),
    ("<td>4.27</td>\n                <td>Multiplicity control (Step 042)</td>",
     "<td>4.1</td>\n                <td>Multiplicity control (Step 042)</td>"),
    ("<td>4.1</td>\n                <td>Cross-station pooling; Grasse concentration</td>",
     "<td>4.2</td>\n                <td>Cross-station pooling; Grasse concentration</td>"),
    ("<td>4.26</td>\n                <td>Station-balance power stress (Step 046)</td>",
     "<td>4.3</td>\n                <td>Station-balance power stress (Step 046)</td>"),
    ("<td>4.2–4.5</td>\n                <td>Core regression and per-station detail</td>",
     "<td>4.7–4.10</td>\n                <td>Core regression and per-station detail</td>"),
    ("<td>4.6</td>\n                <td>Hardware epochs</td>",
     "<td>4.11</td>\n                <td>Hardware epochs</td>"),
    ("<td>4.7, 4.28</td>\n                <td>DE430 comparison; linearized post-fit $\\eta$ extraction (Step 056)</td>",
     "<td>4.12, 4.13</td>\n                <td>DE430 comparison; linearized post-fit $\\eta$ extraction (Step 056)</td>"),
    ("<td>4.8–4.25</td>\n                <td>Systematics, spectral specificity, orthogonality simulations</td>",
     "<td>4.14–4.31</td>\n                <td>Systematics, spectral specificity, orthogonality simulations</td>"),
    ("<td>4.18.*</td>\n                <td>Conditional CMB / heliocentric joint layer</td>",
     "<td>4.24.*</td>\n                <td>Conditional CMB / heliocentric joint layer</td>"),

    # Intro paragraph rewrite
    ("""        A correlation analysis between the LLR O-C residuals and the predicted
        TEP Nordtvedt signal modulation $\\cos(D)$, where $D$ is the Moon-Sun
        elongation angle. Section 4.0 states the independent evidence spine;
        Sections 4.1–4.7 develop pooling, regression, and cross-ephemeris
        validation; later sections treat systematics, spectral specificity,
        environmental modulation, orthogonality, and falsification without
        reopening the headline synodic estimand. Section numbers follow
        pipeline step grouping; the read-order table below is the referee
        path through §4.""",
     """        A correlation analysis between the LLR O-C residuals and the predicted
        TEP Nordtvedt signal modulation $\\cos(D)$, where $D$ is the Moon-Sun
        elongation angle. Section 4.0 states the independent evidence spine;
        the subsections that follow develop pooling, regression, and
        cross-ephemeris validation; later sections treat systematics, spectral
        specificity, environmental modulation, orthogonality, and falsification
        without reopening the headline synodic estimand. The table below gives
        the referee path through §4."""),
]


def apply_regex_replacements(text):
    """Replace standalone section numbers using regex to avoid partial matches."""
    for old, new in SECTION_NUMBER_MAP:
        # Pattern: not preceded by a digit, and not followed by a digit
        pattern = rf"(?<![0-9]){re.escape(old)}(?![0-9])"
        text = re.sub(pattern, new, text)
    return text


def main():
    html_files = sorted(COMPONENTS_DIR.glob("*.html"))
    changed_summary = []

    for fpath in html_files:
        original = fpath.read_text(encoding="utf-8")
        text = original

        # --- Exact header replacements (only relevant for 4_results.html) ---
        for old, new in HEADER_REPLACEMENTS:
            text = text.replace(old, new)

        # --- Exact text replacements (table cells, captions, etc.) ---
        for old, new in EXACT_TEXT_REPLACEMENTS:
            text = text.replace(old, new)

        # --- Regex-based section-number replacements ---
        text = apply_regex_replacements(text)

        if text != original:
            fpath.write_text(text, encoding="utf-8")
            # Count changes roughly
            changed_summary.append(str(fpath.name))

    print("Updated files:", ", ".join(changed_summary))


if __name__ == "__main__":
    main()
