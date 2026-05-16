# TEP-LLR Analysis Scripts

This folder contains the reproducible analysis pipeline for Paper 17 (Lunar Laser Ranging).

## Structure

```
scripts/
├── steps/           # Numbered analysis steps (step_000 through step_073, including 006b and 046b)
│   ├── step_000_llr_data_ingestion.py      # Data acquisition
│   ├── step_001_data_preprocessing.py      # Data processing
│   ├── step_002_de430_preprocessing.py     # DE430 ephemeris
│   ├── step_003_statistical_analysis.py    # Statistical analysis
│   ├── ... (74 canonical steps total)
│   └── run_all_steps.py                     # Run complete pipeline
└── utils/           # Shared utilities
```

## Execution

Steps should be run in numerical order. Each step produces outputs in `results/outputs/`.

For a reviewer-facing audit of the evidence spine after running the pipeline:

```bash
python scripts/utils/pipeline_quality_gate.py
```

The quality gate also regenerates `results/outputs/tep_evidence_ledger.md`,
a compact reviewer-facing summary of the positive evidence and bounded risks.
