# TEP-LLR Analysis Scripts

This folder contains the reproducible analysis pipeline for Paper 17 (Lunar Laser Ranging).

## Structure

```
scripts/
├── steps/           # Numbered analysis steps (step_000 through step_043)
│   ├── step_000_llr_data_ingestion.py      # Data acquisition
│   ├── step_001_data_preprocessing.py      # Data processing
│   ├── step_002_de430_preprocessing.py     # DE430 ephemeris
│   ├── step_003_statistical_analysis.py    # Statistical analysis
│   ├── ... (44 steps total)
│   └── run_all_steps.py                     # Run complete pipeline
└── utils/           # Shared utilities
```

## Execution

Steps should be run in numerical order. Each step produces outputs in `results/outputs/`.
