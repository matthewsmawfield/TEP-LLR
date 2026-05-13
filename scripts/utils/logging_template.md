# Standardized Logging Template for TEP-LLR Pipeline

## Purpose
This document defines the standard logging format for all TEP-LLR analysis steps to ensure research-grade reproducibility and debuggability.

## Standard Log Format

Each log file should include the following sections in order:

### 1. Step Header
```
[HH:MM:SS] ═══ Starting Step XXX: [Step Name]...
[HH:MM:SS] ═══ STEP PURPOSE: [Brief description of what this step tests and why it matters]
[HH:MM:SS] ═══ METHOD: [Algorithm/approach used]
[HH:MM:SS] ═══ PARAMETERS: [Key parameter values]
```

### 2. Data Summary
```
[HH:MM:SS] ═══ DATA SUMMARY
[HH:MM:SS]     Dataset: N = [number] observations
[HH:MM:SS]     Time range: [start] - [end] (if applicable)
[HH:MM:SS]     Stations: [list of stations] (if applicable)
[HH:MM:SS]     Data quality: [summary of data quality checks]
```

### 3. Analysis Trace
```
[HH:MM:SS] ═══ ANALYSIS TRACE
[HH:MM:SS] >>> [Detailed step-by-step execution]
[HH:MM:SS]     [Key intermediate results]
[HH:MM:SS]     [Diagnostics and checks]
```

### 4. Statistical Diagnostics
For regression/fitting steps:
```
[HH:MM:SS] ═══ STATISTICAL DIAGNOSTICS
[HH:MM:SS]     RSS = [value]
[HH:MM:SS]     MSE = [value]
[HH:MM:SS]     χ²_red = [value]
[HH:MM:SS]     Birge Ratio = [value]
[HH:MM:SS]     Condition Number κ(R) = [value]
[HH:MM:SS]     Convergence: [status] (if applicable)
```

### 5. Results Summary
```
[HH:MM:SS] ═══ RESULTS SUMMARY
[HH:MM:SS]     Primary result: [key finding]
[HH:MM:SS]     Significance: [statistical significance]
[HH:MM:SS]     P-value: [if applicable]
[HH:MM:SS]     Confidence interval: [if applicable]
```

### 6. Interpretation
```
[HH:MM:SS] ═══ INTERPRETATION
[HH:MM:SS]     [What the results mean in context]
[HH:MM:SS]     [Whether results meet expectations]
[HH:MM:SS]     [Any issues or anomalies encountered]
[HH:MM:SS]     [Limitations or assumptions]
```

### 7. Reproducibility Information
```
[HH:MM:SS] ═══ REPRODUCIBILITY
[HH:MM:SS]     Output file: [path to output JSON]
[HH:MM:SS]     Random seed: [if applicable]
[HH:MM:SS]     Software version: [if applicable]
[HH:MM:SS]     Execution time: [duration]
```

### 8. Completion Status
```
[HH:MM:SS] Results saved to [output path]
[HH:MM:SS] ✓   Step Complete. [Status summary]
```

## Guidelines

### Verbosity
- **Minimum:** 15-20 lines for simple steps
- **Recommended:** 20-40 lines for analysis steps
- **Maximum:** No hard limit, but be concise and relevant

### Content Requirements
1. **Context:** Always explain what the step is testing and why
2. **Method:** Describe the algorithm or approach used
3. **Parameters:** Report key parameter values
4. **Diagnostics:** Include relevant statistical diagnostics
5. **Interpretation:** Explain what results mean
6. **Reproducibility:** Include output paths and any random seeds

### Formatting
- Use consistent indentation (4 spaces)
- Use section headers with ═══
- Use >>> for major analysis steps
- Use [CALC] for calculations
- Use [DATA] for data information
- Use ✓ for successful completion
- Use ⚠ for warnings
- Use ❌ for errors

### Examples

#### Good Example (step_016_bayesian_analysis.log)
- Includes detailed MCMC configuration
- Shows convergence diagnostics
- Provides interpretation of Bayes factors
- Reports all key parameters

#### Needs Improvement (step_003_statistical_analysis.log)
- Missing regression diagnostics (RSS, MSE, χ²_red, condition number)
- No method description
- No interpretation of results
- Too brief (6 lines)

## Implementation Checklist

For each step script:
- [ ] Add step purpose description
- [ ] Add method description
- [ ] Add parameter reporting
- [ ] Add statistical diagnostics (if applicable)
- [ ] Add result interpretation
- [ ] Add reproducibility information
- [ ] Ensure minimum 15-20 lines
- [ ] Use consistent formatting

## Priority Steps for Improvement

1. **High Priority** (too brief, missing diagnostics):
   - step_003_statistical_analysis.py
   - step_010_ephemeris_independent_analysis.py
   - step_025_solar_cycle_correlation.py
   - step_031_station_power_analysis.py
   - step_032_hardware_epoch_analysis.py

2. **Medium Priority** (good but could be better):
   - step_003_detection_analysis_advanced.py
   - step_004_temporal_drift_analysis.py
   - step_005_multi_ephemeris_comparison.py
   - step_011_systematic_control_analysis.py

3. **Low Priority** (already good):
   - step_016_bayesian_analysis.py (excellent)
   - step_015_null_tests.py (excellent)
   - step_013_subsample_robustness.py (good)
