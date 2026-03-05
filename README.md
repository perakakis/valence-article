# Valence Dynamics and Psychological Well-Being

Code repository for the preprint:

**"Direct Valence Measurement Enables More Predictive Affect Dynamics"**

Preprint: [https://doi.org/10.31234/osf.io/gmkc8_v1](https://doi.org/10.31234/osf.io/gmkc8_v1)

## Overview

Valence is central to affective experience yet remains poorly operationalized in emotion science. Most research infers valence by categorizing discrete emotions as positive or negative, which may conflate distinct aspects of emotional experience. Direct bipolar valence measurement — asking individuals how they feel on a continuum from negative to positive — allows participants to integrate contextual complexity into their affective reports and enables detection and quantification of transitions between positive and negative affective states.

This project tests whether metrics quantifying directional transitions predict psychological well-being more effectively than traditional intensity-based measures. Across three ecological momentary assessment (EMA) datasets (N = 345 participants, >30,000 assessments), the **positive-to-negative affect shift ratio (P2N)** — quantifying the propensity to transition from positive to negative affect — consistently outperformed means and standard deviations of positive and negative affect in predicting well-being outcomes.

## Key Findings

- P2N (positive-to-negative shift ratio) outperformed traditional affect mean/SD metrics in predicting well-being
- Results held across LASSO regression, hierarchical regression, and relative importance analyses
- Advantage was robust even with only 3 daily assessments
- Three datasets: German EMA study, Post-COVID EMA study 1, Post-COVID EMA study 2

## Repository Structure

```
code/
  ASR.r             # Compute P2N and N2P affect shift ratios per subject
  MEAN_SD.r         # Compute mean and SD of positive/negative affect per subject
  RMSSD.r           # Compute Root Mean Square of Successive Differences per subject
  preprocessData.r  # Filter subjects based on compliance and variance thresholds
  Figure1.py        # Generate Figure 1 (predictor convergence across methods)

data/
  00 german.*           # German EMA dataset
  01 postcovid1.*       # Post-COVID EMA dataset 1
  02 postcovid2.*       # Post-COVID EMA dataset 2
```

## Metrics

| Metric | Description |
|--------|-------------|
| P2N | Positive-to-negative affect shift ratio |
| N2P | Negative-to-positive affect shift ratio |
| M_PosA / M_NegA | Mean positive / negative affect |
| SD_PosA / SD_NegA | Standard deviation of positive / negative affect |
| RMSSD | Root mean square of successive differences |

## Dependencies

**R:** `data.table`

**Python:** `pandas`, `numpy`, `matplotlib`

## Data

The datasets are included in this repository and also citable via Zenodo:

- German EMA dataset: [https://doi.org/10.5281/zenodo.11060596](https://doi.org/10.5281/zenodo.11060596)
- Post-COVID EMA datasets: [https://doi.org/10.5281/zenodo.15001651](https://doi.org/10.5281/zenodo.15001651)

## Citation

Goicoechea, C., & Perakakis, P. (2026, February 28). [“How Do You Feel?” Direct Valence Measurement Enables Affect Shift Metrics That Outperform Intensity-Based Predictors of Psychological Well-Being.](https://doi.org/10.31234/osf.io/gmkc8_v1) Preprint.
