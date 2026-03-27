# Global Prevalence Trajectories of Pain across the Lifespan

**A Pooled Analysis of 6.1 Million Individuals from 118 Countries and Territories**

This repository contains the analysis code for:

> [Authors]. Global Prevalence Trajectories of Pain across the Lifespan: A Pooled Analysis of 6.1 Million Individuals from 118 Countries and Territories. *The Lancet*, 2026.

---

## Overview

We harmonised individual-level data from approximately 6.1 million participants across 134 population-based study programmes (894 distinct data sources) spanning 118 countries and territories to construct age-specific reference curves for self-reported pain across eleven anatomical sites, pain intensity, and generalised widespread pain. These trajectories are stratified by sex, world region, national development (HDI), and modifiable risk factors (BMI, smoking, household income).

An interactive benchmarking tool for comparing external cohorts against the global reference curves is available at:
**[Global Pain Benchmarking Tool](https://evppainlab.shinyapps.io/global-pain-benchmark/)**

---

## Repository Structure

```
├── 01_Global_Pain_Trajectories_BySex.R        # Primary lifespan curves (Figures 2 & 3A)
├── 02_Pain_Intensity_And_HighIntensity.R       # Mean intensity & high-intensity pain (Figure 3B, Appendix Fig 18)
├── 03_Generalised_Widespread_Pain.R            # Generalised pain per ACR 2016 criteria (Figure 3C)
├── 04_Regional_Trajectories.R                  # Region-specific prevalence curves (Figure 4)
├── 05_HDI_Development_Trajectories.R           # Trajectories by Human Development Index (Figure 5)
├── 06a_Exposure_Trajectories_BMI.R             # Age-varying prevalence by BMI category (Figure 6A, Appendix Fig 15)
├── 06b_Exposure_Trajectories_Smoking.R         # Age-varying prevalence by smoking status (Figure 6A, Appendix Fig 14)
├── 06c_Exposure_Trajectories_Income.R          # Age-varying prevalence by income quintile (Figure 6A, Appendix Fig 16)
├── 06d_Exposure_AgeAdjusted_RR.R               # Overall age-adjusted risk ratios (Figure 6A body maps, Appendix Table 5)
├── 07_Population_Attributable_Fractions.py     # Country- and region-level PAFs (Figure 6B–C, Appendix Table 6, Appendix Fig 17)
└── README.md
```

## Script Descriptions

### Primary Prevalence Models

**`01_Global_Pain_Trajectories_BySex.R`** — Core analysis script. Fits mixed-effects logistic regression models for 13 pain outcomes (11 anatomical sites plus composite "any bodily pain" and "joint pain") with age modelled via natural cubic splines (df = 3) interacted with sex. Produces sex-stratified and sex-standardised prevalence trajectories, rates of change, risk differences, risk ratios by age band, and deployment reference tables for the benchmarking tool. All inference uses CR2 cluster-robust standard errors.

**`02_Pain_Intensity_And_HighIntensity.R`** — Fits a Gamma(log) GLMM for mean pain intensity (NRS 1–10) among pain sufferers and a binomial(logit) GLMM for high-intensity pain (NRS ≥ 7). Produces sex-stratified trajectories, sex differences (RD and RR) with simulation-based confidence intervals, and deployment packages. Standardised to past-week recall per clinical assessment guidelines.

**`03_Generalised_Widespread_Pain.R`** — Fits a binomial GLMM for generalised widespread pain defined as pain in ≥ 4 of 5 body regions following the 2016 ACR fibromyalgia spatial criteria. Standardised to past-week recall.

### Stratified Analyses

**`04_Regional_Trajectories.R`** — Estimates region-specific prevalence curves for any bodily pain, joint pain, back pain, and headache across eight world regions. Region enters as a fixed effect interacted with age (ns, df = 2), with a study-nested-within-region random intercept. A sex-standardised global reference curve (data-weighted across regions) is included for comparison.

**`05_HDI_Development_Trajectories.R`** — Models the association between country-level Human Development Index and pain prevalence using a flexible surface (ns(Age, df = 2) × ns(HDI, df = 2)). Predicted trajectories are generated at four representative HDI values (0.55, 0.70, 0.80, 0.90). Uses CR2 cluster-robust inference.

### Exposure Analyses

**`06a–c_Exposure_Trajectories_*.R`** — Fit age × exposure interaction models (ns(Age, df = 2) × exposure) for BMI category, smoking status, and household income quintile. Each script produces age-varying prevalence curves per exposure level (standardised to past-month recall and equal sex weights) and age-band risk ratios. Exposure models do not adjust for case definition because the number of chronic-explicit studies with available exposure data was insufficient for stable estimation at many site × exposure combinations.

**`06d_Exposure_AgeAdjusted_RR.R`** — Estimates overall age-adjusted risk ratios (without age × exposure interaction) for all three exposures across 11 pain sites. Uses parametric simulation (multivariate normal draws from the fixed-effects covariance matrix) for confidence intervals.

### Population Attributable Fractions

**`07_Population_Attributable_Fractions.py`** — Combines age- and sex-specific risk ratios from the pooled cohort with country-level exposure prevalence from the Global Burden of Disease study (obesity: GBD 2021; smoking: GBD 2015) and UN population estimates to compute country- and region-level PAFs. Implements the standard binary PAF formula for obesity and smoking, a partial polytomous PAF for income, and the multiplicative independence formula for the combined PAF. Generates choropleth maps (Figure 6B, Appendix Figure 17), regional stacked bar charts (Figure 6C), and Appendix Table 6.

---

## Statistical Methods Summary

All prevalence models use hierarchical generalised linear mixed-effects models (GLMMs) with a binomial distribution and logit link:

```
logit(p) = f(Age) * Exposure + Covariates + (1|Study) + (1|Global_Region)
```

Key specifications:
- **Age**: Natural cubic splines; df = 3 for primary models, df = 2 for exposure and regional interaction models
- **Random effects**: Crossed random intercepts for study programme and world region (primary models); study-nested-within-region for regional analyses
- **Harmonisation covariates**: Recall period (past week / past month / past 3+ months) and case definition (general / chronic-explicit) as fixed effects in primary and regional/HDI models
- **Standardisation**: Predictions standardised to past-month recall and general case definition (past-week for intensity and generalised pain models)
- **Inference**: CR2 cluster-robust variance estimation (primary models) or parametric simulation from the fixed-effects covariance matrix (intensity and age-adjusted RR models)
- **Software**: R 4.4.1 with lme4, marginaleffects, clubSandwich; Python 3.x with pandas, geopandas, matplotlib

---

## Data Availability

Individual-participant data were obtained from 134 study programmes under data-sharing agreements with the respective custodians. Due to the terms of these agreements, the pooled microdata cannot be shared publicly. Investigators wishing to access specific datasets should contact the original study teams listed in Appendix Table 1 of the manuscript.

The following are publicly available:
- GBD exposure prevalence estimates: [Global Health Data Exchange](https://ghdx.healthdata.org/)
- UN Population estimates: [UN Population Division](https://population.un.org/wpp/)
- Human Development Index: [UNDP Human Development Reports](https://hdr.undp.org/data-center)

Model outputs (predicted prevalence trajectories, risk ratios, and deployment reference tables) are provided in the results directories generated by these scripts and are available through the [Global Pain Benchmarking Tool](https://evppainlab.shinyapps.io/global-pain-benchmark/).

---

## Reproducing the Analyses

### Prerequisites

**R packages:**
```r
install.packages(c(
  "lme4", "glmmTMB", "marginaleffects", "clubSandwich",
  "splines", "data.table", "dplyr", "tidyr",
  "ggplot2", "ggthemes", "scales", "cowplot", "patchwork",
  "forcats", "foreach", "doParallel", "MASS", "RColorBrewer"
))
```

**Python packages:**
```bash
pip install pandas numpy geopandas matplotlib seaborn
```

### Configuration

All scripts use hardcoded file paths that must be updated to your local environment before running. Path variables are defined at the top of each script. Update the following:

- `data_path` — path to the harmonised individual-participant data CSV
- `csv_path` / `results_output_path` — output directory for model results
- `save_path` — output directory for figures

### Execution Order

Scripts are numbered to indicate the recommended execution order. Scripts 01–03 (primary models) should be run first as their outputs feed into Script 07 (PAFs). Scripts 04–06 are independent of each other and can be run in parallel. Script 07 requires outputs from both 01 (baseline prevalence) and 06d (risk ratios).

```
01 → 02 → 03                    (primary models)
04, 05, 06a, 06b, 06c, 06d      (stratified analyses, independent)
07                               (PAFs; requires outputs from 01 and 06d)
```

---

## Citation

If you use code or outputs from this repository, please cite:

> [Authors]. Global Prevalence Trajectories of Pain across the Lifespan: A Pooled Analysis of 6.1 Million Individuals from 118 Countries and Territories. *The Lancet*, 2026. DOI: [to be added upon publication]

---

## License

This code is released under the [MIT License](LICENSE).

---

## Contact

For questions about the analysis code, please open an issue in this repository or contact [corresponding author email].
