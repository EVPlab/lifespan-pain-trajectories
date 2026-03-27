# Global Prevalence Trajectories of Pain across the Lifespan

**A Pooled Analysis of 6.1 Million Individuals from 118 Countries and Territories**

This repository contains the analysis code for:

> Global Prevalence Trajectories of Pain across the Lifespan: A Pooled Analysis of 6.1 Million Individuals from 118 Countries and Territories. *In submission*.

---

## Global Pain Benchmarking Tool

An open-access interactive platform for benchmarking external cohort data against the global reference curves derived in this study is available at:

**https://evppainlab.shinyapps.io/global-pain-benchmark/**

The tool enables researchers to upload cohort-level data and compare observed pain prevalence and intensity against age-, sex-, and recall-period-specific global reference trajectories without requiring access to the underlying individual-participant data. See Appendix Section 10 of the manuscript for full documentation of the benchmarking methodology.

---

## Analysis Scripts

Scripts are numbered in the recommended execution order and map to manuscript figures as indicated.

| Script | Description | Figures |
|---|---|---|
| `01_global_pain_trajectories.R` | Lifespan prevalence curves for 11 anatomical sites, any bodily pain, and joint pain, stratified by sex | Figures 2, 3A |
| `02_pain_intensity.R` | Mean pain intensity (Gamma GLMM) and high-intensity pain (NRS ≥ 7) trajectories | Figure 3B, Appendix Fig 18 |
| `03_widespread_pain.R` | Generalised widespread pain (≥ 4/5 body regions, ACR 2016 criteria) | Figure 3C |
| `04_regional_trajectories.R` | Region-specific prevalence curves across eight world regions | Figure 4 |
| `05_hdi_trajectories.R` | Prevalence trajectories by Human Development Index | Figure 5 |
| `06a_smoking_trajectories.R` | Age-varying prevalence by smoking status | Figure 6A, Appendix Fig 14 |
| `06b_bmi_trajectories.R` | Age-varying prevalence by BMI category | Figure 6A, Appendix Fig 15 |
| `06c_income_trajectories.R` | Age-varying prevalence by household income quintile | Figure 6A, Appendix Fig 16 |
| `06d_age_adjusted_rr.R` | Overall age-adjusted risk ratios for BMI, smoking, and income | Figure 6A body maps, Appendix Table 5 |
| `07_population_attributable_fractions.py` | Country- and region-level PAFs combining study-derived RRs with GBD exposure prevalence | Figure 6B–C, Appendix Table 6, Appendix Fig 17 |

---

## Data Availability

Individual-participant data were obtained from 134 study programmes under data-sharing agreements with the respective custodians. Due to the terms of these agreements, the pooled microdata cannot be shared publicly. Investigators wishing to access specific datasets should contact the original study teams listed in Appendix Table 1 of the manuscript.

Model outputs (predicted prevalence trajectories, risk ratios, and deployment reference tables) are available through the [Global Pain Benchmarking Tool](https://evppainlab.shinyapps.io/global-pain-benchmark/).

---

## Contact

Matthew Fillingim — matthew.fillingim@mail.mcgill.ca

For questions about the analysis code, please open an issue in this repository.
