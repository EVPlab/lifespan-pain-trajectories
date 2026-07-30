# Global Prevalence Trajectories of Pain across the Lifespan
### A Pooled Analysis of 6.1 Million Individuals from 118 Countries and Territories

This repository contains the analysis code for:

**Global Prevalence Trajectories of Pain across the Lifespan: A Pooled Analysis of 6.1 Million Individuals from 118 Countries and Territories.

## Global Pain Benchmarking Tool

An open-access interactive platform for benchmarking external cohort data against the global reference curves derived in this study is available at:

https://evppainlab.shinyapps.io/global-pain-benchmark/

The tool enables researchers to upload cohort-level data and compare observed pain prevalence and intensity against age-, sex-, and recall-period-specific global reference trajectories without requiring access to the underlying individual-participant data. See Appendix Section 10 of the manuscript for full documentation of the benchmarking methodology.

## Analysis Scripts

Scripts are in the `scripts/` folder, numbered in the recommended execution order and mapped to manuscript figures as indicated.

| Script | Description | Figures |
|---|---|---|
| `01_global_pain_trajectories.R` | Lifespan prevalence curves for 11 anatomical sites stratified by sex | Figures 2, 3A |
| `02_pain_intensity.R` | Mean pain intensity and high-intensity pain (NRS ≥ 7) trajectories | Figure 3B, Appendix Fig 18 |
| `03_widespread_pain.R` | Generalised widespread pain (≥ 4/5 body regions) | Figure 3C |
| `04_regional_trajectories.R` | Region-specific prevalence curves across eight world regions | Figure 4 |
| `05_hdi_trajectories.R` | Prevalence trajectories by Human Development Index | Figure 5 |
| `06a_smoking_trajectories.R` | Age-varying prevalence by smoking status | Figure 6A, Appendix Fig 14 |
| `06b_bmi_trajectories.R` | Age-varying prevalence by BMI category | Figure 6A, Appendix Fig 15 |
| `06c_income_trajectories.R` | Age-varying prevalence by household income quintile | Figure 6A, Appendix Fig 16 |
| `06d_age_adjusted_rr.R` | Overall age-adjusted risk ratios for BMI, smoking, and income | Figure 6A body maps, Appendix Table 5 |
| `07_population_attributable_fractions.py` | Country- and region-level PAFs | Figure 6B–C, Appendix Table 6, Appendix Fig 17 |

A full prose description of the modelling (the "pseudocode" for the pipeline) is provided in the **Methods → Statistical analysis** section of the manuscript and in the corresponding **Appendix** section.

## System requirements

- **Operating system:** platform-independent. Tested on macOS 14 (Apple Silicon). No non-standard hardware is required; a multi-core CPU speeds up the parallelised scripts but is optional.
- **R:** version 4.4.1. Packages (versions used): `lme4` (1.1-35), `splines` (base), `marginaleffects` (0.20+), `emmeans` (1.10+), `clubSandwich` (0.5+), `data.table` (1.15+), `dplyr` (1.1+), `tidyr` (1.3+), `ggplot2` (3.5+), `foreach` + `doParallel`, `scales`.
- **Python:** version 3.11, for `07_population_attributable_fractions.py`. Packages: `pandas`, `numpy`, `scipy`, `statsmodels`, `matplotlib`, `seaborn`, `geopandas`.

## Installation

```r
install.packages(c("lme4","marginaleffects","emmeans","clubSandwich",
                   "data.table","dplyr","tidyr","ggplot2","foreach",
                   "doParallel","scales"))
```
```bash
pip install pandas numpy scipy statsmodels matplotlib seaborn geopandas
```

## Demo

Because the individual-participant data cannot be shared (see Data Availability), the `demo/` folder provides a **fully synthetic** dataset — randomly generated, containing no real records — so the pipeline can be installed and run end-to-end. Run from the repository root:

```bash
Rscript demo/generate_synthetic_data.R   # writes data/PainData_CS_synthetic.csv
Rscript demo/run_demo.R                   # fits the core GLMM and produces outputs
```

**Expected output** (in `demo/output/`):
- `BackPain_demo_trajectory.csv` — 192 rows (age 5–100 × sex) of predicted prevalence.
- `BackPain_demo_trajectory.png` — a predicted back-pain prevalence curve rising across adulthood, by sex.

**Expected run time:** ~1–3 minutes on a normal desktop. Because the demo data are random, the exact values are not meaningful; the demo confirms the code installs and runs and shows the shape of the expected outputs.

## Running on your own data

The scripts expect a single participant-level CSV (`PainData_CS.csv`) with one row per participant and columns: `Study`, `Sex` (1 = Male, 2 = Female), `Age`, `Global_Region`, `time_window_category` (`past_week`/`past_month`/`past_3_12months`), `case_definition` (`general`/`chronic_explicit`), `HDI`, `Year`, `BMI_Category`, `SmokingStatus`, `Income_class_quintile`, the 0/1 pain-site outcomes (`Pain`, `JointPain`, `Headache`, `FacialPain`, `NeckShoulderPain`, `ChestPain`, `StomachAbdominalPain`, `BackPain`, `HipPain`, `KneePain`, `HandPain`, `ElbowPain`, `FootPain`), and optional per-site overrides `time_window_category_<Site>` / `case_definition_<Site>`. Set the `data_path` variable near the top of each script to your CSV and run `scripts/01_…` through `scripts/07_…` in order.

## Data Availability

Individual-participant data were obtained from 138 study programmes under data-sharing agreements with the respective custodians. Due to the terms of these agreements, the pooled microdata cannot be shared publicly. Investigators wishing to access specific datasets should contact the original study teams listed in the manuscript Appendix.

Model outputs (predicted prevalence trajectories, risk ratios, and deployment reference tables) are available through the Global Pain Benchmarking Tool.

## License

Released under the [MIT License](LICENSE).

## Contact

Matthew Fillingim — matthew.fillingim@mail.mcgill.ca

For questions about the analysis code, please open an issue in this repository.
