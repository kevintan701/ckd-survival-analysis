# CKD Progression & Mortality Risk Factors — NHANES 2017–2023

End-to-end survival analysis of all-cause mortality in U.S. adults with chronic kidney disease (CKD) using NHANES 2017–2023 linked to the NCHS National Death Index. Covers cohort construction, CKD-EPI 2021 eGFR staging, multivariable Cox regression, cross-validated model discrimination, multiple imputation, and sensitivity analyses — implemented entirely in R.

**[View Research Report →](https://kevintan701.github.io/ckd-survival-analysis/)**

---

## Project Structure

```
ckd-survival-analysis/
├── R/
│   └── ckd_analysis.R           # End-to-end pipeline: download → clean → EDA → Cox → CV → sensitivity
├── output/
│   ├── figures/                 # PNG figures (ggplot2 / survminer)
│   └── tables/                  # Model output and summary tables (CSV + HTML)
├── data/
│   ├── raw/                     # NHANES XPT files (J and L cycles) + NCHS mortality linkage
│   └── processed/               # Cleaned analytic dataset and fitted model objects (.rds / .csv)
├── index.html                   # Self-contained research report (GitHub Pages)
└── README.md
```

---

## Methods

### Data Source

- **NHANES J cycle** (2017–2018) and **L cycle** (2021–2023). The K cycle (2019–2020) was suspended due to COVID-19 and never released by CDC; the pipeline handles this automatically.
- Mortality outcome from the **NCHS Public Use Linked Mortality Files** (NDI linkage through December 2019).
- All files are downloaded at runtime via `nhanesA` (XPT) and direct URL (fixed-width `.dat`).

### Cohort

| Step | N |
|---|---|
| NHANES adults ≥18 | 14,009 |
| Has serum creatinine / eGFR | 10,812 |
| Linked mortality record | 5,124 |
| Complete covariates (BMI, HbA1c) | **5,038** |

Median follow-up: **2.1 years**; **102 deaths** (2.0%).

### Key Variables

| Domain | Variables |
|---|---|
| **Outcome** | All-cause mortality (time-to-event) |
| **Renal function** | eGFR (CKD-EPI 2021, race-free), UACR (log-transformed), CKD G-stage (KDIGO 2024) |
| **Metabolic** | Diabetes (HbA1c ≥6.5% or diagnosis), hypertension (questionnaire), BMI category |
| **Lifestyle** | Physical activity (MET-min/week), smoking status |
| **Demographics** | Age (per 10 years), sex, race/ethnicity, poverty-to-income ratio category |

### Statistical Approach

| Stage | Method | Package |
|---|---|---|
| Descriptive statistics | Table 1 stratified by CKD status | `gtsummary` |
| Survival visualization | Kaplan–Meier curves by CKD stage | `survminer` |
| Primary inference | Multivariable Cox proportional hazards (3 nested models) | `survival` |
| PH assumption | Schoenfeld residuals (global + per-variable) | `survival` |
| Functional form | Martingale residuals, restricted cubic splines | `rms` |
| Model discrimination | Harrell's C-statistic + 10-fold stratified CV | `tidymodels` + `censored` |
| Missing data | Multiple imputation by chained equations (m=20, PMM) | `mice` |
| Sensitivity | Cause-specific Cox, subgroup forest plot | `survival` |

---

## Results

### Model Discrimination

| Model | N | Events | C-statistic |
|---|---|---|---|
| Model 1: CKD stage only | 5,038 | 102 | 0.744 |
| Model 2: + Demographics (age, sex, race/ethnicity) | 3,790 | 87 | 0.833 |
| Model 3: Full model (all covariates) | 3,739 | 78 | 0.847 |
| Model 3: 10-fold cross-validated | — | — | **0.810** |

Proportional hazards assumption passed globally (Schoenfeld residuals, p = 0.920).

### Key Findings — Model 3 Adjusted Hazard Ratios

| Predictor | HR (95% CI) | p |
|---|---|---|
| Age (per 10 years) | **1.73 (1.40–2.15)** | <0.001 |
| Diabetes | **1.73 (1.03–2.90)** | 0.038 |
| Below poverty line | **2.75 (1.18–6.37)** | 0.019 |
| Low income | **2.48 (1.18–5.19)** | 0.016 |
| Overweight BMI | **0.32 (0.17–0.58)** | <0.001 |
| Obese BMI | **0.37 (0.21–0.65)** | <0.001 |
| CKD G3a vs G1 | 1.20 (0.53–2.72) | 0.670 |
| CKD G3b vs G1 | 1.75 (0.71–4.34) | 0.230 |
| log(UACR) | 1.19 (1.00–1.41) | 0.057 |

Age, diabetes, and poverty were the dominant predictors. CKD G-stage was not independently significant after full adjustment, consistent with the short follow-up window — eGFR-based mortality gradients require longer observation to manifest. The inverse BMI–mortality association (obesity paradox) is a well-documented phenomenon in CKD populations.

---

## Reproducing the Analysis

### Requirements

```r
install.packages("pacman")  # handles all remaining dependencies

# Packages used across the pipeline:
# nhanesA, tidyverse, haven, here, janitor, lubridate
# survival, survminer, gtsummary, gt, broom, ggpubr, patchwork, scales
# tidymodels, censored, mice, rms
```

### Run

```bash
Rscript R/ckd_analysis.R   # ~13 min (downloads ~25 MB from CDC, then runs all stages)
```

> The K cycle (2019–2020) is unavailable — CDC never released these files. The pipeline handles this gracefully and proceeds with J and L cycles only.

---

## References

1. Inker LA, et al. New Creatinine- and Cystatin C–Based Equations to Estimate GFR without Race. *N Engl J Med.* 2021;385:1737–1749.
2. KDIGO CKD Work Group. 2024 Clinical Practice Guideline for CKD Evaluation and Management. *Kidney Int.* 2024;105(4S):S117–S314.
3. CDC/NCHS. NHANES Public Use Linked Mortality Files. https://www.cdc.gov/nchs/data-linkage/mortality-public.htm
4. Therneau TM, Grambsch PM. *Modeling Survival Data: Extending the Cox Model.* Springer; 2000.
5. van Buuren S, Groothuis-Oudshoorn K. mice: Multivariate Imputation by Chained Equations in R. *J Stat Softw.* 2011;45(3).
6. Kuhn M, Wickham H. Tidymodels. https://www.tidymodels.org/

---

## Author

**Yuntao (Kevin) Tan** · tyuntao@umich.edu
