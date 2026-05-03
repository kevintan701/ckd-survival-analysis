# CKD Progression & Mortality Risk Factors — NHANES 2017–2023

**[View the Research Report →](https://kevintan701.github.io/ckd-survival-analysis/)**

Survival analysis of **all-cause mortality in U.S. adults with chronic kidney disease (CKD)** using NHANES 2017–2023 linked to National Death Index records. Implemented end-to-end in R: data download, cohort construction, Cox regression, cross-validation, multiple imputation, and a self-contained Quarto report.

---

## Project Structure

```
ckd-survival-analysis/
├── R/
│   ├── 01_data_download.R       # NHANES XPT retrieval via nhanesA; mortality .dat parsing
│   ├── 02_data_cleaning.R       # Cohort construction, CKD-EPI 2021 eGFR, variable harmonization
│   ├── 03_eda.R                 # Descriptive statistics, Table 1, survival curves (gtsummary + survminer)
│   ├── 04_survival_analysis.R   # Cox regression (3 models), PH assumption, model discrimination
│   ├── 05_tidymodels.R          # 10-fold stratified CV via tidymodels + censored
│   └── 06_sensitivity.R         # Multiple imputation (MICE), spline dose-response, subgroup, cause-specific
├── report/
│   └── analysis.qmd             # Quarto report → HTML
├── output/
│   ├── figures/                 # ggplot2 / survminer plots (PNG)
│   └── tables/                  # gtsummary tables (HTML + CSV)
├── data/
│   ├── raw/                     # Downloaded NHANES XPT files (J and L cycles)
│   └── processed/               # Cleaned analysis dataset (.rds)
└── README.md
```

---

## Methods

### Data Source

- **NHANES J cycle** (2017–2018) and **L cycle** (2021–2023) — two completed survey cycles with publicly available XPT files. The K cycle (2019–2020) was suspended mid-collection due to COVID-19 and was never released as standalone files by CDC.
- Mortality follow-up linked via the **NHANES Public Use Linked Mortality Files** (NCHS/NDI), providing time-to-death outcomes.

### Cohort

| Step | N |
|---|---|
| NHANES adults ≥18 | 14,009 |
| Has serum creatinine / eGFR | 10,812 |
| Linked mortality record | 5,124 |
| Complete covariates (BMI, HbA1c) | **5,038** |

Follow-up: median **2.1 years**; **102 deaths** (2.0%).

### Key Variables

| Domain | Variables |
|---|---|
| **Outcome** | All-cause mortality (time-to-event, NHANES linked mortality file) |
| **Renal function** | eGFR (CKD-EPI 2021, race-free), UACR (log-transformed), CKD G-stage |
| **Metabolic** | Diabetes (HbA1c ≥6.5% or diagnosis), hypertension (questionnaire), BMI category |
| **Lifestyle** | Physical activity (MET-min/week from PAQ/PAD modules), smoking |
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
| Missing data | Multiple imputation by chained equations (m=20) | `mice` |
| Sensitivity | Cause-specific Cox, subgroup forest plot | `survival` |

---

## Results

### Cohort & Model Discrimination

| Model | N | Events | C-statistic |
|---|---|---|---|
| Model 1: CKD stage only (unadjusted) | 5,038 | 102 | 0.744 |
| Model 2: + Demographics | 3,790 | 87 | 0.833 |
| Model 3: Full (all covariates) | 3,739 | 78 | 0.847 |
| Model 3: 10-fold CV (unbiased) | — | — | **0.810** |

Proportional hazards assumption passed globally (Schoenfeld residuals, p = 0.92).

### Key Findings — Model 3 Adjusted Hazard Ratios

| Predictor | HR (95% CI) | p |
|---|---|---|
| Age (per 10 years) | **1.73 (1.40–2.15)** | <0.001 |
| Diabetes | **1.73 (1.03–2.90)** | 0.038 |
| Below poverty line | **2.75 (1.18–6.37)** | 0.019 |
| Low income | **2.48 (1.18–5.19)** | 0.016 |
| Overweight BMI | **0.32 (0.17–0.58)** | <0.001 |
| Obese BMI | **0.37 (0.21–0.65)** | <0.001 |
| CKD stage G3a vs G1 | 1.20 (0.53–2.72) | 0.67 |
| CKD stage G3b vs G1 | 1.75 (0.71–4.34) | 0.23 |
| log(UACR) | 1.19 (1.00–1.41) | 0.057 |

**Age, diabetes, and poverty** were the dominant predictors. CKD stage was not independently significant after full adjustment, consistent with the short median follow-up (2.1 years) in this cross-sectional sample linked to mortality — eGFR-based staging requires longer observation to separate mortality gradients. The inverse BMI–mortality association (obesity paradox) is a recognized phenomenon in CKD populations and is discussed in the report.

---

## Interactive Report (GitHub Pages)

The `index.html` at the repo root is a self-contained research report with interactive visualizations. To publish it via GitHub Pages:

1. Push this repository to GitHub
2. Go to **Settings → Pages → Source** and set branch to `main`, folder to `/ (root)`
3. The report will be live at `https://<your-username>.github.io/ckd-survival-analysis/`

All figure paths in `index.html` are relative (`output/figures/*.png`), so they resolve correctly under GitHub Pages with no changes needed.

> **Note:** `data/raw/` contains the original NHANES XPT files and NCHS mortality linkage file (publicly available from CDC). `data/processed/` contains the cleaned analytic dataset. Both are included in the repo for reproducibility.

---

## Reproducing the Analysis

### Requirements

```r
install.packages(c(
  "tidyverse", "here", "haven",
  "nhanesA",
  "survival", "survminer",
  "gtsummary", "gt", "broom.helpers",
  "tidymodels", "censored",
  "mice",
  "rms",
  "quarto", "sessioninfo"
))
```

### Run Order

```bash
Rscript R/01_data_download.R      # ~5 min (downloads ~50 MB from CDC)
Rscript R/02_data_cleaning.R
Rscript R/03_eda.R
Rscript R/04_survival_analysis.R
Rscript R/05_tidymodels.R         # ~3 min (10-fold CV)
Rscript R/06_sensitivity.R        # ~5 min (MICE m=20)
quarto render report/analysis.qmd
```

> **Note:** The K cycle (2019–2020) raw folder will be empty — CDC never released these files. Scripts handle this gracefully and proceed with J and L cycles only.

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
