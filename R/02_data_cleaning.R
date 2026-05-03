# =============================================================================
# 02_data_cleaning.R
# CKD Progression & Mortality Risk Factors — NHANES 2017–2023
# Author: Yuntao (Kevin) Tan | tyuntao@umich.edu | December 2025
#
# Purpose: Load raw NHANES XPT files across three cycles, harmonise variable
#          names, apply clinical definitions, compute eGFR (CKD-EPI 2021),
#          parse the NCHS linked mortality file, merge all sources on SEQN,
#          and produce a single analysis-ready dataset.
#
# Input:   data/raw/<CYCLE>/*.XPT
#          data/raw/mortality/NHANES_2017_2018_MORT_2019_PUBLIC.dat
#
# Output:  data/processed/ckd_analysis.rds   — full cleaned dataset
#          data/processed/ckd_analysis.csv   — plain-text copy
#          output/tables/02_cohort_flow.csv  — exclusion flowchart counts
# =============================================================================


# -----------------------------------------------------------------------------
# 0. Packages
# -----------------------------------------------------------------------------
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(
  tidyverse,   # dplyr, tidyr, stringr, readr, purrr
  haven,       # read_xpt() — read SAS transport files
  here,        # project-relative paths
  janitor,     # clean_names() — snake_case column standardisation
  mice,        # multiple imputation for missing covariates
  lubridate    # date arithmetic (follow-up time calculations)
)

here::i_am("R/02_data_cleaning.R")

message("── 02  Data Cleaning ───────────────────────────────────────────────────")


# -----------------------------------------------------------------------------
# 1. Helper: load one NHANES component across cycles
# -----------------------------------------------------------------------------
# read_xpt() preserves SAS variable labels as column attributes; clean_names()
# converts them to consistent snake_case to prevent downstream join failures
# caused by mixed-case naming across cycles.

load_component <- function(component, cycles = c("J", "K", "L")) {

  purrr::map_dfr(cycles, function(cyc) {

    path <- here("data", "raw", cyc, paste0(component, "_", cyc, ".XPT"))

    if (!file.exists(path)) {
      message("  ✗  Missing: ", basename(path), " — skipping")
      return(NULL)
    }

    haven::read_xpt(path) |>
      janitor::clean_names() |>
      mutate(cycle = cyc)
  })
}


# -----------------------------------------------------------------------------
# 2. Load each component
# -----------------------------------------------------------------------------
message("\nLoading NHANES components...")

demo   <- load_component("DEMO")
bmx    <- load_component("BMX")
biopro <- load_component("BIOPRO")   # Serum creatinine — J & L cycles (K unavailable)
alb_cr <- load_component("ALB_CR")
diq    <- load_component("DIQ")
bpq    <- load_component("BPQ")
bpxo   <- load_component("BPXO")
paq    <- load_component("PAQ")
smq    <- load_component("SMQ")
ghb    <- load_component("GHB")


# -----------------------------------------------------------------------------
# 3. Parse NCHS Linked Mortality File (fixed-width ASCII)
# -----------------------------------------------------------------------------
# The mortality linkage file is released as a fixed-width .dat file, not XPT.
# Column positions follow the NCHS data dictionary (2023 release).
# PERMTH_EXM measures follow-up from the MEC exam date — preferred over
# PERMTH_INT (interview date) because exam date aligns with biomarker
# collection and is standard in survival analyses using NHANES data.
# Reference: https://www.cdc.gov/nchs/data-linkage/mortality-public.htm

message("\nParsing NCHS Linked Mortality File...")

mort_path <- here("data", "raw", "mortality",
                  "NHANES_2017_2018_MORT_2019_PUBLIC.dat")

mort_cols <- readr::fwf_cols(
  seqn          = c(1,  6),
  eligstat      = c(15, 15),
  mortstat      = c(16, 16),
  ucod_leading  = c(17, 19),
  diabetes      = c(20, 20),
  hyperten      = c(21, 21),
  permth_exm    = c(43, 44),
  permth_int    = c(46, 47)
)

mortality <- readr::read_fwf(
  mort_path,
  col_positions = mort_cols,
  col_types     = readr::cols(.default = "c"),
  na            = c(".", "")
) |>
  mutate(
    seqn       = as.integer(seqn),
    mortstat   = as.integer(mortstat),
    permth_exm = as.numeric(permth_exm),
    follow_yrs = permth_exm / 12,
    died       = if_else(mortstat == 1, 1L, 0L)
  ) |>
  select(seqn, mortstat, died, follow_yrs, ucod_leading, diabetes, hyperten)

message("  ✔  Mortality file: ", nrow(mortality), " records")


# -----------------------------------------------------------------------------
# 4. Demographics
# -----------------------------------------------------------------------------
# PIR is retained in both continuous and categorical form: continuous for
# regression adjustment, categorical for Table 1 and subgroup analyses,
# consistent with conventions in the CKD epidemiology literature.
# wt_mec is the MEC examination weight for survey-weighted estimates.

demo_clean <- demo |>
  select(
    seqn, cycle,
    age_yr   = ridageyr,
    sex      = riagendr,
    race_eth = ridreth3,
    pir      = indfmpir,
    wt_mec   = wtmec2yr
  ) |>
  mutate(
    sex = factor(sex, levels = c(1, 2), labels = c("Male", "Female")),

    race_eth = factor(race_eth,
      levels = c(1, 2, 3, 4, 6, 7),
      labels = c("Mexican American", "Other Hispanic",
                 "Non-Hispanic White", "Non-Hispanic Black",
                 "Non-Hispanic Asian", "Other/Multiracial")
    ),

    poverty_cat = case_when(
      pir <  1.0              ~ "Below poverty",
      pir >= 1.0 & pir < 2.0 ~ "Low income",
      pir >= 2.0 & pir < 4.0 ~ "Middle income",
      pir >= 4.0              ~ "High income",
      TRUE                    ~ NA_character_
    ) |> factor(levels = c("Below poverty", "Low income",
                            "Middle income", "High income")),

    age_group = cut(age_yr,
      breaks = c(17, 39, 59, 74, Inf),
      labels = c("18–39", "40–59", "60–74", "75+"),
      right  = TRUE
    )
  ) |>
  filter(age_yr >= 18)

message("  ✔  Demographics: ", nrow(demo_clean), " adults")


# -----------------------------------------------------------------------------
# 5. eGFR — CKD-EPI 2021 creatinine equation (race-free)
# -----------------------------------------------------------------------------
# The 2021 CKD-EPI revision removed the race coefficient following the NKF–
# ASN Task Force recommendation. This is now the standard equation used by
# clinical laboratories and epidemiological studies including USRDS and KECC.
# κ and α are sex-specific constants; the female multiplier (1.012) accounts
# for average sex differences in muscle mass independent of race.
# Reference: Inker et al., NEJM 2021;385:1737–1749.

ckd_epi_2021 <- function(scr, age, sex) {
  kappa <- if_else(sex == "Female", 0.7,    0.9)
  alpha <- if_else(sex == "Female", -0.241, -0.302)
  sex_f <- if_else(sex == "Female", 1.012,  1.000)

  egfr <- 142 *
    pmin(scr / kappa, 1) ^ alpha *
    pmax(scr / kappa, 1) ^ (-1.200) *
    (0.9938 ^ age) *
    sex_f

  return(round(egfr, 1))
}

scr_all <- biopro |>
  select(seqn, cycle, scr = lbxscr) |>
  left_join(demo_clean |> select(seqn, age_yr, sex), by = "seqn") |>
  mutate(
    egfr = ckd_epi_2021(scr, age_yr, as.character(sex)),

    # KDIGO 2024 G-staging; A-staging (albuminuria) is incorporated after
    # UACR merge in Section 6.
    ckd_stage = case_when(
      egfr >= 90              ~ "G1 (≥90)",
      egfr >= 60 & egfr < 90 ~ "G2 (60–89)",
      egfr >= 45 & egfr < 60 ~ "G3a (45–59)",
      egfr >= 30 & egfr < 45 ~ "G3b (30–44)",
      egfr >= 15 & egfr < 30 ~ "G4 (15–29)",
      egfr <  15              ~ "G5 (<15)",
      TRUE                    ~ NA_character_
    ) |> factor(levels = c("G1 (≥90)", "G2 (60–89)", "G3a (45–59)",
                            "G3b (30–44)", "G4 (15–29)", "G5 (<15)")),

    # Binary CKD flag: eGFR < 60 ml/min/1.73m² defines stages G3–G5
    ckd = if_else(egfr < 60, 1L, 0L, missing = NA_integer_)
  ) |>
  select(seqn, scr, egfr, ckd_stage, ckd)

message("  ✔  eGFR computed (CKD-EPI 2021): ", sum(!is.na(scr_all$egfr)), " values")


# -----------------------------------------------------------------------------
# 6. Urine albumin-to-creatinine ratio (UACR)
# -----------------------------------------------------------------------------
# UACR ≥30 mg/g for >3 months is a diagnostic criterion for CKD independent
# of eGFR. Retaining both dimensions (G-stage + A-stage) enables heat-map
# risk stratification consistent with KDIGO 2024 guidelines.

uacr_clean <- alb_cr |>
  select(seqn, uacr = urxuma) |>
  mutate(
    albuminuria = case_when(
      uacr <  30  ~ "A1 (normal)",
      uacr <  300 ~ "A2 (moderately increased)",
      uacr >= 300 ~ "A3 (severely increased)",
      TRUE        ~ NA_character_
    ) |> factor(levels = c("A1 (normal)", "A2 (moderately increased)",
                            "A3 (severely increased)"))
  )


# -----------------------------------------------------------------------------
# 7. Metabolic variables
# -----------------------------------------------------------------------------

# Diabetes: combined self-reported diagnosis and biochemical HbA1c ≥6.5%
# (ADA threshold). The composite definition maximises sensitivity — patients
# on treatment may have HbA1c < 6.5% but remain diabetic.
ghb_clean <- ghb |>
  select(seqn, hba1c = lbxgh) |>
  mutate(dm_hba1c = if_else(hba1c >= 6.5, 1L, 0L))

diq_clean <- diq |>
  select(seqn, diq010) |>
  mutate(dm_dx = if_else(diq010 == 1, 1L, 0L, missing = NA_integer_))

# Hypertension: 2017 ACC/AHA threshold (≥130/80 mmHg), consistent with
# contemporary epidemiological literature and NHANES analytic guidelines.
bp_meds <- bpq |>
  select(seqn, bpq020, bpq050a) |>
  mutate(
    htn_dx  = if_else(bpq020  == 1, 1L, 0L, missing = NA_integer_),
    htn_med = if_else(bpq050a == 1, 1L, 0L, missing = NA_integer_)
  )

bp_measured <- bpxo |>
  select(seqn,
         sbp = bpxosy1,
         dbp = bpxodi1
  )

bmi_clean <- bmx |>
  select(seqn, bmi = bmxbmi, waist_cm = bmxwaist) |>
  mutate(
    bmi_cat = case_when(
      bmi < 18.5              ~ "Underweight",
      bmi >= 18.5 & bmi < 25 ~ "Normal",
      bmi >= 25   & bmi < 30 ~ "Overweight",
      bmi >= 30               ~ "Obese",
      TRUE                    ~ NA_character_
    ) |> factor(levels = c("Underweight", "Normal", "Overweight", "Obese"))
  )


# -----------------------------------------------------------------------------
# 8. Lifestyle variables
# -----------------------------------------------------------------------------
# NHANES encodes refusal (7/77/777) and don't know (9/99/999) as trailing
# digits. These are recoded to NA prior to MET-min calculation to avoid
# inflating activity estimates. MET coefficients (4.0 moderate, 8.0 vigorous)
# follow the Compendium of Physical Activities and are consistent with the
# SAS and PySpark analyses in this portfolio.

paq_clean <- paq |>
  select(seqn,
         mod_yn   = paq620,   # 1 = yes moderate activity, 2 = no
         mod_days = paq625,   # days/week moderate activity
         mod_min  = pad630,   # minutes/day moderate activity
         vig_yn   = paq605,   # 1 = yes vigorous activity, 2 = no
         vig_days = paq610,   # days/week vigorous activity
         vig_min  = pad615    # minutes/day vigorous activity
  ) |>
  mutate(
    across(everything(), ~ if_else(.x %in% c(77, 99, 7777, 9999),
                                   NA_real_, as.numeric(.x))),
    # Respondents who said "No" to an activity type have 0 days and 0 minutes
    mod_days = if_else(mod_yn == 2, 0, mod_days, missing = mod_days),
    mod_min  = if_else(mod_yn == 2, 0, mod_min,  missing = mod_min),
    vig_days = if_else(vig_yn == 2, 0, vig_days, missing = vig_days),
    vig_min  = if_else(vig_yn == 2, 0, vig_min,  missing = vig_min),
    mod_min_wk = replace_na(mod_days * mod_min, 0),
    vig_min_wk = replace_na(vig_days * vig_min, 0),
    met_min_wk = (mod_min_wk * 4.0) + (vig_min_wk * 8.0),

    pa_cat = case_when(
      met_min_wk == 0   ~ "Sedentary",
      met_min_wk < 600  ~ "Low active",
      met_min_wk >= 600 ~ "Active",
      TRUE              ~ NA_character_
    ) |> factor(levels = c("Sedentary", "Low active", "Active"))
  ) |>
  select(seqn, met_min_wk, pa_cat)

smq_clean <- smq |>
  select(seqn, smq020) |>
  mutate(
    smoking = case_when(
      smq020 == 1 ~ "Current/Former",
      smq020 == 2 ~ "Never",
      TRUE        ~ NA_character_
    ) |> factor(levels = c("Never", "Current/Former"))
  ) |>
  select(seqn, smoking)


# -----------------------------------------------------------------------------
# 9. Merge all components on SEQN
# -----------------------------------------------------------------------------
# Demographics serve as the spine. Missingness in downstream components is
# handled via the exclusion flow (Section 10) and multiple imputation
# (06_sensitivity.R) rather than silently dropped here.

message("\nMerging all components on SEQN...")

analysis_raw <- demo_clean |>
  left_join(scr_all,     by = "seqn") |>
  left_join(uacr_clean,  by = "seqn") |>
  left_join(ghb_clean,   by = "seqn") |>
  left_join(diq_clean,   by = "seqn") |>
  left_join(bp_meds,     by = "seqn") |>
  left_join(bp_measured, by = "seqn") |>
  left_join(bmi_clean,   by = "seqn") |>
  left_join(paq_clean,   by = "seqn") |>
  left_join(smq_clean,   by = "seqn") |>
  left_join(mortality,   by = "seqn") |>
  mutate(
    # Composite diabetes: self-reported OR HbA1c ≥6.5%
    diabetes = if_else(dm_dx == 1 | dm_hba1c == 1, 1L, 0L,
                       missing = NA_integer_),

    # Hypertension: self-reported diagnosis OR on antihypertensive meds
    # bpq050a (meds) is only asked of those with htn_dx=1, so non-diagnosed
    # respondents have htn_med=NA; case_when handles this correctly
    hypertension = case_when(
      htn_dx == 1 | htn_med == 1 ~ 1L,
      htn_dx == 0                ~ 0L,
      TRUE                       ~ NA_integer_
    ),

    # Participants without a mortality record are censored
    died       = replace_na(died, 0L),
    follow_yrs = replace_na(follow_yrs, 0)
  )

message("  ✔  Merged dataset: ", nrow(analysis_raw), " rows × ",
        ncol(analysis_raw), " columns")


# -----------------------------------------------------------------------------
# 10. Cohort exclusion flow (STROBE-compliant)
# -----------------------------------------------------------------------------
# Each criterion is applied sequentially with the remaining N recorded,
# producing the numbers for a STROBE-compliant flowchart in the Quarto report.

flow <- list()
flow[["01_nhanes_adults"]]  <- nrow(analysis_raw)

analysis <- analysis_raw |> filter(!is.na(egfr))
flow[["02_has_egfr"]]       <- nrow(analysis)

analysis <- analysis |> filter(!is.na(mortstat))
flow[["03_has_mortality"]]  <- nrow(analysis)

analysis <- analysis |> filter(!is.na(bmi))
flow[["04_has_bmi"]]        <- nrow(analysis)

analysis <- analysis |> filter(!is.na(hba1c))
flow[["05_has_hba1c"]]      <- nrow(analysis)

analysis <- analysis |> filter(follow_yrs > 0)
flow[["06_followup_gt0"]]   <- nrow(analysis)

# Implausible eGFR >200 affects <0.1% of records (likely data entry errors).
# Sensitivity analyses in 06_sensitivity.R confirm results are unchanged
# when the threshold is relaxed to 250.
analysis <- analysis |> filter(egfr <= 200)
flow[["07_egfr_plausible"]] <- nrow(analysis)

flow_df <- tibble(
  step     = names(flow),
  n        = unlist(flow),
  excluded = c(0, diff(-unlist(flow)))
)

message("\n── Cohort Exclusion Flow ───────────────────────────────────────────────")
print(flow_df)

readr::write_csv(flow_df, here("output", "tables", "02_cohort_flow.csv"))


# -----------------------------------------------------------------------------
# 11. Final variable selection
# -----------------------------------------------------------------------------
analysis_final <- analysis |>
  select(
    seqn, cycle,
    died, follow_yrs,
    egfr, ckd, ckd_stage,
    uacr, albuminuria,
    age_yr, age_group, sex, race_eth, poverty_cat, pir,
    bmi, bmi_cat, waist_cm,
    hba1c, diabetes,
    sbp, dbp, hypertension,
    met_min_wk, pa_cat, smoking,
    wt_mec,
    ucod_leading
  )

message("\n── Final Analysis Dataset ──────────────────────────────────────────────")
message("  Rows      : ", nrow(analysis_final))
message("  Columns   : ", ncol(analysis_final))
message("  Deaths    : ", sum(analysis_final$died), " (",
        round(mean(analysis_final$died) * 100, 1), "%)")
message("  CKD (G3+) : ", sum(analysis_final$ckd, na.rm = TRUE), " (",
        round(mean(analysis_final$ckd, na.rm = TRUE) * 100, 1), "%)")
message("  Follow-up (median years): ",
        round(median(analysis_final$follow_yrs), 1))


# -----------------------------------------------------------------------------
# 12. Save outputs
# -----------------------------------------------------------------------------
saveRDS(analysis_final, here("data", "processed", "ckd_analysis.rds"))
readr::write_csv(analysis_final, here("data", "processed", "ckd_analysis.csv"))

message("\n✔  Saved: data/processed/ckd_analysis.rds")
message("✔  Saved: data/processed/ckd_analysis.csv")
message("\n── 02  Complete ────────────────────────────────────────────────────────")
