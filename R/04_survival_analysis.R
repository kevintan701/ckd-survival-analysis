# =============================================================================
# 04_survival_analysis.R
# CKD Progression & Mortality Risk Factors — NHANES 2017–2023
# Author: Yuntao (Kevin) Tan | tyuntao@umich.edu | December 2025
#
# Purpose: Primary inference stage — fit unadjusted and multivariable Cox
#          proportional hazards models, verify the proportional hazards (PH)
#          assumption, and produce publication-ready hazard ratio tables and
#          forest plots. This script is the analytical centrepiece of the
#          project; all downstream tidymodels work in 05_tidymodels.R extends
#          and cross-validates these findings.
#
# Models fitted:
#   Model 1 — Unadjusted (eGFR / CKD stage only)
#   Model 2 — Demographically adjusted (+ age, sex, race/ethnicity, income)
#   Model 3 — Fully adjusted (+ diabetes, hypertension, BMI, PA, smoking)
#
# Input:   data/processed/ckd_analysis.rds
#
# Output:  output/tables/04_cox_model1.csv   — unadjusted HR table
#          output/tables/04_cox_model3.csv   — fully adjusted HR table
#          output/tables/04_ph_test.csv      — Schoenfeld PH test results
#          output/figures/04_forest_plot.png — HR forest plot (Model 3)
#          output/figures/04_schoenfeld.png  — PH assumption diagnostics
#          output/figures/04_martingale.png  — functional form check
# =============================================================================


# -----------------------------------------------------------------------------
# 0. Packages
# -----------------------------------------------------------------------------
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(
  tidyverse,    # data manipulation
  here,         # project-relative paths
  survival,     # coxph(), Surv(), cox.zph(), survfit()
  survminer,    # ggcoxdiagnostics(), ggforest()
  gtsummary,    # tbl_regression() — publication-ready HR tables
  gt,           # table export
  broom,        # tidy() — extract model coefficients as tibbles
  patchwork     # multi-panel figure layout
)

here::i_am("R/04_survival_analysis.R")

theme_set(
  theme_bw(base_size = 12) +
    theme(
      plot.title       = element_text(face = "bold", size = 13),
      plot.subtitle    = element_text(colour = "grey40", size = 10),
      panel.grid.minor = element_blank(),
      legend.position  = "bottom"
    )
)

message("── 04  Survival Analysis ────────────────────────────────────────────────")


# -----------------------------------------------------------------------------
# 1. Load data and define analysis sample
# -----------------------------------------------------------------------------
df <- readRDS(here("data", "processed", "ckd_analysis.rds"))

# Continuous eGFR is scaled to a clinically interpretable unit: per 10
# mL/min/1.73m² decline. This makes the HR directly comparable to published
# USRDS and CKD-PC consortium estimates and avoids the near-null HR that
# results from a per-unit (1 mL/min) parameterisation.
df_cox <- df |>
  mutate(
    egfr_10       = egfr / 10,          # HR per 10-unit decline (reverse coded below)
    egfr_10_dec   = (egfr * -1) / 10,   # positive coefficient = higher risk with lower eGFR
    log_uacr      = log(uacr + 1),      # log-transform to reduce leverage of extreme values
    age_10        = age_yr / 10,        # HR per decade — standard in geriatric epidemiology
    met_100       = met_min_wk / 100    # HR per 100 MET-min/week increment
  ) |>
  # Relevel reference categories to match epidemiological conventions
  mutate(
    ckd_stage   = relevel(ckd_stage,   ref = "G1 (≥90)"),
    race_eth    = relevel(race_eth,    ref = "Non-Hispanic White"),
    poverty_cat = relevel(poverty_cat, ref = "High income"),
    bmi_cat     = relevel(bmi_cat,     ref = "Normal"),
    pa_cat      = relevel(pa_cat,      ref = "Active"),
    smoking     = relevel(smoking,     ref = "Never")
  )

message("  Analysis sample: ", nrow(df_cox), " participants | ",
        sum(df_cox$died), " events | ",
        round(sum(df_cox$follow_yrs), 0), " person-years")


# -----------------------------------------------------------------------------
# 2. Define Surv() object
# -----------------------------------------------------------------------------
# Surv(time, event) is the response object for all Cox models. time = years
# of follow-up from MEC exam date; event = 1 if died, 0 if censored.
# Right-censoring is assumed non-informative — standard in linked mortality
# analyses where censoring reflects end of NDI follow-up window, not dropout.

surv_obj <- Surv(time = df_cox$follow_yrs, event = df_cox$died)


# -----------------------------------------------------------------------------
# 3. Model 1 — Unadjusted: eGFR stage and continuous eGFR
# -----------------------------------------------------------------------------
# Unadjusted estimates establish the crude association before confounding
# adjustment. Both stage (categorical) and continuous eGFR are fitted
# separately to report the full clinical picture.

cox_m1_stage <- coxph(
  surv_obj ~ ckd_stage,
  data    = df_cox,
  ties    = "efron"    # Efron approximation — preferred over Breslow for tied
                       # survival times, which are common in annual surveys
)

cox_m1_cont <- coxph(
  surv_obj ~ egfr_10_dec,
  data = df_cox,
  ties = "efron"
)

message("\n── Model 1: Unadjusted ─────────────────────────────────────────────────")
print(summary(cox_m1_stage)$coefficients)


# -----------------------------------------------------------------------------
# 4. Model 2 — Demographically adjusted
# -----------------------------------------------------------------------------
# Age, sex, race/ethnicity, and income are pre-specified confounders
# regardless of their p-values, per established epidemiological practice.
# Forcing these into the model prevents residual confounding that would
# otherwise bias the CKD–mortality association estimate.

cox_m2 <- coxph(
  surv_obj ~ ckd_stage + age_10 + sex + race_eth + poverty_cat,
  data = df_cox,
  ties = "efron"
)

message("\n── Model 2: Demographically adjusted ───────────────────────────────────")
print(summary(cox_m2)$coefficients)


# -----------------------------------------------------------------------------
# 5. Model 3 — Fully adjusted (primary inference model)
# -----------------------------------------------------------------------------
# The fully adjusted model includes all pre-specified covariates from the
# causal model. Metabolic (diabetes, hypertension, BMI) and lifestyle
# (physical activity, smoking) variables are included as they represent
# both confounders and mediators of the CKD–mortality pathway. Continuous
# UACR is included alongside eGFR stage to capture the independent
# prognostic contribution of albuminuria, consistent with KDIGO 2024 guidance.

cox_m3 <- coxph(
  surv_obj ~ ckd_stage + log_uacr +
    age_10 + sex + race_eth + poverty_cat +
    diabetes + hypertension + bmi_cat +
    pa_cat + smoking,
  data = df_cox,
  ties = "efron"
)

message("\n── Model 3: Fully adjusted (primary) ───────────────────────────────────")
print(summary(cox_m3))


# -----------------------------------------------------------------------------
# 6. Proportional hazards assumption — Schoenfeld residual test
# -----------------------------------------------------------------------------
# The PH assumption requires that the hazard ratio for each covariate is
# constant over time. cox.zph() regresses scaled Schoenfeld residuals on
# time; a significant p-value (< 0.05) indicates time-varying hazards.
# Variables that violate PH are handled via time-interaction terms or
# stratification in the sensitivity analysis (06_sensitivity.R).

ph_test <- cox.zph(cox_m3, transform = "km")

message("\n── Proportional Hazards Test (Schoenfeld residuals) ────────────────────")
print(ph_test)

ph_df <- as.data.frame(ph_test$table) |>
  rownames_to_column("term") |>
  rename(chisq = chisq, df = df, p_value = p) |>
  mutate(
    ph_violated = p_value < 0.05,
    p_value     = round(p_value, 4)
  )

readr::write_csv(ph_df, here("output", "tables", "04_ph_test.csv"))

# Schoenfeld residual plots — one panel per covariate
png(here("output", "figures", "04_schoenfeld.png"),
    width = 12, height = 10, units = "in", res = 300)
par(mfrow = c(4, 3), mar = c(4, 4, 2, 1))
plot(ph_test)
dev.off()
message("  ✔  Figure: Schoenfeld residual plots")


# -----------------------------------------------------------------------------
# 7. Martingale residuals — functional form of continuous predictors
# -----------------------------------------------------------------------------
# Martingale residuals from a null Cox model are plotted against each
# continuous predictor to verify that the assumed linear functional form
# is appropriate. A non-linear LOESS trend indicates the need for
# restricted cubic splines or categorisation.

cox_null <- coxph(surv_obj ~ 1, data = df_cox, ties = "efron")
df_cox$mart_resid <- residuals(cox_null, type = "martingale")

p_mart_egfr <- ggplot(df_cox, aes(x = egfr, y = mart_resid)) +
  geom_point(alpha = 0.15, size = 0.8, colour = "grey50") +
  geom_smooth(method = "loess", se = TRUE,
              colour = "#2166AC", fill = "#2166AC", alpha = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "firebrick") +
  labs(x = "eGFR (mL/min/1.73m²)", y = "Martingale residual",
       title = "eGFR")

p_mart_age <- ggplot(df_cox, aes(x = age_yr, y = mart_resid)) +
  geom_point(alpha = 0.15, size = 0.8, colour = "grey50") +
  geom_smooth(method = "loess", se = TRUE,
              colour = "#2166AC", fill = "#2166AC", alpha = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "firebrick") +
  labs(x = "Age (years)", y = "Martingale residual",
       title = "Age")

p_mart_uacr <- ggplot(df_cox |> filter(uacr < 1000),
                      aes(x = log_uacr, y = mart_resid)) +
  geom_point(alpha = 0.15, size = 0.8, colour = "grey50") +
  geom_smooth(method = "loess", se = TRUE,
              colour = "#2166AC", fill = "#2166AC", alpha = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "firebrick") +
  labs(x = "log(UACR + 1)", y = "Martingale residual",
       title = "log(UACR)")

p_mart <- (p_mart_egfr | p_mart_age | p_mart_uacr) +
  plot_annotation(
    title    = "Martingale Residuals — Functional Form Check",
    subtitle = "LOESS smooth with 95% CI | Deviation from 0 suggests non-linearity",
    theme    = theme(plot.title = element_text(face = "bold", size = 13))
  )

ggsave(here("output", "figures", "04_martingale.png"),
       p_mart, width = 12, height = 4, dpi = 300)
message("  ✔  Figure: Martingale residual plots")


# -----------------------------------------------------------------------------
# 8. Publication-ready HR tables via gtsummary
# -----------------------------------------------------------------------------
# tbl_regression() extracts coefficients, exponentiated to HR, with 95% CI
# and p-values. The layout matches Table 2 format in nephrology journals
# (JASN, CJASN, KI). Three-model comparison table allows readers to trace
# how estimates change with progressive confounder adjustment — a standard
# approach in observational epidemiology to assess confounding magnitude.

all_labels <- list(
  ckd_stage    ~ "CKD stage (ref: G1 ≥90)",
  log_uacr     ~ "log(UACR + 1), per unit",
  age_10       ~ "Age, per decade",
  sex          ~ "Sex (ref: Male)",
  race_eth     ~ "Race/ethnicity (ref: Non-Hispanic White)",
  poverty_cat  ~ "Income category (ref: High income)",
  diabetes     ~ "Diabetes (ref: No)",
  hypertension ~ "Hypertension (ref: No)",
  bmi_cat      ~ "BMI category (ref: Normal)",
  pa_cat       ~ "Physical activity (ref: Active)",
  smoking      ~ "Smoking (ref: Never)"
)

fmt_hr <- function(model, col_label) {
  # Keep only labels for variables actually in this model
  model_vars  <- unique(gsub("([A-Za-z_]+).*", "\\1", names(coef(model))))
  valid_labels <- Filter(
    function(x) as.character(x[[2]]) %in% model_vars,
    all_labels
  )

  tbl_regression(
    model,
    exponentiate = TRUE,
    label        = valid_labels,
    conf.int     = TRUE,
    pvalue_fun   = ~ style_pvalue(.x, digits = 3)
  ) |>
    bold_p(t = 0.05) |>
    bold_labels() |>
    modify_header(estimate ~ glue::glue("**{col_label}**"))
}

tbl_m1 <- fmt_hr(cox_m1_stage, "HR (Unadjusted)")
tbl_m2 <- fmt_hr(cox_m2,       "HR (Model 2)")
tbl_m3 <- fmt_hr(cox_m3,       "HR (Model 3)")


# Merge the three models into a single comparison table
tbl_merged <- tbl_merge(
  tbls        = list(tbl_m1, tbl_m2, tbl_m3),
  tab_spanner = c("**Model 1**<br>Unadjusted",
                  "**Model 2**<br>+ Demographics",
                  "**Model 3**<br>Fully adjusted")
) |>
  modify_caption(
    "**Table 2.** Hazard ratios (95% CI) for all-cause mortality from Cox
     proportional hazards models. Model 1: unadjusted. Model 2: adjusted for
     age, sex, race/ethnicity, and income. Model 3: additionally adjusted for
     diabetes, hypertension, BMI, physical activity, and smoking.
     Bold p-values indicate statistical significance at α = 0.05."
  )

tbl_merged |>
  as_gt() |>
  gt::gtsave(here("output", "tables", "04_cox_hr_table.html"))

# Individual model CSV exports for supplementary materials
broom::tidy(cox_m3, exponentiate = TRUE, conf.int = TRUE) |>
  readr::write_csv(here("output", "tables", "04_cox_model3.csv"))

broom::tidy(cox_m1_stage, exponentiate = TRUE, conf.int = TRUE) |>
  readr::write_csv(here("output", "tables", "04_cox_model1.csv"))

message("  ✔  HR tables saved")


# -----------------------------------------------------------------------------
# 9. Forest plot — Model 3 fully adjusted HRs
# -----------------------------------------------------------------------------
# The forest plot is the primary visual for the inference results. Terms are
# ordered by HR magnitude. Reference categories and the global p-value from
# the likelihood ratio test are annotated.

forest_data <- broom::tidy(cox_m3, exponentiate = TRUE, conf.int = TRUE) |>
  filter(!str_detect(term, "^ckd_stageG1")) |>   # drop reference category
  # Drop terms with non-estimable CIs (sparse cells, e.g. G5 n<5 events)
  filter(is.finite(conf.low) & conf.low > 0 & is.finite(conf.high)) |>
  mutate(
    term_label = case_when(
      term == "ckd_stageG2 (60–89)"         ~ "CKD G2 vs G1",
      term == "ckd_stageG3a (45–59)"        ~ "CKD G3a vs G1",
      term == "ckd_stageG3b (30–44)"        ~ "CKD G3b vs G1",
      term == "ckd_stageG4 (15–29)"         ~ "CKD G4 vs G1",
      term == "ckd_stageG5 (<15)"                ~ "CKD G5 vs G1",
      term == "log_uacr"                         ~ "log(UACR), per unit",
      term == "age_10"                           ~ "Age, per decade",
      term == "sexFemale"                        ~ "Female vs Male",
      str_detect(term, "race_eth")               ~ str_remove(term, "race_eth"),
      str_detect(term, "poverty_cat")            ~ str_remove(term, "poverty_cat"),
      term == "diabetes"                         ~ "Diabetes",
      term == "hypertension"                     ~ "Hypertension",
      str_detect(term, "bmi_cat")                ~ str_remove(term, "bmi_cat"),
      term == "pa_catSedentary"                  ~ "Sedentary vs Active",
      term == "pa_catLow active"                 ~ "Low active vs Active",
      term == "smokingCurrent/Former"            ~ "Current/Former smoker",
      TRUE                                       ~ term
    ),
    significant = p.value < 0.05,
    direction   = if_else(estimate >= 1, "Increased risk", "Decreased risk")
  ) |>
  arrange(estimate)

p_forest <- ggplot(forest_data,
                   aes(x = estimate, y = reorder(term_label, estimate),
                       colour = direction)) +
  geom_vline(xintercept = 1, linetype = "dashed",
             colour = "grey40", linewidth = 0.7) +
  geom_errorbar(aes(xmin = conf.low, xmax = conf.high),
                width = 0.25, linewidth = 0.6, orientation = "y") +
  geom_point(aes(size = significant), shape = 18) +
  scale_size_manual(values = c("TRUE" = 4, "FALSE" = 2.5),
                    guide = "none") +
  scale_colour_manual(
    values = c("Increased risk" = "#D6604D", "Decreased risk" = "#2166AC")
  ) +
  scale_x_log10(
    breaks = c(0.25, 0.5, 0.75, 1, 1.5, 2, 3, 5),
    labels = c("0.25", "0.50", "0.75", "1.00", "1.50", "2.00", "3.00", "5.00")
  ) +
  labs(
    title    = "Hazard Ratios for All-Cause Mortality — Model 3 (Fully Adjusted)",
    subtitle = "Cox proportional hazards | NHANES 2017–2023 | Diamond size: p < 0.05",
    x        = "Hazard Ratio (log scale, 95% CI)",
    y        = NULL,
    colour   = NULL
  ) +
  theme(
    legend.position  = "top",
    axis.text.y      = element_text(size = 9),
    panel.grid.major.y = element_line(colour = "grey93")
  )

ggsave(here("output", "figures", "04_forest_plot.png"),
       p_forest, width = 10, height = 9, dpi = 300)
message("  ✔  Figure: Forest plot")


# -----------------------------------------------------------------------------
# 10. Model fit statistics
# -----------------------------------------------------------------------------
# Concordance (C-statistic) is the survival-analysis analogue of AUC-ROC.
# Values are compared across models to quantify the discriminative gain
# from adding demographic and clinical confounders.

extract_fit <- function(model, name) {
  s <- summary(model)
  tibble(
    model         = name,
    n             = s$n,
    events        = s$nevent,
    concordance   = round(s$concordance[["C"]], 3),
    concordance_se= round(s$concordance[["se(C)"]], 4),
    loglik_null   = round(s$loglik[1], 1),
    loglik_model  = round(s$loglik[2], 1),
    lr_chisq      = round(s$logtest[["test"]], 2),
    lr_df         = s$logtest[["df"]],
    lr_p          = signif(s$logtest[["pvalue"]], 3),
    wald_chisq    = round(s$waldtest[["test"]], 2)
  )
}

fit_stats <- bind_rows(
  extract_fit(cox_m1_stage, "Model 1: Unadjusted"),
  extract_fit(cox_m2,       "Model 2: + Demographics"),
  extract_fit(cox_m3,       "Model 3: Fully adjusted")
)

message("\n── Model Fit Statistics ────────────────────────────────────────────────")
print(fit_stats)

readr::write_csv(fit_stats, here("output", "tables", "04_model_fit.csv"))


# -----------------------------------------------------------------------------
# 11. Save model objects for downstream use in 05_tidymodels.R
# -----------------------------------------------------------------------------
saveRDS(cox_m3,   here("data", "processed", "cox_m3.rds"))
saveRDS(df_cox,   here("data", "processed", "df_cox.rds"))
saveRDS(surv_obj, here("data", "processed", "surv_obj.rds"))

message("\n✔  Model objects saved to data/processed/")
message("\n── 04  Complete ────────────────────────────────────────────────────────")
