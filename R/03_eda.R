# =============================================================================
# 03_eda.R
# CKD Progression & Mortality Risk Factors — NHANES 2017–2023
# Author: Yuntao (Kevin) Tan | tyuntao@umich.edu | December 2025
#
# Purpose: Exploratory data analysis — produce a STROBE-compliant Table 1,
#          descriptive plots, and bivariate summaries stratified by CKD status
#          and vital status. All figures are saved to output/figures/ and all
#          tables to output/tables/ for direct inclusion in the Quarto report.
#
# Input:   data/processed/ckd_analysis.rds
#
# Output:  output/tables/03_table1.html         — publication-ready Table 1
#          output/tables/03_table1.csv          — machine-readable copy
#          output/tables/03_missing_summary.csv — missingness audit
#          output/figures/03_egfr_distribution.png
#          output/figures/03_km_ckd_stage.png
#          output/figures/03_uacr_boxplot.png
#          output/figures/03_risk_heatmap.png
# =============================================================================


# -----------------------------------------------------------------------------
# 0. Packages
# -----------------------------------------------------------------------------
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(
  tidyverse,    # core data manipulation and ggplot2
  here,         # project-relative paths
  gtsummary,    # publication-ready Table 1
  gt,           # table formatting and export
  survival,     # Surv() object for KM curves
  survminer,    # ggsurvplot() — styled KM curves
  ggpubr,       # ggarrange() — multi-panel figure layout
  scales,       # axis formatting helpers
  patchwork     # plot composition
)

here::i_am("R/03_eda.R")

# Consistent ggplot2 theme across all figures
theme_set(
  theme_bw(base_size = 12) +
    theme(
      plot.title    = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(colour = "grey40", size = 10),
      strip.background = element_rect(fill = "grey95"),
      panel.grid.minor = element_blank(),
      legend.position  = "bottom"
    )
)

message("── 03  Exploratory Data Analysis ───────────────────────────────────────")


# -----------------------------------------------------------------------------
# 1. Load analysis dataset
# -----------------------------------------------------------------------------
df <- readRDS(here("data", "processed", "ckd_analysis.rds"))

message("  Dataset: ", nrow(df), " participants | ",
        sum(df$died), " deaths | ",
        sum(df$ckd, na.rm = TRUE), " CKD cases")


# -----------------------------------------------------------------------------
# 2. Missingness audit
# -----------------------------------------------------------------------------
# Reported before Table 1 per STROBE guideline item 12(c). Variables with
# >20% missingness are flagged for sensitivity analysis in 06_sensitivity.R.

missing_summary <- df |>
  summarise(across(everything(), ~ sum(is.na(.x)))) |>
  pivot_longer(everything(), names_to = "variable", values_to = "n_missing") |>
  mutate(
    pct_missing = round(n_missing / nrow(df) * 100, 1),
    flag        = if_else(pct_missing > 20, "⚠ >20%", "")
  ) |>
  arrange(desc(pct_missing))

readr::write_csv(missing_summary,
                 here("output", "tables", "03_missing_summary.csv"))

message("\n── Missingness (top 10 variables) ─────────────────────────────────────")
print(head(missing_summary, 10))


# -----------------------------------------------------------------------------
# 3. Table 1 — baseline characteristics stratified by CKD status
# -----------------------------------------------------------------------------
# gtsummary::tbl_summary() produces a publication-ready table with
# continuous variables summarised as median (IQR) and categorical variables
# as n (%). The add_p() call appends appropriate tests (Wilcoxon rank-sum for
# continuous, chi-square for categorical). add_overall() prepends a full-
# cohort column. The layout mirrors standard nephrology journal format.

table1_vars <- df |>
  select(
    # Outcome
    died, follow_yrs,
    # Renal function
    egfr, ckd_stage, uacr, albuminuria,
    # Demographics
    age_yr, sex, race_eth, poverty_cat,
    # Metabolic
    bmi, bmi_cat, hba1c, diabetes, sbp, dbp, hypertension,
    # Lifestyle
    pa_cat, smoking,
    # Stratification variable
    ckd
  ) |>
  mutate(
    ckd = factor(ckd, levels = c(0, 1),
                 labels = c("eGFR ≥60 (No CKD)", "eGFR <60 (CKD G3–G5)"))
  )

tbl1 <- tbl_summary(
  table1_vars,
  by        = ckd,
  missing   = "ifany",
  statistic = list(
    all_continuous()  ~ "{median} ({p25}, {p75})",
    all_categorical() ~ "{n} ({p}%)"
  ),
  label = list(
    died        ~ "All-cause mortality, n (%)",
    follow_yrs  ~ "Follow-up, years",
    egfr        ~ "eGFR, mL/min/1.73m²",
    ckd_stage   ~ "CKD stage (KDIGO G-category)",
    uacr        ~ "UACR, mg/g",
    albuminuria ~ "Albuminuria category",
    age_yr      ~ "Age, years",
    sex         ~ "Sex",
    race_eth    ~ "Race/ethnicity",
    poverty_cat ~ "Income category",
    bmi         ~ "BMI, kg/m²",
    bmi_cat     ~ "BMI category",
    hba1c       ~ "HbA1c, %",
    diabetes    ~ "Diabetes, n (%)",
    sbp         ~ "Systolic BP, mmHg",
    dbp         ~ "Diastolic BP, mmHg",
    hypertension ~ "Hypertension, n (%)",
    pa_cat      ~ "Physical activity category",
    smoking     ~ "Smoking status"
  ),
  digits = list(
    all_continuous()  ~ 1,
    all_categorical() ~ c(0, 1)
  )
) |>
  add_overall(last = FALSE) |>
  add_p(
    test = list(
      all_continuous()  ~ "wilcox.test",
      all_categorical() ~ "chisq.test"
    ),
    pvalue_fun = ~ style_pvalue(.x, digits = 3)
  ) |>
  add_n() |>
  bold_labels() |>
  modify_header(
    label      ~ "**Characteristic**",
    stat_0     ~ "**Overall**  \nN = {N}",
    stat_1     ~ "**eGFR ≥60**  \nN = {n}",
    stat_2     ~ "**CKD G3–G5**  \nN = {n}",
    p.value    ~ "**p-value**"
  ) |>
  modify_caption(
    "**Table 1.** Baseline characteristics of NHANES 2017–2023 participants
     stratified by chronic kidney disease status. Continuous variables are
     presented as median (IQR); categorical variables as n (%).
     p-values from Wilcoxon rank-sum test (continuous) and
     chi-square test (categorical)."
  )

# Export as HTML (for Quarto report embedding)
tbl1 |>
  as_gt() |>
  gt::gtsave(here("output", "tables", "03_table1.html"))

# Export as CSV for supplementary materials
tbl1 |>
  as_tibble() |>
  readr::write_csv(here("output", "tables", "03_table1.csv"))

message("  ✔  Table 1 saved")


# -----------------------------------------------------------------------------
# 4. eGFR distribution by CKD stage
# -----------------------------------------------------------------------------
p_egfr <- df |>
  filter(!is.na(ckd_stage)) |>
  ggplot(aes(x = egfr, fill = ckd_stage)) +
  geom_histogram(binwidth = 5, colour = "white", linewidth = 0.2) +
  geom_vline(xintercept = 60, linetype = "dashed",
             colour = "firebrick", linewidth = 0.8) +
  annotate("text", x = 63, y = Inf, label = "eGFR = 60\n(CKD threshold)",
           hjust = 0, vjust = 1.5, size = 3.2, colour = "firebrick") +
  scale_x_continuous(limits = c(0, 180), breaks = seq(0, 180, 30)) +
  scale_fill_brewer(palette = "Blues", direction = -1) +
  labs(
    title    = "eGFR Distribution by CKD Stage",
    subtitle = "NHANES 2017–2023 | CKD-EPI 2021 creatinine equation (race-free)",
    x        = "eGFR (mL/min/1.73m²)",
    y        = "Count",
    fill     = "KDIGO G-stage"
  )

ggsave(here("output", "figures", "03_egfr_distribution.png"),
       p_egfr, width = 9, height = 5, dpi = 300)
message("  ✔  Figure: eGFR distribution")


# -----------------------------------------------------------------------------
# 5. Kaplan–Meier survival curves by CKD stage
# -----------------------------------------------------------------------------
# KM curves are stratified by KDIGO G-stage (collapsed to 3 groups for visual
# clarity: G1–G2, G3a–G3b, G4–G5). Log-rank test p-value is displayed.
# This figure motivates the Cox regression in 04_survival_analysis.R.

df_km <- df |>
  filter(!is.na(ckd_stage)) |>
  mutate(
    ckd_grp = case_when(
      ckd_stage %in% c("G1 (≥90)", "G2 (60–89)")     ~ "G1–G2 (eGFR ≥60)",
      ckd_stage %in% c("G3a (45–59)", "G3b (30–44)") ~ "G3 (eGFR 30–59)",
      ckd_stage %in% c("G4 (15–29)", "G5 (<15)")     ~ "G4–G5 (eGFR <30)",
      TRUE ~ NA_character_
    ) |> factor(levels = c("G1–G2 (eGFR ≥60)", "G3 (eGFR 30–59)",
                            "G4–G5 (eGFR <30)"))
  )

km_fit <- survfit(Surv(follow_yrs, died) ~ ckd_grp, data = df_km)

km_plot <- ggsurvplot(
  km_fit,
  data          = df_km,
  pval          = TRUE,
  pval.method   = TRUE,
  conf.int      = TRUE,
  risk.table    = TRUE,
  risk.table.height = 0.28,
  palette       = c("#2166AC", "#F4A582", "#D6604D"),
  xlab          = "Follow-up (years)",
  ylab          = "Overall Survival",
  title         = "Kaplan–Meier Survival by CKD Stage",
  subtitle      = "NHANES 2017–2023 | Log-rank test",
  legend.title  = "CKD Stage",
  legend.labs   = levels(df_km$ckd_grp),
  ggtheme       = theme_bw(base_size = 12),
  fontsize      = 3.5,
  tables.theme  = theme_cleantable()
)

png(here("output", "figures", "03_km_ckd_stage.png"),
    width = 9, height = 7, units = "in", res = 300)
print(km_plot)
dev.off()
message("  ✔  Figure: Kaplan–Meier curves")


# -----------------------------------------------------------------------------
# 6. UACR by CKD stage (boxplot)
# -----------------------------------------------------------------------------
# UACR is log-transformed for display due to extreme right skew. The
# log10 scale is labelled with original units for clinical interpretability.

p_uacr <- df |>
  filter(!is.na(uacr), !is.na(ckd_stage), uacr > 0) |>
  ggplot(aes(x = ckd_stage, y = uacr, fill = ckd_stage)) +
  geom_boxplot(outlier.shape = 21, outlier.size = 1,
               outlier.alpha = 0.4, linewidth = 0.4) +
  geom_hline(yintercept = 30,  linetype = "dashed",
             colour = "orange3", linewidth = 0.7) +
  geom_hline(yintercept = 300, linetype = "dashed",
             colour = "firebrick", linewidth = 0.7) +
  annotate("text", x = 0.55, y = 35,  label = "A2 threshold (30 mg/g)",
           hjust = 0, size = 3, colour = "orange3") +
  annotate("text", x = 0.55, y = 330, label = "A3 threshold (300 mg/g)",
           hjust = 0, size = 3, colour = "firebrick") +
  scale_y_log10(
    breaks = c(1, 10, 30, 100, 300, 1000, 5000),
    labels = scales::comma
  ) +
  scale_fill_brewer(palette = "Blues", direction = -1) +
  labs(
    title    = "Urine Albumin-to-Creatinine Ratio by CKD Stage",
    subtitle = "log₁₀ scale | KDIGO albuminuria thresholds shown",
    x        = "KDIGO G-stage",
    y        = "UACR (mg/g, log scale)",
    fill     = "CKD stage"
  ) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 25, hjust = 1))

ggsave(here("output", "figures", "03_uacr_boxplot.png"),
       p_uacr, width = 8, height = 5, dpi = 300)
message("  ✔  Figure: UACR boxplot")


# -----------------------------------------------------------------------------
# 7. KDIGO risk heat map (G-stage × A-stage)
# -----------------------------------------------------------------------------
# The KDIGO 2024 heat map visualises combined eGFR and albuminuria risk —
# the standard clinical tool for CKD risk communication. Green = low risk,
# orange = moderately increased, red = high/very high risk.

heat_data <- df |>
  filter(!is.na(ckd_stage), !is.na(albuminuria)) |>
  count(ckd_stage, albuminuria) |>
  mutate(
    risk_colour = case_when(
      ckd_stage == "G1 (≥90)"    & albuminuria == "A1 (normal)"                ~ "Low",
      ckd_stage == "G1 (≥90)"    & albuminuria == "A2 (moderately increased)"  ~ "Moderately increased",
      ckd_stage == "G1 (≥90)"    & albuminuria == "A3 (severely increased)"    ~ "High",
      ckd_stage == "G2 (60–89)"  & albuminuria == "A1 (normal)"                ~ "Low",
      ckd_stage == "G2 (60–89)"  & albuminuria == "A2 (moderately increased)"  ~ "Moderately increased",
      ckd_stage == "G2 (60–89)"  & albuminuria == "A3 (severely increased)"    ~ "High",
      ckd_stage == "G3a (45–59)" & albuminuria == "A1 (normal)"                ~ "Moderately increased",
      ckd_stage == "G3a (45–59)" & albuminuria == "A2 (moderately increased)"  ~ "High",
      ckd_stage == "G3a (45–59)" & albuminuria == "A3 (severely increased)"    ~ "Very high",
      ckd_stage == "G3b (30–44)" & albuminuria == "A1 (normal)"                ~ "High",
      ckd_stage == "G3b (30–44)" & albuminuria == "A2 (moderately increased)"  ~ "Very high",
      ckd_stage == "G3b (30–44)" & albuminuria == "A3 (severely increased)"    ~ "Very high",
      ckd_stage %in% c("G4 (15–29)", "G5 (<15)")                              ~ "Very high",
      TRUE ~ "Very high"
    ) |> factor(levels = c("Low", "Moderately increased", "High", "Very high"))
  )

p_heat <- heat_data |>
  ggplot(aes(x = albuminuria, y = fct_rev(ckd_stage), fill = risk_colour)) +
  geom_tile(colour = "white", linewidth = 1.2) +
  geom_text(aes(label = scales::comma(n)), size = 3.5, fontface = "bold") +
  scale_fill_manual(
    values = c(
      "Low"                  = "#4DAF4A",
      "Moderately increased" = "#FFFF33",
      "High"                 = "#FF7F00",
      "Very high"            = "#E41A1C"
    )
  ) +
  scale_x_discrete(position = "top") +
  labs(
    title    = "KDIGO 2024 Risk Heat Map — CKD Stage × Albuminuria",
    subtitle = "Cell values = number of NHANES participants | Colour = prognosis category",
    x        = "Albuminuria Category (A-stage)",
    y        = "eGFR Category (G-stage)",
    fill     = "Risk category"
  ) +
  theme(
    axis.text.x  = element_text(angle = 20, hjust = 0),
    panel.grid   = element_blank(),
    panel.border = element_blank()
  )

ggsave(here("output", "figures", "03_risk_heatmap.png"),
       p_heat, width = 9, height = 6, dpi = 300)
message("  ✔  Figure: KDIGO risk heat map")


# -----------------------------------------------------------------------------
# 8. Mortality rate by key risk factors (descriptive)
# -----------------------------------------------------------------------------
# Unadjusted mortality rates per 100 person-years provide the descriptive
# context for the Cox regression. Rates are calculated as events / total
# person-years × 100, the standard epidemiological measure.

calc_rate <- function(data, group_var) {
  data |>
    group_by({{ group_var }}) |>
    summarise(
      n          = n(),
      deaths     = sum(died),
      pyears     = sum(follow_yrs),
      rate_100py = round(deaths / pyears * 100, 2),
      .groups    = "drop"
    ) |>
    filter(!is.na({{ group_var }}))
}

rates_ckd  <- calc_rate(df, ckd_stage)
rates_dm   <- calc_rate(df, diabetes) |>
  mutate(diabetes = factor(diabetes, 0:1, c("No diabetes", "Diabetes")))
rates_htn  <- calc_rate(df, hypertension) |>
  mutate(hypertension = factor(hypertension, 0:1, c("No hypertension", "Hypertension")))
rates_pa   <- calc_rate(df, pa_cat)

rate_plots <- list(
  ckd  = ggplot(rates_ckd,  aes(x = fct_rev(ckd_stage),  y = rate_100py)) +
    labs(x = "CKD Stage", title = "CKD Stage"),
  dm   = ggplot(rates_dm,   aes(x = diabetes,   y = rate_100py)) +
    labs(x = "Diabetes", title = "Diabetes"),
  htn  = ggplot(rates_htn,  aes(x = hypertension, y = rate_100py)) +
    labs(x = "Hypertension", title = "Hypertension"),
  pa   = ggplot(rates_pa,   aes(x = pa_cat,     y = rate_100py)) +
    labs(x = "Physical Activity", title = "Physical Activity")
)

rate_plots <- purrr::map(rate_plots, function(p) {
  p +
    geom_col(fill = "#2166AC", alpha = 0.85, width = 0.6) +
    geom_text(aes(label = rate_100py), vjust = -0.4, size = 3.2) +
    labs(y = "Mortality rate\n(per 100 person-years)") +
    coord_flip() +
    theme(axis.text.y = element_text(size = 9))
})

p_rates <- (rate_plots$ckd | rate_plots$dm) /
           (rate_plots$htn | rate_plots$pa) +
  plot_annotation(
    title    = "Unadjusted All-Cause Mortality Rates by Key Risk Factors",
    subtitle = "NHANES 2017–2023 | Per 100 person-years",
    theme    = theme(plot.title = element_text(face = "bold", size = 13))
  )

ggsave(here("output", "figures", "03_mortality_rates.png"),
       p_rates, width = 11, height = 7, dpi = 300)
message("  ✔  Figure: Unadjusted mortality rates")


# -----------------------------------------------------------------------------
# 9. Summary statistics to console
# -----------------------------------------------------------------------------
message("\n── Cohort Summary ──────────────────────────────────────────────────────")
message("  Median eGFR        : ", median(df$egfr, na.rm = TRUE))
message("  CKD prevalence     : ",
        round(mean(df$ckd, na.rm = TRUE) * 100, 1), "%")
message("  Diabetes prevalence: ",
        round(mean(df$diabetes, na.rm = TRUE) * 100, 1), "%")
message("  HTN prevalence     : ",
        round(mean(df$hypertension, na.rm = TRUE) * 100, 1), "%")
message("  Mortality rate     : ",
        round(sum(df$died) / sum(df$follow_yrs) * 100, 2),
        " per 100 person-years")

message("\n── 03  Complete ────────────────────────────────────────────────────────")
