# =============================================================================
# 05_tidymodels.R
# CKD Progression & Mortality Risk Factors — NHANES 2017–2023
# Author: Yuntao (Kevin) Tan | tyuntao@umich.edu | December 2025
#
# Purpose: Implement the primary Cox model within the tidymodels framework
#          to demonstrate reproducible, cross-validated survival modelling.
#          The tidymodels workflow wraps the same inferential model from
#          04_survival_analysis.R in a standardised pipeline that supports
#          resampling-based performance evaluation, recipe-driven preprocessing,
#          and model comparison — capabilities absent from base coxph().
#
#          This script complements rather than replaces 04_survival_analysis.R:
#          inference (HR, CI, p-values) is reported from the full-data coxph()
#          fit; tidymodels provides the cross-validated C-statistic and
#          preprocessing audit trail that reviewers increasingly expect.
#
# Framework:
#   rsample   — stratified V-fold cross-validation
#   recipes   — reproducible preprocessing pipeline
#   parsnip   — model specification (proportional_hazards engine)
#   workflows — recipe + model bundling
#   tune      — resampling execution and metric collection
#   yardstick — C-statistic (concordance_survival) evaluation
#   censored  — survival model extensions for tidymodels
#
# Input:   data/processed/ckd_analysis.rds
#          data/processed/cox_m3.rds
#
# Output:  output/tables/05_cv_results.csv       — fold-level C-statistics
#          output/tables/05_model_comparison.csv — tidymodels vs base coxph
#          output/figures/05_cv_cstat.png        — cross-validated C-stat plot
#          output/figures/05_calibration.png     — observed vs predicted risk
# =============================================================================


# -----------------------------------------------------------------------------
# 0. Packages
# -----------------------------------------------------------------------------
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(
  tidyverse,    # core
  here,         # paths
  tidymodels,   # rsample, recipes, parsnip, workflows, tune, yardstick
  censored,     # survival model support in tidymodels (proportional_hazards)
  survival,     # Surv(), survfit() — underlying engine
  broom,        # tidy() model outputs
  patchwork     # figure composition
)

here::i_am("R/05_tidymodels.R")

# tidymodels produces verbose output during resampling; suppress for cleaner
# logs while retaining error messages
tidymodels_prefer(quiet = TRUE)

theme_set(
  theme_bw(base_size = 12) +
    theme(
      plot.title       = element_text(face = "bold", size = 13),
      plot.subtitle    = element_text(colour = "grey40", size = 10),
      panel.grid.minor = element_blank(),
      legend.position  = "bottom"
    )
)

message("── 05  tidymodels Workflow ──────────────────────────────────────────────")


# -----------------------------------------------------------------------------
# 1. Load data
# -----------------------------------------------------------------------------
df_cox <- readRDS(here("data", "processed", "df_cox.rds"))
cox_m3 <- readRDS(here("data", "processed", "cox_m3.rds"))

# tidymodels survival models expect the outcome encoded as a Surv object
# within the dataframe, passed via the formula interface. The outcome column
# is named .surv_col by convention in censored.
df_tm <- df_cox |>
  mutate(
    surv_outcome = Surv(follow_yrs, died)
  ) |>
  # Drop rows with NA in any model covariate — tidymodels recipes handle
  # imputation explicitly; ad-hoc NA dropping here ensures fold integrity
  drop_na(ckd_stage, log_uacr, age_10, sex, race_eth, poverty_cat,
          diabetes, hypertension, bmi_cat, pa_cat, smoking)

message("  Modelling sample after complete-case restriction: ", nrow(df_tm))


# -----------------------------------------------------------------------------
# 2. Cross-validation strategy
# -----------------------------------------------------------------------------
# 10-fold stratified cross-validation, stratified on the event indicator
# (died) to ensure each fold has a representative proportion of events.
# A fixed seed guarantees reproducibility across runs.
#
# 10 folds are used rather than 5 because the event rate in this sample is
# moderate; more folds reduce bias in the C-statistic estimate at the cost
# of higher variance, which is acceptable given N > 5,000.

set.seed(2025)

cv_folds <- vfold_cv(
  df_tm,
  v       = 10,
  strata  = died    # stratify on event indicator
)

message("  Cross-validation: ", cv_folds$id |> length(), "-fold, stratified on event")


# -----------------------------------------------------------------------------
# 3. Preprocessing recipe
# -----------------------------------------------------------------------------
# The recipe defines all preprocessing steps in a declarative, reproducible
# manner. Each step is estimated on the training fold and applied to the
# assessment fold — preventing data leakage that would inflate performance
# estimates. This is the key advantage of the recipe paradigm over ad-hoc
# preprocessing in base R.

cox_recipe <- recipe(surv_outcome ~ ckd_stage + log_uacr + age_10 +
                       sex + race_eth + poverty_cat +
                       diabetes + hypertension + bmi_cat +
                       pa_cat + smoking,
                     data = df_tm) |>

  # Collapse rare race/ethnicity levels to prevent zero-cell folds.
  # Levels representing < 5% of the training sample are pooled into "Other".
  step_other(race_eth, poverty_cat, threshold = 0.05) |>

  # One-hot encode all nominal predictors. Cox regression in the survival
  # package requires dummy variables; tidymodels handles this transparently.
  step_dummy(all_nominal_predictors(), one_hot = FALSE) |>

  # Normalise continuous predictors to mean = 0, SD = 1.
  # Standardisation is not required for Cox regression (which is scale-
  # invariant) but facilitates numerical stability in edge-case folds with
  # restricted covariate ranges.
  step_normalize(all_numeric_predictors()) |>

  # Remove zero-variance predictors that may arise in small folds
  step_zv(all_predictors())

message("  Recipe defined: ", length(cox_recipe$steps), " preprocessing steps")


# -----------------------------------------------------------------------------
# 4. Model specification
# -----------------------------------------------------------------------------
# proportional_hazards() is the parsnip model type for Cox regression.
# The "survival" engine calls coxph() internally, ensuring numerical
# equivalence with the results in 04_survival_analysis.R.
# Setting mode = "censored regression" is required by censored to
# distinguish survival outcomes from standard regression.

cox_spec <- proportional_hazards(
  penalty = NULL    # no regularisation — inference model, not penalised
) |>
  set_engine("survival") |>
  set_mode("censored regression")


# -----------------------------------------------------------------------------
# 5. Workflow — bundle recipe and model
# -----------------------------------------------------------------------------
# A workflow object couples the preprocessing recipe and model specification.
# This ensures the same transformations are applied consistently across all
# training and assessment splits without manual repetition.

cox_workflow <- workflow() |>
  add_recipe(cox_recipe) |>
  add_model(cox_spec,
            formula = surv_outcome ~ .)   # "." expands to all recipe outputs


# -----------------------------------------------------------------------------
# 6. Fit across cross-validation folds
# -----------------------------------------------------------------------------
# fit_resamples() trains the workflow on each training fold and evaluates on
# the held-out assessment fold. The C-statistic (concordance_survival) is
# computed on each assessment fold and aggregated.
#
# concordance_survival is the Harrell's C-statistic — the probability that,
# for two randomly selected participants, the one who dies first had the
# higher predicted risk. It is the survival analogue of AUC-ROC.

message("\n  Fitting Cox model across 10 CV folds...")

cox_cv_results <- fit_resamples(
  cox_workflow,
  resamples  = cv_folds,
  metrics    = metric_set(concordance_survival),
  control    = control_resamples(
    save_pred    = TRUE,
    verbose      = FALSE,
    allow_par    = FALSE    # sequential for reproducibility
  )
)

message("  ✔  Cross-validation complete")


# -----------------------------------------------------------------------------
# 7. Extract and summarise CV performance
# -----------------------------------------------------------------------------
cv_metrics <- collect_metrics(cox_cv_results)

message("\n── Cross-validated Performance ─────────────────────────────────────────")
print(cv_metrics)

# Fold-level results for distributional summary and plot
cv_fold_metrics <- collect_metrics(cox_cv_results, summarize = FALSE) |>
  filter(.metric == "concordance_survival")

readr::write_csv(cv_fold_metrics,
                 here("output", "tables", "05_cv_results.csv"))


# -----------------------------------------------------------------------------
# 8. Compare tidymodels CV C-stat vs base coxph full-data C-stat
# -----------------------------------------------------------------------------
# The full-data C-statistic from coxph() is optimistically biased because
# it is evaluated on the same data used for fitting. The CV C-statistic is
# the unbiased estimate. Reporting both demonstrates methodological rigour
# and quantifies the magnitude of overfitting optimism.

base_cstat <- summary(cox_m3)$concordance

comparison <- tibble(
  method          = c("Base coxph (full data, optimistic)",
                      "tidymodels 10-fold CV (unbiased)"),
  c_statistic     = c(round(base_cstat[["C"]], 3),
                      round(cv_metrics$mean, 3)),
  se              = c(round(base_cstat[["se(C)"]], 4),
                      round(cv_metrics$std_err, 4)),
  note            = c("Evaluated on training data — subject to overfitting bias",
                      "Evaluated on held-out folds — preferred for reporting")
)

message("\n── Model Comparison: Full-data vs Cross-validated C-statistic ──────────")
print(comparison)

readr::write_csv(comparison,
                 here("output", "tables", "05_model_comparison.csv"))


# -----------------------------------------------------------------------------
# 9. Figure: fold-level C-statistic distribution
# -----------------------------------------------------------------------------
p_cv <- cv_fold_metrics |>
  mutate(fold = as.integer(str_extract(id, "\\d+"))) |>
  ggplot(aes(x = fold, y = .estimate)) +
  geom_hline(
    yintercept = cv_metrics$mean,
    linetype   = "dashed", colour = "#2166AC", linewidth = 0.8
  ) +
  geom_hline(
    yintercept = base_cstat[["C"]],
    linetype   = "dotted", colour = "firebrick", linewidth = 0.8
  ) +
  geom_col(fill = "#4393C3", alpha = 0.8, width = 0.6) +
  geom_point(size = 3, colour = "#1A5276") +
  annotate("text",
           x = 10.4, y = cv_metrics$mean + 0.003,
           label = sprintf("CV mean = %.3f", cv_metrics$mean),
           colour = "#2166AC", hjust = 1, size = 3.5) +
  annotate("text",
           x = 10.4, y = base_cstat[["C"]] + 0.003,
           label = sprintf("Full-data = %.3f", base_cstat[["C"]]),
           colour = "firebrick", hjust = 1, size = 3.5) +
  scale_x_continuous(breaks = 1:10) +
  scale_y_continuous(limits = c(0.5, 0.85),
                     breaks = seq(0.5, 0.85, 0.05)) +
  labs(
    title    = "10-Fold Cross-Validated C-Statistic by Fold",
    subtitle = "Dashed = CV mean | Dotted = full-data (optimistic) | Cox PH model",
    x        = "Fold",
    y        = "Concordance (C-statistic)"
  )

ggsave(here("output", "figures", "05_cv_cstat.png"),
       p_cv, width = 9, height = 5, dpi = 300)
message("  ✔  Figure: CV C-statistic by fold")


# -----------------------------------------------------------------------------
# 10. Calibration — observed vs predicted survival at 5 years
# -----------------------------------------------------------------------------
# Calibration assesses whether predicted survival probabilities match
# observed event rates. The fitted workflow is retrained on the full dataset
# and predicted survival at 5 years is extracted, then decile-binned against
# observed Kaplan–Meier survival to produce a calibration plot.
# Perfect calibration lies on the 45-degree diagonal.

cox_final_fit <- fit(cox_workflow, data = df_tm)

# Extract 5-year survival predictions
pred_surv <- predict(cox_final_fit,
                     new_data  = df_tm,
                     type      = "survival",
                     eval_time = 5) |>
  bind_cols(df_tm |> select(follow_yrs, died))

# Decile-bin predicted survival and compute observed KM survival per decile
pred_surv <- pred_surv |>
  mutate(
    pred_5yr  = .pred |> purrr::map_dbl(~ .x$.pred_survival),
    decile    = ntile(pred_5yr, 10)
  )

km_by_decile <- pred_surv |>
  group_by(decile) |>
  summarise(
    mean_pred   = mean(pred_5yr),
    km_fit      = list(survfit(Surv(follow_yrs, died) ~ 1,
                               data = cur_data())),
    .groups     = "drop"
  ) |>
  mutate(
    km_surv_5yr = purrr::map_dbl(km_fit, function(fit) {
      # Extract KM survival estimate at 5 years
      idx <- max(which(fit$time <= 5), 1)
      fit$surv[idx]
    })
  ) |>
  select(decile, mean_pred, km_surv_5yr)

p_calib <- ggplot(km_by_decile,
                  aes(x = mean_pred, y = km_surv_5yr, label = decile)) +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", colour = "grey50", linewidth = 0.8) +
  geom_point(size = 4, colour = "#2166AC") +
  geom_text(vjust = -0.8, size = 3.2, colour = "grey30") +
  geom_smooth(method = "lm", se = FALSE,
              colour = "#D6604D", linewidth = 0.8, linetype = "solid") +
  scale_x_continuous(limits = c(0.4, 1), labels = scales::percent_format()) +
  scale_y_continuous(limits = c(0.4, 1), labels = scales::percent_format()) +
  labs(
    title    = "Calibration Plot — 5-Year Survival (Decile Groups)",
    subtitle = "Points = risk deciles | Dashed = perfect calibration | Red = fitted line",
    x        = "Mean predicted 5-year survival",
    y        = "Observed KM 5-year survival"
  )

ggsave(here("output", "figures", "05_calibration.png"),
       p_calib, width = 7, height = 6, dpi = 300)
message("  ✔  Figure: Calibration plot")


# -----------------------------------------------------------------------------
# 11. Save final fitted workflow
# -----------------------------------------------------------------------------
saveRDS(cox_final_fit, here("data", "processed", "cox_final_workflow.rds"))

message("\n── 05  Summary ─────────────────────────────────────────────────────────")
message(sprintf("  Full-data C-statistic : %.3f (SE = %.4f)",
                base_cstat[["C"]], base_cstat[["se(C)"]]))
message(sprintf("  10-fold CV C-statistic: %.3f (SE = %.4f)",
                cv_metrics$mean, cv_metrics$std_err))
message(sprintf("  Optimism bias         : %.3f",
                base_cstat[["C"]] - cv_metrics$mean))

message("\n✔  Saved: data/processed/cox_final_workflow.rds")
message("\n── 05  Complete ────────────────────────────────────────────────────────")
