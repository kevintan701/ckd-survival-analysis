# =============================================================================
# 01_data_download.R
# CKD Progression & Mortality Risk Factors — NHANES 2017–2023
# Author: Yuntao (Kevin) Tan | tyuntao@umich.edu | December 2025
#
# Purpose: Download all required NHANES XPT files and the NCHS Public Use
#          Linked Mortality File across three survey cycles:
#            - 2017–2018  (cycle J)
#            - 2019–2020  (cycle K, pre-pandemic)
#            - 2021–2023  (cycle L, redesigned)
#
# Output:  data/raw/<CYCLE>/<COMPONENT>.XPT
#          data/raw/mortality/NHANES_2019_MORT_2019_PUBLIC.dat  (linked file)
# =============================================================================


# -----------------------------------------------------------------------------
# 0. Package setup
# -----------------------------------------------------------------------------
# nhanesA  — programmatic NHANES XPT download (no manual CDC navigation)
# tidyverse — data wrangling utilities used throughout pipeline
# here     — project-relative file paths (reproducibility best practice)

if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(nhanesA, tidyverse, here, haven)

# here::here() always resolves relative to the project root (.Rproj or .here),
# so paths are portable across machines — important for reproducible research.
here::i_am("R/01_data_download.R")

# Create directory structure if it does not already exist
dirs <- c(
  here("data", "raw", "J"),          # 2017–2018
  here("data", "raw", "K"),          # 2019–2020
  here("data", "raw", "L"),          # 2021–2023
  here("data", "raw", "mortality"),  # NDI linked mortality
  here("data", "processed"),
  here("output", "figures"),
  here("output", "tables")
)
purrr::walk(dirs, dir.create, recursive = TRUE, showWarnings = FALSE)

message("✔  Directory structure initialised")


# -----------------------------------------------------------------------------
# 1. Define components to download
# -----------------------------------------------------------------------------
# Each NHANES component is a separate XPT file identified by a short code.
# We download the same component across all three cycles for harmonisation.
#
# Component reference:
#   DEMO   — Demographics (age, sex, race/ethnicity, poverty-income ratio)
#   BMX    — Body Measures (BMI, weight, height, waist circumference)
#   GFR    — Kidney Health (serum creatinine → used to compute eGFR)
#   ALB_CR — Albumin & Creatinine — Urine (UACR, albuminuria staging)
#   DIQ    — Diabetes Questionnaire (diagnosed diabetes, insulin use)
#   BPQ    — Blood Pressure Questionnaire (hypertension diagnosis, meds)
#   BPXO   — Oscillometric Blood Pressure (measured SBP/DBP)
#   PAQ    — Physical Activity Questionnaire (MET-min/week)
#   SMQ    — Smoking Questionnaire (current/former/never)
#   GHB    — Glycohemoglobin (HbA1c — diabetes biochemical definition)

components <- list(
  DEMO   = "Demographics",
  BMX    = "Examination",
  BIOPRO = "Laboratory",   # Serum creatinine in 2017–2020
  GFR    = "Laboratory",   # eGFR in 2021–2023 (redesigned)
  ALB_CR = "Laboratory",
  DIQ    = "Questionnaire",
  BPQ    = "Questionnaire",
  BPXO   = "Examination",
  PAQ    = "Questionnaire",
  SMQ    = "Questionnaire",
  GHB    = "Laboratory"
)

# NHANES cycle suffix codes used by nhanesA
cycles <- c(J = "2017-2018", K = "2019-2020", L = "2021-2023")


# -----------------------------------------------------------------------------
# 2. Download helper function
# -----------------------------------------------------------------------------
# nhanes() — returns an R data frame directly from the CDC server
# haven::write_xpt() — saves the dataframe as an XPT file on disk
#
# We wrap this in a robust helper that:
#   (a) skips the download if the file already exists (idempotent runs)
#   (b) catches errors gracefully so one missing component doesn't abort all
#   (c) logs progress clearly

download_component <- function(component, cycle_code, cycle_label, out_dir) {

  # nhanesA naming convention: e.g. "DEMO_J", "GHB_L"
  table_name <- paste0(component, "_", cycle_code)
  out_path   <- file.path(out_dir, paste0(table_name, ".XPT"))

  # Idempotency check — skip if file already downloaded
  if (file.exists(out_path)) {
    message("  ↷  ", table_name, " already exists — skipping")
    return(invisible(NULL))
  }

  tryCatch({
    message("  ↓  Downloading ", table_name, " (", cycle_label, ") ...")
    df <- nhanes(table_name)           # fetch from CDC as data frame
    haven::write_xpt(df, out_path)     # save as XPT using haven
    message("  ✔  ", table_name, " saved  [", nrow(df), " rows]")
  }, error = function(e) {
    # Some components changed names across cycles (handled in 02_cleaning.R)
    message("  ✗  ", table_name, " — not available: ", conditionMessage(e))
  })
}


# -----------------------------------------------------------------------------
# 3. Execute downloads across all cycles and components
# -----------------------------------------------------------------------------
message("\n── Downloading NHANES components ──────────────────────────────────────")

for (code in names(cycles)) {
  label   <- cycles[[code]]
  out_dir <- here("data", "raw", code)
  message("\nCycle ", code, " (", label, ")")

  for (comp in names(components)) {
    download_component(comp, code, label, out_dir)
  }
}


# -----------------------------------------------------------------------------
# 4. Download NCHS Public Use Linked Mortality File
# -----------------------------------------------------------------------------
# The mortality linkage file is NOT distributed as an XPT — it is a fixed-
# width ASCII (.dat) file released by NCHS. It links NHANES SEQN identifiers
# to National Death Index (NDI) records, providing:
#   - MORTSTAT  : vital status (0 = assumed alive, 1 = deceased)
#   - PERMTH_EXM: person-months of follow-up from NHANES exam date
#   - UCOD_LEADING: underlying cause-of-death (ICD-10 category)
#
# This is what transforms NHANES from a cross-sectional survey into a
# longitudinal survival dataset — the key methodological feature of this
# project vs. prior cross-sectional analyses.
#
# Source: https://www.cdc.gov/nchs/data-linkage/mortality-public.htm

mort_url  <- paste0(
  "https://ftp.cdc.gov/pub/Health_Statistics/NCHS/",
  "datalinkage/linked_mortality/",
  "NHANES_2017_2018_MORT_2019_PUBLIC.dat"
)
mort_path <- here("data", "raw", "mortality", "NHANES_2017_2018_MORT_2019_PUBLIC.dat")

if (!file.exists(mort_path)) {
  message("\n── Downloading NCHS Linked Mortality File ──────────────────────────────")
  tryCatch({
    download.file(mort_url, destfile = mort_path, mode = "wb")
    message("✔  Mortality file saved: ", mort_path)
  }, error = function(e) {
    message("✗  Mortality download failed — check URL or download manually:")
    message("   ", mort_url)
  })
} else {
  message("\n↷  Mortality file already exists — skipping")
}


# -----------------------------------------------------------------------------
# 5. Verify downloads — inventory report
# -----------------------------------------------------------------------------
message("\n── Download Inventory ──────────────────────────────────────────────────")

all_files <- list.files(here("data", "raw"), recursive = TRUE,
                        full.names = TRUE, pattern = "\\.(XPT|dat)$")

inventory <- tibble(
  file     = basename(all_files),
  cycle    = str_extract(all_files, "(?<=/raw/)[A-Z]+"),
  size_kb  = round(file.size(all_files) / 1024, 1),
  path     = all_files
)

print(inventory, n = Inf)

message(
  "\n✔  Download complete — ",
  nrow(inventory), " files / ",
  round(sum(inventory$size_kb) / 1024, 1), " MB total"
)

# Save inventory to output for documentation
readr::write_csv(inventory, here("output", "tables", "01_download_inventory.csv"))
