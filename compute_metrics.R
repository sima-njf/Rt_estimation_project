# ================================================================================
# ACCURACY & COVERAGE METRICS for Rt estimators vs ABM ground truth
# Days 12-75. True Rt = ABM mean_rt (the black line in the figures).
#
# For each (model, transition, R0, method, scenario) computes:
#   MAE      = mean |estimate_mean - abm_mean|
#   RMSE     = sqrt(mean (estimate_mean - abm_mean)^2)
#   Coverage(band)   = fraction of days where abm_mean is within the estimator's
#                      2.5-97.5% ACROSS-SIMULATION band (what the plots show)
#   Coverage(perSim) = fraction of days where abm_mean is within the estimator's
#                      OWN reported CI columns (ci_lower/ci_upper), if present
#
# Scenarios: Real Data (baseline), 50% Delayed, 70% Reporting -- for both
# Cori and Wallinga-Teunis.
# ================================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
})

# -- CONFIG (same paths as the plotting scripts) ---------------------------------
DATE_MIN     <- 12
DATE_MAX     <- 75

R0_VALUES    <- c(1.5, 2.0, 3.0, 5.0)
MODEL_TYPES  <- c("full", "partial")
TRANSITIONS  <- c("susceptible_to_exposed", "exposed_to_infected")

BASELINE_DIR <- "cori_wallinga_results"
DELAY_DIR    <- "complete_delay_results"
MISSPEC_DIR  <- "complete_misspec_results"
ABM_DIR      <- "saved_data"

OUT_DIR      <- "metrics_results"

# -- LOADERS (mirror the plotting scripts) ---------------------------------------
safe_load <- function(f, ...) tryCatch(f(...), error = function(e) NULL)

load_baseline <- function(model_type, R0_val, transition) {
  fp <- file.path(BASELINE_DIR, paste0(model_type, "_R0_", R0_val, "_", transition, ".csv"))
  if (!file.exists(fp)) return(NULL)
  read_csv(fp, show_col_types = FALSE) %>%
    rename(ci_lower = lower, ci_upper = upper) %>%   # <-- ADD THIS
    mutate(scenario = "Real Data")
}

load_delay <- function(model_type, R0_val, transition) {
  fp <- file.path(DELAY_DIR, paste0("delay_results_", model_type, "_R0_", R0_val, "_", transition, ".csv"))
  if (!file.exists(fp)) return(NULL)
  read_csv(fp, show_col_types = FALSE) %>%
    filter(delay_scenario == "high_delay") %>%
    mutate(scenario = "50% Delayed")
}

load_misspec <- function(model_type, R0_val, transition) {
  fp <- file.path(MISSPEC_DIR, paste0("misspec_results_", model_type, "_R0_", R0_val, "_", transition, ".csv"))
  if (!file.exists(fp)) return(NULL)
  read_csv(fp, show_col_types = FALSE) %>% mutate(scenario = "70% Reporting")
}

load_abm <- function(model_type, R0_val) {
  fp <- file.path(ABM_DIR, paste0(model_type, "_R0_", R0_val, "_n_1e+05_nsim_1000_rt_ci.rds"))
  if (!file.exists(fp)) {
    csvp <- gsub("\\.rds$", ".csv", fp)
    if (file.exists(csvp)) return(read_csv(csvp, show_col_types = FALSE))
    return(NULL)
  }
  readRDS(fp)
}

# -- DIAGNOSTIC: print the columns of one file of each type ----------------------
inspect_columns <- function() {
  cat("=== Column inspection (first available file of each type) ===\n\n")
  b <- safe_load(load_baseline, "full", 1.5, "exposed_to_infected")
  d <- safe_load(load_delay,    "full", 1.5, "exposed_to_infected")
  m <- safe_load(load_misspec,  "full", 1.5, "exposed_to_infected")
  a <- safe_load(load_abm,      "full", 1.5)
  if (!is.null(b)) cat("BASELINE cols:\n  ", paste(names(b), collapse = ", "), "\n\n")
  if (!is.null(d)) cat("DELAY cols:\n  ",    paste(names(d), collapse = ", "), "\n\n")
  if (!is.null(m)) cat("MISSPEC cols:\n  ",  paste(names(m), collapse = ", "), "\n\n")
  if (!is.null(a)) cat("ABM cols:\n  ",      paste(names(a), collapse = ", "), "\n\n")
  invisible(list(baseline = b, delay = d, misspec = m, abm = a))
}

# -- BUILD ABM TRUTH (mean_rt per day) -------------------------------------------
get_abm_truth <- function(model_type, R0_val) {
  abm <- safe_load(load_abm, model_type, R0_val)
  if (is.null(abm)) return(NULL)
  abm %>%
    rename(date = source_exposure_date) %>%
    filter(date >= DATE_MIN, date <= DATE_MAX) %>%
    transmute(date, abm_mean = mean_rt)
}

# -- SUMMARISE ONE METHOD x SCENARIO ---------------------------------------------
# Produces, per day: estimate mean, and across-sim 2.5/97.5 band.
# Also, IF the raw file has per-row ci_lower/ci_upper, the per-sim CI averaged
# to a daily interval (mean of lowers, mean of uppers) is returned for the
# per-sim coverage variant.
summarise_method_scenario <- function(df) {
  has_persim_ci <- all(c("ci_lower", "ci_upper") %in% names(df))
  
  base <- df %>%
    filter(date >= DATE_MIN, date <= DATE_MAX) %>%
    group_by(method, date) %>%
    summarise(
      est_mean  = mean(median_rt, na.rm = TRUE),
      band_low  = quantile(median_rt, 0.025, na.rm = TRUE),
      band_high = quantile(median_rt, 0.975, na.rm = TRUE),
      .groups   = "drop"
    )
  
  if (has_persim_ci) {
    persim <- df %>%
      filter(date >= DATE_MIN, date <= DATE_MAX) %>%
      group_by(method, date) %>%
      summarise(
        persim_low  = mean(ci_lower, na.rm = TRUE),
        persim_high = mean(ci_upper, na.rm = TRUE),
        .groups     = "drop"
      )
    base <- base %>% left_join(persim, by = c("method", "date"))
  } else {
    base <- base %>% mutate(persim_low = NA_real_, persim_high = NA_real_)
  }
  base
}

# -- COMPUTE METRICS FOR ONE (model, transition) ---------------------------------
metrics_for_combo <- function(model_type, transition) {
  
  out <- list()
  
  for (R0_val in R0_VALUES) {
    
    abm_truth <- get_abm_truth(model_type, R0_val)
    if (is.null(abm_truth) || nrow(abm_truth) == 0) next
    
    baseline <- safe_load(load_baseline, model_type, R0_val, transition)
    delay    <- safe_load(load_delay,    model_type, R0_val, transition)
    misspec  <- safe_load(load_misspec,  model_type, R0_val, transition)
    
    scen_list <- list(baseline, delay, misspec)
    
    for (df in scen_list) {
      if (is.null(df) || nrow(df) == 0) next
      scen_name <- df$scenario[1]
      
      summ <- summarise_method_scenario(df) %>%
        left_join(abm_truth, by = "date") %>%
        filter(!is.na(abm_mean))
      
      if (nrow(summ) == 0) next
      
      m <- summ %>%
        group_by(method) %>%
        summarise(
          n_days        = n(),
          MAE           = mean(abs(est_mean - abm_mean), na.rm = TRUE),
          RMSE          = sqrt(mean((est_mean - abm_mean)^2, na.rm = TRUE)),
          mean_bias     = mean(est_mean - abm_mean, na.rm = TRUE),
          # Coverage against across-sim band (what plots show)
          coverage_band = mean(abm_mean >= band_low & abm_mean <= band_high, na.rm = TRUE),
          # Coverage against per-sim reported CI (NA if columns absent)
          coverage_persim = if (all(is.na(persim_low))) NA_real_
          else mean(abm_mean >= persim_low & abm_mean <= persim_high, na.rm = TRUE),
          .groups = "drop"
        ) %>%
        mutate(model = model_type, transition = transition,
               R0 = R0_val, scenario = scen_name, .before = 1)
      
      out[[paste(R0_val, scen_name, sep = "_")]] <- m
    }
  }
  
  bind_rows(out)
}

# -- MAIN ------------------------------------------------------------------------
compute_all_metrics <- function() {
  if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)
  
  cat("\n")
  info <- inspect_columns()
  has_ci <- !is.null(info$baseline) && all(c("ci_lower","ci_upper") %in% names(info$baseline))
  cat("Per-simulation CI columns present in baseline file: ", has_ci, "\n")
  cat("(If FALSE, coverage_persim will be NA; coverage_band is always computed.)\n\n")
  
  all_rows <- list()
  for (mt in MODEL_TYPES) {
    for (tr in TRANSITIONS) {
      cat("Computing:", mt, "-", tr, "\n")
      res <- metrics_for_combo(mt, tr)
      if (!is.null(res) && nrow(res) > 0) all_rows[[paste(mt, tr, sep = "_")]] <- res
    }
  }
  
  metrics <- bind_rows(all_rows) %>%
    arrange(model, transition, R0, method, scenario)
  
  # Round for readability
  metrics_out <- metrics %>%
    mutate(across(c(MAE, RMSE, mean_bias, coverage_band, coverage_persim),
                  ~ round(.x, 4)))
  
  out_csv <- file.path(OUT_DIR, "rt_metrics_days12-75.csv")
  write_csv(metrics_out, out_csv)
  
  cat("\n================ METRICS (days", DATE_MIN, "-", DATE_MAX, ") ================\n")
  print(as.data.frame(metrics_out), row.names = FALSE)
  cat("\nSaved:", out_csv, "\n")
  
  invisible(metrics_out)
}

if (sys.nframe() == 0) {
  compute_all_metrics()
}

# To run interactively:
#   CODE_DIR <- "/uufs/chpc.utah.edu/common/home/u1418987/sima/Rt_estimation_project"
#   setwd("/scratch/general/vast/u1418987")
#   source(file.path(CODE_DIR, "compute_metrics.R"))
#   compute_all_metrics()