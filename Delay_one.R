# delay_one.R
# Runs the SI-based reporting-delays analysis for ONE (model_type, R0)
# combination (both transitions) and saves per-combo results + plots.
# Designed to run as a single Slurm array task.
# ---------------------------------------------------------------------------

# ==== ABSOLUTE PATHS (edit only if your layout changes) ====
CODE_DIR <- "/uufs/chpc.utah.edu/common/home/u1418987/sima/Rt_estimation_project"
BASE     <- "/scratch/general/vast/u1418987"
SAVED    <- file.path(BASE, "saved_data")               # *_nsim_1000 simulation files
RESULTS  <- file.path(BASE, "complete_delay_results")   # per-combo delay results
NSIM     <- 1000

# The file that defines run_complete_delay_analysis(), analyze_delays_all_sims(),
# apply_reporting_delays(), calculate_serial_intervals(), load_seir_data(), etc.
# CONFIRM this filename with:  ls ~/sima/Rt_estimation_project/*.R
DELAY_CODE_FILE <- "delay_analysis005.R"
# Delay scenarios and how many sims to use per combination.
# NOTE: max_sims samples this many of the available sims. Set to NSIM (1000)
# to use all of them, but each task gets MUCH heavier. 60 mirrors your current code.
DELAY_SCENARIOS <- list("medium_delay" = 0.3, "high_delay" = 0.5)
MAX_SIMS        <- 60
START_DAY       <- 0

process_one_delay_combo <- function(combo) {
  
  # --- Attach packages on the node ---
  suppressPackageStartupMessages({
    library(ggplot2); library(dplyr); library(ggpubr)
    library(EpiEstim); library(tidyr); library(purrr); library(readr)
  })
  
  # --- Source the delay-analysis functions on the node (ABSOLUTE path) ---
  # This file defines run_complete_delay_analysis() and all helpers.
  source(file.path(CODE_DIR, DELAY_CODE_FILE))
  
  model_type <- combo$model_type
  R0_val     <- combo$R0
  
  if (!dir.exists(RESULTS)) dir.create(RESULTS, recursive = TRUE, showWarnings = FALSE)
  
  # run_complete_delay_analysis() loops internally over R0/model/transition and
  # reads/writes RELATIVE paths ("saved_data/", "complete_delay_results/").
  # To keep everything on scratch and pointed at nsim_1000, we set the working
  # directory to BASE so those relative paths resolve under scratch.
  old_wd <- getwd()
  setwd(BASE)
  on.exit(setwd(old_wd), add = TRUE)
  
  # Run only this one combination (both transitions), all 1000-sim data.
  res <- tryCatch(
    run_complete_delay_analysis(
      R0_values             = R0_val,                       # single R0
      model_types           = model_type,                   # single model
      transitions           = c("susceptible_to_exposed",
                                "exposed_to_infected"),
      delay_scenarios       = DELAY_SCENARIOS,
      max_sims_per_scenario = MAX_SIMS,
      start_day             = START_DAY
    ),
    error = function(e) {
      message("ERROR in ", model_type, " R0=", R0_val, ": ", e$message)
      NULL
    }
  )
  
  if (is.null(res)) {
    return(paste0("FAILED: ", model_type, " R0=", R0_val))
  }
  paste0("OK: ", model_type, " R0=", R0_val,
         " (", nrow(res$results), " result rows)")
}