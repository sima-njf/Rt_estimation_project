# misspec_one.R
# Runs the SI-based underreporting-misspecification analysis for ONE
# (model_type, R0) combination (both transitions) and saves per-combo
# results + plots. Designed to run as a single Slurm array task.
# ---------------------------------------------------------------------------

# ==== ABSOLUTE PATHS (edit only if your layout changes) ====
CODE_DIR <- "/uufs/chpc.utah.edu/common/home/u1418987/sima/Rt_estimation_project"
BASE     <- "/scratch/general/vast/u1418987"
SAVED    <- file.path(BASE, "saved_data")               # *_nsim_1000 simulation files
RESULTS  <- file.path(BASE, "complete_misspec_results") # per-combo misspec results
NSIM     <- 1000

# The file that defines run_complete_misspec_analysis(), analyze_misspec_all_sims(),
# apply_underreporting(), calculate_serial_intervals(), load_seir_data(), etc.
# CONFIRM this filename with:  ls ~/sima/Rt_estimation_project/*.R
MISSPEC_CODE_FILE <- "misespecification_analysis006.R"   # <-- CONFIRM / EDIT

# Analysis knobs
REPORTING_RATE <- 0.7
MAX_SIMS       <- 30      # samples this many sims; raise toward 1000 to use all
START_DAY      <- 0
RNG_SEED       <- 12345

process_one_misspec_combo <- function(combo) {
  
  # --- Attach packages on the node ---
  suppressPackageStartupMessages({
    library(ggplot2); library(dplyr); library(ggpubr)
    library(EpiEstim); library(tidyr); library(purrr); library(readr)
  })
  
  # --- Source the misspec-analysis functions on the node (ABSOLUTE path) ---
  source(file.path(CODE_DIR, MISSPEC_CODE_FILE))
  
  model_type <- combo$model_type
  R0_val     <- combo$R0
  
  if (!dir.exists(RESULTS)) dir.create(RESULTS, recursive = TRUE, showWarnings = FALSE)
  
  # run_complete_misspec_analysis() reads/writes RELATIVE paths ("saved_data/",
  # "complete_misspec_results/"). Set the working directory to BASE so those
  # relative paths resolve under scratch.
  old_wd <- getwd()
  setwd(BASE)
  on.exit(setwd(old_wd), add = TRUE)
  
  # Run only this one combination (both transitions), all 1000-sim data.
  res <- tryCatch(
    run_complete_misspec_analysis(
      R0_values             = R0_val,                       # single R0
      model_types           = model_type,                   # single model
      transitions           = c("susceptible_to_exposed",
                                "exposed_to_infected"),
      reporting_rate        = REPORTING_RATE,
      max_sims_per_scenario = MAX_SIMS,
      start_day             = START_DAY,
      results_dir           = RESULTS,                       # absolute scratch dir
      rng_seed              = RNG_SEED
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