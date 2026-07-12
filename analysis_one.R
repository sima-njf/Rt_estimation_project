# analysis_one.R
# Processes ONE (model_type, R0) combination through Cori + Wallinga-Teunis
# for both transitions, and saves results. Designed to run as one Slurm array task.
# ---------------------------------------------------------------------------

# ==== ABSOLUTE PATHS (edit only if your layout changes) ====
CODE_DIR <- "/uufs/chpc.utah.edu/common/home/u1418987/sima/Rt_estimation_project"
BASE     <- "/scratch/general/vast/u1418987"
SAVED    <- file.path(BASE, "saved_data")               # where the *_nsim_1000 sim files live
RESULTS  <- file.path(BASE, "cori_wallinga_results")    # where CW results get written
NSIM     <- 1000

# ---------------------------------------------------------------------------
# Process a single combination. Everything the worker needs is loaded INSIDE
# the function, because each Slurm array task runs in a fresh R session on a
# node that has NOT sourced this driver or attached these packages.
# ---------------------------------------------------------------------------
process_one_combo <- function(combo) {
  
  # --- Attach packages on the node ---
  suppressMessages({
    library(EpiEstim)
    library(dplyr)
    library(tidyr)
    library(purrr)
  })
  
  # --- Source dependency functions on the node (ABSOLUTE paths) ---
  # data001.R  -> provides load_seir_data()
  # <cori file> -> provides calculate_serial_intervals(),
  #                apply_both_methods_all_sims(), etc.
  source(file.path(CODE_DIR, "data001.R"))
  source(file.path(CODE_DIR, "cori_wallinga_all_sims003.R"))  # <-- CONFIRM this filename
  
  model_type  <- combo$model_type
  R0_val      <- combo$R0
  transitions <- c("susceptible_to_exposed", "exposed_to_infected")
  
  # Make sure the results directory exists (safe if run concurrently)
  if (!dir.exists(RESULTS)) {
    dir.create(RESULTS, recursive = TRUE, showWarnings = FALSE)
  }
  
  # Build the simulation-file prefix. NOTE on R0 formatting:
  # paste0("R0_", 2.0) -> "R0_2", which matches your existing filenames
  # (full_R0_2_..., full_R0_5_...). Do not "fix" this to 2.0 on one side only.
  file_prefix <- file.path(
    SAVED,
    paste0(model_type, "_R0_", R0_val, "_n_1e+05_nsim_", NSIM)
  )
  
  # Guard: if the sim data isn't there, report cleanly instead of crashing.
  if (!file.exists(paste0(file_prefix, "_metadata.rds"))) {
    return(paste("MISSING:", file_prefix))
  }
  
  # Load the simulation data for this combo
  saved_data <- load_seir_data(file_prefix)
  
  out <- character(0)
  
  for (transition in transitions) {
    
    results <- tryCatch(
      apply_both_methods_all_sims(
        saved_data,
        how_infected = transition,
        start_day    = 0,
        verbose      = TRUE
      ),
      error = function(e) {
        message("ERROR in ", model_type, " R0=", R0_val,
                " ", transition, ": ", e$message)
        NULL
      }
    )
    
    if (!is.null(results) && nrow(results) > 0) {
      result_file <- file.path(
        RESULTS,
        paste0(model_type, "_R0_", R0_val, "_", transition, ".rds")
      )
      saveRDS(results, result_file)
      write.csv(results, sub("\\.rds$", ".csv", result_file), row.names = FALSE)
      
      out <- c(out, paste0(
        "OK: ", basename(result_file),
        " (", length(unique(results$sim_id)), " sims)"
      ))
    } else {
      out <- c(out, paste0("EMPTY: ", model_type, " R0=", R0_val, " ", transition))
    }
  }
  
  paste(out, collapse = " | ")
}