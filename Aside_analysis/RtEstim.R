# =============================================================================
# RtEstim – FULL DATA, ALL SCENARIOS (mclapply + proper seed filtering)
# =============================================================================

library(rtestim)
library(dplyr)
library(parallel)

source("data001.R")   # load_seir_data()

# =============================================================================
# Process ONE simulation
# =============================================================================
process_single_sim_rtestim <- function(
    sim_id,
    transitions_data,
    generation_data,
    transmission_data,
    how_infected,
    verbose_debug = FALSE
) {
  
  tryCatch({
    
    # ---- Filter to this simulation ------------------------------------------
    sim_transitions  <- transitions_data  %>% filter(id == sim_id)
    sim_generation   <- generation_data   %>% filter(sim_num == sim_id)
    sim_transmission <- transmission_data %>% filter(sim_num == sim_id)
    
    # ---- Incidence from $transitions (already aggregated) -------------------
    if (how_infected == "susceptible_to_exposed") {
      incidence <- sim_transitions$susceptible_to_exposed
    } else if (how_infected == "exposed_to_infected") {
      incidence <- sim_transitions$exposed_to_infected
    } else {
      return(NULL)
    }
    
    days <- sim_transitions$date
    
    # For S->E: day 0 includes seeds, so remove seed count from day 0
    if (how_infected == "susceptible_to_exposed") {
      n_seeds <- sum(sim_transmission$source == -1)
      incidence[1] <- incidence[1] - n_seeds
      if (incidence[1] < 0) incidence[1] <- 0
    }
    
    # Trim leading zeros (rtestim needs first value > 0)
    first_nz <- which(incidence > 0)[1]
    if (is.na(first_nz)) return(NULL)
    
    # Trim trailing zeros (zeros at the tail cause NAs in cv_estimate_rt)
    last_nz <- max(which(incidence > 0))
    
    incidence <- incidence[first_nz:last_nz]
    days      <- days[first_nz:last_nz]
    
    if (length(incidence) < 10) return(NULL)
    
    # ---- Generation interval PMF (remove seeds) ----------------------------
    seed_ids <- sim_transmission %>%
      filter(source == -1) %>%
      pull(target)
    
    sim_gen_clean <- sim_generation %>%
      filter(!(source %in% seed_ids))
    
    gen_times <- sim_gen_clean$gentime
    gen_times <- gen_times[!is.na(gen_times) & gen_times > 0]
    
    if (length(gen_times) < 10) return(NULL)
    
    max_gen    <- max(gen_times)
    gen_counts <- table(factor(gen_times, levels = 0:max_gen))
    gen_int    <- as.numeric(gen_counts) / sum(gen_counts)
    
    if (verbose_debug) {
      cat(sprintf(
        "Sim %d | days=%d | total_cases=%d | GI_len=%d | mean_GI=%.1f\n",
        sim_id, length(incidence), sum(incidence),
        length(gen_int), mean(gen_times)
      ))
    }
    
    # ---- Rt estimation (cross-validated) ------------------------------------
    rtestim_fit <- suppressWarnings(cv_estimate_rt(
      observed_counts = incidence,
      delay_distn     = gen_int,
      nsol            = 20
    ))
    
    # ---- Confidence bands ---------------------------------------------------
    # lambda.1se can be -Inf or NA in some cases; fall back to lambda.min
    chosen_lambda <- rtestim_fit$lambda.1se
    if (is.na(chosen_lambda) || is.infinite(chosen_lambda) || chosen_lambda < 0) {
      chosen_lambda <- rtestim_fit$lambda.min
    }
    
    rtestim_cb <- suppressWarnings(confband(rtestim_fit, lambda = chosen_lambda))
    
    # ---- Output (trim early unreliable days) --------------------------------
    # Use 95th percentile of generation time as burn-in, not the max
    # (one outlier gen time of 40 shouldn't wipe out 40 days of estimates)
    burn_in <- ceiling(quantile(gen_times, 0.95))
    
    out <- data.frame(
      sim_id    = sim_id,
      method    = "RtEstim",
      date      = days,
      median_rt = rtestim_cb$fit,
      q025_rt   = rtestim_cb$`2.5%`,
      q975_rt   = rtestim_cb$`97.5%`
    )
    
    out <- out %>% filter(date > burn_in)
    out
    
  }, error = function(e) {
    if (verbose_debug) {
      cat(sprintf("Sim %d ERROR: %s\n", sim_id, e$message))
    }
    NULL
  })
}

# =============================================================================
# Apply to ALL simulations — parallel with mclapply
# =============================================================================
apply_rtestim_parallel <- function(
    saved_data,
    how_infected,
    ncores = ncores,
    verbose = TRUE
) {
  
  sim_nums <- unique(saved_data$transitions$id)
  
  if (verbose) {
    cat("Running RtEstim for", length(sim_nums),
        "simulations on", ncores, "cores\n")
  }
  
  results <- parallel::mclapply(sim_nums, FUN = function(sid) {
    
    process_single_sim_rtestim(
      sim_id            = sid,
      transitions_data  = saved_data$transitions,
      generation_data   = saved_data$generation,
      transmission_data = saved_data$transmission,
      how_infected      = how_infected,
      verbose_debug     = FALSE
    )
    
  }, mc.cores = ncores)
  
  # Remove NULLs and failed results (mclapply returns error objects on failure)
  is_valid <- sapply(results, function(x) is.data.frame(x))
  n_errors <- sum(sapply(results, function(x) inherits(x, "try-error") || inherits(x, "error")))
  
  results <- results[is_valid]
  
  if (verbose) {
    cat("  Successful:", length(results), "/", length(sim_nums), "\n")
    if (n_errors > 0) cat("  Errors:", n_errors, "\n")
  }
  
  if (length(results) == 0) return(NULL)
  bind_rows(results)
}

# =============================================================================
# Run ALL MODEL TYPES × R0 × TRANSITIONS
# =============================================================================
process_rtestim_all <- function(ncores = parallel::detectCores() - 1) {
  
  R0_values   <- c(5)
  model_types <- c("full")
  transitions <- c( "exposed_to_infected")
  
  dir.create("rtestim_results", showWarnings = FALSE)
  
  for (model_type in model_types) {
    
    cat("\n==============================\n")
    cat("MODEL TYPE:", model_type, "\n")
    cat("==============================\n")
    
    for (R0_val in R0_values) {
      
      prefix <- paste0(
        "saved_data/",
        model_type, "_R0_", R0_val, "_n_1e+05_nsim_100"
      )
      
      if (!file.exists(paste0(prefix, "_metadata.rds"))) {
        cat("  skipping:", prefix, "\n")
        next
      }
      
      cat("\n=== R0 =", R0_val, "===\n")
      saved_data <- load_seir_data(prefix)
      
      for (transition in transitions) {
        
        cat("Transition:", transition, "\n")
        
        results <- apply_rtestim_parallel(
          saved_data,
          how_infected = transition,
          ncores       = ncores
        )
        
        if (!is.null(results)) {
          
          out_file <- paste0(
            "rtestim_results/",
            model_type, "_R0_", R0_val, "_", transition, ".rds"
          )
          
          saveRDS(results, out_file)
          write.csv(results, sub(".rds", ".csv", out_file), row.names = FALSE)
          
          cat("  saved:", basename(out_file), "\n")
        }
      }
    }
  }
  
  cat("\n=== RtEstim ALL SCENARIOS COMPLETE ===\n")
}

# =============================================================================
# RUN
# =============================================================================

process_rtestim_all(ncores = 9)