# MODIFIED CODE - Using serial intervals from transmission data instead of generation time
# ===================================================================================

library(EpiEstim)
library(dplyr)
library(tidyr)
library(purrr)

# Function to calculate serial intervals from transmission data
calculate_serial_intervals <- function(transmission_data) {
  transmission_clean <- transmission_data %>%
    filter(source != -1) %>%  # Remove initial infections
    mutate(serial_interval = date - source_exposure_date) %>%
    filter(serial_interval > 0) %>%
    group_by(sim_num) %>%
    summarise(
      mean_si = mean(serial_interval, na.rm = TRUE),
      sd_si = sd(serial_interval, na.rm = TRUE),
      n_transmissions = n(),
      .groups = 'drop'
    ) %>%
    filter(!is.na(mean_si), !is.na(sd_si), sd_si > 0.1, mean_si > 0.5, mean_si < 20)
  
  return(transmission_clean)
}

# Function to apply Cori method to all simulations - USING SERIAL INTERVALS
apply_cori_to_all_sims <- function(saved_data, how_infected = "susceptible_to_exposed", 
                                   start_day = 0, verbose = TRUE) {
  
  # Extract data
  df <- saved_data$transitions
  transmission_data <- saved_data$transmission
  metadata <- saved_data$metadata
  
  # Calculate serial intervals for all simulations
  serial_intervals <- calculate_serial_intervals(transmission_data)
  
  if (verbose) {
    cat("Calculated serial intervals for", nrow(serial_intervals), "simulations\n")
    cat("Mean serial interval across sims:", round(mean(serial_intervals$mean_si), 2), "days\n")
    cat("Mean SD across sims:", round(mean(serial_intervals$sd_si), 2), "days\n")
  }
  
  # Check available columns
  if (verbose) {
    cat("Available transition columns:\n")
    transition_cols <- grep("_to_", colnames(df), value = TRUE)
    cat(paste(transition_cols, collapse = ", "), "\n")
  }
  
  # Check if the requested column exists
  if (!how_infected %in% colnames(df)) {
    cat("ERROR: Column", how_infected, "not found!\n")
    cat("Available columns:", paste(colnames(df), collapse = ", "), "\n")
    return(NULL)
  }
  
  # Get unique simulation IDs that have serial interval data
  sim_ids <- intersect(unique(df$id), serial_intervals$sim_num)
  
  if (verbose) {
    cat("Applying Cori method to", length(sim_ids), "simulations with serial interval data...\n")
    cat("Using transition:", how_infected, "\n")
  }
  
  # Function to apply Cori to a single simulation
  apply_cori_single <- function(sim_id) {
    
    if (verbose && sim_id %% 10 == 0) {
      cat("Processing simulation", sim_id, "of", max(sim_ids), "\n")
    }
    
    tryCatch({
      # Get simulation data
      sim_data <- df[df$id == sim_id & df$date >= start_day, ]
      
      if (nrow(sim_data) == 0) {
        return(NULL)
      }
      
      # Get serial intervals for this simulation
      si_data <- serial_intervals %>% 
        filter(sim_num == sim_id)
      
      if (nrow(si_data) == 0) {
        return(NULL)
      }
      
      SI_mean <- si_data$mean_si
      SI_sd <- si_data$sd_si
      
      # Skip if SD is too small or mean is unrealistic
      if (is.na(SI_sd) || SI_sd < 0.1 || SI_mean < 0.5 || SI_mean > 20) {
        return(NULL)
      }
      
      # Get incidence data
      incidence <- sim_data[[how_infected]]
      
      # Skip if all zeros or too few data points
      if (sum(incidence, na.rm = TRUE) == 0 || length(incidence) < 7) {
        return(NULL)
      }
      
      # Apply Cori method using serial intervals
      result2 <- estimate_R(
        incidence,
        method = "parametric_si",
        config = make_config(mean_si = SI_mean, std_si = SI_sd)
      )
      
      # Extract R_t values
      rt_values2 <- result2$R[, "Median(R)"]
      rt_values2_ul <- result2$R[, "Quantile.0.025(R)"]
      rt_values2_ll <- result2$R[, "Quantile.0.975(R)"]
      date <- result2$R[, "t_start"]
      
      # Create results dataframe
      cori_results <- data.frame(
        sim_id = sim_id,
        date = date,
        median_rt = rt_values2,
        upper = rt_values2_ll,
        lower = rt_values2_ul,
        mean_si = SI_mean,
        std_si = SI_sd,
        n_transmissions = si_data$n_transmissions,
        method = "Cori"
      )
      
      return(cori_results)
      
    }, error = function(e) {
      if (verbose) cat("Error in simulation", sim_id, ":", e$message, "\n")
      return(NULL)
    })
  }
  
  # Apply to all simulations
  all_results <- map_dfr(sim_ids, apply_cori_single)
  
  if (verbose) {
    cat("Cori method completed!\n")
    if (!is.null(all_results) && nrow(all_results) > 0) {
      cat("Successful simulations:", length(unique(all_results$sim_id)), "out of", length(sim_ids), "\n")
    } else {
      cat("No successful results!\n")
    }
  }
  
  return(all_results)
}

# Function to apply Wallinga-Teunis method to all simulations - USING SERIAL INTERVALS
apply_wallinga_to_all_sims <- function(saved_data, how_infected = "susceptible_to_exposed", 
                                       start_day = 10, verbose = TRUE) {
  
  # Extract data
  df <- saved_data$transitions
  transmission_data <- saved_data$transmission
  metadata <- saved_data$metadata
  
  # Calculate serial intervals for all simulations
  serial_intervals <- calculate_serial_intervals(transmission_data)
  
  if (verbose) {
    cat("Calculated serial intervals for", nrow(serial_intervals), "simulations\n")
  }
  
  # Check available columns
  if (verbose) {
    cat("Available transition columns:\n")
    transition_cols <- grep("_to_", colnames(df), value = TRUE)
    cat(paste(transition_cols, collapse = ", "), "\n")
  }
  
  # Check if the requested column exists
  if (!how_infected %in% colnames(df)) {
    cat("ERROR: Column", how_infected, "not found!\n")
    cat("Available columns:", paste(colnames(df), collapse = ", "), "\n")
    return(NULL)
  }
  
  # Get unique simulation IDs that have serial interval data
  sim_ids <- intersect(unique(df$id), serial_intervals$sim_num)
  
  if (verbose) {
    cat("Applying Wallinga-Teunis method to", length(sim_ids), "simulations...\n")
    cat("Using transition:", how_infected, "\n")
  }
  
  # Function to apply Wallinga-Teunis to a single simulation
  apply_wallinga_single <- function(sim_id) {
    
    if (verbose && sim_id %% 10 == 0) {
      cat("Processing simulation", sim_id, "of", max(sim_ids), "\n")
    }
    
    tryCatch({
      # Get simulation data
      sim_data <- df[df$id == sim_id & df$date >= start_day, ]
      
      if (nrow(sim_data) == 0) {
        return(NULL)
      }
      
      # Get serial intervals for this simulation
      si_data <- serial_intervals %>% 
        filter(sim_num == sim_id)
      
      if (nrow(si_data) == 0) {
        return(NULL)
      }
      
      SI_mean <- si_data$mean_si
      SI_sd <- si_data$sd_si
      
      # Skip if SD is too small or mean is unrealistic
      if (is.na(SI_sd) || SI_sd < 0.1 || SI_mean < 0.5 || SI_mean > 20) {
        return(NULL)
      }
      
      # Get incidence data
      incidence <- sim_data[[how_infected]]
      
      # Skip if all zeros or too few data points
      if (sum(incidence, na.rm = TRUE) == 0 || length(incidence) < 7) {
        return(NULL)
      }
      
      # Apply Wallinga-Teunis method using serial intervals
      config <- make_config(
        incid = incidence,
        method = "parametric_si",
        mean_si = SI_mean,
        std_si = SI_sd
      )
      
      result <- wallinga_teunis(
        incid = incidence,
        method = "parametric_si",
        config = config
      )
      
      # Extract results
      rt_values_wt <- result$R[,"Mean(R)"]
      rt_values_wt_ul <- result$R[,"Quantile.0.025(R)"]
      rt_values_wt_ll <- result$R[,"Quantile.0.975(R)"]
      date <- result$R[,"t_start"]
      
      # Create results dataframe
      wt_results <- data.frame(
        sim_id = sim_id,
        date = date,
        median_rt = rt_values_wt,
        upper = rt_values_wt_ll,
        lower = rt_values_wt_ul,
        mean_si = SI_mean,
        std_si = SI_sd,
        n_transmissions = si_data$n_transmissions,
        method = "Wallinga-Teunis"
      )
      
      return(wt_results)
      
    }, error = function(e) {
      if (verbose) cat("Error in simulation", sim_id, ":", e$message, "\n")
      return(NULL)
    })
  }
  
  # Apply to all simulations
  all_results <- map_dfr(sim_ids, apply_wallinga_single)
  
  if (verbose) {
    cat("Wallinga-Teunis method completed!\n")
    if (!is.null(all_results) && nrow(all_results) > 0) {
      cat("Successful simulations:", length(unique(all_results$sim_id)), "out of", length(sim_ids), "\n")
    } else {
      cat("No successful results!\n")
    }
  }
  
  return(all_results)
}

# Function to apply both methods and combine results - USING SERIAL INTERVALS
apply_both_methods_all_sims <- function(saved_data, how_infected = "susceptible_to_exposed", 
                                        start_day =10, verbose = TRUE) {
  
  if (verbose) {
    cat("Applying both Cori and Wallinga-Teunis methods using serial intervals...\n")
    cat("Dataset info:\n")
    cat("- R0:", saved_data$metadata$R0, "\n")
    cat("- Model type:", saved_data$metadata$model_type, "\n")
    cat("- Number of simulations:", saved_data$metadata$nsim, "\n")
    cat("- Requested transition:", how_infected, "\n")
  }
  
  # Apply both methods
  cori_results <- apply_cori_to_all_sims(saved_data, how_infected, start_day, verbose)
  wt_results <- apply_wallinga_to_all_sims(saved_data, how_infected, start_day, verbose)
  
  # Check if we have any results
  if ((is.null(cori_results) || nrow(cori_results) == 0) && 
      (is.null(wt_results) || nrow(wt_results) == 0)) {
    cat("No results from either method!\n")
    return(NULL)
  }
  
  # Combine results
  all_results <- bind_rows(cori_results, wt_results)
  
  if (nrow(all_results) == 0) {
    cat("No successful results to combine!\n")
    return(NULL)
  }
  
  # Add metadata
  all_results$R0_true <- saved_data$metadata$R0
  all_results$model_type <- saved_data$metadata$model_type
  all_results$how_infected <- how_infected
  
  if (verbose) {
    cat("\nCombined results summary:\n")
    cat("- Total rows:", nrow(all_results), "\n")
    cat("- Methods:", paste(unique(all_results$method), collapse = ", "), "\n")
    cat("- Date range:", min(all_results$date), "to", max(all_results$date), "\n")
    cat("- Simulation IDs:", length(unique(all_results$sim_id)), "\n")
    cat("- Serial interval range:", round(range(all_results$mean_si), 2), "days\n")
  }
  
  return(all_results)
}

# Test function - USING SERIAL INTERVALS
test_single_simulation <- function(saved_data, sim_id = 1, how_infected = "susceptible_to_exposed") {
  
  cat("Testing single simulation", sim_id, "with transition", how_infected, "\n")
  
  # Check data
  df <- saved_data$transitions
  sim_data <- df[df$id == sim_id, ]
  
  cat("Simulation data rows:", nrow(sim_data), "\n")
  cat("Date range:", min(sim_data$date), "to", max(sim_data$date), "\n")
  
  # Check incidence
  if (how_infected %in% colnames(sim_data)) {
    incidence <- sim_data[[how_infected]]
    cat("Incidence - Total:", sum(incidence), "Max:", max(incidence), "Non-zero days:", sum(incidence > 0), "\n")
    
    # Show first few days
    cat("First 10 days of incidence:\n")
    print(head(incidence, 10))
  } else {
    cat("ERROR: Column", how_infected, "not found!\n")
    cat("Available columns:", paste(colnames(sim_data), collapse = ", "), "\n")
  }
  
  # Check serial intervals
  serial_intervals <- calculate_serial_intervals(saved_data$transmission)
  si_data <- serial_intervals %>% filter(sim_num == sim_id)
  
  if (nrow(si_data) > 0) {
    cat("Serial intervals - Mean:", round(si_data$mean_si, 2), "SD:", round(si_data$sd_si, 2), 
        "N transmissions:", si_data$n_transmissions, "\n")
  } else {
    cat("No serial interval data found for simulation", sim_id, "\n")
  }
  
  return(invisible(NULL))
}

# Function to view serial interval summary across all simulations
view_serial_interval_summary <- function(saved_data) {
  serial_intervals <- calculate_serial_intervals(saved_data$transmission)
  
  cat("Serial Interval Summary Across All Simulations:\n")
  cat("Number of simulations with SI data:", nrow(serial_intervals), "\n")
  cat("Mean SI across sims:", round(mean(serial_intervals$mean_si), 2), "±", round(sd(serial_intervals$mean_si), 2), "days\n")
  cat("Mean SD across sims:", round(mean(serial_intervals$sd_si), 2), "±", round(sd(serial_intervals$sd_si), 2), "days\n")
  cat("Range of mean SI:", round(range(serial_intervals$mean_si), 2), "days\n")
  cat("Range of SD SI:", round(range(serial_intervals$sd_si), 2), "days\n")
  
  return(serial_intervals)
}
