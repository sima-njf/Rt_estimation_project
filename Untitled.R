# ================================================================================
# COMPLETE REPORTING DELAYS ANALYSIS (SI-BASED) - ALL R0 VALUES
# Uses serial intervals from transmission data (NOT generation times)
# Uniform(1,7) reporting delays; ABM = True baseline, Cori/Wallinga = Delayed reporting
# Plus: Comprehensive plots built FROM the saved delay-results outputs
# FIXED: Proper handling of delayed cases that extend beyond original window
# PLOTS SAVED TO: comparison_delays2
# ================================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(ggpubr)
  library(EpiEstim)
  library(tidyr)
  library(purrr)
  library(readr)
})

# ── HELPER FUNCTIONS ─────────────────────────────────────────────────────────────

load_seir_data <- function(file_prefix) {
  metadata_file     <- paste0(file_prefix, "_metadata.rds")
  transitions_file  <- paste0(file_prefix, "_transitions.rds")
  transmission_file <- paste0(file_prefix, "_transmission.rds")
  
  if (!file.exists(metadata_file))     stop("Metadata file not found: ", metadata_file)
  if (!file.exists(transitions_file))  stop("Transitions file not found: ", transitions_file)
  if (!file.exists(transmission_file)) stop("Transmission file not found: ", transmission_file)
  
  cat("Loading files:\n")
  cat("  • Metadata:",     basename(metadata_file), "\n")
  cat("  • Transitions:",  basename(transitions_file), "\n")
  cat("  • Transmission:", basename(transmission_file), "\n")
  
  metadata     <- readRDS(metadata_file)
  transitions  <- readRDS(transitions_file)
  transmission <- readRDS(transmission_file)
  
  if (!is.data.frame(transitions))  stop("Transitions data is not a data frame")
  if (!is.data.frame(transmission)) stop("Transmission data is not a data frame")
  
  required_cols <- c("id", "date")
  missing_cols  <- setdiff(required_cols, colnames(transitions))
  if (length(missing_cols) > 0) stop("Missing required columns in transitions: ", paste(missing_cols, collapse = ", "))
  
  transition_cols <- grep("_to_", colnames(transitions), value = TRUE)
  if (length(transition_cols) == 0) stop("No transition columns found (looking for '*_to_*' pattern)")
  
  cat("  • Available transitions:", paste(transition_cols, collapse = ", "), "\n")
  cat("  • Simulation IDs:", length(unique(transitions$id)), "\n")
  cat("  • Date range:", min(transitions$date), "to", max(transitions$date), "\n")
  
  return(list(
    transitions  = transitions,
    transmission = transmission,
    metadata     = metadata
  ))
}

inspect_data_directory <- function(data_dir = "saved_data/") {
  if (!dir.exists(data_dir)) {
    cat("Directory not found:", data_dir, "\n")
    return(NULL)
  }
  
  cat("Inspecting directory:", data_dir, "\n")
  cat(rep("=", 50), "\n")
  
  all_files <- list.files(data_dir, full.names = FALSE)
  if (length(all_files) == 0) {
    cat("No files found in directory.\n")
    return(NULL)
  }
  
  rds_files <- all_files[grepl("\\.rds$", all_files)]
  csv_files <- all_files[grepl("\\.csv$", all_files)]
  
  cat("RDS files found:\n"); for (file in rds_files) cat("  •", file, "\n")
  cat("\nCSV files found:\n"); for (file in csv_files) cat("  •", file, "\n")
  
  cat("\nIdentifying file patterns:\n")
  r0_patterns <- unique(gsub(".*R0_([0-9.]+).*", "\\1", rds_files[grepl("R0_", rds_files)]))
  if (length(r0_patterns) > 0) cat("  • R0 values found:", paste(r0_patterns, collapse = ", "), "\n")
  if (any(grepl("full", rds_files)))     cat("  • 'full' model files found\n")
  if (any(grepl("partial", rds_files)))  cat("  • 'partial' model files found\n")
  if (any(grepl("metadata", rds_files))) cat("  • Metadata files found\n")
  if (any(grepl("transitions", rds_files))) cat("  • Transitions files found\n")
  if (any(grepl("transmission", rds_files))) cat("  • Transmission files found\n")
  if (any(grepl("rt_ci", rds_files)))    cat("  • ABM R_t files found\n")
  
  return(list(rds_files = rds_files, csv_files = csv_files))
}

# ── DELAY FUNCTIONS (FIXED) ──────────────────────────────────────────────────────

apply_reporting_delays <- function(incidence_vector, delay_percentage = 0.3, max_delay = 7) {
  n_days <- length(incidence_vector)
  if (sum(incidence_vector) == 0 || delay_percentage == 0) return(incidence_vector)
  
  # Create extended vector to accommodate delays beyond original window
  delayed_incidence <- numeric(n_days + max_delay)
  
  for (day in seq_len(n_days)) {
    daily_cases <- incidence_vector[day]
    if (daily_cases > 0) {
      # Randomly select which cases are delayed
      n_delayed <- rbinom(1, daily_cases, delay_percentage)
      n_on_time <- daily_cases - n_delayed
      
      # Report on-time cases immediately
      delayed_incidence[day] <- delayed_incidence[day] + n_on_time
      
      # Distribute delayed cases across future days
      if (n_delayed > 0) {
        delays <- sample(1:max_delay, n_delayed, replace = TRUE)
        for (d in delays) {
          delay_day <- day + d
          # No need to check bounds since we created extended vector
          delayed_incidence[delay_day] <- delayed_incidence[delay_day] + 1
        }
      }
    }
  }
  
  # FIXED: Return the full extended vector instead of truncating
  # This preserves all delayed cases
  return(delayed_incidence)
}

# Diagnostic function to verify delay preservation
check_delay_preservation <- function(original, delayed, delay_percentage) {
  total_original <- sum(original)
  total_delayed <- sum(delayed)
  
  cat("Delay preservation check:\n")
  cat("  Original cases:", total_original, "\n")
  cat("  Delayed cases:", total_delayed, "\n")
  cat("  Difference:", total_delayed - total_original, "\n")
  cat("  Match:", abs(total_delayed - total_original) < 0.01, "\n")
  
  if (abs(total_delayed - total_original) > 0.01) {
    warning("Case count mismatch! Cases may have been lost or gained.")
  }
  
  invisible(list(original = total_original, delayed = total_delayed))
}

# ── SERIAL INTERVALS (from transmission) ─────────────────────────────────────────

# Expected transmission columns: sim_num, source, date, source_exposure_date
# 'source' == -1 are initial infections; SI = (case onset/exposure date) - (source exposure date)
calculate_serial_intervals <- function(transmission_data) {
  transmission_data %>%
    filter(source != -1) %>%
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
}

# ── CORE ANALYSIS (single sim) using SI ──────────────────────────────────────────

analyze_delays_single_sim <- function(saved_data, sim_id = 1, 
                                      how_infected = "susceptible_to_exposed",
                                      delay_scenarios = list("medium_delay" = 0.3,
                                                             "high_delay"   = 0.5),
                                      start_day = 0,
                                      max_delay = 7,
                                      verbose = FALSE) {
  if (verbose) cat("Analyzing delays for simulation", sim_id, "\n")
  
  original_sim <- saved_data$transitions %>%
    filter(id == sim_id, date >= start_day) %>%
    arrange(date)
  
  if (!how_infected %in% colnames(original_sim)) {
    cat("  ERROR: Column", how_infected, "not found!\n")
    return(NULL)
  }
  original_incidence <- original_sim[[how_infected]]
  if (sum(original_incidence) == 0 || length(original_incidence) < 7) {
    if (verbose) cat("  Insufficient incidence data\n")
    return(NULL)
  }
  
  si_all <- calculate_serial_intervals(saved_data$transmission)
  si_data <- si_all %>% filter(sim_num == sim_id)
  if (nrow(si_data) == 0) {
    if (verbose) cat("  No valid serial interval data for sim", sim_id, "\n")
    return(NULL)
  }
  
  SI_mean <- si_data$mean_si
  SI_sd   <- si_data$sd_si
  if (is.na(SI_sd) || SI_sd < 0.1 || SI_mean < 0.5 || SI_mean > 20) {
    if (verbose) cat("  Invalid SI: mean =", SI_mean, ", sd =", SI_sd, "\n")
    return(NULL)
  }
  
  all_results <- data.frame()
  
  for (scenario_name in names(delay_scenarios)) {
    delay_pct <- delay_scenarios[[scenario_name]]
    if (verbose) cat("  Processing", scenario_name, "- delaying", delay_pct*100, "% of cases\n")
    
    # FIXED: Now returns extended vector
    delayed_incidence <- apply_reporting_delays(original_incidence, delay_pct, max_delay)
    
    # Optional: Check preservation
    if (verbose) {
      check_delay_preservation(original_incidence, delayed_incidence, delay_pct)
    }
    
    if (sum(delayed_incidence) == 0) {
      if (verbose) cat("    No cases after applying delays\n")
      next
    }
    
    # Cori
    tryCatch({
      cori_result <- estimate_R(
        delayed_incidence,
        method = "parametric_si",
        config = make_config(mean_si = SI_mean, std_si = SI_sd)
      )
      
      # Store results with original day indexing (accounting for extension)
      cori_dates <- cori_result$R[, "t_start"]
      
      all_results <- bind_rows(all_results, data.frame(
        sim_id = sim_id,
        date = cori_dates,
        median_rt = cori_result$R[, "Median(R)"],
        ci_lower = cori_result$R[, "Quantile.0.025(R)"],
        ci_upper = cori_result$R[, "Quantile.0.975(R)"],
        method = "Cori",
        delay_scenario = scenario_name,
        delay_percentage = delay_pct,
        mean_si = SI_mean,
        std_si = SI_sd
      ))
      if (verbose) cat("    ✓ Cori method successful\n")
    }, error = function(e) if (verbose) cat("    ✗ Cori method failed:", e$message, "\n"))
    
    # Wallinga-Teunis
    tryCatch({
      wt_result <- wallinga_teunis(
        incid = delayed_incidence,
        method = "parametric_si",
        config = make_config(
          incid = delayed_incidence,
          method = "parametric_si",
          mean_si = SI_mean,
          std_si = SI_sd
        )
      )
      
      wt_dates <- wt_result$R[, "t_start"]
      
      all_results <- bind_rows(all_results, data.frame(
        sim_id = sim_id,
        date = wt_dates,
        median_rt = wt_result$R[, "Mean(R)"],
        ci_lower = wt_result$R[, "Quantile.0.025(R)"],
        ci_upper = wt_result$R[, "Quantile.0.975(R)"],
        method = "Wallinga-Teunis",
        delay_scenario = scenario_name,
        delay_percentage = delay_pct,
        mean_si = SI_mean,
        std_si = SI_sd
      ))
      if (verbose) cat("    ✓ Wallinga-Teunis method successful\n")
    }, error = function(e) if (verbose) cat("    ✗ Wallinga-Teunis method failed:", e$message, "\n"))
  }
  
  all_results
}

# ── MULTI-SIM ANALYSIS (SI-based) ────────────────────────────────────────────────

analyze_delays_all_sims <- function(saved_data, how_infected = "susceptible_to_exposed",
                                    delay_scenarios = list("medium_delay" = 0.3,
                                                           "high_delay"   = 0.5),
                                    max_sims = 60, start_day = 0, max_delay = 7, verbose = TRUE) {
  if (verbose) {
    cat("Analyzing delays across multiple simulations\n")
    cat("Transition:", how_infected, "\n")
    cat("Delay scenarios:", paste(names(delay_scenarios), "=", 
                                  paste0(unlist(delay_scenarios)*100, "%"), collapse = ", "), "\n")
    cat("Max delay:", max_delay, "days\n")
  }
  
  sim_ids <- unique(saved_data$transitions$id)
  if (length(sim_ids) > max_sims) {
    sim_ids <- sample(sim_ids, max_sims)
    if (verbose) cat("Using", max_sims, "random simulations out of", length(unique(saved_data$transitions$id)), "\n")
  }
  
  all_results <- data.frame()
  successful_sims <- 0
  
  for (i in seq_along(sim_ids)) {
    sim_id <- sim_ids[i]
    if (verbose && i %% 10 == 0) {
      cat("Processing simulation", i, "of", length(sim_ids), 
          "(", successful_sims, "successful so far)\n")
    }
    
    sim_results <- analyze_delays_single_sim(
      saved_data, sim_id, how_infected, delay_scenarios, 
      start_day = start_day, max_delay = max_delay, verbose = FALSE
    )
    
    if (!is.null(sim_results) && nrow(sim_results) > 0) {
      all_results <- bind_rows(all_results, sim_results)
      successful_sims <- successful_sims + 1
    }
  }
  
  if (nrow(all_results) > 0) {
    all_results$R0_value        <- saved_data$metadata$R0
    all_results$model_type      <- saved_data$metadata$model_type
    all_results$transition_type <- how_infected
  }
  
  if (verbose) {
    cat("Completed! Successful simulations:", successful_sims, "out of", length(sim_ids), "\n")
    cat("Total results:", nrow(all_results), "rows\n")
  }
  
  all_results
}

# ── PLOTTING (single-combo) ─────────────────────────────────────────────────────

create_combined_delay_plot <- function(delay_results, abm_data, R0_val, model_type, transition,
                                       date_min = 12, date_max = 75) {
  cat("Creating comparison-style plot for", model_type, "R0 =", R0_val, "\n")
  
  # Aggregate Cori/WT across sims; keep delay scenarios
  plot_data <- delay_results %>%
    filter(date >= date_min, date <= date_max) %>%
    group_by(method, date, delay_scenario, delay_percentage) %>%
    summarise(
      mean_R   = mean(median_rt, na.rm = TRUE),
      ci_lower = quantile(median_rt, 0.025, na.rm = TRUE),
      ci_upper = quantile(median_rt, 0.975, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Add ABM (true) baseline, harmonized to "ABM" and same columns
  if (!is.null(abm_data)) {
    abm_plot_data <- abm_data %>%
      rename(date = source_exposure_date) %>%
      filter(date >= date_min, date <= date_max) %>%
      transmute(
        method = "ABM",
        date   = date,
        mean_R = mean_rt,
        ci_lower = ci_lower,
        ci_upper = ci_upper,
        delay_scenario   = "true_baseline",
        delay_percentage = 0
      )
    plot_data <- bind_rows(plot_data, abm_plot_data)
  }
  
  # Ensure consistent method labels and ordering
  plot_data$method <- factor(plot_data$method, levels = c("ABM", "Cori", "Wallinga-Teunis"))
  
  delay_labels <- c(
    "true_baseline" = "ABM: No Delays",
    "medium_delay"  = "DOTTED: Medium Delay (30%)",
    "high_delay"    = "DASHED: High Delay (50%)"
  )
  method_colors <- c("ABM" = "blue", "Cori" = "red", "Wallinga-Teunis" = "green")
  
  ggplot(plot_data, aes(x = date, y = mean_R, color = method, fill = method)) +
    geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.20, color = NA) +
    geom_line(aes(linetype = delay_scenario), linewidth = 1.1) +
    geom_hline(yintercept = R0_val, linetype = "dashed", color = "black", linewidth = 0.9) +
    scale_color_manual(values = method_colors, name = "Method") +
    scale_fill_manual(values = method_colors, guide = "none") +
    scale_linetype_manual(
      values = c("true_baseline" = "solid", "medium_delay" = "dotted", "high_delay" = "dashed"),
      labels = delay_labels,
      name = "Delay Scenario"
    ) +
    scale_x_continuous(limits = c(0, date_max), breaks = seq(0, date_max, by = 25)) +
    labs(
      title = "R_t Estimates: ABM vs Cori vs Wallinga-Teunis",
      subtitle = paste0(
        toupper(model_type), " Model | ",
        gsub("_to_", " → ", toupper(transition)), " Transition | R0 = ", R0_val
      ),
      x = "Day", y = "R_t Estimate"
    ) +
    theme_minimal() +
    theme(
      legend.position   = "bottom",
      plot.title        = element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.subtitle     = element_text(hjust = 0.5, size = 11),
      axis.title        = element_text(size = 12),
      axis.text         = element_text(size = 10),
      panel.grid.minor  = element_blank()
    ) +
    annotate("text", x = 60, y = R0_val + 0.3,
             label = paste0("R0 = ", R0_val),
             color = "black", size = 4)
}

# ── MAIN WRAPPER: RUN COMPLETE DELAY ANALYSIS ────────────────────────────────────

run_complete_delay_analysis <- function(R0_values = c(1.5, 2.0, 3.0, 5.0),
                                        model_types = c("full", "partial"),
                                        transitions = c("susceptible_to_exposed"),
                                        delay_scenarios = list("medium_delay" = 0.3,
                                                               "high_delay"   = 0.5),
                                        max_sims_per_scenario = 60,
                                        start_day = 0,
                                        max_delay = 7,
                                        results_dir = "complete_delay_results",
                                        plots_dir = "comparison_delays2") {
  cat("╔══════════════════════════════════════════════════════════════╗\n")
  cat("║     COMPLETE REPORTING DELAYS ANALYSIS (SI from transmission)║\n")
  cat("║            Processing All R0 Values with Uniform(1,7)        ║\n")
  cat("║                        FIXED VERSION                         ║\n")
  cat("╚══════════════════════════════════════════════════════════════╝\n")
  
  cat("\nParameters:\n")
  cat("• R0 values:", paste(R0_values, collapse = ", "), "\n")
  cat("• Model types:", paste(model_types, collapse = ", "), "\n")
  cat("• Transitions:", paste(transitions, collapse = ", "), "\n")
  cat("• Delay scenarios:", paste(names(delay_scenarios), "=", 
                                  paste0(unlist(delay_scenarios)*100, "%"), collapse = ", "), "\n")
  cat("• Max delay:", max_delay, "days\n")
  cat("• Max simulations per scenario:", max_sims_per_scenario, "\n")
  cat("• Start day for incidence:", start_day, "\n")
  cat("• Results directory:", results_dir, "\n")
  cat("• Plots directory:", plots_dir, "\n")
  
  if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)
  if (!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE)
  
  all_results <- data.frame()
  processing_summary <- data.frame()
  
  for (model_type in model_types) {
    for (transition in transitions) {
      for (R0_val in R0_values) {
        combination_name <- paste0(model_type, "_R0_", R0_val, "_", transition)
        cat("\n", rep("=", 70), "\n")
        cat("Processing:", combination_name, "\n")
        cat(rep("=", 70), "\n")
        
        file_prefix <- paste0("saved_data/", model_type, "_R0_", R0_val, "_n_1e+05_nsim_100")
        
        if (!file.exists(paste0(file_prefix, "_metadata.rds"))) {
          cat("⚠️  Simulation data not found:", file_prefix, "\n")
          processing_summary <- bind_rows(processing_summary, data.frame(
            combination = combination_name,
            R0_value = R0_val,
            model_type = model_type,
            transition = transition,
            n_results = 0,
            status = "❌ Data file not found"
          ))
          next
        }
        
        tryCatch({
          cat("📂 Loading simulation data...\n")
          saved_data <- load_seir_data(file_prefix)
          
          abm_file <- paste0(file_prefix, "_rt_ci.rds")
          abm_data <- if (file.exists(abm_file)) readRDS(abm_file) else { 
            cat("⚠️  ABM data not found, will plot without true baseline\n"); NULL 
          }
          
          cat("🔄 Running delay analysis (SI-based, FIXED)...\n")
          delay_results <- analyze_delays_all_sims(
            saved_data,
            how_infected = transition,
            delay_scenarios = delay_scenarios,
            max_sims = max_sims_per_scenario,
            start_day = start_day,
            max_delay = max_delay,
            verbose = TRUE
          )
          
          if (!is.null(delay_results) && nrow(delay_results) > 0) {
            all_results <- bind_rows(all_results, delay_results)
            
            individual_file <- file.path(results_dir, paste0("delay_results_", combination_name, ".rds"))
            saveRDS(delay_results, individual_file)
            csv_file <- file.path(results_dir, paste0("delay_results_", combination_name, ".csv"))
            write.csv(delay_results, csv_file, row.names = FALSE)
            
            cat("🎨 Creating plot...\n")
            plot_obj <- create_combined_delay_plot(delay_results, abm_data, R0_val, model_type, transition)
            plot_file <- file.path(plots_dir, paste0("combined_plot_", combination_name, ".png"))
            ggsave(plot_file, plot_obj, width = 14, height = 10, dpi = 300)
            
            processing_summary <- bind_rows(processing_summary, data.frame(
              combination = combination_name,
              R0_value = R0_val,
              model_type = model_type,
              transition = transition,
              n_results = nrow(delay_results),
              status = "✅ Success"
            ))
            
            cat("✅ Completed successfully!\n")
            cat("   • Results:", nrow(delay_results), "rows\n")
            cat("   • Simulations:", length(unique(delay_results$sim_id)), "\n")
            cat("   • Files saved:", basename(individual_file), ",", basename(plot_file), "\n")
          } else {
            processing_summary <- bind_rows(processing_summary, data.frame(
              combination = combination_name,
              R0_value = R0_val,
              model_type = model_type,
              transition = transition,
              n_results = 0,
              status = "❌ No results obtained"
            ))
            cat("❌ No results obtained\n")
          }
        }, error = function(e) {
          cat("💥 Error:", e$message, "\n")
          processing_summary <<- bind_rows(processing_summary, data.frame(
            combination = combination_name,
            R0_value = R0_val,
            model_type = model_type,
            transition = transition,
            n_results = 0,
            status = paste("💥 Error:", e$message)
          ))
        })
      }
    }
  }
  
  if (nrow(all_results) > 0) {
    combined_file <- file.path(results_dir, "combined_all_results.rds")
    saveRDS(all_results, combined_file)
    combined_csv <- file.path(results_dir, "combined_all_results.csv")
    write.csv(all_results, combined_csv, row.names = FALSE)
    
    summary_stats <- all_results %>%
      filter(date >= 10 & date <= 75) %>%
      group_by(method, delay_scenario, delay_percentage, R0_value, model_type, transition_type) %>%
      summarise(
        n_observations = n(),
        mean_rt  = mean(median_rt, na.rm = TRUE),
        median_rt = median(median_rt, na.rm = TRUE),
        bias     = mean_rt - R0_value,
        abs_bias = abs(bias),
        rmse     = sqrt(mean((median_rt - R0_value)^2, na.rm = TRUE)),
        .groups = "drop"
      )
    summary_file <- file.path(results_dir, "summary_statistics.csv")
    write.csv(summary_stats, summary_file, row.names = FALSE)
    
    cat("\n💾 Combined results saved:", nrow(all_results), "total results\n")
  }
  
  summary_csv <- file.path(results_dir, "processing_summary.csv")
  write.csv(processing_summary, summary_csv, row.names = FALSE)
  
  cat("\n📊 PROCESSING SUMMARY:\n")
  print(processing_summary)
  
  cat("\n")
  cat("╔══════════════════════════════════════════════════════════════╗\n")
  cat("║                    ANALYSIS COMPLETE!                       ║\n")
  cat("╚══════════════════════════════════════════════════════════════╝\n")
  
  cat("\n📁 Results saved in:", results_dir, "/\n")
  cat("├── 💾 combined_all_results.rds/.csv - All results\n") 
  cat("├── 📊 summary_statistics.csv       - Summary stats by scenario\n")
  cat("├── 📋 processing_summary.csv       - Processing status\n")
  cat("├── 📈 Individual result .rds/.csv  - Per combination\n")
  cat("\n📁 Plots saved in:", plots_dir, "/\n")
  cat("└── 🖼️  Combined plots per combo     - ABM + Cori + Wallinga\n")
  
  return(list(results = all_results, summary = processing_summary))
}

# ── TEST / DIAGNOSTICS (optional) ────────────────────────────────────────────────

test_single_simulation_delays <- function(model_type = "full", R0_val = 2.0, 
                                          sim_id = 1, transition = "susceptible_to_exposed",
                                          start_day = 0, max_delay = 7,
                                          output_dir = "comparison_delays2") {
  cat("╔══════════════════════════════════════════════════════════════╗\n")
  cat("║                TESTING SINGLE SIMULATION (SI) - FIXED       ║\n")
  cat("╚══════════════════════════════════════════════════════════════╝\n")
  
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  cat("🔍 Inspecting data directory...\n")
  file_info <- inspect_data_directory("saved_data/")
  
  possible_prefixes <- c(
    paste0("saved_data/", model_type, "_R0_", R0_val, "_n_1e+05_nsim_100"),
    paste0("saved_data/", model_type, "_R0_", gsub("\\.", "_", R0_val), "_n_1e+05_nsim_100"),
    paste0("saved_data/", model_type, "_R0_", R0_val)
  )
  
  saved_data <- NULL
  file_prefix <- NULL
  
  cat("\n📂 Trying different file name patterns...\n")
  for (prefix in possible_prefixes) {
    cat("Trying:", prefix, "\n")
    if (file.exists(paste0(prefix, "_metadata.rds"))) {
      cat("  ✓ Found matching files!\n")
      file_prefix <- prefix
      break
    } else {
      cat("  ✗ No match\n")
    }
  }
  
  if (is.null(file_prefix)) {
    cat("❌ Could not find data files. Please check:\n")
    cat("1. Files are in 'saved_data/' directory\n")
    cat("2. File naming follows expected pattern\n")
    cat("3. Required files exist: *_metadata.rds, *_transitions.rds, *_transmission.rds\n")
    if (!is.null(file_info)) {
      cat("\nAvailable files suggest trying:\n")
      suggested_files <- file_info$rds_files[grepl("metadata", file_info$rds_files)]
      if (length(suggested_files) > 0) {
        suggested_prefixes <- gsub("_metadata\\.rds$", "", suggested_files)
        cat("Suggested prefixes:\n")
        for (pref in suggested_prefixes[1:min(3, length(suggested_prefixes))]) cat("  •", pref, "\n")
      }
    }
    return(NULL)
  }
  
  tryCatch({
    cat("📂 Loading test data:", file_prefix, "\n")
    saved_data <- load_seir_data(file_prefix)
    
    abm_file <- paste0(file_prefix, "_rt_ci.rds")
    abm_data <- if (file.exists(abm_file)) {
      cat("📊 Loading ABM data:", basename(abm_file), "\n")
      readRDS(abm_file)
    } else { cat("⚠️  ABM data not found, will plot without true baseline\n"); NULL }
    
  }, error = function(e) {
    cat("❌ Error loading data:", e$message, "\n")
    return(NULL)
  })
  
  if (is.null(saved_data)) return(NULL)
  
  cat("📊 Data loaded successfully!\n")
  cat("   • Total simulations:", length(unique(saved_data$transitions$id)), "\n")
  cat("   • ABM data available:", !is.null(abm_data), "\n")
  
  cat("\n🧪 Testing FIXED delay function...\n")
  test_incidence <- c(0, 2, 5, 8, 12, 15, 10, 8, 5, 3, 2, 1, 0, 0, 0)
  cat("Original incidence:", paste(test_incidence, collapse = ", "), "\n")
  cat("Original total:", sum(test_incidence), "cases\n\n")
  
  medium_delayed <- apply_reporting_delays(test_incidence, 0.3, max_delay)
  cat("Medium delayed (30%):", paste(head(medium_delayed, 15), collapse = ", "), "...\n")
  cat("Medium delayed total:", sum(medium_delayed), "cases\n")
  cat("Extended by", length(medium_delayed) - length(test_incidence), "days\n\n")
  
  high_delayed <- apply_reporting_delays(test_incidence, 0.5, max_delay)
  cat("High delayed (50%):", paste(head(high_delayed, 15), collapse = ", "), "...\n")
  cat("High delayed total:", sum(high_delayed), "cases\n")
  cat("Extended by", length(high_delayed) - length(test_incidence), "days\n")
  
  cat("\n🔬 Testing single simulation analysis (SI-based, FIXED)...\n")
  delay_results <- analyze_delays_single_sim(
    saved_data, sim_id = sim_id, how_infected = transition,
    delay_scenarios = list("medium_delay" = 0.3, "high_delay" = 0.5),
    start_day = start_day, max_delay = max_delay, verbose = TRUE
  )
  
  if (!is.null(delay_results) && nrow(delay_results) > 0) {
    cat("✅ Single simulation test successful!\n")
    cat("   • Results:", nrow(delay_results), "rows\n")
    cat("   • Methods:", paste(unique(delay_results$method), collapse = ", "), "\n")
    cat("   • Scenarios:", paste(unique(delay_results$delay_scenario), collapse = ", "), "\n")
    
    cat("\n🎨 Creating test plot...\n")
    plot_obj <- create_combined_delay_plot(delay_results, abm_data, R0_val, model_type, transition)
    test_plot_file <- file.path(output_dir, paste0("test_delay_plot_", model_type, "_R0_", R0_val, "_sim_", sim_id, ".png"))
    ggsave(test_plot_file, plot_obj, width = 12, height = 8, dpi = 300)
    cat("📊 Test plot saved:", test_plot_file, "\n")
    return(delay_results)
  } else {
    cat("❌ Single simulation test failed - no results obtained\n")
    return(NULL)
  }
}

quick_diagnosis <- function() {
  cat("🔍 QUICK DIAGNOSIS (SI-based, FIXED)\n")
  cat("====================================\n")
  
  if (!dir.exists("saved_data/")) {
    cat("❌ 'saved_data/' directory not found\n")
    cat("   Please create it and place your simulation files there\n")
    return(FALSE)
  }
  inspect_data_directory("saved_data/")
  
  cat("\n🧪 Trying to load a sample file...\n")
  test_results <- test_single_simulation_delays(output_dir = "comparison_delays2")
  !is.null(test_results)
}

# ── COMPREHENSIVE PLOTS FROM SAVED DELAY RESULTS ────────────────────────────────

create_comprehensive_plots_from_delay_results <- function(
    results_dir   = "complete_delay_results",
    plots_dir     = "comparison_delays2",
    R0_values     = c(1.5, 2.0, 3.0, 5.0),
    model_types   = c("full","partial"),
    transitions   = c("susceptible_to_exposed","exposed_to_infected"),
    date_min      = 12,
    date_max      = 75
) {
  cat("\n=== CREATING COMPREHENSIVE COMPARISON PLOTS (from delay analysis) ===\n")
  if (!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE)
  
  method_colors <- c("ABM" = "blue", "Cori" = "red", "Wallinga-Teunis" = "green")
  linetypes_map <- c("true_baseline" = "solid", "medium_delay" = "dotted", "high_delay" = "dashed")
  delay_labels  <- c(
    "true_baseline" = "ABM: No Delays",
    "medium_delay"  = "DOTTED: Medium Delay (30%)",
    "high_delay"    = "DASHED: High Delay (50%)"
  )
  
  for (model_type in model_types) {
    for (transition in transitions) {
      cat("\n--- Creating plots for", model_type, "model,", transition, "(from", results_dir, ") ---\n")
      
      plots_list <- vector("list", length(R0_values))
      all_data_combined <- tibble()
      
      for (i in seq_along(R0_values)) {
        R0_val <- R0_values[i]
        combo   <- paste0(model_type, "_R0_", R0_val, "_", transition)
        delay_file <- file.path(results_dir, paste0("delay_results_", combo, ".rds"))
        abm_file   <- file.path("saved_data", paste0(model_type, "_R0_", R0_val, "_n_1e+05_nsim_100_rt_ci.rds"))
        
        if (file.exists(delay_file)) {
          cat("  • Loading delay results for R0 =", R0_val, "\n")
          delay_results <- readRDS(delay_file)
          
          # Aggregate Cori/WT across sims; retain delay scenarios
          cw_plot_data <- delay_results %>%
            filter(date >= date_min, date <= date_max) %>%
            group_by(method, date, delay_scenario, delay_percentage) %>%
            summarise(
              mean_R   = mean(median_rt, na.rm = TRUE),
              ci_lower = quantile(median_rt, 0.025, na.rm = TRUE),
              ci_upper = quantile(median_rt, 0.975, na.rm = TRUE),
              .groups  = "drop"
            ) %>%
            mutate(R0_value = R0_val)
          
          # Add ABM baseline if present
          abm_plot_data <- tibble()
          if (file.exists(abm_file)) {
            abm_plot_data <- readRDS(abm_file) %>%
              rename(date = source_exposure_date) %>%
              filter(date >= date_min, date <= date_max) %>%
              transmute(
                method           = "ABM",
                date,
                mean_R           = mean_rt,
                ci_lower, ci_upper,
                delay_scenario   = "true_baseline",
                delay_percentage = 0,
                R0_value         = R0_val
              )
          } else {
            cat("    ⚠️ ABM baseline not found for", combo, "— plotting CW/WT only\n")
          }
          
          combined_data <- bind_rows(cw_plot_data, abm_plot_data) %>%
            mutate(method = factor(method, levels = c("ABM","Cori","Wallinga-Teunis")))
          
          all_data_combined <- bind_rows(all_data_combined, combined_data)
          
          # Individual small panel (per R0)
          p_individual <- ggplot(
            combined_data,
            aes(x = date, y = mean_R, color = method, fill = method)
          ) +
            geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.20, color = NA) +
            geom_line(aes(linetype = delay_scenario), linewidth = 1.1) +
            geom_hline(yintercept = R0_val, linetype = "dashed", color = "black", linewidth = 0.9) +
            scale_color_manual(values = method_colors) +
            scale_fill_manual(values = method_colors, guide = "none") +
            scale_linetype_manual(values = linetypes_map, labels = delay_labels, name = NULL) +
            scale_x_continuous(limits = c(0, date_max), breaks = seq(0, date_max, by = 25)) +
            labs(title = paste0("R0 = ", R0_val), x = "Day", y = "R_t Estimate") +
            theme_minimal() +
            theme(
              legend.position  = "none",
              plot.title       = element_text(hjust = 0.5, size = 14, face = "bold"),
              axis.title       = element_text(size = 12),
              axis.text        = element_text(size = 10),
              panel.grid.minor = element_blank()
            ) +
            annotate("text", x = date_max - 15, y = R0_val + 0.3,
                     label = paste0("R0 = ", R0_val), color = "black", size = 4)
          
          plots_list[[i]] <- p_individual
          
        } else {
          cat("  ✗ Delay results not found for", combo, "\n")
          plots_list[[i]] <- ggplot() +
            annotate("text", x = 0.5, y = 0.5,
                     label = paste("No delay_results for\nR0 =", R0_val),
                     size = 6) +
            theme_void() +
            labs(title = paste0("R0 = ", R0_val)) +
            theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
        }
      }
      
      # 2×2 grid (panels per R0)
      if (length(plots_list) == 4) {
        combined_grid <- ggarrange(
          plots_list[[1]], plots_list[[2]],
          plots_list[[3]], plots_list[[4]],
          ncol = 2, nrow = 2, common.legend = FALSE
        )
        final_grid <- annotate_figure(
          combined_grid,
          top = text_grob(
            paste0(
              "R_t under Reporting Delays: ABM vs Cori vs Wallinga-Teunis\n",
              toupper(model_type), " Model | ",
              gsub("_to_", " → ", toupper(transition)), " Transition"
            ),
            face = "bold", size = 16
          )
        )
        grid_filename <- file.path(plots_dir, paste0("grid_", model_type, "_", transition, ".png"))
        ggsave(grid_filename, final_grid, width = 16, height = 12, dpi = 300)
        cat("  ✓ 2x2 grid saved:", basename(grid_filename), "\n")
      }
      
      # Faceted plot across R0 (keeps delay scenarios via linetype)
      if (nrow(all_data_combined) > 0) {
        all_data_combined$method <- factor(all_data_combined$method, levels = c("ABM","Cori","Wallinga-Teunis"))
        
        p_facet <- ggplot(
          all_data_combined,
          aes(x = date, y = mean_R, color = method, fill = method, linetype = delay_scenario)
        ) +
          geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.20, color = NA) +
          geom_line(linewidth = 1.05) +
          geom_hline(aes(yintercept = R0_value), linetype = "dashed", color = "black", linewidth = 0.8) +
          scale_color_manual(values = method_colors, name = "Method") +
          scale_fill_manual(values = method_colors, guide = "none") +
          scale_linetype_manual(values = linetypes_map, labels = delay_labels, name = "Delay Scenario") +
          scale_x_continuous(limits = c(0, date_max), breaks = seq(0, date_max, by = 25)) +
          facet_wrap(~paste0("R0 = ", R0_value), ncol = 2, scales = "free_y") +
          labs(
            title = paste0(
              "R_t Across R0 Values under Reporting Delays: ABM vs Cori vs Wallinga-Teunis\n",
              toupper(model_type), " Model | ",
              gsub("_to_", " → ", toupper(transition)), " Transition"
            ),
            x = "Day", y = "R_t Estimate"
          ) +
          theme_minimal() +
          theme(
            legend.position   = "bottom",
            plot.title        = element_text(hjust = 0.5, size = 14, face = "bold"),
            strip.text        = element_text(size = 12, face = "bold"),
            panel.grid.minor  = element_blank()
          )
        
        facet_filename <- file.path(plots_dir, paste0("facet_", model_type, "_", transition, ".png"))
        ggsave(facet_filename, p_facet, width = 14, height = 10, dpi = 300)
        cat("  ✓ Faceted plot saved:", basename(facet_filename), "\n")
        
        # Summary stats by method & delay scenario
        summary_stats <- all_data_combined %>%
          filter(date >= 10, date <= date_max) %>%
          group_by(method, delay_scenario, delay_percentage, R0_value) %>%
          summarise(
            n_obs     = n(),
            mean_rt   = mean(mean_R, na.rm = TRUE),
            median_rt = median(mean_R, na.rm = TRUE),
            sd_rt     = sd(mean_R, na.rm = TRUE),
            bias      = mean_rt - R0_value,
            abs_bias  = abs(bias),
            rmse      = sqrt(mean((mean_R - R0_value)^2, na.rm = TRUE)),
            .groups   = "drop"
          )
        
        summary_filename <- file.path(plots_dir, paste0("summary_", model_type, "_", transition, ".csv"))
        write_csv(summary_stats, summary_filename)
        cat("  ✓ Summary stats saved:", basename(summary_filename), "\n")
        
        # Optional: overall ranking table (averaged over R0) by method & delay scenario
        ranking <- summary_stats %>%
          group_by(method, delay_scenario) %>%
          summarise(
            avg_abs_bias = mean(abs_bias, na.rm = TRUE),
            avg_rmse     = mean(rmse, na.rm = TRUE),
            .groups = "drop"
          ) %>%
          arrange(avg_rmse, avg_abs_bias)
        
        ranking_filename <- file.path(plots_dir, paste0("ranking_", model_type, "_", transition, ".csv"))
        write_csv(ranking, ranking_filename)
        cat("  ✓ Ranking saved:", basename(ranking_filename), "\n")
      }
    }
  }
  
  cat("\n=== DONE: COMPREHENSIVE PLOTS (from delay analysis) ===\n")
  cat("All plots saved to:", plots_dir, "\n")
}

# ── QUICK START (uncomment to run) ───────────────────────────────────────────────

# Test the fix first with diagnostic
# quick_diagnosis()

# Or run a single simulation test
# test_single_simulation_delays(model_type = "full", R0_val = 2.0, sim_id = 1, output_dir = "comparison_delays2")

# 1) Run the full SI-based reporting-delays analysis and save results/plots:
run_complete_delay_analysis(
  R0_values = c(1.5, 2.0, 3.0, 5.0),
  model_types = c("full", "partial"),
  transitions = c("susceptible_to_exposed","exposed_to_infected"),
  delay_scenarios = list("medium_delay" = 0.3, "high_delay" = 0.5),
  max_sims_per_scenario = 60,
  start_day = 0,
  max_delay = 7,
  results_dir = "complete_delay_results",
  plots_dir = "comparison_delays2"
)

# 2) Build the comprehensive multi-R0 plots FROM the saved delay results:
create_comprehensive_plots_from_delay_results(
  results_dir = "complete_delay_results",
  plots_dir   = "comparison_delays2",
  R0_values   = c(1.5, 2.0, 3.0, 5.0),
  model_types = c("full","partial"),
  transitions = c("susceptible_to_exposed","exposed_to_infected")
)