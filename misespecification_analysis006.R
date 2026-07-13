# ================================================================================
# COMPLETE MISSPECIFICATION ANALYSIS (SI-BASED) - ALL R0 VALUES
# 70% reporting rate misspecification (binomial thinning)
# ABM = True baseline, Cori/Wallinga = Underreported data
# Also builds comprehensive multi-R0 plots from the saved results.
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

# ====== CONFIG (edit these if you like) =========================================
DEFAULT_R0_VALUES      <- c(1.5, 2.0, 3.0, 5.0)
DEFAULT_MODEL_TYPES    <- c("full","partial")
DEFAULT_TRANSITIONS    <- c("susceptible_to_exposed","exposed_to_infected")
DEFAULT_REPORTING_RATE <- 0.7
DEFAULT_MAX_SIMS       <- 30
DEFAULT_START_DAY      <- 0
DEFAULT_RESULTS_DIR    <- "complete_misspec_results"
DEFAULT_PLOTS_DIR      <- "misspec_comparison_plots"
DEFAULT_DATE_MIN       <- 12
DEFAULT_DATE_MAX       <- 75
DEFAULT_SEED           <- 12345

# ====== HELPERS =================================================================

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
  
  list(transitions = transitions, transmission = transmission, metadata = metadata)
}

inspect_data_directory <- function(data_dir = "saved_data/") {
  if (!dir.exists(data_dir)) {
    cat("Directory not found:", data_dir, "\n")
    return(NULL)
  }
  cat("Inspecting directory:", data_dir, "\n")
  cat(strrep("=", 50), "\n")
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
  list(rds_files = rds_files, csv_files = csv_files)
}

# ====== MISSPECIFICATION (UNDERREPORTING) =======================================

apply_underreporting <- function(incidence_vector, reporting_rate = 0.7) {
  # Binomial thinning per day
  underreported_incidence <- numeric(length(incidence_vector))
  for (day in seq_along(incidence_vector)) {
    daily_cases <- incidence_vector[day]
    if (daily_cases > 0) {
      underreported_incidence[day] <- rbinom(1, daily_cases, reporting_rate)
    }
  }
  underreported_incidence
}

# ====== SERIAL INTERVALS (from transmission) ====================================
# Expected: transmission has sim_num, source, date, source_exposure_date
# 'source' == -1 are initial infections; SI = (case date) - (source exposure date)
calculate_serial_intervals <- function(transmission_data) {
  transmission_data %>%
    filter(source != -1) %>%
    mutate(serial_interval = date - source_exposure_date) %>%
    filter(serial_interval > 0) %>%
    group_by(sim_num) %>%
    summarise(
      mean_si = mean(serial_interval, na.rm = TRUE),
      sd_si   = sd(serial_interval,   na.rm = TRUE),
      n_transmissions = n(),
      .groups = 'drop'
    ) %>%
    filter(!is.na(mean_si), !is.na(sd_si), sd_si > 0.1, mean_si > 0.5, mean_si < 20)
}

# ====== CORE ANALYSIS (single sim) ==============================================

analyze_misspec_single_sim <- function(saved_data, sim_id = 1, 
                                       how_infected = "susceptible_to_exposed",
                                       reporting_rate = 0.7,
                                       start_day = 0) {
  cat("Analyzing misspecification for simulation", sim_id, "\n")
  
  original_sim <- saved_data$transitions %>%
    filter(id == sim_id, date >= start_day) %>%
    arrange(date)
  
  if (!how_infected %in% colnames(original_sim)) {
    cat("  ERROR: Column", how_infected, "not found!\n")
    return(NULL)
  }
  
  original_incidence <- original_sim[[how_infected]]
  if (sum(original_incidence) == 0 || length(original_incidence) < 7) {
    cat("  Insufficient incidence data\n")
    return(NULL)
  }
  
  si_all  <- calculate_serial_intervals(saved_data$transmission)
  si_data <- si_all %>% filter(sim_num == sim_id)
  if (nrow(si_data) == 0) {
    cat("  No valid serial interval data for sim", sim_id, "\n")
    return(NULL)
  }
  
  SI_mean <- si_data$mean_si
  SI_sd   <- si_data$sd_si
  if (is.na(SI_sd) || SI_sd < 0.1 || SI_mean < 0.5 || SI_mean > 20) {
    cat("  Invalid SI: mean =", SI_mean, ", sd =", SI_sd, "\n")
    return(NULL)
  }
  
  cat("  Applying underreporting at", reporting_rate*100, "%\n")
  underreported_incidence <- apply_underreporting(original_incidence, reporting_rate)
  if (sum(underreported_incidence) == 0) {
    cat("    No cases after applying underreporting\n")
    return(NULL)
  }
  
  all_results <- data.frame()
  
  # Cori
  tryCatch({
    cori_result <- estimate_R(
      underreported_incidence,
      method = "parametric_si",
      config = make_config(mean_si = SI_mean, std_si = SI_sd)
    )
    all_results <- bind_rows(all_results, data.frame(
      sim_id = sim_id,
      date = cori_result$R[, "t_start"],
      median_rt = cori_result$R[, "Median(R)"],
      ci_lower = cori_result$R[, "Quantile.0.025(R)"],
      ci_upper = cori_result$R[, "Quantile.0.975(R)"],
      method = "Cori",
      reporting_rate = reporting_rate,
      mean_si = SI_mean,
      std_si = SI_sd
    ))
    cat("    ✓ Cori method successful\n")
  }, error = function(e) cat("    ✗ Cori method failed:", e$message, "\n"))
  
  # Wallinga-Teunis
  tryCatch({
    wt_result <- wallinga_teunis(
      incid = underreported_incidence,
      method = "parametric_si",
      config = make_config(
        incid = underreported_incidence,
        method = "parametric_si",
        mean_si = SI_mean,
        std_si = SI_sd
      )
    )
    all_results <- bind_rows(all_results, data.frame(
      sim_id = sim_id,
      date = wt_result$R[, "t_start"],
      median_rt = wt_result$R[, "Mean(R)"],
      ci_lower = wt_result$R[, "Quantile.0.025(R)"],
      ci_upper = wt_result$R[, "Quantile.0.975(R)"],
      method = "Wallinga-Teunis",
      reporting_rate = reporting_rate,
      mean_si = SI_mean,
      std_si = SI_sd
    ))
    cat("    ✓ Wallinga-Teunis method successful\n")
  }, error = function(e) cat("    ✗ Wallinga-Teunis method failed:", e$message, "\n"))
  
  all_results
}

# ====== MULTI-SIM ANALYSIS ======================================================

analyze_misspec_all_sims <- function(saved_data, how_infected = "susceptible_to_exposed",
                                     reporting_rate = 0.7,
                                     max_sims = 50, start_day = 0, verbose = TRUE) {
  if (verbose) {
    cat("Analyzing misspecification across multiple simulations\n")
    cat("Transition:", how_infected, "\n")
    cat("Reporting rate:", reporting_rate*100, "%\n")
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
    
    sim_results <- analyze_misspec_single_sim(
      saved_data, sim_id, how_infected, reporting_rate, start_day = start_day
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

# ====== PLOTTING (single combo) ================================================

create_combined_misspec_plot <- function(misspec_results, abm_data, R0_val, model_type, transition) {
  cat("Creating combined plot for", model_type, "R0 =", R0_val, "\n")
  
  plot_data <- misspec_results %>%
    filter(date >= DEFAULT_DATE_MIN, date <= DEFAULT_DATE_MAX) %>%
    group_by(method, date, reporting_rate) %>%
    summarise(
      mean_R   = mean(median_rt, na.rm = TRUE),
      ci_lower = quantile(median_rt, 0.025, na.rm = TRUE),
      ci_upper = quantile(median_rt, 0.975, na.rm = TRUE),
      .groups  = "drop"
    )
  
  if (!is.null(abm_data)) {
    abm_plot_data <- abm_data %>%
      rename(date = source_exposure_date) %>%
      filter(date >= DEFAULT_DATE_MIN, date <= DEFAULT_DATE_MAX) %>%
      transmute(
        method = "ABM",
        date = date,
        mean_R = mean_rt,
        ci_lower, ci_upper,
        reporting_rate = 1.0
      )
    plot_data <- bind_rows(plot_data, abm_plot_data)
  }
  
  plot_data$method <- factor(plot_data$method, levels = c("ABM","Cori","Wallinga-Teunis"))
  plot_data$rate_label <- paste0(round(plot_data$reporting_rate*100), "%")
  
  # map reporting rates to linetypes
  .make_linetypes <- function(rate_labels) {
    base_types <- c("solid","dotted","dashed","twodash","longdash","dotdash")
    rate_labels <- unique(rate_labels)
    if ("100%" %in% rate_labels) {
      others <- setdiff(rate_labels, "100%")
      types <- c("100%" = "solid")
      if (length(others) > 0) {
        types <- c(types, setNames(base_types[2:(length(others)+1)], others))
      }
    } else {
      types <- setNames(base_types[seq_along(rate_labels)], rate_labels)
    }
    types
  }
  ltypes <- .make_linetypes(plot_data$rate_label)
  
  method_colors <- c("ABM" = "blue", "Cori" = "red", "Wallinga-Teunis" = "green")
  
  ggplot(plot_data, aes(x = date, y = mean_R, color = method, fill = method)) +
    geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.20, color = NA) +
    geom_line(aes(linetype = rate_label), linewidth = 1.1) +
    geom_hline(yintercept = R0_val, linetype = "dashed", color = "black", linewidth = 0.9) +
    scale_color_manual(values = method_colors, name = "Method") +
    scale_fill_manual(values = method_colors, guide = "none") +
    scale_linetype_manual(values = ltypes, name = "Reporting rate") +
    scale_x_continuous(limits = c(0, DEFAULT_DATE_MAX), breaks = seq(0, DEFAULT_DATE_MAX, by = 25)) +
    labs(
      title = "R_t under Underreporting Misspecification",
      subtitle = paste0(
        toupper(model_type), " Model | ",
        gsub("_to_", " \u2192 ", toupper(transition)), " Transition | R0 = ", R0_val
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
    )
}

# ====== MAIN MISSPEC ANALYSIS WRAPPER ===========================================

run_complete_misspec_analysis <- function(R0_values = DEFAULT_R0_VALUES,
                                          model_types = DEFAULT_MODEL_TYPES,
                                          transitions = DEFAULT_TRANSITIONS,
                                          reporting_rate = DEFAULT_REPORTING_RATE,
                                          max_sims_per_scenario = DEFAULT_MAX_SIMS,
                                          start_day = DEFAULT_START_DAY,
                                          results_dir = DEFAULT_RESULTS_DIR,
                                          rng_seed = DEFAULT_SEED) {
  cat("╔══════════════════════════════════════════════════════════════╗\n")
  cat("║     COMPLETE MISSPECIFICATION ANALYSIS (SI from transmission)║\n")
  cat("║             Processing All R0 Values with Reporting < 100%   ║\n")
  cat("╚══════════════════════════════════════════════════════════════╝\n")
  
  set.seed(rng_seed)
  
  cat("\nParameters:\n")
  cat("• R0 values:", paste(R0_values, collapse = ", "), "\n")
  cat("• Model types:", paste(model_types, collapse = ", "), "\n")
  cat("• Transitions:", paste(transitions, collapse = ", "), "\n")
  cat("• Reporting rate:", reporting_rate*100, "%\n")
  cat("• Max simulations per scenario:", max_sims_per_scenario, "\n")
  cat("• Start day for incidence:", start_day, "\n")
  cat("• Results dir:", results_dir, "\n")
  
  if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)
  
  all_results <- data.frame()
  processing_summary <- data.frame()
  
  for (model_type in model_types) {
    for (transition in transitions) {
      for (R0_val in R0_values) {
        combination_name <- paste0(model_type, "_R0_", R0_val, "_", transition)
        cat("\n", strrep("=", 70), "\n")
        cat("Processing:", combination_name, "\n")
        cat(strrep("=", 70), "\n")
        
        file_prefix <- paste0("saved_data/", model_type, "_R0_", R0_val, "_n_1e+05_nsim_1000")
        
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
            cat("⚠️  ABM data not found, plots will omit true baseline for this combo\n"); NULL 
          }
          
          cat("🔄 Running misspecification analysis (SI-based)...\n")
          misspec_results <- analyze_misspec_all_sims(
            saved_data,
            how_infected = transition,
            reporting_rate = reporting_rate,
            max_sims = max_sims_per_scenario,
            start_day = start_day,
            verbose = TRUE
          )
          
          if (!is.null(misspec_results) && nrow(misspec_results) > 0) {
            all_results <- bind_rows(all_results, misspec_results)
            
            individual_file <- file.path(results_dir, paste0("misspec_results_", combination_name, ".rds"))
            saveRDS(misspec_results, individual_file)
            csv_file <- file.path(results_dir, paste0("misspec_results_", combination_name, ".csv"))
            write.csv(misspec_results, csv_file, row.names = FALSE)
            
            cat("🎨 Creating per-combination plot...\n")
            plot_obj <- create_combined_misspec_plot(misspec_results, abm_data, R0_val, model_type, transition)
            plot_file <- file.path(results_dir, paste0("combined_plot_", combination_name, ".png"))
            ggsave(plot_file, plot_obj, width = 14, height = 10, dpi = 300)
            
            processing_summary <- bind_rows(processing_summary, data.frame(
              combination = combination_name,
              R0_value = R0_val,
              model_type = model_type,
              transition = transition,
              n_results = nrow(misspec_results),
              status = "✅ Success"
            ))
            
            cat("✅ Completed successfully!\n")
            cat("   • Results:", nrow(misspec_results), "rows\n")
            cat("   • Simulations:", length(unique(misspec_results$sim_id)), "\n")
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
      filter(date >= 10 & date <= DEFAULT_DATE_MAX) %>%
      group_by(method, reporting_rate, R0_value, model_type, transition_type) %>%
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
    
    cat("\n💾 Combined results saved:", nrow(all_results), "total rows\n")
  }
  
  summary_csv <- file.path(results_dir, "processing_summary.csv")
  write.csv(processing_summary, summary_csv, row.names = FALSE)
  
  cat("\n📊 PROCESSING SUMMARY:\n")
  print(processing_summary)
  
  cat("\n")
  cat("╔══════════════════════════════════════════════════════════════╗\n")
  cat("║                    MISSPEC ANALYSIS COMPLETE!                ║\n")
  cat("╚══════════════════════════════════════════════════════════════╝\n")
  
  cat("\n📁 Results saved in:", results_dir, "/\n")
  cat("├── 💾 combined_all_results.rds/.csv - All results\n") 
  cat("├── 📊 summary_statistics.csv        - Summary stats\n")
  cat("├── 📋 processing_summary.csv        - Processing status\n")
  cat("├── 📈 Individual result files per combo (.rds/.csv)\n")
  cat("└── 🖼️  combined_plot_<combo>.png     - ABM + Cori + Wallinga\n")
  
  invisible(list(results = all_results, summary = processing_summary))
}

# ====== COMPREHENSIVE PLOTS FROM SAVED MISSPEC RESULTS ==========================

.rate_label <- function(x) paste0(round(100 * x), "%")
.make_linetypes <- function(rate_labels) {
  base_types <- c("solid","dotted","dashed","twodash","longdash","dotdash")
  rate_labels <- unique(rate_labels)
  if ("100%" %in% rate_labels) {
    others <- setdiff(rate_labels, "100%")
    types <- c("100%" = "solid")
    if (length(others) > 0) types <- c(types, setNames(base_types[2:(length(others)+1)], others))
  } else {
    types <- setNames(base_types[seq_along(rate_labels)], rate_labels)
  }
  types
}

create_comprehensive_plots_from_misspec_results <- function(
    results_dir   = DEFAULT_RESULTS_DIR,
    plots_dir     = DEFAULT_PLOTS_DIR,
    R0_values     = DEFAULT_R0_VALUES,
    model_types   = DEFAULT_MODEL_TYPES,
    transitions   = DEFAULT_TRANSITIONS,
    date_min      = DEFAULT_DATE_MIN,
    date_max      = DEFAULT_DATE_MAX
) {
  cat("\n=== COMPREHENSIVE PLOTS (MISSPECIFICATION / UNDERREPORTING) ===\n")
  if (!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE)
  
  method_colors <- c("ABM" = "blue", "Cori" = "red", "Wallinga-Teunis" = "green")
  
  for (model_type in model_types) {
    for (transition in transitions) {
      cat("\n--- Creating plots for", model_type, "model,", transition, "(from", results_dir, ") ---\n")
      
      plots_list <- vector("list", length(R0_values))
      all_data_combined <- tibble()
      
      for (i in seq_along(R0_values)) {
        R0_val <- R0_values[i]
        combo   <- paste0(model_type, "_R0_", R0_val, "_", transition)
        misspec_file <- file.path(results_dir, paste0("misspec_results_", combo, ".rds"))
        abm_file <- file.path("saved_data", paste0(model_type, "_R0_", R0_val, "_n_1e+05_nsim_1000_rt_ci.rds"))        
        if (file.exists(misspec_file)) {
          cat("  • Loading misspec results for R0 =", R0_val, "\n")
          misspec_results <- readRDS(misspec_file)
          
          cw_plot_data <- misspec_results %>%
            filter(date >= date_min, date <= date_max) %>%
            group_by(method, date, reporting_rate) %>%
            summarise(
              mean_R   = mean(median_rt, na.rm = TRUE),
              ci_lower = quantile(median_rt, 0.025, na.rm = TRUE),
              ci_upper = quantile(median_rt, 0.975, na.rm = TRUE),
              .groups  = "drop"
            ) %>%
            mutate(
              R0_value   = R0_val,
              rate_label = .rate_label(reporting_rate)
            )
          
          abm_plot_data <- tibble()
          if (file.exists(abm_file)) {
            abm_plot_data <- readRDS(abm_file) %>%
              rename(date = source_exposure_date) %>%
              filter(date >= date_min, date <= date_max) %>%
              transmute(
                method         = "ABM",
                date,
                mean_R         = mean_rt,
                ci_lower, ci_upper,
                reporting_rate = 1.0,
                rate_label     = "100%",
                R0_value       = R0_val
              )
          } else {
            cat("    ⚠️ ABM baseline not found for", combo, "— plotting methods only\n")
          }
          
          combined_data <- bind_rows(cw_plot_data, abm_plot_data) %>%
            mutate(method = factor(method, levels = c("ABM","Cori","Wallinga-Teunis")))
          
          all_data_combined <- bind_rows(all_data_combined, combined_data)
          
          ltypes <- .make_linetypes(unique(combined_data$rate_label))
          
          p_individual <- ggplot(
            combined_data,
            aes(x = date, y = mean_R, color = method, fill = method)
          ) +
            geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.20, color = NA) +
            geom_line(aes(linetype = rate_label), linewidth = 1.1) +
            geom_hline(yintercept = R0_val, linetype = "dashed", color = "black", linewidth = 0.9) +
            scale_color_manual(values = method_colors) +
            scale_fill_manual(values = method_colors, guide = "none") +
            scale_linetype_manual(values = ltypes, name = "Reporting rate") +
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
          cat("  ✗ Misspec results not found for", combo, "\n")
          plots_list[[i]] <- ggplot() +
            annotate("text", x = 0.5, y = 0.5,
                     label = paste("No misspec_results for\nR0 =", R0_val),
                     size = 6) +
            theme_void() +
            labs(title = paste0("R0 = ", R0_val)) +
            theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
        }
      }
      
      # 2×2 grid (panels per R0)
      if (length(plots_list) == 4) {
        combined_grid <- ggpubr::ggarrange(
          plots_list[[1]], plots_list[[2]],
          plots_list[[3]], plots_list[[4]],
          ncol = 2, nrow = 2, common.legend = FALSE
        )
        final_grid <- ggpubr::annotate_figure(
          combined_grid,
          top = ggpubr::text_grob(
            paste0(
              "R_t under Underreporting Misspecification: ABM vs Cori vs Wallinga-Teunis\n",
              toupper(model_type), " Model | ",
              gsub("_to_", " \u2192 ", toupper(transition)), " Transition"
            ),
            face = "bold", size = 16
          )
        )
        grid_filename <- file.path(plots_dir, paste0("grid_", model_type, "_", transition, ".png"))
        ggsave(grid_filename, final_grid, width = 16, height = 12, dpi = 300)
        cat("  ✓ 2x2 grid saved:", basename(grid_filename), "\n")
      }
      
      # Faceted plot across R0 (keeps reporting-rate linetypes)
      if (nrow(all_data_combined) > 0) {
        ltypes_all <- .make_linetypes(unique(all_data_combined$rate_label))
        all_data_combined$method <- factor(all_data_combined$method, levels = c("ABM","Cori","Wallinga-Teunis"))
        
        p_facet <- ggplot(
          all_data_combined,
          aes(x = date, y = mean_R, color = method, fill = method, linetype = rate_label)
        ) +
          geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.20, color = NA) +
          geom_line(linewidth = 1.05) +
          geom_hline(aes(yintercept = R0_value), linetype = "dashed", color = "black", linewidth = 0.8) +
          scale_color_manual(values = method_colors, name = "Method") +
          scale_fill_manual(values = method_colors, guide = "none") +
          scale_linetype_manual(values = ltypes_all, name = "Reporting rate") +
          scale_x_continuous(limits = c(0, date_max), breaks = seq(0, date_max, by = 25)) +
          facet_wrap(~paste0("R0 = ", R0_value), ncol = 2, scales = "free_y") +
          labs(
            title = paste0(
              "R_t Across R0 Values under Underreporting Misspecification\n",
              toupper(model_type), " Model | ",
              gsub("_to_", " \u2192 ", toupper(transition)), " Transition"
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
        
        # Summary stats by method & reporting rate
        summary_stats <- all_data_combined %>%
          filter(date >= 10, date <= date_max) %>%
          group_by(method, reporting_rate, rate_label, R0_value) %>%
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
        
        # Overall ranking averaged over R0 by method & reporting rate
        ranking <- summary_stats %>%
          group_by(method, rate_label) %>%
          summarise(
            avg_abs_bias = mean(abs_bias, na.rm = TRUE),
            avg_rmse     = mean(rmse,     na.rm = TRUE),
            .groups = "drop"
          ) %>%
          arrange(avg_rmse, avg_abs_bias)
        
        ranking_filename <- file.path(plots_dir, paste0("ranking_", model_type, "_", transition, ".csv"))
        write_csv(ranking, ranking_filename)
        cat("  ✓ Ranking saved:", basename(ranking_filename), "\n")
      }
    }
  }
  
  cat("\n=== DONE: COMPREHENSIVE PLOTS (MISSPECIFICATION) ===\n")
}

# ====== QUICK DIAGNOSTICS (optional) ============================================

test_single_simulation_misspec <- function(model_type = "full", R0_val = 2.0, 
                                           sim_id = 1, transition = "susceptible_to_exposed",
                                           start_day = 0) {
  cat("╔══════════════════════════════════════════════════════════════╗\n")
  cat("║            TESTING SINGLE SIMULATION (Misspec + SI)         ║\n")
  cat("╚══════════════════════════════════════════════════════════════╝\n")
  inspect_data_directory("saved_data/")
  
  possible_prefixes <- c(
    paste0("saved_data/", model_type, "_R0_", R0_val, "_n_1e+05_nsim_100"),
    paste0("saved_data/", model_type, "_R0_", gsub("\\.", "_", R0_val), "_n_1e+05_nsim_100"),
    paste0("saved_data/", model_type, "_R0_", R0_val)
  )
  
  saved_data <- NULL; file_prefix <- NULL
  cat("\n📂 Trying different file name patterns...\n")
  for (prefix in possible_prefixes) {
    cat("Trying:", prefix, "\n")
    if (file.exists(paste0(prefix, "_metadata.rds"))) { file_prefix <- prefix; cat("  ✓ Found!\n"); break } else { cat("  ✗ No match\n") }
  }
  if (is.null(file_prefix)) { cat("❌ Could not find data files.\n"); return(NULL) }
  
  cat("📂 Loading test data:", file_prefix, "\n")
  saved_data <- load_seir_data(file_prefix)
  
  abm_file <- paste0(file_prefix, "_rt_ci.rds")
  abm_data <- if (file.exists(abm_file)) { readRDS(abm_file) } else { cat("⚠️  ABM data not found\n"); NULL }
  
  test_incidence <- c(0, 2, 5, 8, 12, 15, 10, 8, 5, 3, 2, 1, 0, 0, 0)
  cat("\n🧪 Underreporting sample:\n")
  cat("Original:", paste(test_incidence, collapse = ", "), "\n")
  set.seed(DEFAULT_SEED)
  underreported <- apply_underreporting(test_incidence, 0.7)
  cat("70% rate:", paste(underreported, collapse = ", "), "\n")
  
  misspec_results <- analyze_misspec_single_sim(
    saved_data, sim_id = sim_id, how_infected = transition,
    reporting_rate = 0.7, start_day = start_day
  )
  if (!is.null(misspec_results) && nrow(misspec_results) > 0) {
    cat("✅ Single simulation OK (", nrow(misspec_results), " rows )\n", sep = "")
    plot_obj <- create_combined_misspec_plot(misspec_results, abm_data, R0_val, model_type, transition)
    out <- paste0("test_misspec_plot_", model_type, "_R0_", R0_val, "_sim_", sim_id, ".png")
    ggsave(out, plot_obj, width = 12, height = 8, dpi = 300)
    cat("📊 Test plot saved:", out, "\n")
  } else {
    cat("❌ Single simulation failed - no results\n")
  }
  invisible(misspec_results)
}

# ====== AUTO-RUN (so `Rscript misspec_complete.R` just works) ====================
if (sys.nframe() == 0) {
  # 1) Run the full misspec analysis and write results + per-combo plots
  run_complete_misspec_analysis(
    R0_values = DEFAULT_R0_VALUES,
    model_types = DEFAULT_MODEL_TYPES,
    transitions = DEFAULT_TRANSITIONS,
    reporting_rate = DEFAULT_REPORTING_RATE,
    max_sims_per_scenario = DEFAULT_MAX_SIMS,
    start_day = DEFAULT_START_DAY,
    results_dir = DEFAULT_RESULTS_DIR,
    rng_seed = DEFAULT_SEED
  )
  
  # 2) Build comprehensive multi-R0 plots from the saved results
  create_comprehensive_plots_from_misspec_results(
    results_dir = DEFAULT_RESULTS_DIR,
    plots_dir   = DEFAULT_PLOTS_DIR,
    R0_values   = DEFAULT_R0_VALUES,
    model_types = DEFAULT_MODEL_TYPES,
    transitions = DEFAULT_TRANSITIONS,
    date_min    = DEFAULT_DATE_MIN,
    date_max    = DEFAULT_DATE_MAX
  )
  
  cat("\nAll done! See:\n  -", DEFAULT_RESULTS_DIR, " for results & per-combo plots\n  -", DEFAULT_PLOTS_DIR, " for 2×2 grid, faceted plots, and summaries\n\n")
}
