# ================================================================================
# DEGREE COMPARISON ANALYSIS - PARTIAL CONTACT MODEL
# Comparing network degree effects on R_t estimation
# Focus: R0 = 2.0 and 3.0 with degree = 10 and 20
# ================================================================================

suppressPackageStartupMessages({
  library(epiworldR)
  library(data.table)
  library(dplyr)
  library(tidyverse)
  library(igraph)
  library(EpiEstim)
  library(ggplot2)
  library(ggpubr)
  library(purrr)
})

# ── CONFIGURATION ───────────────────────────────────────────────────────────────
R0_VALUES <- c(2.0, 3.0)
DEGREES <- c(10, 20)
N <- 1e5
PREVALENCE <- 100/1e5
RECOVERY_RATE <- 1.0/7.0
INCUBATION_DAYS <- 4
NDAYS <- 150
NSIM <- 100
SEED <- 1234
SAVE_DIR <- "degree_comparison_data/"
PLOTS_DIR <- "degree_comparison_plots/"

# Create directories
if (!dir.exists(SAVE_DIR)) dir.create(SAVE_DIR, recursive = TRUE)
if (!dir.exists(PLOTS_DIR)) dir.create(PLOTS_DIR, recursive = TRUE)

# ── STEP 1: GENERATE SIMULATION DATA ───────────────────────────────────────────

generate_partial_data <- function(R0, deg, n = N, prevalence = PREVALENCE,
                                  recovery_rate = RECOVERY_RATE,
                                  incubation_days = INCUBATION_DAYS,
                                  ndays = NDAYS, nsim = NSIM, seed = SEED) {
  
  cat("\n╔══════════════════════════════════════════════════════════════╗\n")
  cat("║  Generating PARTIAL contact data: R0 =", R0, ", Degree =", deg, "  ║\n")
  cat("╚══════════════════════════════════════════════════════════════╝\n")
  
  set.seed(seed)
  
  # Generate Erdős–Rényi graph
  p_er <- deg / (n - 1)
  g <- erdos.renyi.game(n = n, p.or.m = p_er, type = "gnp", directed = FALSE)
  net <- as_edgelist(g, names = FALSE)
  net <- apply(net, 2, as.integer) - 1L  # convert to 0-indexed integers
  
  cat("Network created:\n")
  cat("  • Nodes:", n, "\n")
  cat("  • Edges:", nrow(net), "\n")
  cat("  • Target degree:", deg, "\n")
  cat("  • Actual avg degree:", round(2 * nrow(net) / n, 2), "\n")
  
  # Calculate transmission rate
  transmission_rate <- R0 * recovery_rate / (deg * 2)
  cat("  • Transmission rate:", transmission_rate, "\n")
  
  # Create SEIR model
  model <- ModelSEIR(
    name = paste("Covid_ER_deg", deg),
    prevalence = prevalence,
    transmission_rate = transmission_rate,
    incubation_days = incubation_days,
    recovery_rate = recovery_rate
  )
  
  # Create agents from network
  agents_from_edgelist(
    model = model,
    source = net[, 1],
    target = net[, 2],
    size = as.integer(n),
    directed = FALSE
  )
  
  # Create saver
  saver <- make_saver(
    "total_hist",
    "transmission",
    "transition",
    "reproductive",
    "generation"
  )
  
  # Run simulations
  cat("\nRunning", nsim, "simulations for", ndays, "days...\n")
  run_multiple(model, ndays = ndays, nsim = nsim, seed = seed, saver = saver, nthreads = 18)
  
  # Get results
  cat("Extracting results...\n")
  results <- run_multiple_get_results(model, nthreads = 18)
  
  # Prepare transition data
  df_transitions <- results$transition %>%
    rename(id = sim_num) %>%
    unite("transition", from, to, sep = " to ") %>%
    pivot_wider(names_from = transition, values_from = counts, values_fill = list(counts = 0))
  
  # Clean column names
  colnames(df_transitions) <- gsub(" ", "_", colnames(df_transitions))
  colnames(df_transitions) <- tolower(colnames(df_transitions))
  
  # Floor the date and transition counts
  df_transitions$date <- floor(df_transitions$date)
  transition_cols <- grep("_to_", colnames(df_transitions), value = TRUE)
  for (col in transition_cols) {
    if (col %in% colnames(df_transitions)) {
      df_transitions[[col]] <- floor(df_transitions[[col]])
    }
  }
  
  # Get reproductive numbers
  reproductive_data <- results$reproductive %>%
    filter(source >= 0)
  
  # Calculate Rt confidence intervals
  ci_results_rt <- reproductive_data %>%
    group_by(sim_num, source_exposure_date) %>%
    summarise(rt = mean(rt), .groups = "drop") %>%
    group_by(source_exposure_date) %>%
    summarise(
      mean_rt = mean(rt),
      ci_lower = quantile(rt, probs = 0.025),
      ci_upper = quantile(rt, probs = 0.975),
      n = n(),
      .groups = "drop"
    )
  
  # Get generation interval data
  generation_data <- results$generation %>%
    filter(gentime > 0)
  
  # Calculate generation interval statistics
  gi_stats <- generation_data %>%
    group_by(sim_num) %>%
    summarise(
      mean_gi = mean(gentime),
      sd_gi = sd(gentime),
      median_gi = median(gentime),
      n_gi = n(),
      .groups = "drop"
    )
  
  # Overall GI stats
  overall_gi_stats <- list(
    mean_gi = mean(generation_data$gentime),
    sd_gi = sd(generation_data$gentime),
    median_gi = median(generation_data$gentime),
    n_total = nrow(generation_data)
  )
  
  # Metadata
  metadata <- list(
    model_type = "partial",
    R0 = R0,
    degree = deg,
    n = n,
    prevalence = prevalence,
    transmission_rate = transmission_rate,
    recovery_rate = recovery_rate,
    incubation_days = incubation_days,
    ndays = ndays,
    nsim = nsim,
    seed = seed,
    run_date = Sys.time(),
    overall_gi_stats = overall_gi_stats,
    network_type = "Erdos-Renyi",
    avg_degree = deg
  )
  
  # Save files
  file_prefix <- paste0(SAVE_DIR, "partial_R0_", R0, "_deg_", deg, "_nsim_", nsim)
  
  cat("\nSaving results to:", file_prefix, "\n")
  saveRDS(df_transitions, paste0(file_prefix, "_transitions.rds"))
  saveRDS(reproductive_data, paste0(file_prefix, "_reproductive.rds"))
  saveRDS(ci_results_rt, paste0(file_prefix, "_rt_ci.rds"))
  saveRDS(generation_data, paste0(file_prefix, "_generation.rds"))
  saveRDS(gi_stats, paste0(file_prefix, "_gi_stats.rds"))
  saveRDS(results$transmission, paste0(file_prefix, "_transmission.rds"))
  saveRDS(results$total_hist, paste0(file_prefix, "_total_hist.rds"))
  saveRDS(metadata, paste0(file_prefix, "_metadata.rds"))
  
  cat("✅ Data saved successfully!\n")
  
  list(
    transitions = df_transitions,
    transmission = results$transmission,
    reproductive = reproductive_data,
    rt_ci = ci_results_rt,
    generation = generation_data,
    gi_stats = gi_stats,
    metadata = metadata
  )
}

# ── STEP 2: SERIAL INTERVAL CALCULATION ────────────────────────────────────────

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

# ── STEP 3: APPLY CORI AND WALLINGA-TEUNIS ─────────────────────────────────────

apply_cori_single_sim <- function(saved_data, sim_id, how_infected = "susceptible_to_exposed", 
                                  start_day = 0) {
  
  sim_data <- saved_data$transitions %>%
    filter(id == sim_id, date >= start_day) %>%
    arrange(date)
  
  if (!how_infected %in% colnames(sim_data)) return(NULL)
  
  incidence <- sim_data[[how_infected]]
  if (sum(incidence, na.rm = TRUE) == 0 || length(incidence) < 7) return(NULL)
  
  # Get serial intervals
  si_all <- calculate_serial_intervals(saved_data$transmission)
  si_data <- si_all %>% filter(sim_num == sim_id)
  if (nrow(si_data) == 0) return(NULL)
  
  SI_mean <- si_data$mean_si
  SI_sd <- si_data$sd_si
  if (is.na(SI_sd) || SI_sd < 0.1 || SI_mean < 0.5 || SI_mean > 20) return(NULL)
  
  tryCatch({
    result <- estimate_R(
      incidence,
      method = "parametric_si",
      config = make_config(mean_si = SI_mean, std_si = SI_sd)
    )
    
    data.frame(
      sim_id = sim_id,
      date = result$R[, "t_start"],
      median_rt = result$R[, "Median(R)"],
      ci_lower = result$R[, "Quantile.0.025(R)"],
      ci_upper = result$R[, "Quantile.0.975(R)"],
      mean_si = SI_mean,
      std_si = SI_sd,
      method = "Cori"
    )
  }, error = function(e) NULL)
}

apply_wallinga_single_sim <- function(saved_data, sim_id, how_infected = "susceptible_to_exposed",
                                      start_day = 0) {
  
  sim_data <- saved_data$transitions %>%
    filter(id == sim_id, date >= start_day) %>%
    arrange(date)
  
  if (!how_infected %in% colnames(sim_data)) return(NULL)
  
  incidence <- sim_data[[how_infected]]
  if (sum(incidence, na.rm = TRUE) == 0 || length(incidence) < 7) return(NULL)
  
  # Get serial intervals
  si_all <- calculate_serial_intervals(saved_data$transmission)
  si_data <- si_all %>% filter(sim_num == sim_id)
  if (nrow(si_data) == 0) return(NULL)
  
  SI_mean <- si_data$mean_si
  SI_sd <- si_data$sd_si
  if (is.na(SI_sd) || SI_sd < 0.1 || SI_mean < 0.5 || SI_mean > 20) return(NULL)
  
  tryCatch({
    result <- wallinga_teunis(
      incid = incidence,
      method = "parametric_si",
      config = make_config(
        incid = incidence,
        method = "parametric_si",
        mean_si = SI_mean,
        std_si = SI_sd
      )
    )
    
    data.frame(
      sim_id = sim_id,
      date = result$R[, "t_start"],
      median_rt = result$R[, "Mean(R)"],
      ci_lower = result$R[, "Quantile.0.025(R)"],
      ci_upper = result$R[, "Quantile.0.975(R)"],
      mean_si = SI_mean,
      std_si = SI_sd,
      method = "Wallinga-Teunis"
    )
  }, error = function(e) NULL)
}

apply_methods_all_sims <- function(saved_data, how_infected = "susceptible_to_exposed",
                                   start_day = 0, verbose = TRUE) {
  
  sim_ids <- unique(saved_data$transitions$id)
  
  if (verbose) {
    cat("Applying Cori and Wallinga-Teunis to", length(sim_ids), "simulations\n")
  }
  
  all_results <- data.frame()
  
  for (i in seq_along(sim_ids)) {
    sim_id <- sim_ids[i]
    
    if (verbose && i %% 10 == 0) {
      cat("  Processing simulation", i, "of", length(sim_ids), "\n")
    }
    
    cori_res <- apply_cori_single_sim(saved_data, sim_id, how_infected, start_day)
    wt_res <- apply_wallinga_single_sim(saved_data, sim_id, how_infected, start_day)
    
    all_results <- bind_rows(all_results, cori_res, wt_res)
  }
  
  if (verbose) {
    cat("✅ Completed! Total results:", nrow(all_results), "rows\n")
  }
  
  all_results
}

# ── STEP 4: CREATE COMPARISON PLOTS ─────────────────────────────────────────────

create_degree_comparison_plot <- function(R0_val, transition = "susceptible_to_exposed",
                                          date_min = 12, date_max = 75) {
  
  cat("\n📊 Creating comparison plot for R0 =", R0_val, "\n")
  
  all_data <- data.frame()
  
  for (deg in DEGREES) {
    cat("  Loading degree =", deg, "...\n")
    
    file_prefix <- paste0(SAVE_DIR, "partial_R0_", R0_val, "_deg_", deg, "_nsim_", NSIM)
    
    # Check if files exist
    if (!file.exists(paste0(file_prefix, "_metadata.rds"))) {
      cat("    ⚠️ Data not found, skipping\n")
      next
    }
    
    # Load saved data
    saved_data <- list(
      transitions = readRDS(paste0(file_prefix, "_transitions.rds")),
      transmission = readRDS(paste0(file_prefix, "_transmission.rds")),
      rt_ci = readRDS(paste0(file_prefix, "_rt_ci.rds")),
      metadata = readRDS(paste0(file_prefix, "_metadata.rds"))
    )
    
    # Apply Cori and Wallinga-Teunis
    cat("    Applying estimation methods...\n")
    cw_results <- apply_methods_all_sims(saved_data, how_infected = transition, 
                                         start_day = 0, verbose = FALSE)
    
    if (!is.null(cw_results) && nrow(cw_results) > 0) {
      # Aggregate across simulations
      cw_plot_data <- cw_results %>%
        filter(date >= date_min, date <= date_max) %>%
        group_by(method, date) %>%
        summarise(
          mean_rt = mean(median_rt, na.rm = TRUE),
          ci_lower = quantile(median_rt, 0.025, na.rm = TRUE),
          ci_upper = quantile(median_rt, 0.975, na.rm = TRUE),
          .groups = "drop"
        ) %>%
        mutate(degree = deg, data_source = "Estimation")
      
      # Add ABM data
      abm_plot_data <- saved_data$rt_ci %>%
        rename(date = source_exposure_date) %>%
        filter(date >= date_min, date <= date_max) %>%
        transmute(
          method = "ABM",
          date,
          mean_rt,
          ci_lower,
          ci_upper,
          degree = deg,
          data_source = "ABM"
        )
      
      all_data <- bind_rows(all_data, cw_plot_data, abm_plot_data)
    }
  }
  
  if (nrow(all_data) == 0) {
    cat("  ❌ No data available for plotting\n")
    return(NULL)
  }
  
  # Create plot
  all_data$method <- factor(all_data$method, levels = c("ABM", "Cori", "Wallinga-Teunis"))
  all_data$degree_label <- paste0("Degree = ", all_data$degree)
  
  method_colors <- c("ABM" = "blue", "Cori" = "red", "Wallinga-Teunis" = "green")
  
  p <- ggplot(all_data, aes(x = date, y = mean_rt, color = method, fill = method)) +
    geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.2, color = NA) +
    geom_line(size = 1.2) +
    geom_hline(yintercept = R0_val, linetype = "dashed", color = "black", size = 1) +
    facet_wrap(~degree_label, ncol = 2) +
    scale_color_manual(values = method_colors) +
    scale_fill_manual(values = method_colors) +
    scale_x_continuous(limits = c(0, date_max), breaks = seq(0, date_max, by = 25)) +
    labs(
      title = paste0("R_t Estimates by Network Degree: R0 = ", R0_val),
      subtitle = paste0(
        "PARTIAL Model (Erdős-Rényi) | ",
        gsub("_to_", " → ", toupper(transition)), " Transition"
      ),
      x = "Day",
      y = "R_t Estimate",
      color = "Method",
      fill = "Method"
    ) +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 12),
      strip.text = element_text(size = 12, face = "bold"),
      strip.background = element_rect(fill = "lightblue"),
      panel.grid.minor = element_blank()
    ) +
    annotate("text", x = date_max - 15, y = R0_val + 0.3,
             label = paste0("R0 = ", R0_val), color = "black", size = 4)
  
  # Save plot
  filename <- paste0(PLOTS_DIR, "degree_comparison_R0_", R0_val, "_", transition, ".png")
  ggsave(filename, p, width = 14, height = 8, dpi = 300)
  cat("  ✅ Plot saved:", filename, "\n")
  
  # Calculate summary statistics
  summary_stats <- all_data %>%
    filter(date >= 10, date <= date_max) %>%
    group_by(method, degree) %>%
    summarise(
      mean_rt_overall = mean(mean_rt, na.rm = TRUE),
      median_rt_overall = median(mean_rt, na.rm = TRUE),
      sd_rt = sd(mean_rt, na.rm = TRUE),
      bias = mean_rt_overall - R0_val,
      abs_bias = abs(bias),
      rmse = sqrt(mean((mean_rt - R0_val)^2, na.rm = TRUE)),
      .groups = "drop"
    )
  
  summary_filename <- paste0(PLOTS_DIR, "summary_R0_", R0_val, "_", transition, ".csv")
  write.csv(summary_stats, summary_filename, row.names = FALSE)
  cat("  ✅ Summary stats saved:", summary_filename, "\n")
  
  print(summary_stats)
  
  invisible(list(plot = p, summary = summary_stats, data = all_data))
}

# ── STEP 5: MAIN EXECUTION ──────────────────────────────────────────────────────

run_degree_comparison_analysis <- function() {
  
  cat("╔══════════════════════════════════════════════════════════════╗\n")
  cat("║         DEGREE COMPARISON ANALYSIS - PARTIAL MODEL          ║\n")
  cat("║            R0 = 2.0, 3.0 | Degree = 10, 20                  ║\n")
  cat("╚══════════════════════════════════════════════════════════════╝\n")
  
  # Generate all simulation data
  cat("\n=== STEP 1: GENERATING SIMULATION DATA ===\n")
  
  for (R0_val in R0_VALUES) {
    for (deg in DEGREES) {
      generate_partial_data(R0 = R0_val, deg = deg)
    }
  }
  
  # Create comparison plots
  cat("\n=== STEP 2: CREATING COMPARISON PLOTS ===\n")
  
  for (R0_val in R0_VALUES) {
    create_degree_comparison_plot(R0_val, transition = "susceptible_to_exposed")
  }
  
  cat("\n╔══════════════════════════════════════════════════════════════╗\n")
  cat("║                    ANALYSIS COMPLETE!                       ║\n")
  cat("╚══════════════════════════════════════════════════════════════╝\n")
  
  cat("\n📁 Results saved in:\n")
  cat("  • Data:", SAVE_DIR, "\n")
  cat("  • Plots:", PLOTS_DIR, "\n\n")
}

# ── RUN THE ANALYSIS ────────────────────────────────────────────────────────────

# Uncomment to run:
run_degree_comparison_analysis()
