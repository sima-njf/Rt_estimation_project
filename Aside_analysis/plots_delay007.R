library(ggplot2)
library(dplyr)
library(patchwork)
library(tools)

# File directory for delay results
results_dir <- "complete_delay_results/"
files <- list.files(results_dir, pattern = "^delay_results_.*\\.csv$", full.names = TRUE)

# Initialize plot list
plot_list <- list()

# Loop through all files
for (f in sort(files)) {
  df <- read.csv(f)
  
  R0_val <- as.numeric(gsub(".*_R0_([0-9.]+)_.*", "\\1", f))
  model_type_raw <- if (grepl("full", f)) "full" else "partial"
  
  # Load ABM data
  abm_file <- paste0("saved_data/", model_type_raw, "_R0_", R0_val, "_n_1e+05_nsim_100_rt_ci.csv")
  abm_plot_data <- NULL
  if (file.exists(abm_file)) {
    abm_data <- read.csv(abm_file)
    abm_plot_data <- abm_data %>%
      rename(date = source_exposure_date) %>%
      filter(date <= 75) %>%
      mutate(
        mean_R = mean_rt,
        method = "ABM (True)",
        delay_scenario = "true_baseline",
        delay_percentage = 0,
        n_sims = 1
      ) %>%
      select(method, date, mean_R, ci_lower, ci_upper, delay_scenario, delay_percentage, n_sims)
  }
  
  # Compute summary across sims
  plot_data <- df %>%
    filter(date <= 75) %>%
    group_by(method, date, delay_scenario, delay_percentage) %>%
    summarise(
      mean_R = mean(median_rt, na.rm = TRUE),
      ci_lower = quantile(median_rt, 0.025, na.rm = TRUE),
      ci_upper = quantile(median_rt, 0.975, na.rm = TRUE),
      n_sims = n(),
      .groups = "drop"
    )
  
  # Combine with ABM data
  if (!is.null(abm_plot_data)) {
    plot_data <- bind_rows(plot_data, abm_plot_data)
  }
  
  plot_data$method <- factor(plot_data$method, levels = c("ABM (True)", "Cori", "Wallinga-Teunis"))
  
  # Line type settings
  plot_data$delay_scenario <- factor(plot_data$delay_scenario,
                                     levels = c("true_baseline", "medium_delay", "high_delay"))
  
  p <- ggplot(plot_data, aes(x = date, y = mean_R, color = method, fill = method)) +
    geom_line(aes(linetype = delay_scenario), size = 1.2) +
    geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.12, linetype = 0) +
    geom_hline(yintercept = R0_val, linetype = "dashed", color = "gray30", size = 1) +
    scale_color_manual(values = c("ABM (True)" = "black", "Cori" = "red", "Wallinga-Teunis" = "blue")) +
    scale_fill_manual(values = c("ABM (True)" = "black", "Cori" = "red", "Wallinga-Teunis" = "blue")) +
    scale_linetype_manual(
      values = c("true_baseline" = "solid", "medium_delay" = "dotted", "high_delay" = "dashed"),
      labels = c("ABM (True)", "Medium Delay (30%)", "High Delay (50%)")
    ) +
    labs(
      title = bquote(R[0] ~ "=" ~ .(R0_val)),
      x = "Day", y = expression(R[t]),
      linetype = "Delay Scenario"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold", hjust = 0.5, size = 13),
      axis.title.x = element_text(size = 12),
      axis.title.y = element_text(size = 12),
      axis.text = element_text(size = 11)
    )
  
  plot_list[[length(plot_list) + 1]] <- p
}

# Separate full and partial panels
full_plots <- plot_list[1:4]
partial_plots <- plot_list[5:8]

# ────────────────────────────────
# 🟦 Full Contact Panel
# ────────────────────────────────
full_panel <- (full_plots[[1]] | full_plots[[2]]) /
  (full_plots[[3]] | full_plots[[4]]) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

full_panel <- full_panel + plot_annotation(
  title = "Rt Estimation Under Reporting Delays — Full Contact Model",
  subtitle = "ABM (True), Cori, Wallinga-Teunis | Delays: None (solid), 30% (dotted), 50% (dashed)",
  theme = theme(
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 14, hjust = 0.5, margin = margin(b = 10))
  )
)

ggsave("Rt_delay_full_contact_panel.png", full_panel, width = 14, height = 11, dpi = 300)
cat("✅ Saved: Rt_delay_full_contact_panel.png\n")

# ────────────────────────────────
# 🟥 Partial Contact Panel
# ────────────────────────────────
partial_panel <- (partial_plots[[1]] | partial_plots[[2]]) /
  (partial_plots[[3]] | partial_plots[[4]]) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

partial_panel <- partial_panel + plot_annotation(
  title = "Rt Estimation Under Reporting Delays — Partial Contact Model",
  subtitle = "ABM (True), Cori, Wallinga-Teunis | Delays: None (solid), 30% (dotted), 50% (dashed)",
  theme = theme(
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 14, hjust = 0.5, margin = margin(b = 10))
  )
)

ggsave("Rt_delay_partial_contact_panel.png", partial_panel, width = 14, height = 11, dpi = 300)
cat("✅ Saved: Rt_delay_partial_contact_panel.png\n")
