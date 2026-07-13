# ================================================================================
# PUBLICATION FIGURES - "Estimating the Reproduction Number" (static, for LaTeX)
# ggplot2 + faceting. Confidence bands for ALL series.
# Warm/cool method-family palette: Cori = warm, Wallinga-Teunis = cool, ABM = black.
# Color + linetype so series stay distinct in grayscale.
# Outputs vector PDF and 300-dpi PNG per (model, transition).
# ================================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
})

# -- CONFIG ----------------------------------------------------------------------
OUT_DIR      <- "paper_figures"
DATE_MIN     <- 12
DATE_MAX     <- 75

R0_VALUES    <- c(1.5, 2.0, 3.0, 5.0)
MODEL_TYPES  <- c("full", "partial")
TRANSITIONS  <- c("susceptible_to_exposed", "exposed_to_infected")

BASELINE_DIR <- "cori_wallinga_results"
DELAY_DIR    <- "complete_delay_results"
MISSPEC_DIR  <- "complete_misspec_results"
ABM_DIR      <- "saved_data"

# -- SERIES ORDER, COLORS, LINETYPES ---------------------------------------------
SERIES_LEVELS <- c(
  "ABM (True Rt)",
  "Cori: Real Data",
  "Wallinga-Teunis: Real Data",
  "Cori: 50% Delayed",
  "Wallinga-Teunis: 50% Delayed",
  "Cori: 70% Reporting",
  "Wallinga-Teunis: 70% Reporting"
)

# Warm/cool palette: Cori = warm (reds/oranges), WT = cool (blues/greens), ABM = near-black
SERIES_COLORS <- c(
  "ABM (True Rt)"                  = "#1A1A1A",  # near-black
  "Cori: Real Data"               = "#C0392B",  # deep red
  "Cori: 50% Delayed"             = "#E67E22",  # orange
  "Cori: 70% Reporting"           = "#F1C40F",  # amber
  "Wallinga-Teunis: Real Data"    = "#1F6FB2",  # deep blue
  "Wallinga-Teunis: 50% Delayed"  = "#16A085",  # teal green
  "Wallinga-Teunis: 70% Reporting"= "#8E44AD"   # violet
)

SERIES_LINETYPES <- c(
  "ABM (True Rt)"                  = "solid",
  "Cori: Real Data"               = "solid",
  "Wallinga-Teunis: Real Data"    = "solid",
  "Cori: 50% Delayed"             = "dashed",
  "Wallinga-Teunis: 50% Delayed"  = "dashed",
  "Cori: 70% Reporting"           = "dotted",
  "Wallinga-Teunis: 70% Reporting"= "dotted"
)

# -- LOADERS ---------------------------------------------------------------------
safe_load <- function(f, ...) tryCatch(f(...), error = function(e) NULL)

load_baseline <- function(model_type, R0_val, transition) {
  fp <- file.path(BASELINE_DIR, paste0(model_type, "_R0_", R0_val, "_", transition, ".csv"))
  if (!file.exists(fp)) return(NULL)
  read_csv(fp, show_col_types = FALSE) %>% mutate(scenario = "Real Data")
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

# -- BUILD PLOT DATA -------------------------------------------------------------
build_grid_data <- function(model_type, transition) {
  
  rows <- list()
  
  for (R0_val in R0_VALUES) {
    baseline <- safe_load(load_baseline, model_type, R0_val, transition)
    delay    <- safe_load(load_delay,    model_type, R0_val, transition)
    misspec  <- safe_load(load_misspec,  model_type, R0_val, transition)
    abm      <- safe_load(load_abm,      model_type, R0_val)
    
    cw <- bind_rows(baseline, delay, misspec)
    if (!is.null(cw) && nrow(cw) > 0) {
      s <- cw %>%
        filter(date >= DATE_MIN, date <= DATE_MAX) %>%
        group_by(method, date, scenario) %>%
        summarise(
          mean_rt  = mean(median_rt, na.rm = TRUE),
          ci_lower = quantile(median_rt, 0.025, na.rm = TRUE),
          ci_upper = quantile(median_rt, 0.975, na.rm = TRUE),
          .groups  = "drop"
        ) %>%
        mutate(R0 = R0_val, series = paste0(method, ": ", scenario)) %>%
        select(date, mean_rt, ci_lower, ci_upper, R0, series)
      rows[[paste0("cw_", R0_val)]] <- s
    }
    
    if (!is.null(abm)) {
      a <- abm %>%
        rename(date = source_exposure_date) %>%
        filter(date >= DATE_MIN, date <= DATE_MAX) %>%
        transmute(date, mean_rt, ci_lower, ci_upper,
                  R0 = R0_val, series = "ABM (True Rt)")
      rows[[paste0("abm_", R0_val)]] <- a
    }
  }
  
  df <- bind_rows(rows)
  if (nrow(df) == 0) return(df)
  
  present <- SERIES_LEVELS[SERIES_LEVELS %in% unique(df$series)]
  df$series <- factor(df$series, levels = present)
  attr(df, "present") <- present
  df
}

# -- MAKE ONE FIGURE -------------------------------------------------------------
make_figure <- function(model_type, transition) {
  
  df <- build_grid_data(model_type, transition)
  if (nrow(df) == 0) {
    cat("  No data for", model_type, transition, "- skipping\n")
    return(invisible(NULL))
  }
  present <- attr(df, "present")
  
  facet_labeller <- as_labeller(function(x) paste0("R[0] == ", x), label_parsed)
  
  trans_label <- transition
  trans_label <- gsub("susceptible_to_exposed", "S \u2192 E", trans_label)
  trans_label <- gsub("exposed_to_infected",    "E \u2192 I", trans_label)
  model_label <- tools::toTitleCase(model_type)
  
  cols <- SERIES_COLORS[present]
  ltys <- SERIES_LINETYPES[present]
  
  p <- ggplot(df, aes(x = date)) +
    geom_ribbon(
      aes(ymin = ci_lower, ymax = ci_upper, fill = series),
      alpha = 0.12, color = NA
    ) +
    geom_line(
      aes(y = mean_rt, color = series, linetype = series),
      linewidth = 0.75
    ) +
    geom_hline(
      data = data.frame(R0 = R0_VALUES),
      aes(yintercept = R0),
      linetype = "dashed", color = "grey40", linewidth = 0.4,
      inherit.aes = FALSE
    ) +
    facet_wrap(~ R0, ncol = 2, scales = "free_y", labeller = facet_labeller) +
    scale_color_manual(values = cols, name = NULL) +
    scale_fill_manual(values = cols, name = NULL) +
    scale_linetype_manual(values = ltys, name = NULL) +
    scale_x_continuous(limits = c(0, 75), breaks = seq(0, 75, 25),
                       expand = expansion(mult = c(0.01, 0.02))) +
    labs(
      title    = "Estimating the Reproduction Number",
      subtitle = paste0(model_label, " model  |  ", trans_label, " transition"),
      x = "Day",
      y = expression(R[t])
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title       = element_text(face = "bold", size = 14, hjust = 0),
      plot.subtitle    = element_text(size = 11, color = "grey30", hjust = 0,
                                      margin = margin(b = 6)),
      strip.background = element_rect(fill = "grey92", color = NA),
      strip.text       = element_text(face = "bold", size = 11),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
      legend.position  = "bottom",
      legend.key.width = unit(1.8, "lines"),
      legend.text      = element_text(size = 9),
      legend.margin    = margin(t = 4)
    ) +
    guides(
      color    = guide_legend(nrow = 3, byrow = TRUE),
      fill     = guide_legend(nrow = 3, byrow = TRUE),
      linetype = guide_legend(nrow = 3, byrow = TRUE)
    )
  
  base <- file.path(OUT_DIR, paste0("reproduction_number_", model_type, "_", transition))
  ggsave(paste0(base, ".pdf"), p, width = 9, height = 8.5, device = cairo_pdf)
  ggsave(paste0(base, ".png"), p, width = 9, height = 8.5, dpi = 300)
  cat("  Saved:", basename(paste0(base, ".pdf")), "and .png\n")
  
  invisible(p)
}

# -- RUN ALL ---------------------------------------------------------------------
make_all_figures <- function() {
  if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)
  for (mt in MODEL_TYPES) {
    for (tr in TRANSITIONS) {
      cat("Building:", mt, "-", tr, "\n")
      make_figure(mt, tr)
    }
  }
  cat("\nDone. Figures in:", OUT_DIR, "/\n")
}

if (sys.nframe() == 0) {
  make_all_figures()
}

# To run interactively:
#   CODE_DIR <- "/uufs/chpc.utah.edu/common/home/u1418987/sima/Rt_estimation_project"
#   setwd("/scratch/general/vast/u1418987")
#   source(file.path(CODE_DIR, "paper_figures.R"))
#   make_all_figures()