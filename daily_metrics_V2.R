# ================================================================================
# DAILY per-simulation metrics vs ABM ground truth  (VERSION 2)
# All outputs prefixed "v2_" so they never collide with earlier results.
#
# True Rt on day t = ABM mean_rt (the black line).
#
# For each (model, transition, R0, method, scenario) and each DAY t, using all
# 1000 simulations:
#
#   MAE_t   = mean over sims of |median_rt_{i,t} - abm_mean_t|
#   RMSE_t  = sqrt( mean over sims of (median_rt_{i,t} - abm_mean_t)^2 )
#
#   COVERAGE -- computed BOTH ways:
#   (1) coverage_persim_t = fraction of the 1000 sims whose OWN CI
#                           [ci_lower_{i,t}, ci_upper_{i,t}] contains abm_mean_t.
#                           (روز اول 80%، روز دوم 20%، ... -- your definition)
#   (2) coverage_band_t   = 1 if abm_mean_t lies within the 2.5-97.5% band of the
#                           sim point estimates on day t, else 0.
#
# Plots: one 2x2 grid per (model, transition, metric), faceted by R0,
#        x-axis = Day, 6 coloured lines (method x scenario), main-figure palette.
# ================================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(ggplot2)
})

DATE_MIN <- 12; DATE_MAX <- 75
R0_VALUES   <- c(1.5, 2.0, 3.0, 5.0)
MODEL_TYPES <- c("full", "partial")
TRANSITIONS <- c("susceptible_to_exposed", "exposed_to_infected")

BASELINE_DIR <- "cori_wallinga_results"
DELAY_DIR    <- "complete_delay_results"
MISSPEC_DIR  <- "complete_misspec_results"
ABM_DIR      <- "saved_data"
OUT_DIR      <- "metrics_results"

SERIES_LEVELS <- c(
  "Cori: Real Data", "Wallinga-Teunis: Real Data",
  "Cori: 50% Delayed", "Wallinga-Teunis: 50% Delayed",
  "Cori: 70% Reporting", "Wallinga-Teunis: 70% Reporting"
)
SERIES_COLORS <- c(
  "Cori: Real Data"                = "#C0392B",
  "Cori: 50% Delayed"              = "#E67E22",
  "Cori: 70% Reporting"            = "#F1C40F",
  "Wallinga-Teunis: Real Data"     = "#1F6FB2",
  "Wallinga-Teunis: 50% Delayed"   = "#16A085",
  "Wallinga-Teunis: 70% Reporting" = "#8E44AD"
)
SERIES_LTY <- c(
  "Cori: Real Data" = "solid",  "Wallinga-Teunis: Real Data" = "solid",
  "Cori: 50% Delayed" = "dashed","Wallinga-Teunis: 50% Delayed" = "dashed",
  "Cori: 70% Reporting" = "dotted","Wallinga-Teunis: 70% Reporting" = "dotted"
)

safe_load <- function(f, ...) tryCatch(f(...), error = function(e) NULL)

load_baseline <- function(model_type, R0_val, transition) {
  fp <- file.path(BASELINE_DIR, paste0(model_type, "_R0_", R0_val, "_", transition, ".csv"))
  if (!file.exists(fp)) return(NULL)
  df <- read_csv(fp, show_col_types = FALSE)
  if (all(c("lower","upper") %in% names(df)) && !all(c("ci_lower","ci_upper") %in% names(df)))
    df <- df %>% rename(ci_lower = lower, ci_upper = upper)
  df %>% mutate(scenario = "Real Data")
}
load_delay <- function(model_type, R0_val, transition) {
  fp <- file.path(DELAY_DIR, paste0("delay_results_", model_type, "_R0_", R0_val, "_", transition, ".csv"))
  if (!file.exists(fp)) return(NULL)
  read_csv(fp, show_col_types = FALSE) %>% filter(delay_scenario == "high_delay") %>%
    mutate(scenario = "50% Delayed")
}
load_misspec <- function(model_type, R0_val, transition) {
  fp <- file.path(MISSPEC_DIR, paste0("misspec_results_", model_type, "_R0_", R0_val, "_", transition, ".csv"))
  if (!file.exists(fp)) return(NULL)
  read_csv(fp, show_col_types = FALSE) %>% mutate(scenario = "70% Reporting")
}
load_abm <- function(model_type, R0_val) {
  fp <- file.path(ABM_DIR, paste0(model_type, "_R0_", R0_val, "_n_1e+05_nsim_1000_rt_ci.rds"))
  if (!file.exists(fp)) { csvp <- gsub("\\.rds$",".csv",fp)
  if (file.exists(csvp)) return(read_csv(csvp, show_col_types = FALSE)); return(NULL) }
  readRDS(fp)
}

abm_daily_mean <- function(model_type, R0_val) {
  abm <- safe_load(load_abm, model_type, R0_val)
  if (is.null(abm)) return(NULL)
  abm %>% rename(date = source_exposure_date) %>%
    filter(date >= DATE_MIN, date <= DATE_MAX) %>%
    group_by(date) %>% summarise(abm_mean = mean(mean_rt, na.rm = TRUE), .groups = "drop")
}

# One-time check: do the per-scenario files carry per-row CI columns?
inspect_ci <- function() {
  cat("=== CI column check ===\n")
  for (nm in c("baseline","delay","misspec")) {
    df <- switch(nm,
                 baseline = safe_load(load_baseline, "full", 1.5, "exposed_to_infected"),
                 delay    = safe_load(load_delay,    "full", 1.5, "exposed_to_infected"),
                 misspec  = safe_load(load_misspec,  "full", 1.5, "exposed_to_infected"))
    if (is.null(df)) { cat(sprintf("  %-9s : file not found\n", nm)); next }
    has <- all(c("ci_lower","ci_upper") %in% names(df))
    cat(sprintf("  %-9s : per-row CI present = %s   (cols: %s)\n",
                nm, has, paste(names(df), collapse=", ")))
  }
  cat("\n")
}

# Daily metrics for one (model, transition). Per-sim MAE/RMSE + both coverages.
daily_metrics_for_combo <- function(model_type, transition) {
  out <- list()
  for (R0_val in R0_VALUES) {
    truth <- abm_daily_mean(model_type, R0_val)
    if (is.null(truth)) next
    for (df in list(load_baseline(model_type,R0_val,transition),
                    load_delay(model_type,R0_val,transition),
                    load_misspec(model_type,R0_val,transition))) {
      if (is.null(df) || nrow(df) == 0) next
      scen <- df$scenario[1]
      has_ci <- all(c("ci_lower","ci_upper") %in% names(df))
      
      dd <- df %>% filter(date >= DATE_MIN, date <= DATE_MAX) %>%
        left_join(truth, by = "date") %>% filter(!is.na(abm_mean))
      if (nrow(dd) == 0) next
      
      daily <- dd %>% group_by(method, date) %>%
        summarise(
          n_sims = n(),
          # per-sim error metrics vs ABM mean
          MAE  = mean(abs(median_rt - abm_mean), na.rm = TRUE),
          RMSE = sqrt(mean((median_rt - abm_mean)^2, na.rm = TRUE)),
          # coverage (1): fraction of sims whose OWN CI contains ABM truth
          coverage_persim = if (has_ci)
            mean(abm_mean >= ci_lower & abm_mean <= ci_upper, na.rm = TRUE)
          else NA_real_,
          # coverage (2): ABM inside the 2.5-97.5% band of sim point estimates
          band_low  = quantile(median_rt, 0.025, na.rm = TRUE),
          band_high = quantile(median_rt, 0.975, na.rm = TRUE),
          abm_mean  = first(abm_mean),
          .groups   = "drop"
        ) %>%
        mutate(
          coverage_band = as.integer(abm_mean >= band_low & abm_mean <= band_high),
          R0 = R0_val, scenario = scen,
          series = paste0(method, ": ", scenario)
        )
      out[[paste(R0_val, scen)]] <- daily
    }
  }
  res <- bind_rows(out)
  if (nrow(res)) { res$model <- model_type; res$transition <- transition }
  res
}

# One 2x2 grid figure (facet by R0) for a (model, transition, metric)
make_grid_figure <- function(dat, metric, ylab, model_type, transition, add_095 = FALSE) {
  dat <- dat %>% filter(!is.na(.data[[metric]]))
  if (nrow(dat) == 0) { cat("  (no data for", metric, model_type, transition, ")\n"); return() }
  present <- SERIES_LEVELS[SERIES_LEVELS %in% unique(dat$series)]
  dat$series <- factor(dat$series, levels = present)
  
  trans_lab <- gsub("susceptible_to_exposed","S \u2192 E",transition)
  trans_lab <- gsub("exposed_to_infected","E \u2192 I",trans_lab)
  facet_lab <- as_labeller(function(x) paste0("R[0] == ", x), label_parsed)
  
  p <- ggplot(dat, aes(x = date, y = .data[[metric]],
                       color = series, linetype = series)) +
    { if (add_095) geom_hline(yintercept = 0.95, color="grey60",
                              linetype="dashed", linewidth=0.4) } +
    geom_line(linewidth = 0.7) +
    facet_wrap(~ R0, ncol = 2, scales = "free_y", labeller = facet_lab) +
    scale_color_manual(values = SERIES_COLORS[present], name = NULL) +
    scale_linetype_manual(values = SERIES_LTY[present], name = NULL) +
    scale_x_continuous(limits = c(DATE_MIN, DATE_MAX), breaks = seq(0,75,25)) +
    labs(title = paste0("Daily ", ylab, " vs ABM ground truth"),
         subtitle = paste0(tools::toTitleCase(model_type), " model  |  ",
                           trans_lab, " transition"),
         x = "Day", y = ylab) +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(face="bold", size=13),
          plot.subtitle = element_text(color="grey30"),
          strip.background = element_rect(fill="grey92", color=NA),
          strip.text = element_text(face="bold"),
          panel.grid.minor = element_blank(),
          legend.position = "bottom", legend.text = element_text(size=9)) +
    guides(color = guide_legend(nrow = 3, byrow = TRUE),
           linetype = guide_legend(nrow = 3, byrow = TRUE))
  
  base <- file.path(OUT_DIR, paste0("v2_daily_", metric, "_", model_type, "_", transition))
  ggsave(paste0(base,".png"), p, width = 9, height = 8.5, dpi = 300)
  ggsave(paste0(base,".pdf"), p, width = 9, height = 8.5, device = cairo_pdf)
  cat("  Saved:", basename(base), "(.png/.pdf)\n")
}

build_daily_v2 <- function() {
  if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)
  inspect_ci()
  
  all_rows <- list()
  for (mt in MODEL_TYPES) for (tr in TRANSITIONS) {
    cat("Computing:", mt, "-", tr, "\n")
    d <- daily_metrics_for_combo(mt, tr)
    if (nrow(d) == 0) next
    all_rows[[paste(mt,tr)]] <- d
    make_grid_figure(d, "MAE",             "MAE",                     mt, tr)
    make_grid_figure(d, "RMSE",            "RMSE",                    mt, tr)
    make_grid_figure(d, "coverage_persim", "Coverage (per-sim CI)",   mt, tr, add_095 = TRUE)
    make_grid_figure(d, "coverage_band",   "Coverage (band, 0/1)",    mt, tr)
  }
  daily <- bind_rows(all_rows)
  out_csv <- file.path(OUT_DIR, "v2_rt_daily_metrics_days12-75.csv")
  write_csv(daily, out_csv)
  cat("\nSaved CSV:", out_csv, "\n")
  cat("All v2 figures in", OUT_DIR, "/\n")
  invisible(daily)
}

if (sys.nframe() == 0) build_daily_v2()

# Interactive:
#   CODE_DIR <- "/uufs/chpc.utah.edu/common/home/u1418987/sima/Rt_estimation_project"
#   setwd("/scratch/general/vast/u1418987")
#   source(file.path(CODE_DIR, "daily_metrics_v2.R")); build_daily_v2()