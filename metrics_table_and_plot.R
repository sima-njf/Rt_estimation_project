# ================================================================================
# ADD-ON to compute_metrics.R : LaTeX summary table + metrics plots.
# Source this AFTER compute_metrics.R (it reuses the same OUT_DIR / metrics CSV),
# or paste these functions into compute_metrics.R before the sys.nframe() block.
#
# Produces:
#   metrics_results/rt_metrics_table.tex   -- Cori vs WT side-by-side, Real Data
#   metrics_results/metrics_MAE_RMSE.png/pdf
#   metrics_results/metrics_coverage.png/pdf
# ================================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
})

# Nice facet labels
.trans_lab <- function(x) {
  x <- gsub("susceptible_to_exposed", "S\u2192E", x)
  x <- gsub("exposed_to_infected",    "E\u2192I", x)
  x
}
.model_lab <- function(x) tools::toTitleCase(x)

# -- 1. LaTeX TABLE (Real Data scenario, Cori vs WT side by side) ----------------
# One row per (model, transition, R0); columns: Cori MAE/RMSE/Cov, WT MAE/RMSE/Cov.
make_latex_table <- function(metrics, scenario_pick = "Real Data",
                             out_tex = file.path(OUT_DIR, "rt_metrics_table.tex")) {
  
  d <- metrics %>%
    filter(scenario == scenario_pick) %>%
    mutate(cov = ifelse(is.na(coverage_band), NA, coverage_band)) %>%
    select(model, transition, R0, method, MAE, RMSE, cov) %>%
    pivot_wider(
      names_from  = method,
      values_from = c(MAE, RMSE, cov),
      names_sep   = "_"
    ) %>%
    arrange(model, transition, R0)
  
  fmt <- function(x) ifelse(is.na(x), "--", formatC(x, format = "f", digits = 3))
  
  lines <- c(
    "% Auto-generated metrics table (Real Data scenario). Requires \\usepackage{booktabs}.",
    "\\begin{table}[ht]",
    "\\centering",
    "\\small",
    paste0("\\caption{Accuracy and empirical band coverage of the Cori and ",
           "Wallinga--Teunis estimators against the agent-based ground truth ",
           "($R_t^{\\text{ABM}}$) over days 12--75, perfect-reporting scenario. ",
           "MAE and RMSE are in units of $R_t$; coverage is the fraction of days ",
           "on which the true value lies within the estimator's 2.5--97.5\\% ",
           "across-replicate band.}"),
    "\\label{tab:metrics}",
    "\\begin{tabular}{ll r rrr rrr}",
    "\\toprule",
    " & & & \\multicolumn{3}{c}{\\textbf{Cori}} & \\multicolumn{3}{c}{\\textbf{Wallinga--Teunis}} \\\\",
    "\\cmidrule(lr){4-6}\\cmidrule(lr){7-9}",
    "Model & Incidence & $R_0$ & MAE & RMSE & Cov. & MAE & RMSE & Cov. \\\\",
    "\\midrule"
  )
  
  for (i in seq_len(nrow(d))) {
    r <- d[i, ]
    lines <- c(lines, paste(
      .model_lab(r$model),
      .trans_lab(r$transition),
      formatC(r$R0, format = "g"),
      fmt(r$`MAE_Cori`),  fmt(r$`RMSE_Cori`),  fmt(r$`cov_Cori`),
      fmt(r$`MAE_Wallinga-Teunis`), fmt(r$`RMSE_Wallinga-Teunis`), fmt(r$`cov_Wallinga-Teunis`),
      sep = " & "
    ))
    lines[length(lines)] <- paste0(lines[length(lines)], " \\\\")
  }
  
  lines <- c(lines, "\\bottomrule", "\\end{tabular}", "\\end{table}", "")
  writeLines(lines, out_tex)
  cat("Wrote LaTeX table:", out_tex, "\n")
  invisible(out_tex)
}

# -- 2. PLOTS --------------------------------------------------------------------
make_metric_plots <- function(metrics) {
  
  pd <- metrics %>%
    mutate(
      transition = .trans_lab(transition),
      model      = .model_lab(model),
      panel      = paste0(model, " | ", transition)
    )
  
  # (a) RMSE vs R0, coloured by method, line style by scenario
  p1 <- ggplot(pd, aes(x = R0, y = RMSE, colour = method,
                       linetype = scenario, group = interaction(method, scenario))) +
    geom_line(linewidth = 0.7) +
    geom_point(size = 1.6) +
    facet_wrap(~ panel, ncol = 2, scales = "free_y") +
    scale_colour_manual(values = c("Cori" = "#C0392B",
                                   "Wallinga-Teunis" = "#1F6FB2"), name = NULL) +
    scale_linetype_manual(values = c("Real Data" = "solid",
                                     "50% Delayed" = "dashed",
                                     "70% Reporting" = "dotted"), name = NULL) +
    labs(title = "Estimator error grows with R0",
         subtitle = "RMSE vs ABM ground truth, days 12-75",
         x = expression(R[0]), y = "RMSE") +
    theme_bw(base_size = 11) +
    theme(legend.position = "bottom", panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold"))
  
  ggsave(file.path(OUT_DIR, "metrics_MAE_RMSE.png"), p1, width = 9, height = 7, dpi = 300)
  ggsave(file.path(OUT_DIR, "metrics_MAE_RMSE.pdf"), p1, width = 9, height = 7, device = cairo_pdf)
  
  # (b) Coverage vs R0
  p2 <- ggplot(pd, aes(x = R0, y = coverage_band, colour = method,
                       linetype = scenario, group = interaction(method, scenario))) +
    geom_hline(yintercept = 0.95, colour = "grey60", linetype = "dashed", linewidth = 0.4) +
    geom_line(linewidth = 0.7) +
    geom_point(size = 1.6) +
    facet_wrap(~ panel, ncol = 2) +
    scale_colour_manual(values = c("Cori" = "#C0392B",
                                   "Wallinga-Teunis" = "#1F6FB2"), name = NULL) +
    scale_linetype_manual(values = c("Real Data" = "solid",
                                     "50% Delayed" = "dashed",
                                     "70% Reporting" = "dotted"), name = NULL) +
    scale_y_continuous(limits = c(0, 1)) +
    labs(title = "Band coverage of the ground truth falls as R0 rises",
         subtitle = "Fraction of days (12-75) ABM truth is inside the 2.5-97.5% band; grey line = 0.95",
         x = expression(R[0]), y = "Coverage") +
    theme_bw(base_size = 11) +
    theme(legend.position = "bottom", panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold"))
  
  ggsave(file.path(OUT_DIR, "metrics_coverage.png"), p2, width = 9, height = 7, dpi = 300)
  ggsave(file.path(OUT_DIR, "metrics_coverage.pdf"), p2, width = 9, height = 7, device = cairo_pdf)
  
  cat("Wrote plots: metrics_MAE_RMSE.(png/pdf), metrics_coverage.(png/pdf)\n")
}

# -- CONVENIENCE: run everything from a saved metrics CSV ------------------------
build_table_and_plots <- function() {
  csv <- file.path(OUT_DIR, "rt_metrics_days12-75.csv")
  if (!file.exists(csv)) stop("Run compute_all_metrics() first to create ", csv)
  metrics <- read_csv(csv, show_col_types = FALSE)
  make_latex_table(metrics)
  make_metric_plots(metrics)
  cat("\nDone. Table + plots in", OUT_DIR, "/\n")
}

if (sys.nframe() == 0) {
  build_table_and_plots()
}

# To run interactively (AFTER compute_all_metrics() has produced the CSV):
#   CODE_DIR <- "/uufs/chpc.utah.edu/common/home/u1418987/sima/Rt_estimation_project"
#   setwd("/scratch/general/vast/u1418987")
#   source(file.path(CODE_DIR, "compute_metrics.R")); compute_all_metrics()
#   source(file.path(CODE_DIR, "metrics_table_and_plot.R")); build_table_and_plots()