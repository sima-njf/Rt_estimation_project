# ================================================================================
# BEAUTIFUL 2×2 INTERACTIVE GRID - REPORTING DELAYS
# ABM (True) + Cori + WT + RtEstim + EpiLPS — Baseline (solid) vs Delays (dotted/dashed)
# ================================================================================

suppressPackageStartupMessages({
  library(plotly)
  library(dplyr)
  library(readr)
  library(htmlwidgets)
})

# ── CONFIG ──────────────────────────────────────────────────────────────────────
PLOTS_DIR       <- "beautiful_grid_plots_delay"
DATE_MIN        <- 12
DATE_MAX        <- 75

R0_VALUES       <- c(1.5, 2.0, 3.0, 5.0)
MODEL_TYPES     <- c("full", "partial")
TRANSITIONS     <- c("susceptible_to_exposed", "exposed_to_infected")

# Data directories
BASELINE_CW_DIR      <- "cori_wallinga_results"
DELAY_CW_DIR         <- "complete_delay_results"
BASELINE_RTESTIM_DIR <- "rtestim_results"
DELAY_RTESTIM_DIR    <- "rtestim_delay_results"
BASELINE_EPILPS_DIR  <- "epilps_results"
DELAY_EPILPS_DIR     <- "epilps_delay_results"
ABM_DIR              <- "saved_data"

# Color palette
COLORS <- list(
  ABM     = list(line = "#2C3E50", ribbon = "rgba(44, 62, 80, 0.12)"),
  Cori    = list(line = "#E74C3C", ribbon = "rgba(231, 76, 60, 0.15)"),
  WT      = list(line = "#3498DB", ribbon = "rgba(52, 152, 219, 0.15)"),
  RtEstim = list(line = "#27AE60", ribbon = "rgba(39, 174, 96, 0.15)"),
  EpiLPS  = list(line = "#8E44AD", ribbon = "rgba(142, 68, 173, 0.15)")
)

# ── LOAD HELPERS ────────────────────────────────────────────────────────────────

safe_load <- function(loader_func, ...) {
  tryCatch(loader_func(...), error = function(e) NULL)
}

load_csv_or_rds <- function(dir, filename_base) {
  csv_path <- file.path(dir, paste0(filename_base, ".csv"))
  rds_path <- file.path(dir, paste0(filename_base, ".rds"))
  if (file.exists(csv_path)) return(read_csv(csv_path, show_col_types = FALSE))
  if (file.exists(rds_path)) return(readRDS(rds_path))
  return(NULL)
}

load_abm <- function(model_type, R0_val) {
  file_path <- file.path(ABM_DIR, paste0(model_type, "_R0_", R0_val, "_n_1e+05_nsim_100_rt_ci.rds"))
  if (!file.exists(file_path)) {
    csv_path <- gsub("\\.rds$", ".csv", file_path)
    if (file.exists(csv_path)) return(read_csv(csv_path, show_col_types = FALSE))
    return(NULL)
  }
  readRDS(file_path)
}

# ── SUMMARIZE HELPERS ───────────────────────────────────────────────────────────

summarize_sims <- function(df, date_min, date_max) {
  if (is.null(df) || nrow(df) == 0) return(NULL)
  df %>%
    filter(date >= date_min, date <= date_max) %>%
    group_by(date) %>%
    summarise(
      mean_rt  = mean(median_rt, na.rm = TRUE),
      ci_lower = quantile(median_rt, 0.025, na.rm = TRUE),
      ci_upper = quantile(median_rt, 0.975, na.rm = TRUE),
      .groups  = "drop"
    )
}

summarize_sims_by_method <- function(df, date_min, date_max) {
  if (is.null(df) || nrow(df) == 0) return(NULL)
  df %>%
    filter(date >= date_min, date <= date_max) %>%
    group_by(method, date) %>%
    summarise(
      mean_rt  = mean(median_rt, na.rm = TRUE),
      ci_lower = quantile(median_rt, 0.025, na.rm = TRUE),
      ci_upper = quantile(median_rt, 0.975, na.rm = TRUE),
      .groups  = "drop"
    )
}

# Filter delay results by scenario
summarize_delay_sims <- function(df, scenario, date_min, date_max) {
  if (is.null(df) || nrow(df) == 0) return(NULL)
  df %>%
    filter(delay_scenario == scenario, date >= date_min, date <= date_max) %>%
    group_by(date) %>%
    summarise(
      mean_rt  = mean(median_rt, na.rm = TRUE),
      ci_lower = quantile(median_rt, 0.025, na.rm = TRUE),
      ci_upper = quantile(median_rt, 0.975, na.rm = TRUE),
      .groups  = "drop"
    )
}

summarize_delay_sims_by_method <- function(df, scenario, date_min, date_max) {
  if (is.null(df) || nrow(df) == 0) return(NULL)
  df %>%
    filter(delay_scenario == scenario, date >= date_min, date <= date_max) %>%
    group_by(method, date) %>%
    summarise(
      mean_rt  = mean(median_rt, na.rm = TRUE),
      ci_lower = quantile(median_rt, 0.025, na.rm = TRUE),
      ci_upper = quantile(median_rt, 0.975, na.rm = TRUE),
      .groups  = "drop"
    )
}

# ── PREPARE DATA ────────────────────────────────────────────────────────────────

prepare_panel_data <- function(model_type, R0_val, transition) {
  
  fname <- paste0(model_type, "_R0_", R0_val, "_", transition)
  delay_cw_fname <- paste0("delay_results_", fname)
  
  # --- Baseline ---
  baseline_cw      <- safe_load(load_csv_or_rds, BASELINE_CW_DIR, fname)
  baseline_rtestim <- safe_load(load_csv_or_rds, BASELINE_RTESTIM_DIR, fname)
  baseline_epilps  <- safe_load(load_csv_or_rds, BASELINE_EPILPS_DIR, fname)
  abm              <- safe_load(load_abm, model_type, R0_val)
  
  # --- Delay results ---
  delay_cw      <- safe_load(load_csv_or_rds, DELAY_CW_DIR, delay_cw_fname)
  delay_rtestim <- safe_load(load_csv_or_rds, DELAY_RTESTIM_DIR, fname)
  delay_epilps  <- safe_load(load_csv_or_rds, DELAY_EPILPS_DIR, fname)
  
  # Summarize baseline
  base_cw_sum      <- summarize_sims_by_method(baseline_cw, DATE_MIN, DATE_MAX)
  base_rtestim_sum <- summarize_sims(baseline_rtestim, DATE_MIN, DATE_MAX)
  base_epilps_sum  <- summarize_sims(baseline_epilps, DATE_MIN, DATE_MAX)
  
  # Summarize medium delay (30%)
  med_cw_sum      <- summarize_delay_sims_by_method(delay_cw, "medium_delay", DATE_MIN, DATE_MAX)
  med_rtestim_sum <- summarize_delay_sims(delay_rtestim, "medium_delay", DATE_MIN, DATE_MAX)
  med_epilps_sum  <- summarize_delay_sims(delay_epilps, "medium_delay", DATE_MIN, DATE_MAX)
  
  # Summarize high delay (50%)
  high_cw_sum      <- summarize_delay_sims_by_method(delay_cw, "high_delay", DATE_MIN, DATE_MAX)
  high_rtestim_sum <- summarize_delay_sims(delay_rtestim, "high_delay", DATE_MIN, DATE_MAX)
  high_epilps_sum  <- summarize_delay_sims(delay_epilps, "high_delay", DATE_MIN, DATE_MAX)
  
  # ABM
  abm_summary <- NULL
  if (!is.null(abm)) {
    abm_summary <- abm %>%
      rename(date = source_exposure_date) %>%
      filter(date >= DATE_MIN, date <= DATE_MAX) %>%
      transmute(method = "ABM", date, mean_rt, ci_lower, ci_upper)
  }
  
  list(
    abm = abm_summary,
    # Baseline
    base_cw      = base_cw_sum,
    base_rtestim = base_rtestim_sum,
    base_epilps  = base_epilps_sum,
    # Medium delay (30%)
    med_cw      = med_cw_sum,
    med_rtestim = med_rtestim_sum,
    med_epilps  = med_epilps_sum,
    # High delay (50%)
    high_cw      = high_cw_sum,
    high_rtestim = high_rtestim_sum,
    high_epilps  = high_epilps_sum,
    R0 = R0_val
  )
}

# ── ADD TRACE TRIPLET (baseline solid + medium dotted + high dashed) ────────────

add_method_traces <- function(p, base_data, med_data, high_data,
                              method_name, color_key, is_first_panel = FALSE) {
  col <- COLORS[[color_key]]
  
  # Baseline (solid)
  if (!is.null(base_data) && nrow(base_data) > 0) {
    p <- p %>%
      add_ribbons(
        data = base_data,
        x = ~date, ymin = ~ci_lower, ymax = ~ci_upper,
        fillcolor = col$ribbon, line = list(width = 0),
        showlegend = FALSE, legendgroup = paste0(method_name, "_base"), hoverinfo = "skip"
      ) %>%
      add_lines(
        data = base_data,
        x = ~date, y = ~mean_rt,
        line = list(color = col$line, width = 2.8),
        name = paste0(method_name, " (baseline)"),
        legendgroup = paste0(method_name, "_base"),
        showlegend = is_first_panel,
        hovertemplate = paste0("<b>", method_name, " (baseline)</b><br>Day: %{x}<br>Rt: %{y:.3f}<extra></extra>")
      )
  }
  
  # Medium delay (dotted)
  if (!is.null(med_data) && nrow(med_data) > 0) {
    p <- p %>%
      add_ribbons(
        data = med_data,
        x = ~date, ymin = ~ci_lower, ymax = ~ci_upper,
        fillcolor = col$ribbon, line = list(width = 0),
        showlegend = FALSE, legendgroup = paste0(method_name, "_med"), hoverinfo = "skip"
      ) %>%
      add_lines(
        data = med_data,
        x = ~date, y = ~mean_rt,
        line = list(color = col$line, width = 2.2, dash = "dot"),
        name = paste0(method_name, " (30% delay)"),
        legendgroup = paste0(method_name, "_med"),
        showlegend = is_first_panel,
        hovertemplate = paste0("<b>", method_name, " (30% delayed)</b><br>Day: %{x}<br>Rt: %{y:.3f}<extra></extra>")
      )
  }
  
  # High delay (dashed)
  if (!is.null(high_data) && nrow(high_data) > 0) {
    p <- p %>%
      add_ribbons(
        data = high_data,
        x = ~date, ymin = ~ci_lower, ymax = ~ci_upper,
        fillcolor = col$ribbon, line = list(width = 0),
        showlegend = FALSE, legendgroup = paste0(method_name, "_high"), hoverinfo = "skip"
      ) %>%
      add_lines(
        data = high_data,
        x = ~date, y = ~mean_rt,
        line = list(color = col$line, width = 2.2, dash = "dash"),
        name = paste0(method_name, " (50% delay)"),
        legendgroup = paste0(method_name, "_high"),
        showlegend = is_first_panel,
        hovertemplate = paste0("<b>", method_name, " (50% delayed)</b><br>Day: %{x}<br>Rt: %{y:.3f}<extra></extra>")
      )
  }
  
  p
}

# ── CREATE PANEL ────────────────────────────────────────────────────────────────

create_panel <- function(data, R0_val, is_first_panel = FALSE) {
  
  if (is.null(data)) {
    return(plot_ly() %>%
             add_annotations(
               x = 0.5, y = 0.5, text = "No data available",
               showarrow = FALSE, font = list(size = 14, color = "#95A5A6")
             ) %>%
             layout(xaxis = list(visible = FALSE), yaxis = list(visible = FALSE)))
  }
  
  p <- plot_ly()
  
  # --- ABM (True Rt) — always solid, no delay version ---
  if (!is.null(data$abm)) {
    p <- p %>%
      add_ribbons(
        data = data$abm,
        x = ~date, ymin = ~ci_lower, ymax = ~ci_upper,
        fillcolor = COLORS$ABM$ribbon, line = list(width = 0),
        showlegend = FALSE, legendgroup = "ABM", hoverinfo = "skip"
      ) %>%
      add_lines(
        data = data$abm,
        x = ~date, y = ~mean_rt,
        line = list(color = COLORS$ABM$line, width = 3.5),
        name = "ABM (True Rt)", legendgroup = "ABM",
        showlegend = is_first_panel,
        hovertemplate = "<b>ABM (True Rt)</b><br>Day: %{x}<br>Rt: %{y:.3f}<extra></extra>"
      )
  }
  
  # --- Cori ---
  base_cori <- if (!is.null(data$base_cw)) data$base_cw %>% filter(method == "Cori") else NULL
  med_cori  <- if (!is.null(data$med_cw))  data$med_cw  %>% filter(method == "Cori") else NULL
  high_cori <- if (!is.null(data$high_cw)) data$high_cw %>% filter(method == "Cori") else NULL
  p <- add_method_traces(p, base_cori, med_cori, high_cori, "Cori", "Cori", is_first_panel)
  
  # --- Wallinga-Teunis ---
  base_wt <- if (!is.null(data$base_cw)) data$base_cw %>% filter(method == "Wallinga-Teunis") else NULL
  med_wt  <- if (!is.null(data$med_cw))  data$med_cw  %>% filter(method == "Wallinga-Teunis") else NULL
  high_wt <- if (!is.null(data$high_cw)) data$high_cw %>% filter(method == "Wallinga-Teunis") else NULL
  p <- add_method_traces(p, base_wt, med_wt, high_wt, "Wallinga-Teunis", "WT", is_first_panel)
  
  # --- RtEstim ---
  p <- add_method_traces(p, data$base_rtestim, data$med_rtestim, data$high_rtestim,
                         "RtEstim", "RtEstim", is_first_panel)
  
  # --- EpiLPS ---
  p <- add_method_traces(p, data$base_epilps, data$med_epilps, data$high_epilps,
                         "EpiLPS", "EpiLPS", is_first_panel)
  
  # --- Rt = 1 reference line ---
  p <- p %>%
    add_segments(
      x = DATE_MIN, xend = DATE_MAX, y = 1, yend = 1,
      line = list(color = "#95A5A6", width = 1.5, dash = "dash"),
      showlegend = FALSE, hoverinfo = "skip"
    )
  
  # --- Layout ---
  p <- p %>%
    layout(
      xaxis = list(
        title = list(text = "<b>Day</b>", font = list(size = 13)),
        gridcolor = "#BDC3C7", gridwidth = 1.5, showgrid = TRUE,
        range = c(DATE_MIN, DATE_MAX),
        showline = TRUE, linewidth = 2, linecolor = "#2C3E50", mirror = TRUE
      ),
      yaxis = list(
        title = list(text = "<b>R<sub>t</sub></b>", font = list(size = 13)),
        gridcolor = "#BDC3C7", gridwidth = 1.5, showgrid = TRUE,
        showline = TRUE, linewidth = 2, linecolor = "#2C3E50", mirror = TRUE
      ),
      plot_bgcolor = "#FAFAFA",
      paper_bgcolor = "white",
      hovermode = "closest",
      hoverlabel = list(bgcolor = "white", bordercolor = "#BDC3C7", font = list(size = 12)),
      annotations = list(
        list(
          x = 0.5, y = 1.15,
          text = paste0("<b style='font-size:16px; color:#2C3E50'>R<sub>0</sub> = ", R0_val, "</b>"),
          xref = "paper", yref = "paper",
          showarrow = FALSE, xanchor = "center"
        )
      ),
      margin = list(t = 60, b = 40, l = 50, r = 10)
    )
  
  p
}

# ── CREATE 2×2 GRID ─────────────────────────────────────────────────────────────

create_grid <- function(model_type, transition) {
  
  cat("\n📊 Creating delay grid:", model_type, "-", transition, "\n")
  
  panel_data <- list()
  for (R0_val in R0_VALUES) {
    cat("    R0 =", R0_val, "...")
    panel_data[[as.character(R0_val)]] <- prepare_panel_data(model_type, R0_val, transition)
    cat(" ✓\n")
  }
  
  panels <- list()
  for (i in seq_along(R0_VALUES)) {
    R0_val <- R0_VALUES[i]
    panels[[as.character(R0_val)]] <- create_panel(
      panel_data[[as.character(R0_val)]],
      R0_val,
      is_first_panel = (i == 1)
    )
  }
  
  fig <- subplot(
    panels[["1.5"]], panels[["2"]],
    panels[["3"]], panels[["5"]],
    nrows = 2, margin = 0.05,
    shareX = FALSE, shareY = FALSE
  )
  
  trans_label <- gsub("susceptible_to_exposed", "S → E", transition)
  trans_label <- gsub("exposed_to_infected", "E → I", trans_label)
  
  fig <- fig %>%
    layout(
      title = list(
        text = paste0(
          "<b style='font-size:26px; color:#2C3E50'>",
          "Rt Estimation — Reporting Delays: ", toupper(model_type), " Model</b><br>",
          "<span style='font-size:15px; color:#7F8C8D'>",
          trans_label, " Transition | Solid = Baseline, Dotted = 30% Delayed, Dashed = 50% Delayed",
          "</span>"
        ),
        y = 0.98, x = 0.5, xanchor = "center"
      ),
      showlegend = TRUE,
      legend = list(
        title = list(
          text = "<b style='font-size:13px'>Methods</b>",
          font = list(color = "#2C3E50")
        ),
        orientation = "v",
        x = 1.01, y = 0.5,
        bgcolor = "rgba(255,255,255,0.95)",
        bordercolor = "#BDC3C7", borderwidth = 2,
        font = list(size = 11),
        tracegroupgap = 5
      ),
      paper_bgcolor = "white",
      margin = list(t = 120, r = 280, b = 60, l = 60)
    ) %>%
    config(
      displayModeBar = TRUE, displaylogo = FALSE,
      modeBarButtonsToRemove = c("select2d", "lasso2d"),
      toImageButtonOptions = list(
        format = "png",
        filename = paste0("rt_delay_", model_type, "_", transition),
        width = 1600, height = 1200, scale = 2
      )
    )
  
  fig
}

# ── MAIN ────────────────────────────────────────────────────────────────────────

create_all_grids <- function() {
  cat("╔══════════════════════════════════════════════════════════════╗\n")
  cat("║  DELAY Rt GRID: Baseline vs 30%/50% Delays (All Methods)  ║\n")
  cat("╚══════════════════════════════════════════════════════════════╝\n")
  
  if (!dir.exists(PLOTS_DIR)) dir.create(PLOTS_DIR, recursive = TRUE)
  
  for (model_type in MODEL_TYPES) {
    for (transition in TRANSITIONS) {
      
      cat("\n", rep("═", 60), "\n", sep = "")
      
      tryCatch({
        fig <- create_grid(model_type, transition)
        
        filename <- file.path(
          PLOTS_DIR,
          paste0("delay_grid_", model_type, "_", transition, ".html")
        )
        
        cat("  Saving...\n")
        saveWidget(fig, filename, selfcontained = TRUE)
        cat("✅ SAVED:", filename, "\n")
        
      }, error = function(e) {
        cat("❌ ERROR:", e$message, "\n")
      })
    }
  }
  
  cat("\n📄 Creating index page...\n")
  create_index()
  
  cat("\n╔══════════════════════════════════════════════════════════════╗\n")
  cat("║                    ✨ COMPLETE! ✨                          ║\n")
  cat("╚══════════════════════════════════════════════════════════════╝\n")
  cat("\n📁 Open", file.path(PLOTS_DIR, "index.html"), "in your browser\n\n")
}

# ── INDEX PAGE ──────────────────────────────────────────────────────────────────

create_index <- function() {
  
  files <- list.files(PLOTS_DIR, pattern = "^delay_grid.*\\.html$", full.names = FALSE)
  
  cards_html <- ""
  icons <- c("📈", "📉", "📊", "📋")
  
  for (i in seq_along(files)) {
    file <- files[i]
    parts <- strsplit(gsub("delay_grid_|\\.html", "", file), "_")[[1]]
    model <- toupper(parts[1])
    transition <- paste(parts[-1], collapse = " ")
    transition <- gsub("susceptible to exposed", "S → E", transition)
    transition <- gsub("exposed to infected", "E → I", transition)
    
    cards_html <- paste0(cards_html, '
                <div class="plot-card">
                    <a href="', file, '" target="_blank">
                        <span class="icon">', icons[((i - 1) %% length(icons)) + 1], '</span>
                        <h3>', model, ' Model</h3>
                        <p>', transition, ' Transition</p>
                    </a>
                </div>')
  }
  
  html <- paste0('<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <title>Reporting Delays Rt Estimation Grids</title>
    <style>
        * { margin: 0; padding: 0; box-sizing: border-box; }
        body {
            font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
            background: linear-gradient(135deg, #E67E22 0%, #2C3E50 100%);
            min-height: 100vh; padding: 40px 20px;
        }
        .container {
            max-width: 900px; margin: 0 auto; background: white;
            border-radius: 20px; box-shadow: 0 20px 60px rgba(0,0,0,0.3); overflow: hidden;
        }
        .header {
            background: linear-gradient(135deg, #2C3E50 0%, #E67E22 100%);
            color: white; padding: 50px 40px; text-align: center;
        }
        .header h1 { font-size: 32px; margin-bottom: 10px; }
        .header p { font-size: 16px; opacity: 0.9; }
        .info {
            padding: 30px 40px; background: #F8F9FA;
            border-bottom: 1px solid #E9ECEF; font-size: 15px; color: #555;
        }
        .info strong { color: #2C3E50; }
        .methods { display: flex; gap: 15px; margin-top: 15px; flex-wrap: wrap; }
        .method-tag {
            padding: 8px 16px; border-radius: 20px; font-weight: 600;
            font-size: 13px; color: white;
        }
        .legend-note {
            margin-top: 15px; padding: 12px 16px; background: #FFF3CD;
            border-radius: 8px; font-size: 14px; color: #856404;
        }
        .plots { padding: 40px; }
        .plots h2 { font-size: 24px; color: #2C3E50; margin-bottom: 25px; text-align: center; }
        .plot-grid {
            display: grid; grid-template-columns: repeat(auto-fit, minmax(350px, 1fr)); gap: 20px;
        }
        .plot-card {
            background: linear-gradient(135deg, #E67E22 0%, #2C3E50 100%);
            border-radius: 15px; box-shadow: 0 8px 20px rgba(230, 126, 34, 0.3);
            transition: all 0.3s ease; cursor: pointer;
        }
        .plot-card:hover { transform: translateY(-5px); box-shadow: 0 12px 30px rgba(230, 126, 34, 0.5); }
        .plot-card a { display: block; padding: 30px; color: white; text-decoration: none; }
        .plot-card h3 { font-size: 20px; margin-bottom: 5px; }
        .plot-card p { font-size: 14px; opacity: 0.9; }
        .plot-card .icon { font-size: 36px; margin-bottom: 12px; display: block; }
        .footer { text-align: center; padding: 25px; color: #7F8C8D; font-size: 13px; border-top: 1px solid #E9ECEF; }
    </style>
</head>
<body>
<div class="container">
    <div class="header">
        <h1>⏱️ Reporting Delays Analysis</h1>
        <p>Baseline vs 30% / 50% Delayed — All 5 Methods</p>
    </div>
    <div class="info">
        <strong>Methods compared:</strong>
        <div class="methods">
            <span class="method-tag" style="background:#2C3E50">ABM (True Rt)</span>
            <span class="method-tag" style="background:#E74C3C">Cori</span>
            <span class="method-tag" style="background:#3498DB">Wallinga-Teunis</span>
            <span class="method-tag" style="background:#27AE60">RtEstim</span>
            <span class="method-tag" style="background:#8E44AD">EpiLPS</span>
        </div>
        <div class="legend-note">
            <b>Solid lines</b> = Baseline (no delays) &nbsp;|&nbsp;
            <b>Dotted lines</b> = 30% delayed (Uniform 1–7 days) &nbsp;|&nbsp;
            <b>Dashed lines</b> = 50% delayed (Uniform 1–7 days)
        </div>
    </div>
    <div class="plots">
        <h2>Select a Plot</h2>
        <div class="plot-grid">', cards_html, '
        </div>
    </div>
    <div class="footer">
        <p>Interactive 2×2 grids | R₀ = 1.5, 2.0, 3.0, 5.0</p>
    </div>
</div>
</body>
</html>')
  
  writeLines(html, file.path(PLOTS_DIR, "index.html"))
  cat("✅ Created index page\n")
}

# ── RUN ─────────────────────────────────────────────────────────────────────────
if (sys.nframe() == 0) {
  create_all_grids()
}