# ================================================================================
# BEAUTIFUL 2×2 INTERACTIVE GRID - ABM + Cori + WT + RtEstim (Baseline Only)
# ================================================================================

suppressPackageStartupMessages({
  library(plotly)
  library(dplyr)
  library(readr)
  library(htmlwidgets)
})

# ── CONFIG ──────────────────────────────────────────────────────────────────────
PLOTS_DIR       <- "beautiful_grid_plots_combined"
DATE_MIN        <- 12
DATE_MAX        <- 75

R0_VALUES       <- c(1.5, 2.0, 3.0, 5.0)
MODEL_TYPES     <- c("full", "partial")
TRANSITIONS     <- c("susceptible_to_exposed", "exposed_to_infected")

# Data directories
BASELINE_DIR    <- "cori_wallinga_results"
RTESTIM_DIR     <- "rtestim_results"
ABM_DIR         <- "saved_data"

# Color palette
COLORS <- list(
  ABM = list(line = "#2C3E50", ribbon = "rgba(44, 62, 80, 0.12)"),
  Cori = list(line = "#E74C3C", ribbon = "rgba(231, 76, 60, 0.18)"),
  WT = list(line = "#3498DB", ribbon = "rgba(52, 152, 219, 0.18)"),
  RtEstim = list(line = "#27AE60", ribbon = "rgba(39, 174, 96, 0.18)")
)

# ── LOAD FUNCTIONS ──────────────────────────────────────────────────────────────

safe_load <- function(loader_func, ...) {
  tryCatch(loader_func(...), error = function(e) NULL)
}

load_baseline <- function(model_type, R0_val, transition) {
  file_path <- file.path(BASELINE_DIR, paste0(model_type, "_R0_", R0_val, "_", transition, ".csv"))
  if (!file.exists(file_path)) return(NULL)
  read_csv(file_path, show_col_types = FALSE)
}

load_rtestim <- function(model_type, R0_val, transition) {
  # Try .csv first, then .rds
  csv_path <- file.path(RTESTIM_DIR, paste0(model_type, "_R0_", R0_val, "_", transition, ".csv"))
  rds_path <- file.path(RTESTIM_DIR, paste0(model_type, "_R0_", R0_val, "_", transition, ".rds"))
  
  if (file.exists(csv_path)) {
    return(read_csv(csv_path, show_col_types = FALSE))
  } else if (file.exists(rds_path)) {
    return(readRDS(rds_path))
  }
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

# ── PREPARE DATA ────────────────────────────────────────────────────────────────

prepare_panel_data <- function(model_type, R0_val, transition) {
  
  baseline <- safe_load(load_baseline, model_type, R0_val, transition)
  rtestim  <- safe_load(load_rtestim, model_type, R0_val, transition)
  abm      <- safe_load(load_abm, model_type, R0_val)
  
  # Summarize Cori & WT from baseline
  cw_summary <- NULL
  if (!is.null(baseline) && nrow(baseline) > 0) {
    cw_summary <- baseline %>%
      filter(date >= DATE_MIN, date <= DATE_MAX) %>%
      group_by(method, date) %>%
      summarise(
        mean_rt  = mean(median_rt, na.rm = TRUE),
        ci_lower = quantile(median_rt, 0.025, na.rm = TRUE),
        ci_upper = quantile(median_rt, 0.975, na.rm = TRUE),
        .groups  = "drop"
      )
  }
  
  # Summarize RtEstim
  rt_summary <- NULL
  if (!is.null(rtestim) && nrow(rtestim) > 0) {
    rt_summary <- rtestim %>%
      filter(date >= DATE_MIN, date <= DATE_MAX) %>%
      group_by(date) %>%
      summarise(
        mean_rt  = mean(median_rt, na.rm = TRUE),
        ci_lower = quantile(median_rt, 0.025, na.rm = TRUE),
        ci_upper = quantile(median_rt, 0.975, na.rm = TRUE),
        .groups  = "drop"
      ) %>%
      mutate(method = "RtEstim")
  }
  
  # ABM
  abm_summary <- NULL
  if (!is.null(abm)) {
    abm_summary <- abm %>%
      rename(date = source_exposure_date) %>%
      filter(date >= DATE_MIN, date <= DATE_MAX) %>%
      transmute(method = "ABM", date, mean_rt, ci_lower, ci_upper)
  }
  
  list(cw = cw_summary, rtestim = rt_summary, abm = abm_summary, R0 = R0_val)
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
  
  # --- ABM (True Rt) ---
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
  if (!is.null(data$cw)) {
    cori_data <- data$cw %>% filter(method == "Cori")
    if (nrow(cori_data) > 0) {
      p <- p %>%
        add_ribbons(
          data = cori_data,
          x = ~date, ymin = ~ci_lower, ymax = ~ci_upper,
          fillcolor = COLORS$Cori$ribbon, line = list(width = 0),
          showlegend = FALSE, legendgroup = "Cori", hoverinfo = "skip"
        ) %>%
        add_lines(
          data = cori_data,
          x = ~date, y = ~mean_rt,
          line = list(color = COLORS$Cori$line, width = 2.8),
          name = "Cori", legendgroup = "Cori",
          showlegend = is_first_panel,
          hovertemplate = "<b>Cori</b><br>Day: %{x}<br>Rt: %{y:.3f}<extra></extra>"
        )
    }
  }
  
  # --- Wallinga-Teunis ---
  if (!is.null(data$cw)) {
    wt_data <- data$cw %>% filter(method == "Wallinga-Teunis")
    if (nrow(wt_data) > 0) {
      p <- p %>%
        add_ribbons(
          data = wt_data,
          x = ~date, ymin = ~ci_lower, ymax = ~ci_upper,
          fillcolor = COLORS$WT$ribbon, line = list(width = 0),
          showlegend = FALSE, legendgroup = "WT", hoverinfo = "skip"
        ) %>%
        add_lines(
          data = wt_data,
          x = ~date, y = ~mean_rt,
          line = list(color = COLORS$WT$line, width = 2.8),
          name = "Wallinga-Teunis", legendgroup = "WT",
          showlegend = is_first_panel,
          hovertemplate = "<b>Wallinga-Teunis</b><br>Day: %{x}<br>Rt: %{y:.3f}<extra></extra>"
        )
    }
  }
  
  # --- RtEstim ---
  if (!is.null(data$rtestim)) {
    p <- p %>%
      add_ribbons(
        data = data$rtestim,
        x = ~date, ymin = ~ci_lower, ymax = ~ci_upper,
        fillcolor = COLORS$RtEstim$ribbon, line = list(width = 0),
        showlegend = FALSE, legendgroup = "RtEstim", hoverinfo = "skip"
      ) %>%
      add_lines(
        data = data$rtestim,
        x = ~date, y = ~mean_rt,
        line = list(color = COLORS$RtEstim$line, width = 2.8),
        name = "RtEstim", legendgroup = "RtEstim",
        showlegend = is_first_panel,
        hovertemplate = "<b>RtEstim</b><br>Day: %{x}<br>Rt: %{y:.3f}<extra></extra>"
      )
  }
  
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
  
  cat("\n📊 Creating grid:", model_type, "-", transition, "\n")
  
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
          "<b style='font-size:28px; color:#2C3E50'>",
          "Rt Estimation: ", toupper(model_type), " Model</b><br>",
          "<span style='font-size:16px; color:#7F8C8D'>",
          trans_label, " Transition | ABM + Cori + Wallinga-Teunis + RtEstim",
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
        font = list(size = 12),
        tracegroupgap = 15
      ),
      paper_bgcolor = "white",
      margin = list(t = 120, r = 220, b = 60, l = 60)
    ) %>%
    config(
      displayModeBar = TRUE, displaylogo = FALSE,
      modeBarButtonsToRemove = c("select2d", "lasso2d"),
      toImageButtonOptions = list(
        format = "png",
        filename = paste0("rt_combined_", model_type, "_", transition),
        width = 1600, height = 1200, scale = 2
      )
    )
  
  fig
}

# ── MAIN ────────────────────────────────────────────────────────────────────────

create_all_grids <- function() {
  cat("╔══════════════════════════════════════════════════════════════╗\n")
  cat("║     COMBINED Rt GRID: ABM + Cori + WT + RtEstim           ║\n")
  cat("╚══════════════════════════════════════════════════════════════╝\n")
  
  if (!dir.exists(PLOTS_DIR)) dir.create(PLOTS_DIR, recursive = TRUE)
  
  for (model_type in MODEL_TYPES) {
    for (transition in TRANSITIONS) {
      
      cat("\n", rep("═", 60), "\n", sep = "")
      
      tryCatch({
        fig <- create_grid(model_type, transition)
        
        filename <- file.path(
          PLOTS_DIR,
          paste0("combined_grid_", model_type, "_", transition, ".html")
        )
        
        cat("  Saving...\n")
        saveWidget(fig, filename, selfcontained = TRUE)
        cat("✅ SAVED:", filename, "\n")
        
      }, error = function(e) {
        cat("❌ ERROR:", e$message, "\n")
      })
    }
  }
  
  # Index page
  cat("\n📄 Creating index page...\n")
  create_index()
  
  cat("\n╔══════════════════════════════════════════════════════════════╗\n")
  cat("║                    ✨ COMPLETE! ✨                          ║\n")
  cat("╚══════════════════════════════════════════════════════════════╝\n")
  cat("\n📁 Open", file.path(PLOTS_DIR, "index.html"), "in your browser\n\n")
}

# ── INDEX PAGE ──────────────────────────────────────────────────────────────────

create_index <- function() {
  
  files <- list.files(PLOTS_DIR, pattern = "^combined_grid.*\\.html$", full.names = FALSE)
  
  cards_html <- ""
  icons <- c("📈", "📉", "📊", "📋")
  
  for (i in seq_along(files)) {
    file <- files[i]
    parts <- strsplit(gsub("combined_grid_|\\.html", "", file), "_")[[1]]
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
    <title>Combined Rt Estimation Grids</title>
    <style>
        * { margin: 0; padding: 0; box-sizing: border-box; }
        body {
            font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
            background: linear-gradient(135deg, #27AE60 0%, #2C3E50 100%);
            min-height: 100vh; padding: 40px 20px;
        }
        .container {
            max-width: 900px; margin: 0 auto; background: white;
            border-radius: 20px; box-shadow: 0 20px 60px rgba(0,0,0,0.3); overflow: hidden;
        }
        .header {
            background: linear-gradient(135deg, #2C3E50 0%, #27AE60 100%);
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
        .plots { padding: 40px; }
        .plots h2 { font-size: 24px; color: #2C3E50; margin-bottom: 25px; text-align: center; }
        .plot-grid {
            display: grid; grid-template-columns: repeat(auto-fit, minmax(350px, 1fr)); gap: 20px;
        }
        .plot-card {
            background: linear-gradient(135deg, #27AE60 0%, #2C3E50 100%);
            border-radius: 15px; box-shadow: 0 8px 20px rgba(39, 174, 96, 0.3);
            transition: all 0.3s ease; cursor: pointer;
        }
        .plot-card:hover { transform: translateY(-5px); box-shadow: 0 12px 30px rgba(39, 174, 96, 0.5); }
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
        <h1>🔬 Combined Rt Estimation</h1>
        <p>ABM + Cori + Wallinga-Teunis + RtEstim</p>
    </div>
    <div class="info">
        <strong>Methods compared (baseline data only, no delay or misspecification):</strong>
        <div class="methods">
            <span class="method-tag" style="background:#2C3E50">ABM (True Rt)</span>
            <span class="method-tag" style="background:#E74C3C">Cori</span>
            <span class="method-tag" style="background:#3498DB">Wallinga-Teunis</span>
            <span class="method-tag" style="background:#27AE60">RtEstim</span>
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
