# ================================================================================
# BEAUTIFUL 2×2 INTERACTIVE GRID - PUBLICATION QUALITY
# Polished design with professional styling
# ================================================================================

suppressPackageStartupMessages({
  library(plotly)
  library(dplyr)
  library(readr)
  library(htmlwidgets)
})

# ── CONFIG ──────────────────────────────────────────────────────────────────────
PLOTS_DIR       <- "beautiful_grid_plots"
DATE_MIN        <- 0
DATE_MAX        <- 75

R0_VALUES       <- c(1.5, 2.0, 3.0, 5.0)
MODEL_TYPES     <- c("full", "partial")
TRANSITIONS     <- c("susceptible_to_exposed", "exposed_to_infected")

# Data directories
BASELINE_DIR    <- "cori_wallinga_results"
DELAY_DIR       <- "complete_delay_results"
MISSPEC_DIR     <- "complete_misspec_results"
ABM_DIR         <- "saved_data"

# Beautiful color palette
COLORS <- list(
  ABM = list(
    line = "#2C3E50",      # Dark blue-gray
    ribbon = "#34495E"
  ),
  Cori = list(
    baseline = "#E74C3C",  # Vibrant red
    delay = "#F39C12",     # Warm orange
    misspec = "#C0392B"    # Deep red
  ),
  WT = list(
    baseline = "#3498DB",  # Bright blue
    delay = "#1ABC9C",     # Turquoise
    misspec = "#2980B9"    # Deep blue
  )
)

# ── LOAD FUNCTIONS ──────────────────────────────────────────────────────────────

safe_load <- function(loader_func, ...) {
  tryCatch(loader_func(...), error = function(e) NULL)
}

load_baseline <- function(model_type, R0_val, transition) {
  file_path <- file.path(BASELINE_DIR, paste0(model_type, "_R0_", R0_val, "_", transition, ".csv"))
  if (!file.exists(file_path)) return(NULL)
  read_csv(file_path, show_col_types = FALSE) %>% mutate(scenario = "Baseline")
}

load_delay <- function(model_type, R0_val, transition) {
  file_path <- file.path(DELAY_DIR, paste0("delay_results_", model_type, "_R0_", R0_val, "_", transition, ".csv"))
  if (!file.exists(file_path)) return(NULL)
  read_csv(file_path, show_col_types = FALSE) %>%
    filter(delay_scenario == "high_delay") %>%
    mutate(scenario = "Delay")
}

load_misspec <- function(model_type, R0_val, transition) {
  file_path <- file.path(MISSPEC_DIR, paste0("misspec_results_", model_type, "_R0_", R0_val, "_", transition, ".csv"))
  if (!file.exists(file_path)) return(NULL)
  read_csv(file_path, show_col_types = FALSE) %>% mutate(scenario = "Misspec")
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
  delay <- safe_load(load_delay, model_type, R0_val, transition)
  misspec <- safe_load(load_misspec, model_type, R0_val, transition)
  abm <- safe_load(load_abm, model_type, R0_val)
  
  cw_data <- bind_rows(baseline, delay, misspec)
  if (is.null(cw_data) || nrow(cw_data) == 0) return(NULL)
  
  cw_summary <- cw_data %>%
    filter(date >= DATE_MIN, date <= DATE_MAX) %>%
    group_by(method, date, scenario) %>%
    summarise(
      mean_rt = mean(median_rt, na.rm = TRUE),
      ci_lower = quantile(median_rt, 0.025, na.rm = TRUE),
      ci_upper = quantile(median_rt, 0.975, na.rm = TRUE),
      .groups = "drop"
    )
  
  abm_summary <- NULL
  if (!is.null(abm)) {
    abm_summary <- abm %>%
      rename(date = source_exposure_date) %>%
      filter(date >= DATE_MIN, date <= DATE_MAX) %>%
      transmute(method = "ABM", date, mean_rt, ci_lower, ci_upper, scenario = "ABM")
  }
  
  list(cw = cw_summary, abm = abm_summary, R0 = R0_val)
}

# ── CREATE BEAUTIFUL PANEL ──────────────────────────────────────────────────────

create_beautiful_panel <- function(data, R0_val) {
  
  if (is.null(data)) {
    return(plot_ly() %>%
             add_annotations(
               x = 0.5, y = 0.5,
               text = "No data available",
               showarrow = FALSE,
               font = list(size = 14, color = "#95A5A6")
             ) %>%
             layout(
               xaxis = list(visible = FALSE),
               yaxis = list(visible = FALSE),
               paper_bgcolor = "#FAFAFA"
             ))
  }
  
  p <- plot_ly()
  
  # Helper function for colors
  get_color <- function(method, scenario) {
    if (method == "ABM") return(COLORS$ABM$line)
    if (method == "Cori") {
      if (scenario == "Baseline") return(COLORS$Cori$baseline)
      if (scenario == "Delay") return(COLORS$Cori$delay)
      if (scenario == "Misspec") return(COLORS$Cori$misspec)
    }
    if (method == "Wallinga-Teunis") {
      if (scenario == "Baseline") return(COLORS$WT$baseline)
      if (scenario == "Delay") return(COLORS$WT$delay)
      if (scenario == "Misspec") return(COLORS$WT$misspec)
    }
    return("#95A5A6")
  }
  
  get_dash <- function(scenario) {
    if (scenario %in% c("Baseline", "ABM")) return("solid")
    if (scenario == "Delay") return("dash")
    if (scenario == "Misspec") return("dot")
    return("solid")
  }
  
  # Add ABM (True Rt)
  if (!is.null(data$abm)) {
    abm_df <- data$abm
    
    # ABM ribbon - subtle gray
    p <- p %>%
      add_ribbons(
        data = abm_df,
        x = ~date, ymin = ~ci_lower, ymax = ~ci_upper,
        fillcolor = "rgba(44, 62, 80, 0.12)",
        line = list(width = 0),
        showlegend = FALSE,
        legendgroup = "ABM",
        hoverinfo = "skip"
      )
    
    # ABM line - bold and dark
    p <- p %>%
      add_lines(
        data = abm_df,
        x = ~date, y = ~mean_rt,
        line = list(color = COLORS$ABM$line, width = 3.5),
        name = "ABM (True Rt)",
        legendgroup = "ABM",
        hovertemplate = paste0(
          "<b style='font-size:13px'>ABM (True Rt)</b><br>",
          "<b>Day:</b> %{x}<br>",
          "<b>R<sub>t</sub>:</b> %{y:.3f}",
          "<extra></extra>"
        )
      )
  }
  
  # Add estimation methods
  for (scen in c("Baseline", "Delay", "Misspec")) {
    scen_data <- data$cw %>% filter(scenario == scen)
    if (nrow(scen_data) == 0) next
    
    scenario_labels <- c(
      "Baseline" = "No Intervention",
      "Delay" = "50% Delayed",
      "Misspec" = "70% Reporting"
    )
    
    for (meth in unique(scen_data$method)) {
      meth_data <- scen_data %>% filter(method == meth)
      
      color <- get_color(meth, scen)
      dash <- get_dash(scen)
      label <- paste0(meth, ": ", scenario_labels[scen])
      
      # Ribbon with nice transparency
      rgb_vals <- col2rgb(color)
      p <- p %>%
        add_ribbons(
          data = meth_data,
          x = ~date, ymin = ~ci_lower, ymax = ~ci_upper,
          fillcolor = sprintf("rgba(%d,%d,%d,0.18)", rgb_vals[1], rgb_vals[2], rgb_vals[3]),
          line = list(width = 0),
          showlegend = FALSE,
          legendgroup = scen,
          hoverinfo = "skip"
        )
      
      # Bold line
      p <- p %>%
        add_lines(
          data = meth_data,
          x = ~date, y = ~mean_rt,
          line = list(color = color, width = 2.8, dash = dash),
          name = label,
          legendgroup = scen,
          hovertemplate = paste0(
            "<b style='font-size:13px'>", label, "</b><br>",
            "<b>Day:</b> %{x}<br>",
            "<b>R<sub>t</sub>:</b> %{y:.3f}",
            "<extra></extra>"
          )
        )
    }
  }
  
  # Reference line - elegant dashed
  p <- p %>%
    add_lines(
      x = c(0, 75), y = c(R0_val, R0_val),
      line = list(color = "#7F8C8D", width = 2, dash = "dash"),
      showlegend = FALSE,
      hoverinfo = "skip"
    )
  
  # Beautiful layout
  p <- p %>%
    layout(
      xaxis = list(
        title = list(text = "<b>Day</b>", font = list(size = 13)),
        gridcolor = "#ECF0F1",
        gridwidth = 1,
        range = c(0, 75),
        zeroline = FALSE
      ),
      yaxis = list(
        title = list(text = "<b>R<sub>t</sub></b>", font = list(size = 13)),
        gridcolor = "#ECF0F1",
        gridwidth = 1,
        zeroline = FALSE
      ),
      plot_bgcolor = "#FAFAFA",
      paper_bgcolor = "white",
      hovermode = "closest",
      hoverlabel = list(
        bgcolor = "white",
        bordercolor = "#BDC3C7",
        font = list(size = 12)
      ),
      annotations = list(
        list(
          x = 0.5, y = 1.08,
          text = paste0("<b style='font-size:16px; color:#2C3E50'>R<sub>0</sub> = ", R0_val, "</b>"),
          xref = "paper", yref = "paper",
          showarrow = FALSE,
          xanchor = "center"
        )
      ),
      margin = list(t = 40, b = 40, l = 50, r = 10)
    )
  
  p
}

# ── CREATE BEAUTIFUL 2×2 GRID ───────────────────────────────────────────────────

create_beautiful_grid <- function(model_type, transition) {
  
  cat("\n📊 Creating beautiful 2×2 grid:", model_type, "-", transition, "\n")
  
  # Load data for all panels
  cat("  Loading data...\n")
  panel_data <- list()
  for (R0_val in R0_VALUES) {
    cat("    R0 =", R0_val, "...")
    panel_data[[as.character(R0_val)]] <- prepare_panel_data(model_type, R0_val, transition)
    cat(" ✓\n")
  }
  
  # Create individual panels
  cat("  Creating panels...\n")
  panels <- list()
  for (R0_val in R0_VALUES) {
    panels[[as.character(R0_val)]] <- create_beautiful_panel(
      panel_data[[as.character(R0_val)]], 
      R0_val
    )
  }
  
  # Combine into 2×2 grid
  cat("  Combining into grid...\n")
  fig <- subplot(
    panels[["1.5"]], panels[["2"]],
    panels[["3"]], panels[["5"]],
    nrows = 2,
    margin = 0.08,
    shareX = FALSE,
    shareY = FALSE
  )
  
  # Global styling
  trans_label <- gsub("susceptible_to_exposed", "S → E", transition)
  trans_label <- gsub("exposed_to_infected", "E → I", trans_label)
  
  fig <- fig %>%
    layout(
      title = list(
        text = paste0(
          "<b style='font-size:22px; color:#2C3E50'>",
          "Rt Estimation: ", toupper(model_type), " Model</b><br>",
          "<span style='font-size:14px; color:#7F8C8D'>",
          trans_label, " Transition | ",
          "Click legend to show/hide scenarios across all panels",
          "</span>"
        ),
        y = 0.98,
        x = 0.5,
        xanchor = "center"
      ),
      showlegend = TRUE,
      legend = list(
        title = list(
          text = "<b style='font-size:13px'>Scenarios</b>",
          font = list(color = "#2C3E50")
        ),
        orientation = "v",
        x = 1.01,
        y = 0.5,
        bgcolor = "rgba(255,255,255,0.95)",
        bordercolor = "#BDC3C7",
        borderwidth = 2,
        font = list(size = 11),
        tracegroupgap = 15
      ),
      paper_bgcolor = "white",
      margin = list(t = 120, r = 250, b = 60, l = 60)
    ) %>%
    config(
      displayModeBar = TRUE,
      displaylogo = FALSE,
      modeBarButtonsToRemove = c("select2d", "lasso2d"),
      toImageButtonOptions = list(
        format = "png",
        filename = paste0("rt_", model_type, "_", transition),
        width = 1600,
        height = 1200,
        scale = 2
      )
    )
  
  fig
}

# ── MAIN FUNCTION ───────────────────────────────────────────────────────────────

create_all_beautiful_grids <- function() {
  cat("╔══════════════════════════════════════════════════════════════╗\n")
  cat("║       CREATING BEAUTIFUL 2×2 INTERACTIVE GRIDS             ║\n")
  cat("║              Publication Quality Design                    ║\n")
  cat("╚══════════════════════════════════════════════════════════════╝\n")
  
  if (!dir.exists(PLOTS_DIR)) {
    dir.create(PLOTS_DIR, recursive = TRUE)
    cat("✅ Created directory:", PLOTS_DIR, "\n")
  }
  
  for (model_type in MODEL_TYPES) {
    for (transition in TRANSITIONS) {
      
      cat("\n", rep("═", 60), "\n", sep = "")
      
      tryCatch({
        fig <- create_beautiful_grid(model_type, transition)
        
        filename <- file.path(
          PLOTS_DIR,
          paste0("beautiful_grid_", model_type, "_", transition, ".html")
        )
        
        cat("  Saving...\n")
        saveWidget(fig, filename, selfcontained = TRUE)
        cat("✅ SAVED:", filename, "\n")
        
      }, error = function(e) {
        cat("❌ ERROR:", e$message, "\n")
      })
    }
  }
  
  # Create beautiful index
  cat("\n📄 Creating index page...\n")
  create_beautiful_index()
  
  cat("\n╔══════════════════════════════════════════════════════════════╗\n")
  cat("║                    ✨ COMPLETE! ✨                          ║\n")
  cat("╚══════════════════════════════════════════════════════════════╝\n")
  cat("\n📁 Open", file.path(PLOTS_DIR, "index.html"), "in your browser\n\n")
}

# ── BEAUTIFUL INDEX PAGE ────────────────────────────────────────────────────────

create_beautiful_index <- function() {
  html <- '<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Rt Estimation - Interactive Grids</title>
    <link href="https://fonts.googleapis.com/css2?family=Inter:wght@400;600;700&display=swap" rel="stylesheet">
    <style>
        * { margin: 0; padding: 0; box-sizing: border-box; }
        
        body {
            font-family: "Inter", -apple-system, BlinkMacSystemFont, sans-serif;
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            min-height: 100vh;
            padding: 40px 20px;
        }
        
        .container {
            max-width: 1000px;
            margin: 0 auto;
            background: white;
            border-radius: 20px;
            box-shadow: 0 20px 60px rgba(0,0,0,0.3);
            overflow: hidden;
        }
        
        .header {
            background: linear-gradient(135deg, #2C3E50 0%, #34495E 100%);
            color: white;
            padding: 50px 40px;
            text-align: center;
        }
        
        .header h1 {
            font-size: 36px;
            font-weight: 700;
            margin-bottom: 10px;
        }
        
        .header p {
            font-size: 18px;
            opacity: 0.9;
        }
        
        .info-section {
            padding: 40px;
            background: #F8F9FA;
            border-bottom: 1px solid #E9ECEF;
        }
        
        .info-section h2 {
            font-size: 24px;
            color: #2C3E50;
            margin-bottom: 20px;
            display: flex;
            align-items: center;
        }
        
        .info-section h2::before {
            content: "💡";
            margin-right: 10px;
            font-size: 28px;
        }
        
        .info-grid {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(250px, 1fr));
            gap: 20px;
            margin-top: 20px;
        }
        
        .info-card {
            background: white;
            padding: 20px;
            border-radius: 12px;
            border-left: 4px solid #667eea;
            box-shadow: 0 2px 8px rgba(0,0,0,0.08);
        }
        
        .info-card h3 {
            font-size: 16px;
            color: #2C3E50;
            margin-bottom: 8px;
        }
        
        .info-card p {
            font-size: 14px;
            color: #7F8C8D;
            line-height: 1.6;
        }
        
        .plots-section {
            padding: 40px;
        }
        
        .plots-section h2 {
            font-size: 28px;
            color: #2C3E50;
            margin-bottom: 30px;
            text-align: center;
        }
        
        .plot-grid {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(400px, 1fr));
            gap: 25px;
        }
        
        .plot-card {
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            border-radius: 15px;
            overflow: hidden;
            box-shadow: 0 8px 20px rgba(102, 126, 234, 0.3);
            transition: all 0.3s ease;
            cursor: pointer;
        }
        
        .plot-card:hover {
            transform: translateY(-8px);
            box-shadow: 0 12px 30px rgba(102, 126, 234, 0.5);
        }
        
        .plot-card a {
            display: block;
            padding: 35px;
            color: white;
            text-decoration: none;
        }
        
        .plot-card h3 {
            font-size: 22px;
            font-weight: 700;
            margin-bottom: 8px;
        }
        
        .plot-card p {
            font-size: 15px;
            opacity: 0.9;
        }
        
        .plot-card .icon {
            font-size: 40px;
            margin-bottom: 15px;
            display: block;
        }
        
        .footer {
            text-align: center;
            padding: 30px;
            color: #7F8C8D;
            font-size: 14px;
            border-top: 1px solid #E9ECEF;
        }
    </style>
</head>
<body>
    <div class="container">
        <div class="header">
            <h1>🔬 Rt Estimation Interactive Grids</h1>
            <p>Publication-Quality 2×2 Grid Visualizations</p>
        </div>
        
        <div class="info-section">
            <h2>How to Use</h2>
            <div class="info-grid">
                <div class="info-card">
                    <h3>🖱️ Interactive Legend</h3>
                    <p>Click legend items to show/hide scenarios across all four panels simultaneously</p>
                </div>
                <div class="info-card">
                    <h3>⚫ ABM Reference</h3>
                    <p>True Rt from ABM simulations is always visible in dark gray</p>
                </div>
                <div class="info-card">
                    <h3>🔍 Zoom & Pan</h3>
                    <p>Each panel can be zoomed independently. Double-click to reset</p>
                </div>
                <div class="info-card">
                    <h3>💾 Export</h3>
                    <p>Use the camera icon to download high-resolution PNG images</p>
                </div>
            </div>
        </div>
        
        <div class="plots-section">
            <h2>📊 Select a Plot</h2>
            <div class="plot-grid">
'
  
  files <- list.files(PLOTS_DIR, pattern = "^beautiful_grid.*\\.html$", full.names = FALSE)
  
  icons <- c("📈", "📉", "📊", "📋")
  icon_idx <- 1
  
  for (file in files) {
    parts <- strsplit(gsub("beautiful_grid_|\\.html", "", file), "_")[[1]]
    model <- toupper(parts[1])
    transition <- paste(parts[-1], collapse = " ")
    transition <- gsub("susceptible to exposed", "S → E", transition)
    transition <- gsub("exposed to infected", "E → I", transition)
    transition <- toupper(transition)
    
    html <- paste0(html, '
                <div class="plot-card">
                    <a href="', file, '" target="_blank">
                        <span class="icon">', icons[icon_idx], '</span>
                        <h3>', model, ' Model</h3>
                        <p>', transition, ' Transition</p>
                    </a>
                </div>
')
    icon_idx <- (icon_idx %% length(icons)) + 1
  }
  
  html <- paste0(html, '
            </div>
        </div>
        
        <div class="footer">
            <p><strong>Created with R + plotly</strong></p>
            <p>Interactive visualizations for epidemiological Rt estimation analysis</p>
        </div>
    </div>
</body>
</html>')
  
  writeLines(html, file.path(PLOTS_DIR, "index.html"))
  cat("✅ Created beautiful index page\n")
}

# ── RUN ─────────────────────────────────────────────────────────────────────────
if (sys.nframe() == 0) {
  create_all_beautiful_grids()
}