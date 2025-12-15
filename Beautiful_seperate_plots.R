# ================================================================================
# INTERACTIVE PLOTS - Click to Show Scenario + ABM
# Click "Delay" → Shows only Delay lines + ABM
# Click "Baseline" → Shows only Baseline lines + ABM
# Click "Misspec" → Shows only Misspec lines + ABM
# ================================================================================

suppressPackageStartupMessages({
  library(plotly)
  library(dplyr)
  library(readr)
  library(htmlwidgets)
})

# ── CONFIG ──────────────────────────────────────────────────────────────────────
PLOTS_DIR       <- "interactive_plots"
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

# ── LOAD FUNCTIONS ──────────────────────────────────────────────────────────────

load_baseline_data <- function(model_type, R0_val, transition) {
  file_name <- paste0(model_type, "_R0_", R0_val, "_", transition, ".csv")
  file_path <- file.path(BASELINE_DIR, file_name)
  if (!file.exists(file_path)) return(NULL)
  read_csv(file_path, show_col_types = FALSE) %>%
    mutate(scenario = "Baseline")
}

load_delay_data <- function(model_type, R0_val, transition) {
  file_name <- paste0("delay_results_", model_type, "_R0_", R0_val, "_", transition, ".csv")
  file_path <- file.path(DELAY_DIR, file_name)
  if (!file.exists(file_path)) return(NULL)
  read_csv(file_path, show_col_types = FALSE) %>%
    filter(delay_scenario == "high_delay") %>%
    mutate(scenario = "Delay")
}

load_misspec_data <- function(model_type, R0_val, transition) {
  file_name <- paste0("misspec_results_", model_type, "_R0_", R0_val, "_", transition, ".csv")
  file_path <- file.path(MISSPEC_DIR, file_name)
  if (!file.exists(file_path)) return(NULL)
  read_csv(file_path, show_col_types = FALSE) %>%
    mutate(scenario = "Misspec")
}

load_abm_data <- function(model_type, R0_val) {
  file_name <- paste0(model_type, "_R0_", R0_val, "_n_1e+05_nsim_100_rt_ci.rds")
  file_path <- file.path(ABM_DIR, file_name)
  if (!file.exists(file_path)) {
    csv_path <- gsub("\\.rds$", ".csv", file_path)
    if (file.exists(csv_path)) return(read_csv(csv_path, show_col_types = FALSE))
    return(NULL)
  }
  readRDS(file_path)
}

# ── PREPARE DATA ────────────────────────────────────────────────────────────────

prepare_interactive_data <- function(model_type, R0_val, transition) {
  
  # Load all
  baseline_data <- load_baseline_data(model_type, R0_val, transition)
  delay_data <- load_delay_data(model_type, R0_val, transition)
  misspec_data <- load_misspec_data(model_type, R0_val, transition)
  abm_data <- load_abm_data(model_type, R0_val)
  
  # Combine and aggregate
  cw_data <- bind_rows(baseline_data, delay_data, misspec_data)
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
  
  # ABM
  abm_summary <- NULL
  if (!is.null(abm_data)) {
    abm_summary <- abm_data %>%
      rename(date = source_exposure_date) %>%
      filter(date >= DATE_MIN, date <= DATE_MAX) %>%
      transmute(
        method = "ABM",
        date,
        mean_rt = mean_rt,
        ci_lower,
        ci_upper,
        scenario = "ABM"
      )
  }
  
  list(cw = cw_summary, abm = abm_summary, R0_value = R0_val)
}

# ── CREATE INTERACTIVE PLOT ─────────────────────────────────────────────────────

create_interactive_plot <- function(model_type, R0_val, transition) {
  
  plot_data <- prepare_interactive_data(model_type, R0_val, transition)
  if (is.null(plot_data)) return(NULL)
  
  all_data <- bind_rows(plot_data$cw, plot_data$abm)
  
  # Initialize figure
  fig <- plot_ly()
  
  # Color mapping
  get_color <- function(method, scenario) {
    if (method == "ABM") return("black")
    if (method == "Cori") {
      if (scenario == "Baseline") return("#D62728")
      if (scenario == "Delay") return("#FF7F0E")
      if (scenario == "Misspec") return("#8C564B")
    }
    if (method == "Wallinga-Teunis") {
      if (scenario == "Baseline") return("#1F77B4")
      if (scenario == "Delay") return("#2CA02C")
      if (scenario == "Misspec") return("#9467BD")
    }
    return("gray")
  }
  
  # Dash style
  get_dash <- function(scenario) {
    if (scenario == "Baseline" || scenario == "ABM") return("solid")
    if (scenario == "Delay") return("dash")
    if (scenario == "Misspec") return("dot")
    return("solid")
  }
  
  # Group name for legend
  get_group_name <- function(scenario) {
    if (scenario == "ABM") return("ABM")
    if (scenario == "Baseline") return("Baseline")
    if (scenario == "Delay") return("Delay")
    if (scenario == "Misspec") return("Misspec")
    return(scenario)
  }
  
  # Add ABM trace (always visible with any scenario)
  if (!is.null(plot_data$abm)) {
    abm_df <- plot_data$abm
    
    # ABM ribbon
    fig <- fig %>%
      add_ribbons(
        data = abm_df,
        x = ~date,
        ymin = ~ci_lower,
        ymax = ~ci_upper,
        fillcolor = "rgba(0,0,0,0.2)",
        line = list(width = 0),
        showlegend = FALSE,
        legendgroup = "ABM",
        hoverinfo = "skip"
      )
    
    # ABM line
    fig <- fig %>%
      add_lines(
        data = abm_df,
        x = ~date,
        y = ~mean_rt,
        line = list(color = "black", width = 3),
        name = "ABM (True Rt)",
        legendgroup = "ABM",
        hovertemplate = "<b>ABM (True)</b><br>Day: %{x}<br>Rt: %{y:.3f}<extra></extra>"
      )
  }
  
  # Add Cori/Wallinga traces grouped by scenario
  for (scen in c("Baseline", "Delay", "Misspec")) {
    scen_data <- plot_data$cw %>% filter(scenario == scen)
    
    if (nrow(scen_data) == 0) next
    
    for (meth in unique(scen_data$method)) {
      method_data <- scen_data %>% filter(method == meth)
      
      color <- get_color(meth, scen)
      dash <- get_dash(scen)
      group <- get_group_name(scen)
      
      scenario_label <- case_when(
        scen == "Baseline" ~ "No Intervention",
        scen == "Delay" ~ "50% Delayed",
        scen == "Misspec" ~ "70% Reporting"
      )
      
      trace_name <- paste0(meth, " - ", scenario_label)
      
      # Ribbon
      fig <- fig %>%
        add_ribbons(
          data = method_data,
          x = ~date,
          ymin = ~ci_lower,
          ymax = ~ci_upper,
          fillcolor = paste0("rgba(", 
                             paste(col2rgb(color), collapse = ","), 
                             ",0.25)"),
          line = list(width = 0),
          showlegend = FALSE,
          legendgroup = group,
          hoverinfo = "skip"
        )
      
      # Line
      fig <- fig %>%
        add_lines(
          data = method_data,
          x = ~date,
          y = ~mean_rt,
          line = list(color = color, width = 2.5, dash = dash),
          name = trace_name,
          legendgroup = group,
          hovertemplate = paste0(
            "<b>", trace_name, "</b><br>",
            "Day: %{x}<br>",
            "Rt: %{y:.3f}<br>",
            "<extra></extra>"
          )
        )
    }
  }
  
  # Add reference line
  fig <- fig %>%
    add_lines(
      x = c(DATE_MIN, DATE_MAX),
      y = c(R0_val, R0_val),
      line = list(color = "gray", width = 2, dash = "dash"),
      name = paste0("R0 = ", R0_val),
      showlegend = FALSE,
      hoverinfo = "skip"
    )
  
  # Layout with grouped legend
  trans_label <- gsub("susceptible_to_exposed", "S → E", transition)
  trans_label <- gsub("exposed_to_infected", "E → I", trans_label)
  
  fig <- fig %>%
    layout(
      title = list(
        text = paste0(
          "<b>Rt Estimation: ", toupper(model_type), " Model | ", 
          trans_label, " | R<sub>0</sub> = ", R0_val, "</b>"
        ),
        font = list(size = 18)
      ),
      xaxis = list(
        title = "<b>Day</b>",
        gridcolor = "#f0f0f0",
        zeroline = FALSE,
        range = c(0, 75)
      ),
      yaxis = list(
        title = "<b>R<sub>t</sub></b>",
        gridcolor = "#f0f0f0",
        zeroline = FALSE
      ),
      plot_bgcolor = "white",
      paper_bgcolor = "white",
      hovermode = "closest",
      legend = list(
        title = list(text = "<b>Click to Show/Hide<br>(ABM always shown)</b>"),
        orientation = "v",
        x = 1.02,
        y = 1,
        bgcolor = "rgba(255,255,255,0.9)",
        bordercolor = "gray",
        borderwidth = 1,
        tracegroupgap = 10
      )
    ) %>%
    config(
      displayModeBar = TRUE,
      displaylogo = FALSE,
      modeBarButtonsToRemove = c("pan2d", "select2d", "lasso2d", "autoScale2d")
    )
  
  fig
}

# ── CREATE ALL INTERACTIVE PLOTS ────────────────────────────────────────────────

create_all_interactive <- function() {
  cat("╔══════════════════════════════════════════════════════════════╗\n")
  cat("║          CREATING INTERACTIVE CLICKABLE PLOTS               ║\n")
  cat("║      Click scenario name to show/hide with ABM              ║\n")
  cat("╚══════════════════════════════════════════════════════════════╝\n")
  
  if (!dir.exists(PLOTS_DIR)) dir.create(PLOTS_DIR, recursive = TRUE)
  
  for (model_type in MODEL_TYPES) {
    for (transition in TRANSITIONS) {
      
      cat("\n📊", toupper(model_type), "-", transition, "\n")
      
      for (R0_val in R0_VALUES) {
        cat("  • R0 =", R0_val, "\n")
        
        tryCatch({
          fig <- create_interactive_plot(model_type, R0_val, transition)
          
          if (!is.null(fig)) {
            filename <- file.path(
              PLOTS_DIR,
              paste0("interactive_", model_type, "_R0_", R0_val, "_", transition, ".html")
            )
            
            saveWidget(fig, filename, selfcontained = TRUE)
            cat("    ✅ Saved:", basename(filename), "\n")
          }
          
        }, error = function(e) {
          cat("    ❌ Error:", e$message, "\n")
        })
      }
    }
  }
  
  # Create index page
  create_index_page()
  
  cat("\n╔══════════════════════════════════════════════════════════════╗\n")
  cat("║                     COMPLETE!                               ║\n")
  cat("╚══════════════════════════════════════════════════════════════╝\n")
  cat("\n📁 Interactive plots saved in:", PLOTS_DIR, "/\n")
  cat("   Open index.html in your browser to navigate all plots\n\n")
  cat("💡 How to use:\n")
  cat("   1. Click legend items to show/hide scenarios\n")
  cat("   2. ABM (True Rt) is always shown\n")
  cat("   3. Double-click a legend group to isolate it\n")
  cat("   4. Hover over lines for exact values\n\n")
}

# ── CREATE INDEX PAGE ───────────────────────────────────────────────────────────

create_index_page <- function() {
  
  html_content <- '
<!DOCTYPE html>
<html>
<head>
    <title>Rt Estimation Interactive Plots</title>
    <style>
        body {
            font-family: Arial, sans-serif;
            max-width: 1200px;
            margin: 0 auto;
            padding: 20px;
            background-color: #f5f5f5;
        }
        h1 {
            color: #333;
            text-align: center;
            border-bottom: 3px solid #4CAF50;
            padding-bottom: 10px;
        }
        h2 {
            color: #555;
            margin-top: 30px;
            border-left: 4px solid #2196F3;
            padding-left: 10px;
        }
        .model-section {
            background: white;
            padding: 20px;
            margin: 20px 0;
            border-radius: 8px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }
        .plot-grid {
            display: grid;
            grid-template-columns: repeat(auto-fill, minmax(250px, 1fr));
            gap: 15px;
            margin-top: 15px;
        }
        .plot-link {
            display: block;
            padding: 15px;
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            text-decoration: none;
            border-radius: 6px;
            text-align: center;
            font-weight: bold;
            transition: transform 0.2s, box-shadow 0.2s;
        }
        .plot-link:hover {
            transform: translateY(-2px);
            box-shadow: 0 4px 8px rgba(0,0,0,0.2);
        }
        .instructions {
            background: #e3f2fd;
            padding: 15px;
            border-radius: 6px;
            border-left: 4px solid #2196F3;
            margin: 20px 0;
        }
        .instructions h3 {
            margin-top: 0;
            color: #1976D2;
        }
        .instructions ul {
            margin: 10px 0;
        }
        .instructions li {
            margin: 5px 0;
        }
    </style>
</head>
<body>
    <h1>🔬 Rt Estimation Interactive Plots</h1>
    
    <div class="instructions">
        <h3>💡 How to Use</h3>
        <ul>
            <li><b>Click legend items</b> to show/hide specific scenarios (Baseline, Delay, Misspec)</li>
            <li><b>ABM (True Rt)</b> is always shown in black for reference</li>
            <li><b>Double-click</b> a legend group to isolate only that scenario + ABM</li>
            <li><b>Hover</b> over lines to see exact Rt values</li>
            <li><b>Zoom</b> by clicking and dragging, reset with double-click</li>
        </ul>
    </div>
'
  
  for (model_type in MODEL_TYPES) {
    html_content <- paste0(html_content, '
    <div class="model-section">
        <h2>📊 ', toupper(model_type), ' Contact Model</h2>
')
    
    for (transition in TRANSITIONS) {
      trans_label <- gsub("_", " ", transition)
      trans_label <- gsub("susceptible to exposed", "S → E", trans_label)
      trans_label <- gsub("exposed to infected", "E → I", trans_label)
      
      html_content <- paste0(html_content, '
        <h3>', trans_label, ' Transition</h3>
        <div class="plot-grid">
')
      
      for (R0_val in R0_VALUES) {
        filename <- paste0("interactive_", model_type, "_R0_", R0_val, "_", transition, ".html")
        html_content <- paste0(html_content, '
            <a href="', filename, '" class="plot-link" target="_blank">
                R₀ = ', R0_val, '
            </a>
')
      }
      
      html_content <- paste0(html_content, '
        </div>
')
    }
    
    html_content <- paste0(html_content, '
    </div>
')
  }
  
  html_content <- paste0(html_content, '
    
    <div style="text-align: center; margin-top: 40px; color: #888; font-size: 12px;">
        <p>Created with R + plotly | Click any plot to open in new tab</p>
    </div>
</body>
</html>
')
  
  writeLines(html_content, file.path(PLOTS_DIR, "index.html"))
  cat("✅ Created index.html navigation page\n")
}

# ── RUN ─────────────────────────────────────────────────────────────────────────
if (sys.nframe() == 0) {
  create_all_interactive()
}