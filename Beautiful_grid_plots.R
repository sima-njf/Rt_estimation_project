# ================================================================================
# INTERACTIVE 2x2 GRID - "Estimating the Reproduction Number"
# All-series confidence bands, warm/cool method-family colors, clean layout.
# Cori = warm (reds/oranges), Wallinga-Teunis = cool (blues/greens), ABM = black.
# ================================================================================

suppressPackageStartupMessages({
  library(plotly)
  library(dplyr)
  library(readr)
  library(htmlwidgets)
})

# -- CONFIG ----------------------------------------------------------------------
PLOTS_DIR    <- "beautiful_grid_plots"
DATE_MIN     <- 12
DATE_MAX     <- 75

R0_VALUES    <- c(1.5, 2.0, 3.0, 5.0)
MODEL_TYPES  <- c("full", "partial")
TRANSITIONS  <- c("susceptible_to_exposed", "exposed_to_infected")

BASELINE_DIR <- "cori_wallinga_results"
DELAY_DIR    <- "complete_delay_results"
MISSPEC_DIR  <- "complete_misspec_results"
ABM_DIR      <- "saved_data"

# Warm/cool palette: Cori = warm, WT = cool, ABM = near-black
COLORS <- list(
  ABM = list(line = "#1A1A1A"),
  Cori = list(
    baseline = "#C0392B",   # deep red        (Real Data)
    delay    = "#E67E22",   # orange          (50% Delayed)
    misspec  = "#F1C40F"    # amber           (70% Reporting)
  ),
  WT = list(
    baseline = "#1F6FB2",   # deep blue       (Real Data)
    delay    = "#16A085",   # teal green      (50% Delayed)
    misspec  = "#8E44AD"    # violet          (70% Reporting)
  )
)

# -- LOAD FUNCTIONS --------------------------------------------------------------
safe_load <- function(loader_func, ...) tryCatch(loader_func(...), error = function(e) NULL)

load_baseline <- function(model_type, R0_val, transition) {
  fp <- file.path(BASELINE_DIR, paste0(model_type, "_R0_", R0_val, "_", transition, ".csv"))
  if (!file.exists(fp)) return(NULL)
  read_csv(fp, show_col_types = FALSE) %>% mutate(scenario = "Baseline")
}

load_delay <- function(model_type, R0_val, transition) {
  fp <- file.path(DELAY_DIR, paste0("delay_results_", model_type, "_R0_", R0_val, "_", transition, ".csv"))
  if (!file.exists(fp)) return(NULL)
  read_csv(fp, show_col_types = FALSE) %>%
    filter(delay_scenario == "high_delay") %>%
    mutate(scenario = "Delay")
}

load_misspec <- function(model_type, R0_val, transition) {
  fp <- file.path(MISSPEC_DIR, paste0("misspec_results_", model_type, "_R0_", R0_val, "_", transition, ".csv"))
  if (!file.exists(fp)) return(NULL)
  read_csv(fp, show_col_types = FALSE) %>% mutate(scenario = "Misspec")
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

# -- PREPARE DATA ----------------------------------------------------------------
prepare_panel_data <- function(model_type, R0_val, transition) {
  baseline <- safe_load(load_baseline, model_type, R0_val, transition)
  delay    <- safe_load(load_delay,    model_type, R0_val, transition)
  misspec  <- safe_load(load_misspec,  model_type, R0_val, transition)
  abm      <- safe_load(load_abm,      model_type, R0_val)
  
  cw_data <- bind_rows(baseline, delay, misspec)
  if (is.null(cw_data) || nrow(cw_data) == 0) return(NULL)
  
  cw_summary <- cw_data %>%
    filter(date >= DATE_MIN, date <= DATE_MAX) %>%
    group_by(method, date, scenario) %>%
    summarise(
      mean_rt  = mean(median_rt, na.rm = TRUE),
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

# -- CREATE PANEL ----------------------------------------------------------------
create_beautiful_panel <- function(data, R0_val, is_first_panel = FALSE) {
  
  if (is.null(data)) {
    return(plot_ly() %>%
             add_annotations(x = 0.5, y = 0.5, text = "No data available",
                             showarrow = FALSE,
                             font = list(size = 14, color = "#95A5A6")) %>%
             layout(xaxis = list(visible = FALSE),
                    yaxis = list(visible = FALSE),
                    paper_bgcolor = "white"))
  }
  
  p <- plot_ly()
  
  get_color <- function(method, scenario) {
    if (method == "ABM") return(COLORS$ABM$line)
    if (method == "Cori") {
      if (scenario == "Baseline") return(COLORS$Cori$baseline)
      if (scenario == "Delay")    return(COLORS$Cori$delay)
      if (scenario == "Misspec")  return(COLORS$Cori$misspec)
    }
    if (method == "Wallinga-Teunis") {
      if (scenario == "Baseline") return(COLORS$WT$baseline)
      if (scenario == "Delay")    return(COLORS$WT$delay)
      if (scenario == "Misspec")  return(COLORS$WT$misspec)
    }
    return("#95A5A6")
  }
  
  get_dash <- function(scenario) {
    if (scenario %in% c("Baseline", "ABM")) return("solid")
    if (scenario == "Delay")   return("dash")
    if (scenario == "Misspec") return("dot")
    return("solid")
  }
  
  add_band_and_line <- function(p, df, color, dash, label, show, width = 3, band_alpha = 0.12) {
    rgb_vals <- col2rgb(color)
    p <- p %>%
      add_ribbons(
        data = df, x = ~date, ymin = ~ci_lower, ymax = ~ci_upper,
        fillcolor = sprintf("rgba(%d,%d,%d,%.2f)", rgb_vals[1], rgb_vals[2], rgb_vals[3], band_alpha),
        line = list(width = 0),
        showlegend = FALSE, legendgroup = label, hoverinfo = "skip"
      )
    p %>%
      add_lines(
        data = df, x = ~date, y = ~mean_rt,
        line = list(color = color, width = width, dash = dash),
        name = label, legendgroup = label, showlegend = show,
        hovertemplate = paste0(
          "<b style='font-size:13px'>", label, "</b><br>",
          "<b>Day:</b> %{x}<br><b>R<sub>t</sub>:</b> %{y:.3f}<extra></extra>"
        )
      )
  }
  
  # ABM first (bold anchor, slightly stronger band)
  if (!is.null(data$abm)) {
    p <- add_band_and_line(p, data$abm, COLORS$ABM$line, "solid",
                           "ABM (True Rt)", is_first_panel, width = 4, band_alpha = 0.14)
  }
  
  # Estimation methods, each with its own CI band
  scenario_labels <- c("Baseline" = "Real Data", "Delay" = "50% Delayed", "Misspec" = "70% Reporting")
  for (scen in c("Baseline", "Delay", "Misspec")) {
    scen_data <- data$cw %>% filter(scenario == scen)
    if (nrow(scen_data) == 0) next
    for (meth in unique(scen_data$method)) {
      md    <- scen_data %>% filter(method == meth)
      color <- get_color(meth, scen)
      dash  <- get_dash(scen)
      label <- paste0(meth, ": ", scenario_labels[scen])
      p <- add_band_and_line(p, md, color, dash, label, is_first_panel, width = 3, band_alpha = 0.12)
    }
  }
  
  p %>%
    layout(
      xaxis = list(
        title = list(text = "<b>Day</b>", font = list(size = 14)),
        gridcolor = "#ECF0F1", gridwidth = 1, showgrid = TRUE,
        range = c(0, 75), zeroline = FALSE,
        showline = TRUE, linewidth = 1.5, linecolor = "#95A5A6",
        ticks = "outside", tickfont = list(size = 11)
      ),
      yaxis = list(
        title = list(text = "<b>R<sub>t</sub></b>", font = list(size = 14)),
        gridcolor = "#ECF0F1", gridwidth = 1, showgrid = TRUE, zeroline = FALSE,
        showline = TRUE, linewidth = 1.5, linecolor = "#95A5A6",
        ticks = "outside", tickfont = list(size = 11)
      ),
      plot_bgcolor = "white", paper_bgcolor = "white",
      hovermode = "closest",
      hoverlabel = list(bgcolor = "white", bordercolor = "#BDC3C7", font = list(size = 12)),
      annotations = list(list(
        x = 0.5, y = 1.12,
        text = paste0("<b style='font-size:17px; color:#2C3E50'>R<sub>0</sub> = ", R0_val, "</b>"),
        xref = "paper", yref = "paper", showarrow = FALSE, xanchor = "center"
      )),
      margin = list(t = 55, b = 45, l = 55, r = 15)
    )
}

# -- CREATE 2x2 GRID -------------------------------------------------------------
create_beautiful_grid <- function(model_type, transition) {
  
  cat("\nCreating grid:", model_type, "-", transition, "\n")
  
  panel_data <- list()
  for (R0_val in R0_VALUES) {
    panel_data[[as.character(R0_val)]] <- prepare_panel_data(model_type, R0_val, transition)
  }
  
  panels <- list(); idx <- 1
  for (R0_val in R0_VALUES) {
    panels[[as.character(R0_val)]] <- create_beautiful_panel(
      panel_data[[as.character(R0_val)]], R0_val, is_first_panel = (idx == 1)
    )
    idx <- idx + 1
  }
  
  fig <- subplot(
    panels[["1.5"]], panels[["2"]],
    panels[["3"]], panels[["5"]],
    nrows = 2, margin = 0.07, shareX = FALSE, shareY = FALSE
  )
  
  trans_label <- gsub("susceptible_to_exposed", "S \u2192 E", transition)
  trans_label <- gsub("exposed_to_infected",    "E \u2192 I", trans_label)
  model_label <- tools::toTitleCase(model_type)
  
  fig %>%
    layout(
      title = list(
        text = paste0(
          "<b style='font-size:26px; color:#2C3E50'>Estimating the Reproduction Number</b><br>",
          "<span style='font-size:15px; color:#7F8C8D'>",
          model_label, " model | ", trans_label, " transition | ",
          "click legend to toggle series</span>"
        ),
        y = 0.98, x = 0.5, xanchor = "center"
      ),
      showlegend = TRUE,
      legend = list(
        title = list(text = "<b style='font-size:13px'>Method &amp; scenario</b>",
                     font = list(color = "#2C3E50")),
        orientation = "v", x = 1.01, y = 0.5,
        bgcolor = "rgba(255,255,255,0.95)", bordercolor = "#BDC3C7",
        borderwidth = 1.5, font = list(size = 12), tracegroupgap = 10
      ),
      paper_bgcolor = "white",
      margin = list(t = 120, r = 270, b = 60, l = 60)
    ) %>%
    config(
      displayModeBar = TRUE, displaylogo = FALSE,
      modeBarButtonsToRemove = c("select2d", "lasso2d"),
      toImageButtonOptions = list(
        format = "png",
        filename = paste0("reproduction_number_", model_type, "_", transition),
        width = 1600, height = 1200, scale = 2
      )
    )
}

# -- MAIN ------------------------------------------------------------------------
create_all_beautiful_grids <- function() {
  cat("========================================\n")
  cat("  Estimating the Reproduction Number - interactive grids\n")
  cat("========================================\n")
  
  if (!dir.exists(PLOTS_DIR)) dir.create(PLOTS_DIR, recursive = TRUE)
  
  for (model_type in MODEL_TYPES) {
    for (transition in TRANSITIONS) {
      cat("\n------------------------------------------------------------\n")
      tryCatch({
        fig <- create_beautiful_grid(model_type, transition)
        filename <- file.path(PLOTS_DIR, paste0("grid_", model_type, "_", transition, ".html"))
        saveWidget(fig, filename, selfcontained = TRUE)
        cat("SAVED:", filename, "\n")
      }, error = function(e) cat("ERROR:", e$message, "\n"))
    }
  }
  
  create_beautiful_index()
  cat("\nDone. Open", file.path(PLOTS_DIR, "index.html"), "\n\n")
}

# -- INDEX PAGE ------------------------------------------------------------------
create_beautiful_index <- function() {
  html <- '<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Estimating the Reproduction Number</title>
    <link href="https://fonts.googleapis.com/css2?family=Inter:wght@400;600;700&display=swap" rel="stylesheet">
    <style>
        * { margin: 0; padding: 0; box-sizing: border-box; }
        body { font-family: "Inter", -apple-system, BlinkMacSystemFont, sans-serif;
               background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
               min-height: 100vh; padding: 40px 20px; }
        .container { max-width: 1000px; margin: 0 auto; background: white;
                     border-radius: 20px; box-shadow: 0 20px 60px rgba(0,0,0,0.3); overflow: hidden; }
        .header { background: linear-gradient(135deg, #2C3E50 0%, #34495E 100%);
                  color: white; padding: 50px 40px; text-align: center; }
        .header h1 { font-size: 34px; font-weight: 700; margin-bottom: 10px; }
        .header p { font-size: 18px; opacity: 0.9; }
        .info-section { padding: 40px; background: #F8F9FA; border-bottom: 1px solid #E9ECEF; }
        .info-section h2 { font-size: 24px; color: #2C3E50; margin-bottom: 20px; }
        .info-grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(250px, 1fr)); gap: 20px; }
        .info-card { background: white; padding: 20px; border-radius: 12px;
                     border-left: 4px solid #667eea; box-shadow: 0 2px 8px rgba(0,0,0,0.08); }
        .info-card h3 { font-size: 16px; color: #2C3E50; margin-bottom: 8px; }
        .info-card p { font-size: 14px; color: #7F8C8D; line-height: 1.6; }
        .plots-section { padding: 40px; }
        .plots-section h2 { font-size: 28px; color: #2C3E50; margin-bottom: 30px; text-align: center; }
        .plot-grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(400px, 1fr)); gap: 25px; }
        .plot-card { background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
                     border-radius: 15px; overflow: hidden; box-shadow: 0 8px 20px rgba(102,126,234,0.3);
                     transition: all 0.3s ease; }
        .plot-card:hover { transform: translateY(-8px); box-shadow: 0 12px 30px rgba(102,126,234,0.5); }
        .plot-card a { display: block; padding: 35px; color: white; text-decoration: none; }
        .plot-card h3 { font-size: 22px; font-weight: 700; margin-bottom: 8px; }
        .plot-card p { font-size: 15px; opacity: 0.9; }
        .footer { text-align: center; padding: 30px; color: #7F8C8D; font-size: 14px; border-top: 1px solid #E9ECEF; }
    </style>
</head>
<body>
    <div class="container">
        <div class="header">
            <h1>Estimating the Reproduction Number</h1>
            <p>Comparing Cori and Wallinga-Teunis estimators against ABM ground truth</p>
        </div>
        <div class="info-section">
            <h2>Simulation Parameters</h2>
            <div class="info-grid">
                <div class="info-card"><h3>Population</h3>
                    <p><strong>Size:</strong> 100,000<br><strong>Initial infected:</strong> 100 (0.1%)<br><strong>Simulations:</strong> 1000 per scenario</p></div>
                <div class="info-card"><h3>Disease Model</h3>
                    <p><strong>Model:</strong> SEIR<br><strong>Incubation:</strong> 4 days<br><strong>Recovery:</strong> 7 days (rate = 1/7)</p></div>
                <div class="info-card"><h3>Contact Structure</h3>
                    <p><strong>Full:</strong> contact rate = 20<br><strong>Partial:</strong> Erdos-Renyi, avg degree = 20<br><strong>Days:</strong> 150</p></div>
                <div class="info-card"><h3>R0 Values</h3>
                    <p>R0 = 1.5, 2.0, 3.0, 5.0<br><strong>Transmission rate:</strong> R0 x recovery / contact</p></div>
                <div class="info-card"><h3>Serial Intervals</h3>
                    <p>Estimated from transmission data; used by both Cori and Wallinga-Teunis</p></div>
                <div class="info-card"><h3>Transitions</h3>
                    <p><strong>S -&gt; E</strong> and <strong>E -&gt; I</strong>, analyzed separately</p></div>
                <div class="info-card"><h3>Real Data</h3>
                    <p>100% reporting, no delays - the validation baseline</p></div>
                <div class="info-card"><h3>50% Delayed</h3>
                    <p>Half of cases delayed by Uniform(1,7) days</p></div>
                <div class="info-card"><h3>70% Reporting</h3>
                    <p>Only 70% of cases reported (binomial thinning)</p></div>
            </div>
        </div>
        <div class="plots-section">
            <h2>Select a Figure</h2>
            <div class="plot-grid">
'
  files <- list.files(PLOTS_DIR, pattern = "^grid_.*\\.html$", full.names = FALSE)
  for (file in files) {
    parts <- strsplit(gsub("grid_|\\.html", "", file), "_")[[1]]
    model <- tools::toTitleCase(parts[1])
    transition <- paste(parts[-1], collapse = " ")
    transition <- gsub("susceptible to exposed", "S -> E", transition)
    transition <- gsub("exposed to infected", "E -> I", transition)
    html <- paste0(html, '
                <div class="plot-card">
                    <a href="', file, '" target="_blank">
                        <h3>', model, ' Model</h3>
                        <p>', transition, ' transition</p>
                    </a>
                </div>
')
  }
  html <- paste0(html, '
            </div>
        </div>
        <div class="footer">
            <p><strong>Created with R + plotly</strong></p>
            <p>Interactive comparison of reproduction-number estimators</p>
        </div>
    </div>
</body>
</html>')
  writeLines(html, file.path(PLOTS_DIR, "index.html"))
  cat("Created index page\n")
}

# -- RUN -------------------------------------------------------------------------
if (sys.nframe() == 0) {
  create_all_beautiful_grids()
}

# To run interactively:
#   CODE_DIR <- "/uufs/chpc.utah.edu/common/home/u1418987/sima/Rt_estimation_project"
#   setwd("/scratch/general/vast/u1418987")
#   source(file.path(CODE_DIR, "Beautiful_grid_plots.R"))
#   create_all_beautiful_grids()