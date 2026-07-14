# =============================================================================
# Rt Estimation: Renewal Equation + Compartmental vs Cori + Wallinga-Teunis
# All R0 values (1.5, 2, 3, 5) × model types (full, partial)
# × incidence signals (susceptible_to_exposed, exposed_to_infected)
# =============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)

# =============================================================================
# CONFIGURATION
# =============================================================================

sim_data_dir <- "~/Rt_project/saved_data/"
cori_wt_dir  <- "~/Rt_project/cori_wallinga_results/"
out_dir      <- "~/Rt_project/figures/"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

R0_values  <- c(1.5, 2, 3, 5)
mod_types  <- c("full", "partial")
signals    <- c("susceptible_to_exposed", "exposed_to_infected")

# File naming: R0 = 1.5 keeps the dot; 2, 3, 5 use integer (no decimal)
r0_label <- function(r0) ifelse(r0 == as.integer(r0), as.integer(r0), r0)

# =============================================================================
# HELPER: BUILD GI DISTRIBUTION (gentime > 0 only, pooled across all sims)
# =============================================================================

build_gi_dist <- function(gen_data) {
  gen_clean <- gen_data[gen_data$gentime > 0, ]
  cat("    sims:", n_distinct(gen_clean$sim_num),
      "| rows used:", nrow(gen_clean),
      "| dropped (gentime=0):", nrow(gen_data) - nrow(gen_clean), "\n")
  gi_vals  <- gen_clean$gentime
  max_gi   <- max(gi_vals)
  gi_dist  <- tabulate(gi_vals, nbins = max_gi) / length(gi_vals)
  cat("    mean GI:", round(sum(seq_along(gi_dist) * gi_dist), 3), "\n")
  gi_dist
}

# =============================================================================
# HELPER: RENEWAL EQUATION DIRECT FIT
# =============================================================================

renewal_rt_sim <- function(incidence, gi_dist, min_denom = 1) {
  n      <- length(incidence)
  max_gi <- length(gi_dist)
  Rt     <- rep(NA_real_, n)
  for (t in 2:n) {
    lb     <- min(t - 1, max_gi)
    denom  <- sum(incidence[(t - 1):(t - lb)] * gi_dist[1:lb])
    if (denom >= min_denom) Rt[t] <- incidence[t] / denom
  }
  Rt
}

# =============================================================================
# HELPER: SUMMARISE Rt ACROSS SIMS → one row per date
# =============================================================================

summarise_rt <- function(df, rt_col, sim_col = "sim_num") {
  df %>%
    filter(!is.na(.data[[rt_col]]), is.finite(.data[[rt_col]])) %>%
    group_by(date) %>%
    summarise(
      n_sims    = n_distinct(.data[[sim_col]]),
      mean_Rt   = mean(.data[[rt_col]]),
      median_Rt = median(.data[[rt_col]]),
      lo95      = quantile(.data[[rt_col]], 0.025),
      hi95      = quantile(.data[[rt_col]], 0.975),
      lo50      = quantile(.data[[rt_col]], 0.25),
      hi50      = quantile(.data[[rt_col]], 0.75),
      .groups   = "drop"
    )
}

# =============================================================================
# MAIN LOOP
# =============================================================================

all_results <- list()

for (mt in mod_types) {
  for (r0 in R0_values) {
    
    cat("\n--- Model:", mt, "| R0:", r0, "---\n")
    
    sim_prefix <- file.path(sim_data_dir,
                            paste0(mt, "_R0_", r0, "_n_1e+05_nsim_100"))
    
    # ---- Load sim data ----
    result <- tryCatch({
      
      meta     <- readRDS(paste0(sim_prefix, "_metadata.rds"))
      gen_data <- readRDS(paste0(sim_prefix, "_generation.rds"))
      hist_df  <- readRDS(paste0(sim_prefix, "_total_hist.rds"))
      trans_df <- read.csv(paste0(sim_prefix, "_transitions.csv"))
      
      cat("  Loaded. N =", meta$n, "| nsim =", meta$nsim, "\n")
      
      # ---- GI distribution ----
      cat("  GI distribution:\n")
      gi_dist <- build_gi_dist(gen_data)
      
      # ----------------------------------------------------------
      # METHOD 1: COMPARTMENTAL  Rt = R0 * S_t / N
      # (signal-independent — add to both signals for faceting)
      # ----------------------------------------------------------
      comp_all <- hist_df %>%
        filter(state == "Susceptible") %>%
        select(sim_num, date, S_t = counts) %>%
        mutate(Rt = meta$R0 * S_t / meta$n)
      
      comp_sum <- summarise_rt(comp_all, "Rt", "sim_num")
      cat("  Compartmental sims/date:", paste(range(comp_sum$n_sims), collapse = "–"), "\n")
      
      # store once per signal so it appears in every facet
      for (sig in signals) {
        key <- paste(mt, r0, "Compartmental", sig, sep = "|")
        all_results[[key]] <- comp_sum %>%
          mutate(method = "Compartmental\n(R₀·S/N)", incidence_signal = sig,
                 R0_true = r0, model_type = mt)
      }
      
      # ----------------------------------------------------------
      # METHOD 2: RENEWAL EQUATION — both incidence signals
      # ----------------------------------------------------------
      for (sig in signals) {
        if (!sig %in% colnames(trans_df)) {
          cat("  SKIP renewal [", sig, "] — column not found\n"); next
        }
        
        ren_all <- trans_df %>%
          arrange(id, date) %>%
          group_by(id) %>%
          group_modify(~ tibble(
            date    = .x$date,
            sim_num = .y$id,
            Rt      = renewal_rt_sim(.x[[sig]], gi_dist)
          )) %>%
          ungroup()
        
        ren_sum <- summarise_rt(ren_all, "Rt", "sim_num")
        cat("  Renewal [", sig, "] sims/date:", paste(range(ren_sum$n_sims), collapse = "–"), "\n")
        
        key <- paste(mt, r0, "Renewal", sig, sep = "|")
        all_results[[key]] <- ren_sum %>%
          mutate(method = "Renewal Eq.\n(Direct)", incidence_signal = sig,
                 R0_true = r0, model_type = mt)
      }
      
      # ----------------------------------------------------------
      # METHOD 3 & 4: CORI + WALLINGA-TEUNIS from your .rds files
      # ----------------------------------------------------------
      for (sig in signals) {
        fname <- file.path(cori_wt_dir,
                           paste0(mt, "_R0_", r0_label(r0), "_", sig, ".rds"))
        
        if (!file.exists(fname)) {
          cat("  SKIP Cori/WT [", sig, "] — file not found:", fname, "\n"); next
        }
        
        cw_raw <- readRDS(fname)
        
        # Detect sim ID column
        sim_col <- if ("sim_id" %in% colnames(cw_raw)) "sim_id" else "sim_num"
        cat("  Cori/WT [", sig, "] sims:", n_distinct(cw_raw[[sim_col]]),
            "| methods:", paste(unique(cw_raw$method), collapse = ", "), "\n")
        
        for (est_meth in unique(cw_raw$method)) {
          sub <- cw_raw %>% filter(method == est_meth)
          
          cw_sum <- sub %>%
            group_by(date) %>%
            summarise(
              n_sims    = n_distinct(.data[[sim_col]]),
              mean_Rt   = mean(median_rt),
              median_Rt = median(median_rt),
              lo95      = quantile(median_rt, 0.025),
              hi95      = quantile(median_rt, 0.975),
              lo50      = quantile(median_rt, 0.25),
              hi50      = quantile(median_rt, 0.75),
              .groups   = "drop"
            ) %>%
            mutate(
              method           = est_meth,
              incidence_signal = sig,
              R0_true          = r0,
              model_type       = mt
            )
          
          key <- paste(mt, r0, est_meth, sig, sep = "|")
          all_results[[key]] <- cw_sum
        }
      }
      
      # ----------------------------------------------------------
      # METHOD 5: ABM Rt (epiworldR rt_ci — already summarised across sims)
      # Columns: source_exposure_date, mean_rt, ci_lower, ci_upper, n
      # ----------------------------------------------------------
      abm_file <- file.path(sim_data_dir,
                            paste0(mt, "_R0_", r0, "_n_1e+05_nsim_100_rt_ci.rds"))
      
      if (!file.exists(abm_file)) {
        cat("  SKIP ABM Rt — file not found:", abm_file, "\n")
      } else {
        abm_raw <- readRDS(abm_file)
        cat("  ABM Rt: rows =", nrow(abm_raw),
            "| date range:", min(abm_raw$source_exposure_date),
            "–", max(abm_raw$source_exposure_date), "\n")
        
        abm_sum <- abm_raw %>%
          transmute(
            date      = source_exposure_date,
            n_sims    = n,
            mean_Rt   = mean_rt,
            median_Rt = mean_rt,   # already aggregated — use mean as point estimate
            lo95      = ci_lower,
            hi95      = ci_upper,
            lo50      = NA_real_,  # not available in rt_ci output
            hi50      = NA_real_
          )
        
        # ABM Rt is signal-independent (based on transmission chains, not incidence)
        # Add to both signal facets so it appears as reference in every panel
        for (sig in signals) {
          key <- paste(mt, r0, "ABM", sig, sep = "|")
          all_results[[key]] <- abm_sum %>%
            mutate(method = "ABM\n(epiworldR)", incidence_signal = sig,
                   R0_true = r0, model_type = mt)
        }
      }
      
      "ok"
    }, error = function(e) {
      cat("  ERROR:", conditionMessage(e), "\n"); "error"
    })
  }
}

# =============================================================================
# COMBINE
# =============================================================================

combined <- bind_rows(all_results)

# Clean up method labels for plotting
combined <- combined %>%
  mutate(
    method = case_when(
      method == "Wallinga-Teunis" ~ "Wallinga-\nTeunis",
      TRUE ~ method
    ),
    # Nice facet labels
    signal_label = case_when(
      incidence_signal == "susceptible_to_exposed" ~ "Incidence: S→E",
      incidence_signal == "exposed_to_infected"    ~ "Incidence: E→I"
    ),
    r0_label = paste0("R₀ = ", R0_true),
    model_label = paste0(tools::toTitleCase(model_type), " model")
  )

cat("\nCombined rows:", nrow(combined), "\n")
cat("Methods:      ", paste(unique(combined$method), collapse = " | "), "\n")
cat("R0 values:    ", paste(unique(combined$R0_true), collapse = " | "), "\n")
cat("Model types:  ", paste(unique(combined$model_type), collapse = " | "), "\n")

# =============================================================================
# COLOUR + THEME SYSTEM
# =============================================================================

method_colours <- c(
  "Compartmental\n(R₀·S/N)"  = "#1B998B",   # teal
  "Renewal Eq.\n(Direct)"     = "#E84855",   # red
  "Cori"                      = "#F4A261",   # amber
  "Wallinga-\nTeunis"         = "#6A4C93",   # purple
  "ABM\n(epiworldR)"          = "#2D6A4F"    # dark green
)

method_lty <- c(
  "Compartmental\n(R₀·S/N)"  = "solid",
  "Renewal Eq.\n(Direct)"     = "solid",
  "Cori"                      = "dashed",
  "Wallinga-\nTeunis"         = "dotdash",
  "ABM\n(epiworldR)"          = "longdash"
)

base_theme <- theme_minimal(base_size = 11) +
  theme(
    text             = element_text(family = "sans"),
    plot.title       = element_text(size = 13, face = "bold", margin = margin(b = 4)),
    plot.subtitle    = element_text(size = 9.5, colour = "grey40", margin = margin(b = 10)),
    plot.caption     = element_text(size = 8, colour = "grey55", hjust = 0),
    panel.grid.major = element_line(colour = "grey92", linewidth = 0.4),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "grey99", colour = NA),
    strip.background = element_rect(fill = "grey15", colour = NA),
    strip.text       = element_text(colour = "white", face = "bold", size = 9.5,
                                    margin = margin(4, 6, 4, 6)),
    legend.position  = "bottom",
    legend.title     = element_blank(),
    legend.key.width = unit(1.8, "cm"),
    legend.text      = element_text(size = 9),
    axis.title       = element_text(size = 9.5),
    axis.text        = element_text(size = 8.5),
    plot.margin      = margin(10, 12, 8, 10)
  )

# Annotated true R0 reference lines per facet
r0_refs <- combined %>%
  distinct(r0_label, R0_true, signal_label, model_label)

# =============================================================================
# PLOT 1: FULL MODEL — 2×4 grid (signal × R0), all 4 methods
# =============================================================================

make_main_plot <- function(data, model_filter, title) {
  pd <- data %>% filter(model_type == model_filter)
  
  ggplot(pd, aes(x = date, colour = method, fill = method, linetype = method)) +
    
    # 95% band
    geom_ribbon(aes(ymin = lo95, ymax = hi95), alpha = 0.08, colour = NA) +
    # 50% band
    geom_ribbon(aes(ymin = lo50, ymax = hi50), alpha = 0.18, colour = NA) +
    # Mean line
    geom_line(aes(y = mean_Rt), linewidth = 0.75, na.rm = TRUE) +
    
    # Rt = 1 reference
    geom_hline(yintercept = 1, linetype = "solid", colour = "grey60",
               linewidth = 0.35) +
    
    # True R0 reference per R0 facet
    geom_hline(
      data = r0_refs %>% filter(model_label == paste0(tools::toTitleCase(model_filter), " model")),
      aes(yintercept = R0_true),
      linetype = "dotted", colour = "grey50", linewidth = 0.45,
      inherit.aes = FALSE
    ) +
    
    facet_grid(signal_label ~ r0_label, scales = "free_y") +
    
    scale_colour_manual(values = method_colours, drop = FALSE) +
    scale_fill_manual(values   = method_colours, drop = FALSE) +
    scale_linetype_manual(values = method_lty,   drop = FALSE) +
    scale_x_continuous(breaks = seq(0, 150, 50)) +
    scale_y_continuous(labels = number_format(accuracy = 0.1)) +
    coord_cartesian(ylim = c(0, NA)) +
    
    guides(
      colour   = guide_legend(nrow = 1, override.aes = list(linewidth = 1.4, alpha = 1)),
      fill     = guide_legend(nrow = 1),
      linetype = guide_legend(nrow = 1)
    ) +
    
    labs(
      title    = title,
      subtitle = paste0("Solid line = mean across 100 simulations",
                        "  |  Dark ribbon = 50% band  |  Light ribbon = 95% band",
                        "\nDotted horizontal = true R₀  |  Grey line = Rₜ = 1"),
      x        = "Day",
      y        = expression(R[t]),
      caption  = "Compartmental: R₀·Sₜ/N  |  Renewal Eq: direct fit using empirical GI distribution (gentime > 0, pooled across sims)"
    ) +
    base_theme
}

p_full    <- make_main_plot(combined, "full",    "Rₜ Comparison — Full model (SEIRCONN)")
p_partial <- make_main_plot(combined, "partial", "Rₜ Comparison — Partial model (SEIR small-world)")

ggsave(file.path(out_dir, "rt_full_model.png"),    p_full,    width = 14, height = 7,  dpi = 180)
ggsave(file.path(out_dir, "rt_partial_model.png"), p_partial, width = 14, height = 7,  dpi = 180)
cat("Saved: rt_full_model.png, rt_partial_model.png\n")

# =============================================================================
# PLOT 2: FULL vs PARTIAL SIDE-BY-SIDE — one signal at a time
# =============================================================================

make_model_comparison <- function(data, signal_filter, signal_title) {
  pd <- data %>% filter(incidence_signal == signal_filter)
  
  ggplot(pd, aes(x = date, colour = method, fill = method, linetype = method)) +
    geom_ribbon(aes(ymin = lo95, ymax = hi95), alpha = 0.08, colour = NA) +
    geom_ribbon(aes(ymin = lo50, ymax = hi50), alpha = 0.18, colour = NA) +
    geom_line(aes(y = mean_Rt), linewidth = 0.75, na.rm = TRUE) +
    geom_hline(yintercept = 1, linetype = "solid", colour = "grey60", linewidth = 0.35) +
    geom_hline(
      data = r0_refs,
      aes(yintercept = R0_true),
      linetype = "dotted", colour = "grey50", linewidth = 0.45,
      inherit.aes = FALSE
    ) +
    
    # rows = model type, cols = R0
    facet_grid(model_label ~ r0_label, scales = "free_y") +
    
    scale_colour_manual(values = method_colours, drop = FALSE) +
    scale_fill_manual(values   = method_colours, drop = FALSE) +
    scale_linetype_manual(values = method_lty,   drop = FALSE) +
    scale_x_continuous(breaks = seq(0, 150, 50)) +
    scale_y_continuous(labels = number_format(accuracy = 0.1)) +
    coord_cartesian(ylim = c(0, NA)) +
    guides(
      colour   = guide_legend(nrow = 1, override.aes = list(linewidth = 1.4, alpha = 1)),
      fill     = guide_legend(nrow = 1),
      linetype = guide_legend(nrow = 1)
    ) +
    labs(
      title    = paste0("Full vs Partial model  |  Incidence signal: ", signal_title),
      subtitle = "Solid line = mean  |  Dark ribbon = 50%  |  Light ribbon = 95%  |  Dotted = true R₀",
      x = "Day", y = expression(R[t])
    ) +
    base_theme
}

p_ste <- make_model_comparison(combined, "susceptible_to_exposed", "S→E")
p_eti <- make_model_comparison(combined, "exposed_to_infected",    "E→I")

ggsave(file.path(out_dir, "rt_full_vs_partial_StoE.png"), p_ste, width = 14, height = 8, dpi = 180)
ggsave(file.path(out_dir, "rt_full_vs_partial_EtoI.png"), p_eti, width = 14, height = 8, dpi = 180)
cat("Saved: rt_full_vs_partial_StoE.png, rt_full_vs_partial_EtoI.png\n")

# =============================================================================
# PLOT 3: BIAS vs COMPARTMENTAL TRUTH
# (how much does each statistical estimator deviate from R0·S/N?)
# =============================================================================

comp_ref <- combined %>%
  filter(method == "Compartmental\n(R₀·S/N)") %>%
  select(date, R0_true, model_type, incidence_signal, Rt_true = mean_Rt)

# Bias vs compartmental truth
# ABM is included here intentionally — it shows how epiworldR's chain-based
# Rt compares to the theoretical R0*S/N, which is a meaningful comparison
bias_data <- combined %>%
  filter(method != "Compartmental\n(R₀·S/N)") %>%
  left_join(comp_ref, by = c("date", "R0_true", "model_type", "incidence_signal")) %>%
  mutate(bias = mean_Rt - Rt_true) %>%
  filter(!is.na(bias))

p_bias <- ggplot(bias_data, aes(x = date, y = bias, colour = method, linetype = method)) +
  geom_hline(yintercept = 0, colour = "grey50", linewidth = 0.4) +
  geom_line(linewidth = 0.7, na.rm = TRUE) +
  
  facet_grid(model_label + signal_label ~ r0_label, scales = "free_y") +
  
  scale_colour_manual(values = method_colours) +
  scale_linetype_manual(values = method_lty) +
  scale_x_continuous(breaks = seq(0, 150, 50)) +
  guides(
    colour   = guide_legend(nrow = 1, override.aes = list(linewidth = 1.5)),
    linetype = guide_legend(nrow = 1)
  ) +
  labs(
    title    = "Estimation bias vs. Compartmental truth (R₀·Sₜ/N)",
    subtitle = "Bias = method mean Rₜ − compartmental Rₜ  |  Positive = overestimate",
    x = "Day", y = "Bias"
  ) +
  base_theme +
  theme(strip.text = element_text(size = 8))

ggsave(file.path(out_dir, "rt_bias_all.png"), p_bias, width = 14, height = 12, dpi = 180)
cat("Saved: rt_bias_all.png\n")

# =============================================================================
# PLOT 4: CORI vs WALLINGA-TEUNIS head-to-head (same signal, side by side)
# =============================================================================

cw_data <- combined %>%
  filter(method %in% c("Cori", "Wallinga-\nTeunis"))

# Note: ABM Rt is signal-independent (transmission chains, not incidence),
# so it is duplicated across both signal facets — this is intentional.
# The 50% ribbon is NA for ABM (rt_ci only provides overall 95% CI).

p_cw <- ggplot(cw_data, aes(x = date, colour = method, fill = method, linetype = method)) +
  geom_ribbon(aes(ymin = lo95, ymax = hi95), alpha = 0.10, colour = NA) +
  geom_ribbon(aes(ymin = lo50, ymax = hi50), alpha = 0.20, colour = NA) +
  geom_line(aes(y = mean_Rt), linewidth = 0.8, na.rm = TRUE) +
  geom_hline(yintercept = 1, colour = "grey60", linewidth = 0.35) +
  geom_hline(
    data = r0_refs,
    aes(yintercept = R0_true),
    linetype = "dotted", colour = "grey50", linewidth = 0.45,
    inherit.aes = FALSE
  ) +
  
  facet_grid(model_label + signal_label ~ r0_label, scales = "free_y") +
  
  scale_colour_manual(values = method_colours) +
  scale_fill_manual(values   = method_colours) +
  scale_linetype_manual(values = method_lty) +
  scale_x_continuous(breaks = seq(0, 150, 50)) +
  guides(
    colour   = guide_legend(nrow = 1, override.aes = list(linewidth = 1.5, alpha = 1)),
    fill     = guide_legend(nrow = 1),
    linetype = guide_legend(nrow = 1)
  ) +
  labs(
    title    = "Cori vs Wallinga-Teunis — head-to-head",
    subtitle = "Both use same per-sim SI (mean/SD from empirical GI)  |  Dotted = true R₀",
    x = "Day", y = expression(R[t])
  ) +
  base_theme +
  theme(strip.text = element_text(size = 8))

ggsave(file.path(out_dir, "rt_cori_vs_wt.png"), p_cw, width = 14, height = 12, dpi = 180)
cat("Saved: rt_cori_vs_wt.png\n")

# =============================================================================
# SAVE COMBINED DATA
# =============================================================================

saveRDS(combined, file.path(out_dir, "rt_all_methods_combined.rds"))
write.csv(combined, file.path(out_dir, "rt_all_methods_combined.csv"), row.names = FALSE)
cat("\nAll done. Figures saved to:", out_dir, "\n")