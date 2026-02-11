# =============================================================================
# RtEstim step-by-step with YOUR epiworldR ABM data
# =============================================================================

library(rtestim)
library(dplyr)
library(ggplot2)

source("data001.R")

# =============================================================================
# STEP 1: Load one dataset
# =============================================================================

model_type <- "partial"
R0_val     <- 2.0

prefix <- paste0(
  "saved_data/",
  model_type, "_R0_", R0_val, "_n_1e+05_nsim_100"
)
cat("Loading from:", prefix, "\n")

saved_data <- load_seir_data(prefix)

names(saved_data)
# transitions, reproductive, rt_ci, generation, gi_stats,
# transmission, total_hist, metadata

# =============================================================================
# STEP 2: Pick ONE simulation
# =============================================================================

sim_id <- 1

# =============================================================================
# STEP 3: Build incidence — TWO APPROACHES
# =============================================================================
#
# You have TWO ways to get daily incidence:
#
#   A) $transitions — already has aggregated daily counts per sim
#      columns: susceptible_to_exposed, exposed_to_infected, etc.
#      BUT: day 0 S->E count INCLUDES the 100 seeds!
#
#   B) $transmission — individual-level events
#      has source == -1 for seeds, so you can filter them out
#      date = day the target became Exposed (S->E)
#      source_exposure_date = day the SOURCE was exposed (not the target's E->I)
#

# ---- APPROACH A: $transitions (already aggregated) --------------------------

sim_transitions <- saved_data$transitions %>% filter(id == sim_id)

head(sim_transitions, 10)

# S -> E counts per day
incidence_s2e_A <- sim_transitions$susceptible_to_exposed
# E -> I counts per day
incidence_e2i_A <- sim_transitions$exposed_to_infected
days_A          <- sim_transitions$date

cat("S->E incidence (first 15 days):\n")
print(data.frame(day = days_A[1:15],
                 s2e = incidence_s2e_A[1:15],
                 e2i = incidence_e2i_A[1:15]))

# NOTE: day 0 S->E = 100 includes seeds!

par(mfrow = c(1, 2))
plot(days_A, incidence_s2e_A, type = "h", col = "steelblue",
     main = "S -> E (daily)", xlab = "Day", ylab = "Count")
plot(days_A, incidence_e2i_A, type = "h", col = "coral",
     main = "E -> I (daily)", xlab = "Day", ylab = "Count")
par(mfrow = c(1, 1))


# ---- APPROACH B: $transmission (individual-level, remove seeds) -------------

sim_transmission <- saved_data$transmission %>% filter(sim_num == sim_id)

cat("\nTotal transmission rows:", nrow(sim_transmission), "\n")
cat("Seeds (source == -1):", sum(sim_transmission$source == -1), "\n")
cat("Real infections:", sum(sim_transmission$source != -1), "\n")

# REMOVE SEEDS
sim_trans_clean <- sim_transmission %>% filter(source != -1)

# S -> E: use `date` (day the target became exposed)
dates_s2e_B     <- sim_trans_clean$date
date_range_s2e  <- min(dates_s2e_B):max(dates_s2e_B)
incidence_s2e_B <- as.numeric(table(factor(dates_s2e_B, levels = date_range_s2e)))

cat("\nS->E from $transmission (seeds removed), first 15 days:\n")
print(data.frame(day = date_range_s2e[1:15], count = incidence_s2e_B[1:15]))

# NOTE about E -> I from $transmission:
# source_exposure_date = when the SOURCE was exposed, NOT the target's E->I date.
# You don't have target E->I date in $transmission.
# For E->I incidence, use $transitions$exposed_to_infected instead.

# =============================================================================
# STEP 4: Build generation interval PMF — REMOVE SEEDS
# =============================================================================

sim_generation <- saved_data$generation %>% filter(sim_num == sim_id)

# Find seed agent IDs (those with source == -1 in transmission)
seed_ids <- sim_transmission %>%
  filter(source == -1) %>%
  pull(target)

cat("\nSeed agents:", length(seed_ids), "\n")

# Remove generation times where the source is a seed
sim_gen_clean <- sim_generation %>%
  filter(!(source %in% seed_ids))

cat("Generation rows before:", nrow(sim_generation), "\n")
cat("Generation rows after removing seeds:", nrow(sim_gen_clean), "\n")

gen_times <- sim_gen_clean$gentime
gen_times <- gen_times[!is.na(gen_times)]

cat("\nGeneration time summary (no seeds):\n")
cat("  n:", length(gen_times), "\n")
cat("  mean:", round(mean(gen_times), 2), "\n")
cat("  median:", median(gen_times), "\n")
cat("  range:", min(gen_times), "to", max(gen_times), "\n")

# Build PMF
max_gen    <- max(gen_times)
gen_counts <- table(factor(gen_times, levels = 0:max_gen))
gen_int    <- as.numeric(gen_counts) / sum(gen_counts)

cat("\nGeneration interval PMF:\n")
print(round(gen_int[1:min(15, length(gen_int))], 4))
cat("Sums to:", sum(gen_int), "\n")

barplot(gen_int, names.arg = 0:max_gen,
        main = "Generation interval PMF (seeds removed)",
        xlab = "Days", ylab = "Probability", col = "coral")

# =============================================================================
# STEP 5: Fit rtestim with S -> E incidence
# =============================================================================

# Using approach B (seeds removed from $transmission)
# rtestim needs first value > 0
cat("\nFirst value of S->E incidence:", incidence_s2e_B[1], "\n")

# If first value is 0, trim leading zeros
if (incidence_s2e_B[1] == 0) {
  first_nz <- which(incidence_s2e_B > 0)[1]
  incidence_s2e_B <- incidence_s2e_B[first_nz:length(incidence_s2e_B)]
  date_range_s2e  <- date_range_s2e[first_nz:length(date_range_s2e)]
  cat("Trimmed leading zeros, starts at day", date_range_s2e[1], "\n")
}

# --- 5a: Basic fit (no CV) ---
fit_basic <- estimate_rt(
  observed_counts = incidence_s2e_B,
  x               = seq_along(incidence_s2e_B),
  delay_distn     = gen_int,
  nsol            = 20
)

cat("\nBasic fit: ", length(fit_basic$lambda), "lambda values\n")
plot(fit_basic, main = "All 20 lambda solutions (S->E)")

# --- 5b: Cross-validated fit ---
fit_cv <- cv_estimate_rt(
  observed_counts = incidence_s2e_B,
  x               = seq_along(incidence_s2e_B),
  delay_distn     = gen_int,
  nsol            = 20
)

cat("lambda.min:", fit_cv$lambda.min, "\n")
cat("lambda.1se:", fit_cv$lambda.1se, "\n")

plot(fit_cv)  # CV error curve

# --- 5c: Confidence bands ---
cb_s2e <- confband(fit_cv, lambda = "lambda.1se")

str(cb_s2e)

# --- 5d: Build output & trim ---
out_s2e <- data.frame(
  day   = date_range_s2e,
  rt    = cb_s2e$fit,
  lower = cb_s2e$`2.5%`,
  upper = cb_s2e$`97.5%`
)

cat("Before trim:", nrow(out_s2e), "\n")
out_s2e <- out_s2e %>% filter(day > max_gen)
cat("After trim (day >", max_gen, "):", nrow(out_s2e), "\n")

# --- 5e: Plot ---
ggplot(out_s2e, aes(x = day, y = rt)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "steelblue", alpha = 0.3) +
  geom_line(color = "steelblue", linewidth = 1) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  labs(
    title = paste0("Rt (S->E, seeds removed) — Sim ", sim_id,
                   " (R0=", R0_val, ")"),
    x = "Day", y = "Rt"
  ) +
  theme_bw()

# =============================================================================
# STEP 6: Fit rtestim with E -> I incidence
# =============================================================================

# E -> I comes from $transitions (no individual-level source to filter)
# Day 0 E->I is likely 0 (nobody has transitioned E->I yet on day 0)
incidence_e2i_use <- sim_transitions$exposed_to_infected
days_e2i          <- sim_transitions$date

cat("\nE->I first 10 days:\n")
print(data.frame(day = days_e2i[1:10], e2i = incidence_e2i_use[1:10]))

# Trim leading zeros (rtestim needs first > 0)
first_nz_e2i <- which(incidence_e2i_use > 0)[1]
cat("First non-zero E->I on day:", days_e2i[first_nz_e2i], "\n")

incidence_e2i_use <- incidence_e2i_use[first_nz_e2i:length(incidence_e2i_use)]
days_e2i          <- days_e2i[first_nz_e2i:length(days_e2i)]

fit_cv_e2i <- cv_estimate_rt(
  observed_counts = incidence_e2i_use,
  x               = seq_along(incidence_e2i_use),
  delay_distn     = gen_int,
  nsol            = 20
)

cb_e2i <- confband(fit_cv_e2i, lambda = "lambda.1se")

out_e2i <- data.frame(
  day   = days_e2i,
  rt    = cb_e2i$fit,
  lower = cb_e2i$`2.5%`,
  upper = cb_e2i$`97.5%`
) %>% filter(day > max_gen)

# =============================================================================
# STEP 7: Compare S->E vs E->I
# =============================================================================

out_s2e$transition <- "S -> E"
out_e2i$transition <- "E -> I"
both <- bind_rows(out_s2e, out_e2i)

ggplot(both, aes(x = day, y = rt, color = transition, fill = transition)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.15, color = NA) +
  geom_line(linewidth = 0.8) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  labs(
    title = paste0("Rt: S->E vs E->I — Sim ", sim_id, " (R0=", R0_val, ")"),
    x = "Day", y = "Rt", color = NULL, fill = NULL
  ) +
  theme_bw()

cat("\n=== DONE ===\n")
