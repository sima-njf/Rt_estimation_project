# ==================================================
# Load libraries
# ==================================================
library(epiworldR)    # For ABM simulations
library(ggplot2)      # For plotting
library(dplyr)        # For data manipulation
library(tidyr)        # For pivot_wider
library(data.table)   # Efficient data handling
library(tidyverse)    # Includes dplyr, ggplot2, tidyr, etc.
library(netdiffuseR)  # For small-world network generation

# ==================================================
# Shared Parameters
# ==================================================
R0 <- 2
name <- "Covid"
n <- 1e5
prevalence <- 100 / n
contact_rate <- 20.0
recovery_rate <- 1.0 / 7.0
incubation_days <- 4
ndays <- 150
nsim <- 100
seed <- 1234

# Transmission rate based on contact rate and R0
transmission_rate <- R0 * recovery_rate / contact_rate

# Saver for both models
saver <- make_saver("total_hist", "transmission", "transition", "reproductive", "generation")

# ==================================================
# 1. Fully Connected Network Model (SEIRCONN)
# ==================================================
model_conn <- ModelSEIRCONN(
  name = name,
  n = n,
  prevalence = prevalence,
  contact_rate = contact_rate,
  transmission_rate = transmission_rate,
  recovery_rate = recovery_rate,
  incubation_days = incubation_days
)

set.seed(seed)
cat("Running fully connected model...\n")
run_multiple(model_conn, ndays = ndays, nsim = nsim, seed = seed, saver = saver, nthreads = 18)
results_conn <- run_multiple_get_results(model_conn, nthreads = 18)

# Reproductive number for fully connected
reproductive_conn <- results_conn$reproductive %>% filter(source >= 0)
ci_rt_conn <- reproductive_conn %>%
  group_by(sim_num, source_exposure_date) %>%
  summarise(rt = mean(rt), .groups = "drop") %>%
  group_by(source_exposure_date) %>%
  summarise(
    mean_rt = mean(rt),
    ci_lower = quantile(rt, 0.025),
    ci_upper = quantile(rt, 0.975),
    .groups = "drop"
  ) %>%
  mutate(model = "Fully Connected")

# ==================================================
# 2. Small-World Network Model (netdiffuseR)
# ==================================================
set.seed(seed)

deg <- 50 # Average degree in small-world network
cat("Generating small-world network using netdiffuseR...\n")
net <- rgraph_ws(n = n, p = 1, k = deg, undirected = TRUE)

# Convert adjacency matrix to edgelist (0-based indexing)
net2 <- which(as.matrix(net) != 0, arr.ind = TRUE)
net2 <- net2[net2[, 1] < net2[, 2], ]
net2 <- net2 - 1
net2 <- matrix(as.integer(net2), ncol = 2)

# Transmission rate scaled by degree to match R0
transmission_rate_sw <- R0 * recovery_rate / (deg)

# Build SEIR model
model_sw <- ModelSEIR(
  name = paste(name, "smallworld_netdiffuseR"),
  prevalence = prevalence,
  transmission_rate = transmission_rate_sw,
  incubation_days = incubation_days,
  recovery_rate = recovery_rate
)

# Add the small-world network
agents_from_edgelist(
  model = model_sw,
  source = net2[, 1],
  target = net2[, 2],
  size = as.integer(n),
  directed = FALSE
)

# Run simulation
cat("Running small-world model...\n")
run_multiple(model_sw, ndays = ndays, nsim = nsim, seed = seed, saver = saver, nthreads = 18)
results_sw <- run_multiple_get_results(model_sw, nthreads = 18)

# Reproductive number for small-world
reproductive_sw <- results_sw$reproductive %>% filter(source >= 0)
ci_rt_sw <- reproductive_sw %>%
  group_by(sim_num, source_exposure_date) %>%
  summarise(rt = mean(rt), .groups = "drop") %>%
  group_by(source_exposure_date) %>%
  summarise(
    mean_rt = mean(rt),
    ci_lower = quantile(rt, 0.025),
    ci_upper = quantile(rt, 0.975),
    .groups = "drop"
  ) %>%
  mutate(model = "Small-World")

# ==================================================
# 3. Plotting Both Models Together
# ==================================================
rt_combined <- bind_rows(ci_rt_conn, ci_rt_sw)

ggplot(rt_combined, aes(x = source_exposure_date, y = mean_rt, color = model, fill = model)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.2, color = NA) +
  geom_hline(yintercept = R0, linetype = "dashed", color = "black") +
  annotate("text",
           x = max(rt_combined$source_exposure_date) * 0.9,
           y = R0 + 0.1,
           label = paste("R₀ =", R0),
           color = "black",
           size = 4) +
  labs(
    title = "Estimated Reproductive Number (Rt) Over Time",
    subtitle = "Comparison: Fully Connected vs Small-World Network (netdiffuseR)",
    x = "Days Since Start",
    y = "Rt",
    caption = paste("Simulations:", nsim, "| Population:", format(n, big.mark = ","))
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "top") +
  coord_cartesian(ylim = c(0, max(rt_combined$ci_upper) * 1.1)) +
  scale_color_manual(values = c("Fully Connected" = "blue", "Small-World" = "darkred")) +
  scale_fill_manual(values = c("Fully Connected" = "blue", "Small-World" = "darkred"))
