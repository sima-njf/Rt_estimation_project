# run_one.R
run_one_combo <- function(combo) {
  suppressMessages({
    library(epiworldR); library(data.table); library(dplyr)
    library(tidyverse); library(igraph)
  })
  source("/uufs/chpc.utah.edu/common/home/u1418987/sima/Rt_estimation_project/data001.R")
  
  n <- 1e5
  save_seir_simulation_data(
    model_type      = combo$model_type,
    name            = "Covid",
    n               = n,
    prevalence      = 100 / n,
    contact_rate    = 20.0,
    R0              = combo$R0,
    recovery_rate   = 1.0 / 7.0,
    incubation_days = 4,
    ndays           = 150,
    nsim            = 1000,
    seed            = 1234,
    save_dir        = "/scratch/general/vast/u1418987/saved_data/",  # absolute scratch
    deg             = 20
  )
}