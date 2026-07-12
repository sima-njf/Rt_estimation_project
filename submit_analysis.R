# submit_analysis.R
library(slurmR)
source("analysis_one.R")

R0_values   <- c(1.5, 2.0, 3.0, 5.0)
model_types <- c("full", "partial")

combos <- list()
for (mt in model_types)
  for (r0 in R0_values)
    combos[[length(combos) + 1]] <- list(model_type = mt, R0 = r0)

job <- Slurm_lapply(
  X          = combos,
  FUN        = process_one_combo,
  njobs      = 60,          # 8 array tasks
  mc.cores   = 1,
  job_name   = "cw_analysis",
  plan       = "submit",
  sbatch_opt = list(
    account         = "vegayon-np",
    partition       = "vegayon-np",
    `cpus-per-task` = 8,                 # analysis is lighter than simulation
    `mem-per-cpu`   = "8G",              # transmission tables are the memory driver
    time            = "8:00:00"
  ),
  export = c("process_one_combo", "CODE_DIR", "BASE", "SAVED", "RESULTS", "NSIM"),
  tmp_path   = "/scratch/general/vast/u1418987/slurmR"
)

res <- Slurm_collect(job)
print(res)   # shows OK / EMPTY / MISSING per combo
