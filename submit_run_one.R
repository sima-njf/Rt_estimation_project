# submit_sim.R
library(slurmR)
source("run_one.R")

R0_values   <- c(1.5, 2.0, 3.0, 5.0)
model_types <- c("full", "partial")
combos <- list()
for (mt in model_types)
  for (r0 in R0_values)
    combos[[length(combos) + 1]] <- list(model_type = mt, R0 = r0)

sim_job <- Slurm_lapply(
  X          = combos,
  FUN        = run_one_combo,
  njobs      = length(combos),          # 8 array tasks
  mc.cores   = 1,
  job_name   = "seir_sim",
  plan       = "submit",
  sbatch_opt = list(
    account         = "vegayon-np",
    partition       = "vegayon-np",
    `cpus-per-task` = 18,
    `mem-per-cpu`   = "8G",             # ~144 GB/task; node has 508 GB
    time            = "1-00:00:00"
  ),
  export     = c("run_one_combo"),
  tmp_path   = "/scratch/general/vast/u1418987/slurmR"
)

cat("Submitted simulation job. ID:\n")
print(sim_job)
res <- Slurm_collect(sim_job)
