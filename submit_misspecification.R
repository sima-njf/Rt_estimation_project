# submit_misspec.R
# Fans out the underreporting-misspecification analysis:
# one array task per (model_type, R0). Run AFTER the 1000-sim files exist.
# ---------------------------------------------------------------------------

library(slurmR)
source("misspecific_one.R")

R0_values   <- c(1.5, 2.0, 3.0, 5.0)
model_types <- c("full", "partial")
combos <- list()
for (mt in model_types)
  for (r0 in R0_values)
    combos[[length(combos) + 1]] <- list(model_type = mt, R0 = r0)

misspec_job <- Slurm_lapply(
  X          = combos,
  FUN        = process_one_misspec_combo,
  njobs      = 60,          # 8 array tasks
  mc.cores   = 1,
  job_name   = "misspec_analysis",
  plan       = "submit",
  sbatch_opt = list(
    account         = "vegayon-np",
    partition       = "vegayon-np",
    `cpus-per-task` = 4,
    `mem-per-cpu`   = "8G",
    time            = "8:00:00"
  ),
  # Ship the function AND all the constants it references
  export = c("process_one_misspec_combo", "CODE_DIR", "BASE", "SAVED", "RESULTS", "NSIM",
             "MISSPEC_CODE_FILE", "REPORTING_RATE", "MAX_SIMS", "START_DAY", "RNG_SEED"),
  tmp_path   = "/scratch/general/vast/u1418987/slurmR"
)

res <- Slurm_collect(misspec_job)
print(res)   # OK / FAILED per combo