rm(list = ls())

script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- args[grepl("^--file=", args)]
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]), mustWork = TRUE)))
  }
  getwd()
}
source(file.path(script_dir(), "utils_diagnostic_core_parallel.R"))

# Empirical active-set diagnostic for FD-RIGHT, linear regression.
# This uses the repository solver_right() with T_max = 1 repeatedly to record the trajectory.
# The default signal is intentionally weaker than theta_j = 5, because theta_j = 5 made
# the true support selected almost immediately and left little trajectory information.

cfg <- list(
  model = "linear",
  n = 2000,
  p = 600,
  s_star = 10,
  s_alg = 20,
  theta_seed = 20260425L,
  theta_lower = 1.0,
  theta_upper = 2.0,
  theta_magnitude = NULL,
  random_signs = TRUE,
  df_X = 2.5,
  scale_X = 1.0,
  df_eps = 1.5,
  scale_eps = 1.0,
  eta = 0.01,
  T_max = 500,
  record_every = 5,
  record_dense_until = 40,
  window_records = 5,
  m_min = 10,
  K1_multiplier = 1.0,
  K1_formula = "s_log_p",
  reps = 100,
  n_cores = max(1L, parallel::detectCores(logical = FALSE) - 1L),
  radius = Inf,
  noise_freq_cap = NA_real_,
  seed_base = 2026050101L
)

run_support_localization(cfg)
