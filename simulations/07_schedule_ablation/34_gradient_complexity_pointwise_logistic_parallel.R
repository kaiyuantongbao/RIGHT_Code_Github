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

# Pointwise gradient-estimation diagnostic, logistic regression.
# The population gradient is approximated by Monte Carlo. Increase M_ref for final plots.

cfg <- list(
  model = "logistic",
  n_grad = 1000,
  p = 600,
  s_star = 5,
  s_alg = 10,
  theta_seed = 20260425L,
  theta_lower = 0.5,
  theta_upper = 1.5,
  theta_magnitude = NULL,
  random_signs = TRUE,
  df_X = 2.5,
  scale_X = 1.0,
  df_eps = 1.5,       # unused for logistic; kept for common config structure
  scale_eps = 1.0,    # unused for logistic
  r_grid = c(1.0, 0.5, 0.25, 0.10, 0.05),
  n_directions = 8,
  direction_sparsity = 5,
  avoid_true_support_in_directions = FALSE,
  reps = 20,
  m_min = 10,
  K1_multiplier = 1.0,
  K2_multipliers = c(1, 2, 4, 8),
  K1_formula = "s_log_p",
  n_cores = max(1L, parallel::detectCores(logical = FALSE) - 1L),
  seed_base = 2026050401L,
  M_ref = 30000L,
  ref_chunk_size = 2000L
)

run_gradient_complexity_pointwise(cfg)
