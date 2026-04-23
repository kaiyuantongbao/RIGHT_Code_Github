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

# Pointwise gradient-estimation diagnostic, linear regression.
# This is the corrected version of the previous K1 vs K2 plot: K2 is now calibrated
# over constant multipliers c_K in K2 = ceil(c_K log p). The linear population
# gradient is computed in closed form for the t design, requiring df_X > 2.

cfg <- list(
  model = "linear",
  n_grad = 2000,
  p = 600,
  s_star = 25,
  s_alg = 50,
  theta_seed = 20260425L,
  theta_lower = 1.0,
  theta_upper = 2.0,
  theta_magnitude = NULL,
  random_signs = TRUE,
  df_X = 2.5,
  scale_X = 1.0,
  df_eps = 1.5,
  scale_eps = 1.0,
  r_grid = c(4.0, 2.0, 1.0, 0.5, 0.25, 0.10),
  n_directions = 12,
  direction_sparsity = 10,
  avoid_true_support_in_directions = FALSE,
  reps = 30,
  m_min = 10,
  K1_multiplier = 1.0,
  K2_multipliers = c(10,12,14),
  K1_formula = "s_log_p",
  n_cores = max(1L, parallel::detectCores(logical = FALSE) - 1L),
  seed_base = 2026050301L,
  M_ref = 0L,              # unused for linear
  ref_chunk_size = 0L      # unused for linear
)

run_gradient_complexity_pointwise(cfg)
