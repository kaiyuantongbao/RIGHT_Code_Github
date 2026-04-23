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

# Fixed-iterate versus uniform-class diagnostic, linear regression.
# M = 1 mimics a fixed iterate. Larger M mimics needing simultaneous control over
# many candidate sparse iterates/directions. This is closer to the actual global-to-local
# complexity distinction than the pointwise K1-vs-K2 plot alone.

cfg <- list(
  model = "linear",
  n_grad = 2000,
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
  r_grid = c(2.0, 1.0, 0.5, 0.25),
  candidate_grid = c(1, 5, 20, 50, 100),
  direction_sparsity = 10,
  avoid_true_support_in_directions = FALSE,
  reps = 20,
  m_min = 10,
  K1_multiplier = 1.0,
  K2_multipliers = c(1, 2, 4, 8),
  K1_formula = "s_log_p",
  n_cores = max(1L, parallel::detectCores(logical = FALSE) - 1L),
  seed_base = 2026050501L,
  M_ref = 0L,              # unused for linear
  ref_chunk_size = 0L      # unused for linear
)

run_gradient_complexity_uniform(cfg)
