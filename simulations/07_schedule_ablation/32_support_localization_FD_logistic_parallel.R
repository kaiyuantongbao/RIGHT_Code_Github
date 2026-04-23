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

# Empirical active-set diagnostic for FD-RIGHT, logistic regression.
# The default tail/scale setting is moderate. The older df_X = 2.1, scale_X = 3
# can be used as a stress setting, but it may introduce near-separation.

cfg <- list(
  model = "logistic",
  n = 1000,
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
  eta = 0.02,
  T_max = 800,
  record_every = 5,
  record_dense_until = 60,
  window_records = 5,
  m_min = 10,
  K1_multiplier = 1.0,
  K1_formula = "s_log_p",
  reps = 100,
  n_cores = max(1L, parallel::detectCores(logical = FALSE) - 1L),
  radius = Inf,
  noise_freq_cap = NA_real_,
  seed_base = 2026050201L
)

run_support_localization(cfg)
