# ======================================================================
# TS-RIGHT only: logistic-regression comparison across c_K1 = 1.5, 1.8
# No file writing. Prints summary and plot only.
# Run from the repository root.
# ======================================================================

suppressPackageStartupMessages({
  library(mvtnorm)
  library(MASS)
  library(dplyr)
  library(ggplot2)
  library(here)
  library(foreach)
  library(doParallel)
  library(lhs)
})

source(here("R", "data_generator.R"))
source(here("R", "solvers.R"))
source(here("simulations", "07_schedule_ablation", "utils_schedule_right.R"))

generate_X_t_LHS <- function(n, p, df_X, scale_X = 1, ...) {
  hypercube <- lhs::randomLHS(n, p)
  X_raw <- stats::qt(hypercube, df = df_X)
  X_scaled <- X_raw * scale_X
  return(X_scaled)
}

P_GRID <- c(600)
N_GRID <- c(800,1100,1400,1700,2000)
N_REPS <- 100

C_K1_GRID <- c(0.75)
C_K2 <- 0.50

S_STAR <- 5
S_ALG <- 8
T1 <- 2000
T2 <- 2
Q_TS <- 0.20
ETA_TS <- 0.02
M_MIN <- 10
DF_X <- 2.05
SCALE_X <- 4.0

USE_PARALLEL <- TRUE
N_CORES <- max(1L, parallel::detectCores() - 1L)

env_int <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) return(default)
  as.integer(value)
}

env_num_vec <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) return(default)
  as.numeric(strsplit(value, ",", fixed = TRUE)[[1]])
}

env_bool <- function(name, default) {
  value <- tolower(Sys.getenv(name, unset = ""))
  if (!nzchar(value)) return(default)
  value %in% c("1", "true", "t", "yes", "y")
}

P_GRID <- as.integer(env_num_vec("REC_P_GRID", P_GRID))
N_GRID <- as.integer(env_num_vec("REC_N_GRID", N_GRID))
N_REPS <- env_int("REC_N_REPS", N_REPS)
C_K1_GRID <- env_num_vec("REC_C_K1_GRID", C_K1_GRID)
C_K2 <- env_num_vec("REC_C_K2", C_K2)[1]
USE_PARALLEL <- env_bool("REC_USE_PARALLEL", USE_PARALLEL)
N_CORES <- env_int("REC_N_CORES", N_CORES)

make_theta_recommended <- function(p) {
  theta <- rep(0, p)
  theta[seq_len(S_STAR)] <- c(1, -1, 1, -1, 1)
  theta
}

fit_ts_right_logistic_ck <- function(
    X,
    y,
    theta_star,
    s,
    eta,
    q,
    T1,
    T2,
    m_min,
    c_K1,
    c_K2,
    s_ref,
    seed
) {
  n <- nrow(X)
  p <- ncol(X)
  
  budget <- compute_ts_budget(
    n = n,
    p = p,
    s = s,
    q = q,
    T1 = T1,
    T2 = T2,
    s_ref = s_ref,
    m_min = m_min,
    c_K1 = c_K1,
    c_K2 = c_K2
  )
  
  meta <- data.frame(
    eligible = budget$eligible,
    n1 = budget$n1,
    n2 = budget$n2,
    b2 = budget$b2,
    K1 = budget$K1,
    K2 = budget$K2,
    K1_block_size = budget$n1 / budget$K1,
    K2_block_size = budget$b2 / budget$K2,
    T_fd_budget = budget$T_fd_budget
  )
  
  if (!isTRUE(budget$eligible)) {
    return(list(theta = rep(NA_real_, p), l2_error = NA_real_, meta = meta))
  }
  
  rng_state <- with_local_seed(seed)
  on.exit(restore_local_seed(rng_state), add = TRUE)
  
  idx <- sample.int(n)
  id1 <- idx[seq_len(budget$n1)]
  id2 <- idx[(budget$n1 + 1L):n]
  
  theta_cur <- solver_right(
    X = X[id1, , drop = FALSE],
    y = y[id1],
    s = s,
    eta = eta,
    T_max = T1,
    K = budget$K1,
    theta_init = rep(0, p),
    grad_func_samplewise = grad_logistic_regression_samplewise,
    record_trace = FALSE
  )
  
  use_n2 <- budget$b2 * T2
  id2 <- sample(id2, size = length(id2), replace = FALSE)[seq_len(use_n2)]
  
  for (tt in seq_len(T2)) {
    batch_ids <- id2[((tt - 1L) * budget$b2 + 1L):(tt * budget$b2)]
    theta_cur <- solver_right(
      X = X[batch_ids, , drop = FALSE],
      y = y[batch_ids],
      s = s,
      eta = eta,
      T_max = 1L,
      K = budget$K2,
      theta_init = theta_cur,
      grad_func_samplewise = grad_logistic_regression_samplewise,
      record_trace = FALSE
    )
  }
  
  list(
    theta = as.numeric(theta_cur),
    l2_error = sqrt(sum((theta_cur - theta_star)^2)),
    meta = meta
  )
}

one_run <- function(p, n, rep_id, c_K1) {
  theta_star <- make_theta_recommended(p)
  Sigma_t <- diag(p)
  
  set.seed(1000000 + 1000 * n + rep_id)
  
  sim_data <- generate_logistic_data(
    n = n,
    p = p,
    theta_star = theta_star,
    generator_X = generate_X_t_LHS,
    df_X = DF_X,
    Sigma_X = Sigma_t,
    scale_X = SCALE_X
  )
  
  out_ts <- fit_ts_right_logistic_ck(
    X = sim_data$X,
    y = sim_data$y,
    theta_star = theta_star,
    s = S_ALG,
    eta = ETA_TS,
    q = Q_TS,
    T1 = T1,
    T2 = T2,
    m_min = M_MIN,
    c_K1 = c_K1,
    c_K2 = C_K2,
    s_ref = S_STAR,
    seed = 2000000 + 1000 * n + rep_id
  )
  
  ts_meta <- out_ts$meta
  
  data.frame(
    p = p,
    n = n,
    rep_id = rep_id,
    method = "TS-RIGHT",
    l2_error = out_ts$l2_error,
    df_x = DF_X,
    scale_x = SCALE_X,
    q = Q_TS,
    T1 = T1,
    T2 = T2,
    c_K1 = c_K1,
    c_K2 = C_K2,
    K1 = if ("K1" %in% names(ts_meta)) ts_meta$K1[1] else NA_integer_,
    K2 = if ("K2" %in% names(ts_meta)) ts_meta$K2[1] else NA_integer_,
    K1_block_size = if ("K1_block_size" %in% names(ts_meta)) ts_meta$K1_block_size[1] else NA_real_,
    K2_block_size = if ("K2_block_size" %in% names(ts_meta)) ts_meta$K2_block_size[1] else NA_real_,
    ts_eligible = if ("eligible" %in% names(ts_meta)) ts_meta$eligible[1] else NA,
    seed_data = 1000000 + 1000 * n + rep_id,
    seed_ts = 2000000 + 1000 * n + rep_id
  )
}

param_grid <- expand.grid(
  p = P_GRID,
  n = N_GRID,
  rep_id = seq_len(N_REPS),
  c_K1 = C_K1_GRID,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

message("Running TS-RIGHT only logistic comparison.")
message("Tasks: ", nrow(param_grid))
message("P_GRID = ", paste(P_GRID, collapse = ","))
message("N_GRID = ", paste(N_GRID, collapse = ","))
message("N_REPS = ", N_REPS)
message("C_K1_GRID = ", paste(C_K1_GRID, collapse = ","))
message("C_K2 = ", C_K2)
message("Parallel = ", USE_PARALLEL, "; cores = ", N_CORES)

if (USE_PARALLEL && N_CORES > 1L) {
  cl <- parallel::makeCluster(N_CORES)
  doParallel::registerDoParallel(cl)
  on.exit({
    try(parallel::stopCluster(cl), silent = TRUE)
  }, add = TRUE)
  
  results <- foreach(
    ii = seq_len(nrow(param_grid)),
    .combine = dplyr::bind_rows,
    .packages = c("mvtnorm", "MASS", "dplyr", "here", "lhs")
  ) %dopar% {
    source(here("R", "data_generator.R"))
    source(here("R", "solvers.R"))
    source(here("simulations", "07_schedule_ablation", "utils_schedule_right.R"))
    
    generate_X_t_LHS <- function(n, p, df_X, scale_X = 1, ...) {
      hypercube <- lhs::randomLHS(n, p)
      X_raw <- stats::qt(hypercube, df = df_X)
      X_scaled <- X_raw * scale_X
      return(X_scaled)
    }
    
    one_run(
      p = param_grid$p[ii],
      n = param_grid$n[ii],
      rep_id = param_grid$rep_id[ii],
      c_K1 = param_grid$c_K1[ii]
    )
  }
} else {
  results <- bind_rows(lapply(seq_len(nrow(param_grid)), function(ii) {
    one_run(
      p = param_grid$p[ii],
      n = param_grid$n[ii],
      rep_id = param_grid$rep_id[ii],
      c_K1 = param_grid$c_K1[ii]
    )
  }))
}

summary_tbl <- results %>%
  group_by(p, n, c_K1) %>%
  summarize(
    mean_l2 = mean(l2_error, na.rm = TRUE),
    median_l2 = median(l2_error, na.rm = TRUE),
    sd_l2 = sd(l2_error, na.rm = TRUE),
    q25_l2 = quantile(l2_error, 0.25, na.rm = TRUE),
    q75_l2 = quantile(l2_error, 0.75, na.rm = TRUE),
    eligible_rate = mean(ts_eligible, na.rm = TRUE),
    mean_K1 = mean(K1, na.rm = TRUE),
    mean_K2 = mean(K2, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(p, n, c_K1)

print(summary_tbl)

plot_ts_ck1 <- results %>%
  filter(is.finite(l2_error), !is.na(l2_error), l2_error > 0) %>%
  mutate(c_K1 = factor(c_K1)) %>%
  ggplot(aes(x = factor(n), y = l2_error, color = c_K1, group = c_K1)) +
  stat_summary(fun = median, geom = "line", linewidth = 1.0) +
  stat_summary(fun = median, geom = "point", size = 2.6) +
  stat_summary(
    fun.min = function(x) quantile(x, 0.25, na.rm = TRUE),
    fun.max = function(x) quantile(x, 0.75, na.rm = TRUE),
    geom = "errorbar",
    width = 0.14,
    linewidth = 0.7
  ) +
  scale_y_log10(labels = function(x) formatC(x, format = "fg", digits = 3)) +
  labs(
    x = "Sample Size, n",
    y = "TS-RIGHT L2 Error (log scale)",
    color = "c_K1"
  ) +
  theme_bw(base_size = 18) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

print(plot_ts_ck1)
