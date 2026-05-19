# ======================================================================
# Recommended TS-RIGHT logistic-regression comparison experiment.
# Parallel version with per-task seeds matching the sequential code.
#
# Run from the repository root.
# Required packages: glmnet, mvtnorm, MASS, dplyr, ggplot2, here,
# foreach, doParallel, lhs.
# ======================================================================

suppressPackageStartupMessages({
  library(glmnet)
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

# -----------------------------
# Final recommended parameters
# -----------------------------
P_GRID <- c(600)              # main experiment; optional sensitivity: c(400, 600)
N_GRID <- c(600,800,1200,1600,2000)  # avoid n=400 for TS with this schedule
N_REPS <- 100               # use 8--20 for pilot, 100--200 for paper figures

S_STAR <- 5
S_ALG <- 8
T1 <- 2000
T2 <- 2
Q_TS <- 0.20
ETA_TS <- 0.02
ETA_IHT <- 0.02
T_IHT <- 1500
M_MIN <- 10
C_K1 <- 0.75
C_K2 <- 0.50
DF_X <- 2.05
SCALE_X <- 4.0
TAU_Q_SHRINK <- 0.99

USE_PARALLEL <- TRUE
N_CORES <- max(1L, parallel::detectCores() - 1L)

# Optional smoke-test overrides from PowerShell, for example:
#   $env:REC_N_REPS='1'
#   $env:REC_N_GRID='800'
#   $env:REC_C_K1='0.25'
#   $env:REC_C_K2='0.5'
#   $env:REC_USE_PARALLEL='false'
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
C_K1 <- env_num_vec("REC_C_K1", C_K1)[1]
C_K2 <- env_num_vec("REC_C_K2", C_K2)[1]
USE_PARALLEL <- env_bool("REC_USE_PARALLEL", USE_PARALLEL)
N_CORES <- env_int("REC_N_CORES", N_CORES)

out_tag <- Sys.getenv("REC_OUT_TAG", unset = "recommended_parallel")
out_dir <- here("results", "09_logistic_ts_signal_pilot", out_tag)
raw_dir <- file.path(out_dir, "raw")
summary_dir <- file.path(out_dir, "summary")
fig_dir <- file.path(out_dir, "figures")

dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

make_theta_recommended <- function(p) {
  theta <- rep(0, p)
  # np.random.default_rng(1).uniform(0.5, 1.5, 5); fixed for reproducibility.
  theta[seq_len(S_STAR)] <- c(
    1, -1, 1, -1, 1
  )
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

one_run <- function(p, n, rep_id) {
  theta_star <- make_theta_recommended(p)
  Sigma_t <- diag(p)
  
  # Keep the original data/split RNG mechanism:
  # the train split and cv.glmnet folds inherit the RNG stream after data generation.
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
  X_full <- sim_data$X
  y_full <- sim_data$y
  
  # Unified split for lasso/shrinkage tuning only.
  train_size <- floor(0.8 * n)
  train_idx <- sample(seq_len(n), size = train_size)
  X_train <- X_full[train_idx, , drop = FALSE]
  y_train <- y_full[train_idx]
  
  # TS-RIGHT. The seed formula matches the sequential code.
  out_ts <- fit_ts_right_logistic_ck(
    X = X_full,
    y = y_full,
    theta_star = theta_star,
    s = S_ALG,
    eta = ETA_TS,
    q = Q_TS,
    T1 = T1,
    T2 = T2,
    m_min = M_MIN,
    c_K1 = C_K1,
    c_K2 = C_K2,
    s_ref = S_STAR,
    seed = 2000000 + 1000 * n + rep_id
  )
  err_ts <- out_ts$l2_error
  ts_meta <- out_ts$meta
  
  # IHT baseline.
  theta_iht <- solver_iht(
    X_full,
    y_full,
    s = S_ALG,
    eta = ETA_IHT,
    T_max = T_IHT,
    grad_func = grad_logistic_regression
  )
  err_iht <- sqrt(sum((theta_iht - theta_star)^2))
  
  # Raw logistic lasso.
  fit_lasso_cv <- cv.glmnet(
    X_train,
    y_train,
    family = "binomial",
    standardize = FALSE,
    intercept = FALSE
  )
  lambda_lasso <- fit_lasso_cv$lambda.min
  fit_lasso <- glmnet(
    X_full,
    y_full,
    lambda = lambda_lasso,
    family = "binomial",
    standardize = FALSE,
    intercept = FALSE
  )
  theta_lasso <- as.vector(fit_lasso$beta)
  err_lasso <- sqrt(sum((theta_lasso - theta_star)^2))
  
  # Shrinkage + logistic lasso.
  tau_x <- as.numeric(quantile(abs(X_train), TAU_Q_SHRINK))
  X_train_shrunk <- truncate_operator(X_train, tau_x)
  fit_shrink_cv <- cv.glmnet(
    X_train_shrunk,
    y_train,
    family = "binomial",
    standardize = FALSE,
    intercept = FALSE
  )
  lambda_shrink <- fit_shrink_cv$lambda.min
  theta_shrink <- solver_shrinkage_method(
    X_full,
    y_full,
    tau_x = tau_x,
    tau_y = 999,
    lambda = lambda_shrink,
    family = "binomial"
  )
  err_shrink <- sqrt(sum((theta_shrink - theta_star)^2))
  
  data.frame(
    p = p,
    n = n,
    rep_id = rep_id,
    method = c("TS-RIGHT", "IHT", "Lasso", "Shrinkage"),
    l2_error = c(err_ts, err_iht, err_lasso, err_shrink),
    df_x = DF_X,
    scale_x = SCALE_X,
    q = Q_TS,
    T1 = T1,
    T2 = T2,
    c_K1 = C_K1,
    c_K2 = C_K2,
    K1 = c(if ("K1" %in% names(ts_meta)) ts_meta$K1[1] else NA_integer_, NA, NA, NA),
    K2 = c(if ("K2" %in% names(ts_meta)) ts_meta$K2[1] else NA_integer_, NA, NA, NA),
    K1_block_size = c(if ("K1_block_size" %in% names(ts_meta)) ts_meta$K1_block_size[1] else NA_real_, NA, NA, NA),
    K2_block_size = c(if ("K2_block_size" %in% names(ts_meta)) ts_meta$K2_block_size[1] else NA_real_, NA, NA, NA),
    ts_eligible = c(if ("eligible" %in% names(ts_meta)) ts_meta$eligible[1] else NA, NA, NA, NA),
    seed_data = 1000000 + 1000 * n + rep_id,
    seed_ts = 2000000 + 1000 * n + rep_id
  )
}

param_grid <- expand.grid(
  p = P_GRID,
  n = N_GRID,
  rep_id = seq_len(N_REPS),
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

message("Running recommended logistic TS-RIGHT comparison.")
message("Tasks: ", nrow(param_grid), "; methods per task: 4")
message("P_GRID = ", paste(P_GRID, collapse = ","))
message("N_GRID = ", paste(N_GRID, collapse = ","))
message("N_REPS = ", N_REPS)
message("C_K1 = ", C_K1, "; C_K2 = ", C_K2)
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
    .packages = c("glmnet", "mvtnorm", "MASS", "dplyr", "here", "lhs")
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
      rep_id = param_grid$rep_id[ii]
    )
  }
} else {
  results <- bind_rows(lapply(seq_len(nrow(param_grid)), function(ii) {
    one_run(
      p = param_grid$p[ii],
      n = param_grid$n[ii],
      rep_id = param_grid$rep_id[ii]
    )
  }))
}

summary_tbl <- results %>%
  group_by(p, n, method) %>%
  summarize(
    mean_l2 = mean(l2_error, na.rm = TRUE),
    median_l2 = median(l2_error, na.rm = TRUE),
    sd_l2 = sd(l2_error, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(p, n, median_l2)

print(summary_tbl)

write.csv(results, file.path(raw_dir, "recommended_logistic_ts_right_raw.csv"), row.names = FALSE)
write.csv(summary_tbl, file.path(summary_dir, "recommended_logistic_ts_right_summary.csv"), row.names = FALSE)
saveRDS(results, file.path(raw_dir, "recommended_logistic_ts_right_raw.rds"))
saveRDS(summary_tbl, file.path(summary_dir, "recommended_logistic_ts_right_summary.rds"))

# ----------------------------------------------------------------------
# Grouped quantile-line plot.
# ----------------------------------------------------------------------
raw_results_saved <- results
method_palette <- c(
  "TS-RIGHT" = "#332288",
  "IHT" = "#CC6677",
  #"Adaptive Huber" = "#3B82A1",
  "Shrinkage" = "#44AA99",
  "Lasso" = "#CAB54B"
 
)
method_levels_display <- names(method_palette)
n_levels_display <- sort(unique(raw_results_saved$n))

quantile_plot_raw <- raw_results_saved %>%
  filter(is.finite(l2_error), !is.na(l2_error), l2_error > 0) %>%
  mutate(
    method_display = dplyr::recode(method, "TS-RIGHT" = "TS-RIGHT"),
    method_display = factor(method_display, levels = method_levels_display),
    n_factor = factor(n, levels = n_levels_display)
  )

non_iht_l2 <- quantile_plot_raw %>%
  filter(method_display != "IHT") %>%
  pull(l2_error)

display_cap_l2 <- as.numeric(quantile(non_iht_l2, 0.99, na.rm = TRUE))
display_floor_l2 <- as.numeric(quantile(quantile_plot_raw$l2_error, 0.005, na.rm = TRUE))
display_floor_l2 <- max(display_floor_l2, 1e-8)

dodge_step <- 0.13
tick_half_width <- 0.045
method_mid <- (length(method_levels_display) + 1) / 2

grouped_quantiles <- quantile_plot_raw %>%
  group_by(n, n_factor, method_display) %>%
  summarise(
    q25_l2 = quantile(l2_error, 0.25, na.rm = TRUE),
    q50_l2 = quantile(l2_error, 0.50, na.rm = TRUE),
    q75_l2 = quantile(l2_error, 0.75, na.rm = TRUE),
    n_obs = dplyr::n(),
    .groups = "drop"
  ) %>%
  mutate(
    x_base = as.numeric(n_factor),
    method_index = as.numeric(method_display),
    x_pos = x_base + (method_index - method_mid) * dodge_step
  )

iht_extreme_marks <- quantile_plot_raw %>%
  filter(method_display == "IHT") %>%
  group_by(n, n_factor, method_display) %>%
  summarise(
    n_obs = dplyr::n(),
    n_above_cap = sum(l2_error > display_cap_l2, na.rm = TRUE),
    max_l2 = max(l2_error, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(n_above_cap > 0) %>%
  mutate(
    x_base = as.numeric(n_factor),
    method_index = as.numeric(method_display),
    x_pos = x_base + (method_index - method_mid) * dodge_step,
    y_pos = display_cap_l2,
    label = paste0(n_above_cap, "/", n_obs)
  )

write.csv(
  grouped_quantiles,
  file.path(summary_dir, "grouped_quantile_lines_recommended_logistic.csv"),
  row.names = FALSE
)
write.csv(
  iht_extreme_marks,
  file.path(summary_dir, "iht_extreme_counts_recommended_logistic.csv"),
  row.names = FALSE
)

plot_grouped_quantile_lines <- ggplot(grouped_quantiles, aes(color = method_display)) +
  geom_linerange(
    aes(x = x_pos, ymin = q25_l2, ymax = q75_l2),
    linewidth = 0.95
  ) +
  geom_segment(
    aes(x = x_pos - tick_half_width, xend = x_pos + tick_half_width, y = q25_l2, yend = q25_l2),
    linewidth = 0.65
  ) +
  geom_segment(
    aes(x = x_pos - tick_half_width, xend = x_pos + tick_half_width, y = q50_l2, yend = q50_l2),
    linewidth = 1.25
  ) +
  geom_segment(
    aes(x = x_pos - tick_half_width, xend = x_pos + tick_half_width, y = q75_l2, yend = q75_l2),
    linewidth = 0.65
  ) +
  geom_point(
    data = iht_extreme_marks,
    aes(x = x_pos, y = y_pos),
    inherit.aes = FALSE,
    shape = 24,
    size = 3.1,
    stroke = 0.4,
    color = method_palette["IHT"],
    fill = method_palette["IHT"]
  ) +
  geom_text(
    data = iht_extreme_marks,
    aes(x = x_pos, y = y_pos, label = label),
    inherit.aes = FALSE,
    vjust = 1.5,
    size = 4.8,
    color = method_palette["IHT"]
  ) +
  scale_x_continuous(
    breaks = seq_along(n_levels_display),
    labels = n_levels_display,
    expand = expansion(mult = c(0.05, 0.07))
  ) +
  scale_y_log10(labels = function(x) formatC(x, format = "fg", digits = 3)) +
  coord_cartesian(ylim = c(display_floor_l2, display_cap_l2), clip = "off") +
  scale_color_manual(
    values = method_palette,
    breaks = method_levels_display,
    drop = FALSE
  ) +
  labs(
    x = "Sample Size, n",
    y = "L2 Error (log scale)",
    color = "Method"
  ) +
  theme_bw(base_size = 18) +
  theme(
    axis.title = element_text(size = 19),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 15),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    strip.background = element_rect(fill = "gray90", color = "black"),
    strip.text = element_text(size = 15, face = "bold"),
    plot.margin = margin(t = 18, r = 18, b = 8, l = 8)
  )

print(plot_grouped_quantile_lines)

ggsave(
  file.path(fig_dir, "recommended_logistic_grouped_quantile_lines.pdf"),
  plot_grouped_quantile_lines,
  width = 11,
  height = 5.5
)
ggsave(
  file.path(fig_dir, "recommended_logistic_grouped_quantile_lines.png"),
  plot_grouped_quantile_lines,
  width = 11,
  height = 5.5,
  dpi = 300
)
saveRDS(
  plot_grouped_quantile_lines,
  file.path(fig_dir, "recommended_logistic_grouped_quantile_lines.rds")
)

message("Done. Results saved under: ", out_dir)
