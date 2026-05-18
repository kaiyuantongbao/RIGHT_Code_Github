# ======================================================================
# Logistic TS-RIGHT trajectory pilot at b = 0.55
#
# Purpose:
#   Inspect whether T1 = 200 is sensible before the full signal scan finishes.
#
# Default setup:
#   n = 1000, p = 800, s_star = 10, s = 20
#   df_X = 2.1, scale_X = sqrt((df_X - 2) / df_X)
#   theta_j^* = +/- 0.55 on the first s_star coordinates
#   q = 0.5, T1_max = 300, T2 = 8
#
# Main figure:
#   Median L2 trajectory across replications, with q25-q75 ribbon and
#   individual replication traces in the background.
# ======================================================================

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(ggplot2)
  library(foreach)
  library(doParallel)
})

source(here("simulations", "07_schedule_ablation", "utils_schedule_right.R"))

cfg <- list(
  experiment_id = "logistic_b055_ts_trajectory",

  n = 800,
  p = 600,
  s_star = 5,
  s = 10,
  m_min = 10,
  c_K1 = 0.25,
  c_K2 = 0.50,

  reps = 30,
  seed_base = 20260508L,

  df_X = 2.05,
  scale_X = 4,
  theta_magnitude = 1,

  q = 0.20,
  T1_max = 2200,
  T1_reference = 1500,
  T2 = 4,
  eta_ts = 0.02,
  record_every = 5L,

  use_parallel = TRUE,
  n_cores = max(1L, parallel::detectCores() - 1L)
)

env_int <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) return(default)
  as.integer(value)
}

env_num <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) return(default)
  as.numeric(value)
}

env_bool <- function(name, default) {
  value <- tolower(Sys.getenv(name, unset = ""))
  if (!nzchar(value)) return(default)
  value %in% c("1", "true", "t", "yes", "y")
}

# Optional smoke-test overrides.
cfg$reps <- env_int("TRAJ_REPS", cfg$reps)
cfg$n <- env_int("TRAJ_N", cfg$n)
cfg$theta_magnitude <- env_num("TRAJ_B", cfg$theta_magnitude)
cfg$m_min <- env_int("TRAJ_M_MIN", cfg$m_min)
cfg$c_K1 <- env_num("TRAJ_C_K1", cfg$c_K1)
cfg$c_K2 <- env_num("TRAJ_C_K2", cfg$c_K2)
cfg$T1_max <- env_int("TRAJ_T1_MAX", cfg$T1_max)
cfg$T1_reference <- env_int("TRAJ_T1_REFERENCE", cfg$T1_reference)
cfg$T2 <- env_int("TRAJ_T2", cfg$T2)
cfg$eta_ts <- env_num("TRAJ_ETA_TS", cfg$eta_ts)
cfg$record_every <- env_int("TRAJ_RECORD_EVERY", cfg$record_every)
cfg$use_parallel <- env_bool("TRAJ_USE_PARALLEL", cfg$use_parallel)
cfg$n_cores <- env_int("TRAJ_N_CORES", cfg$n_cores)

out_dir <- here("results", "09_logistic_ts_signal_pilot", "trajectory_b055")
raw_dir <- file.path(out_dir, "raw")
summary_dir <- file.path(out_dir, "summary")
fig_dir <- file.path(out_dir, "figures")

dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

make_theta_star_signed <- function(p, s_star, magnitude) {
  theta_star <- rep(0, p)
  signs <- rep(c(1, -1), length.out = s_star)
  theta_star[seq_len(s_star)] <- signs * magnitude
  theta_star
}

seed_for <- function(base, rep_id, offset = 0L) {
  as.integer(base + offset + rep_id)
}

eta_diagnostics <- function(X, theta_star) {
  eta <- as.numeric(X %*% theta_star)
  prob <- plogis(eta)
  data.frame(
    eta_sd = sd(eta),
    eta_median_abs = median(abs(eta)),
    eta_q90_abs = as.numeric(quantile(abs(eta), 0.90, na.rm = TRUE)),
    eta_gt3_rate = mean(abs(eta) > 3),
    mean_bernoulli_var = mean(prob * (1 - prob))
  )
}

run_ts_right_logistic_ck <- function(
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
  record_every = 1L,
  record_initial = TRUE,
  seed = NULL
) {
  n <- nrow(X)
  p <- ncol(X)
  record_every <- as.integer(record_every)

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

  base <- data.frame(
    arm = "TS",
    model = "logistic",
    q = q,
    frac_stage1 = budget$frac_stage1,
    frac_stage2 = budget$frac_stage2,
    T1 = T1,
    T2 = T2,
    n1 = budget$n1,
    n2 = budget$n2,
    b2 = budget$b2,
    s_ref = s_ref,
    m_min = m_min,
    c_K1 = c_K1,
    c_K2 = c_K2,
    K1_tar = budget$K1_tar,
    K1 = budget$K1,
    K2_tar = budget$K2_tar,
    K2 = budget$K2,
    K1_block_size = budget$n1 / budget$K1,
    K2_block_size = budget$b2 / budget$K2,
    cost_proxy = budget$cost_proxy,
    T_fd_budget = budget$T_fd_budget,
    eligible = budget$eligible
  )

  if (!isTRUE(budget$eligible)) {
    base$t <- NA_integer_
    base$stage <- NA_character_
    base$sample_budget_used <- NA_real_
    base$fd_equiv_iter <- NA_real_
    base$l2_error <- NA_real_
    return(base)
  }

  rng <- with_local_seed(seed)
  on.exit(restore_local_seed(rng), add = TRUE)

  idx <- sample.int(n)
  id1 <- idx[seq_len(budget$n1)]
  id2 <- idx[(budget$n1 + 1L):n]

  trace_t <- integer(0)
  trace_l2 <- numeric(0)

  if (record_initial) {
    trace_t <- c(trace_t, 0L)
    trace_l2 <- c(trace_l2, safe_l2_error(rep(0, p), theta_star))
  }

  X1 <- X[id1, , drop = FALSE]
  y1 <- y[id1]
  theta_cur <- rep(0, p)

  for (tt in seq_len(T1)) {
    theta_cur <- solver_right(
      X = X1,
      y = y1,
      s = s,
      eta = eta,
      T_max = 1L,
      K = budget$K1,
      theta_init = theta_cur,
      grad_func_samplewise = grad_logistic_regression_samplewise,
      record_trace = FALSE
    )

    if (tt %% record_every == 0L || tt == T1) {
      trace_t <- c(trace_t, tt)
      trace_l2 <- c(trace_l2, safe_l2_error(theta_cur, theta_star))
    }
  }

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

    trace_t <- c(trace_t, T1 + tt)
    trace_l2 <- c(trace_l2, safe_l2_error(theta_cur, theta_star))
  }

  budget_axis <- add_budget_axis_to_trace(
    trace_t = trace_t,
    n = n,
    arm = "TS",
    n1 = budget$n1,
    T1 = T1,
    b2 = budget$b2
  )

  out <- base[rep(1, length(trace_t)), ]
  out$t <- trace_t
  out$stage <- budget_axis$stage
  out$sample_budget_used <- budget_axis$sample_budget_used
  out$fd_equiv_iter <- budget_axis$fd_equiv_iter
  out$l2_error <- trace_l2
  out
}

fit_one_trajectory <- function(rep_id, cfg, theta_star) {
  seed_data <- seed_for(cfg$seed_base, rep_id, offset = 1000000L)
  seed_algo <- seed_for(cfg$seed_base, rep_id, offset = 2000000L)

  dat <- gen_one_logistic_dataset(
    seed = seed_data,
    n = cfg$n,
    p = cfg$p,
    theta_star = theta_star,
    df_X = cfg$df_X,
    Sigma_X = diag(cfg$p),
    scale_X = cfg$scale_X
  )

  X <- dat$X
  y <- dat$y

  traj <- run_ts_right_logistic_ck(
    X = X,
    y = y,
    theta_star = theta_star,
    s = cfg$s,
    eta = cfg$eta_ts,
    q = cfg$q,
    T1 = cfg$T1_max,
    T2 = cfg$T2,
    m_min = cfg$m_min,
    c_K1 = cfg$c_K1,
    c_K2 = cfg$c_K2,
    s_ref = cfg$s_star,
    record_initial = TRUE,
    record_every = cfg$record_every,
    seed = seed_algo
  )

  diag <- eta_diagnostics(X, theta_star)

  dplyr::bind_cols(
    data.frame(
      rep_id = rep_id,
      seed_data = seed_data,
      seed_algo = seed_algo
    ),
    traj,
    diag[rep(1, nrow(traj)), ]
  )
}

theta_star <- make_theta_star_signed(
  p = cfg$p,
  s_star = cfg$s_star,
  magnitude = cfg$theta_magnitude
)

message("Running logistic TS-RIGHT trajectory pilot.")
message("n = ", cfg$n, ", p = ", cfg$p, ", b = ", cfg$theta_magnitude)
message("q = ", cfg$q, ", T1_max = ", cfg$T1_max, ", T2 = ", cfg$T2)
message("m_min = ", cfg$m_min, ", c_K1 = ", cfg$c_K1, ", c_K2 = ", cfg$c_K2)
message("reps = ", cfg$reps, ", record_every = ", cfg$record_every)

param_grid <- data.frame(rep_id = seq_len(cfg$reps))
helper_exports <- c(
  "cfg", "theta_star", "make_theta_star_signed", "seed_for",
  "eta_diagnostics", "run_ts_right_logistic_ck", "fit_one_trajectory"
)

if (cfg$use_parallel && cfg$n_cores > 1L) {
  cl <- parallel::makeCluster(cfg$n_cores)
  doParallel::registerDoParallel(cl)
  on.exit({
    try(parallel::stopCluster(cl), silent = TRUE)
  }, add = TRUE)

  raw_trace <- foreach(
    ii = seq_len(nrow(param_grid)),
    .combine = dplyr::bind_rows,
    .packages = c("dplyr", "here"),
    .export = helper_exports
  ) %dopar% {
    source(here("simulations", "07_schedule_ablation", "utils_schedule_right.R"))
    fit_one_trajectory(rep_id = param_grid$rep_id[ii], cfg = cfg, theta_star = theta_star)
  }
} else {
  raw_trace <- lapply(param_grid$rep_id, function(rep_id) {
    fit_one_trajectory(rep_id = rep_id, cfg = cfg, theta_star = theta_star)
  }) %>% dplyr::bind_rows()
}

saveRDS(raw_trace, file.path(raw_dir, "raw_logistic_b055_ts_trajectory.rds"))
write.csv(raw_trace, file.path(raw_dir, "raw_logistic_b055_ts_trajectory.csv"), row.names = FALSE)

summary_trace <- raw_trace %>%
  filter(eligible, is.finite(l2_error), !is.na(t)) %>%
  group_by(t, stage) %>%
  summarise(
    median_l2 = median(l2_error, na.rm = TRUE),
    mean_l2 = mean(l2_error, na.rm = TRUE),
    q25_l2 = quantile(l2_error, 0.25, na.rm = TRUE),
    q75_l2 = quantile(l2_error, 0.75, na.rm = TRUE),
    median_fd_equiv_iter = median(fd_equiv_iter, na.rm = TRUE),
    K1 = first(K1),
    K2 = first(K2),
    K1_block_size = first(K1_block_size),
    K2_block_size = first(K2_block_size),
    .groups = "drop"
  )

write.csv(summary_trace, file.path(summary_dir, "summary_logistic_b055_ts_trajectory.csv"), row.names = FALSE)
saveRDS(summary_trace, file.path(summary_dir, "summary_logistic_b055_ts_trajectory.rds"))

eligible_rate <- raw_trace %>%
  group_by(rep_id) %>%
  summarise(eligible = any(eligible), .groups = "drop") %>%
  summarise(eligible_rate = mean(eligible)) %>%
  pull(eligible_rate)

diagnostics <- raw_trace %>%
  group_by(rep_id) %>%
  summarise(
    eta_sd = first(eta_sd),
    eta_median_abs = first(eta_median_abs),
    eta_gt3_rate = first(eta_gt3_rate),
    mean_bernoulli_var = first(mean_bernoulli_var),
    K1 = first(K1),
    K2 = first(K2),
    K1_block_size = first(K1_block_size),
    K2_block_size = first(K2_block_size),
    .groups = "drop"
  ) %>%
  summarise(
    reps = dplyr::n(),
    eligible_rate = eligible_rate,
    median_eta_sd = median(eta_sd, na.rm = TRUE),
    median_eta_abs = median(eta_median_abs, na.rm = TRUE),
    median_eta_gt3_rate = median(eta_gt3_rate, na.rm = TRUE),
    median_bernoulli_var = median(mean_bernoulli_var, na.rm = TRUE),
    K1 = first(K1),
    K2 = first(K2),
    K1_block_size = first(K1_block_size),
    K2_block_size = first(K2_block_size)
  )

write.csv(diagnostics, file.path(summary_dir, "diagnostics_logistic_b055_ts_trajectory.csv"), row.names = FALSE)

plot_trace <- ggplot() +
  geom_line(
    data = raw_trace %>% filter(eligible, is.finite(l2_error), !is.na(t)),
    aes(x = t, y = l2_error, group = rep_id),
    alpha = 0.18,
    linewidth = 0.35,
    color = "gray45"
  ) +
  geom_ribbon(
    data = summary_trace,
    aes(x = t, ymin = q25_l2, ymax = q75_l2),
    fill = "#7AA6C2",
    alpha = 0.22
  ) +
  geom_line(
    data = summary_trace,
    aes(x = t, y = median_l2),
    linewidth = 1.0,
    color = "#1F5A7A"
  ) +
  geom_vline(xintercept = cfg$T1_reference, linetype = "dashed", color = "#C44E52", linewidth = 0.7) +
  geom_vline(xintercept = cfg$T1_max, linetype = "dotted", color = "gray30", linewidth = 0.7) +
  annotate(
    "text",
    x = cfg$T1_reference,
    y = max(summary_trace$q75_l2, na.rm = TRUE),
    label = paste0("T1 = ", cfg$T1_reference),
    hjust = -0.05,
    vjust = 1.2,
    size = 3.4,
    color = "#C44E52"
  ) +
  annotate(
    "text",
    x = cfg$T1_max,
    y = max(summary_trace$q75_l2, na.rm = TRUE),
    label = "Stage II starts",
    hjust = 1.05,
    vjust = 2.6,
    size = 3.4,
    color = "gray25"
  ) +
  scale_y_log10() +
  labs(
    x = "TS iteration",
    y = "L2 Error (log scale)",
    title = "Logistic TS-RIGHT Trajectory",
    subtitle = sprintf(
      "n=%d, p=%d, b=%.2f, q=%.2f, T2=%d, m_min=%d, cK1=%.2f, cK2=%.2f; median/IQR over %d reps",
      cfg$n, cfg$p, cfg$theta_magnitude, cfg$q, cfg$T2,
      cfg$m_min, cfg$c_K1, cfg$c_K2, cfg$reps
    )
  ) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  )

print(plot_trace)
ggsave(file.path(fig_dir, "logistic_b055_ts_trajectory_l2.pdf"), plot_trace, width = 9, height = 5.5)
ggsave(file.path(fig_dir, "logistic_b055_ts_trajectory_l2.png"), plot_trace, width = 9, height = 5.5, dpi = 300)

message("Done. Results saved under: ", out_dir)
