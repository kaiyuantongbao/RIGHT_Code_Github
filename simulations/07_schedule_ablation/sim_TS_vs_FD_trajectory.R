rm(list = ls())

library(here)
library(dplyr)
library(ggplot2)
library(doParallel)
library(foreach)
library(grid)

# Use the seed-enabled two-stage schedule wrappers.
# This file assumes utils_schedule_right.R is the updated version with q-split and seed control.
source(here("simulations", "07_schedule_ablation", "utils_schedule_right.R"))

out_dir <- "results/07_schedule_ablation/TS_vs_FD_trajectory"
dir.create(file.path(out_dir, "raw"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "summary"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "figures"), recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# 1. Experiment configuration
# ==============================================================================

cfg <- list(
  n = 3200,
  p = 600,
  s_star = 10,
  s = 20,
  eta = 0.02,
  
  # Fixed-fraction TS split: n2 = floor(q * n), n1 = n - n2.
  # For this trajectory figure we vary T1 and keep q,T2 fixed.
  # Later tuning scripts will use q_grid = c(0.1, 0.2, 0.25, 0.5).
  q_grid = c(0.25),
  T1_grid = c(128, 150, 200),
  T2_grid = c(8),
  m_min = 10,
  reps = 100,
  seed_base = 20260417L,
  
  df_X = 2.5,
  scale_X = 1,
  df_eps = 1.5,
  theta_magnitude = 5
)

theta_star <- make_theta_star(cfg$p, cfg$s_star, magnitude = cfg$theta_magnitude)

schedule_grid <- expand.grid(
  q = cfg$q_grid,
  T1 = cfg$T1_grid,
  T2 = cfg$T2_grid,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
) %>%
  arrange(q, T1, T2) %>%
  mutate(
    config_id = row_number(),
    Config = sprintf("TS (q=%.2f, T1=%d, T2=%d)", q, T1, T2)
  )

write.csv(
  schedule_grid,
  file.path(out_dir, "summary", "trajectory_schedule_grid.csv"),
  row.names = FALSE
)

num_cores <- max(1, parallel::detectCores() - 1)
cl <- parallel::makeCluster(num_cores)
doParallel::registerDoParallel(cl)

# ==============================================================================
# 2. Parallel Monte Carlo loop
# ==============================================================================

raw_results <- foreach(
  rep_id = 1:cfg$reps,
  .combine = dplyr::bind_rows,
  .export = c("cfg", "theta_star", "schedule_grid"),
  .packages = c("dplyr", "MASS", "mvtnorm", "here")
) %dopar% {
  source(here("simulations", "07_schedule_ablation", "utils_schedule_right.R"))
  
  seed_data <- cfg$seed_base + rep_id
  seed_algo_ts <- cfg$seed_base + 100000L + rep_id
  
  # Same generated data are used by all methods in this replication.
  dat <- gen_one_linear_dataset(
    seed = seed_data,
    n = cfg$n,
    p = cfg$p,
    theta_star = theta_star,
    df_X = cfg$df_X,
    Sigma_X = diag(cfg$p),
    scale_X = cfg$scale_X,
    df_eps = cfg$df_eps
  )
  
  X <- dat$X
  y <- dat$y
  
  # TS configurations.
  # The same seed_algo_ts is intentionally used across the TS configs in the
  # same replication. When q is the same, this gives the same data split and
  # isolates the effect of T1. If q differs, it still uses the same base
  # permutation before taking different split sizes.
  ans_ts <- lapply(seq_len(nrow(schedule_grid)), function(ii) {
    sc <- schedule_grid[ii, ]
    
    out <- run_ts_right_linear(
      X = X,
      y = y,
      theta_star = theta_star,
      s = cfg$s,
      eta = cfg$eta,
      q = sc$q,
      T1 = sc$T1,
      T2 = sc$T2,
      m_min = cfg$m_min,
      record_trace = TRUE,
      record_l2 = TRUE,
      record_initial = TRUE,
      seed = seed_algo_ts
    )
    
    out$Config <- sc$Config
    out$arm_type <- "TS"
    out$config_id <- sc$config_id
    out$seed_data <- seed_data
    out$seed_algo <- seed_algo_ts
    out
  })
  
  df_ts <- bind_rows(ans_ts)
  
  # FD baseline. It is run long enough to cover both the raw TS iteration axis
  # and the budget-matched FD-equivalent axis.
  max_budget <- max(df_ts$T_fd_budget, na.rm = TRUE)
  max_ts_iter <- max(schedule_grid$T1 + schedule_grid$T2)
  T_fd_run <- ceiling(max(max_ts_iter, max_budget) * 1.1)
  
  out_fd <- run_fd_right_linear(
    X = X,
    y = y,
    theta_star = theta_star,
    s = cfg$s,
    eta = cfg$eta,
    T_fd = T_fd_run,
    m_min = cfg$m_min,
    record_trace = TRUE,
    record_l2 = TRUE,
    record_initial = TRUE
  )
  
  out_fd$Config <- "FD (Baseline)"
  out_fd$arm_type <- "FD"
  out_fd$config_id <- 0L
  out_fd$seed_data <- seed_data
  out_fd$seed_algo <- NA_integer_
  
  bind_rows(df_ts, out_fd) %>%
    mutate(rep_id = rep_id)
}

parallel::stopCluster(cl)

saveRDS(
  raw_results,
  file.path(out_dir, "raw", "raw_TS_vs_FD_trajectory.rds")
)

# ==============================================================================
# 3. Summary and visualization
# ==============================================================================

trace_summary <- raw_results %>%
  filter((arm_type == "FD") | (arm_type == "TS" & eligible == TRUE)) %>%
  group_by(Config, arm_type, t) %>%
  summarise(
    median_l2 = median(l2_error, na.rm = TRUE),
    q25_l2 = quantile(l2_error, 0.25, na.rm = TRUE),
    q75_l2 = quantile(l2_error, 0.75, na.rm = TRUE),
    n_rep = dplyr::n_distinct(rep_id),
    .groups = "drop"
  )

budget_lines <- raw_results %>%
  filter(arm_type == "TS", eligible == TRUE) %>%
  group_by(Config) %>%
  summarise(
    T_budget = first(T_fd_budget),
    fd_equiv_final = max(fd_equiv_iter, na.rm = TRUE),
    q = first(q),
    T1 = first(T1),
    T2 = first(T2),
    n1 = first(n1),
    n2 = first(n2),
    b2 = first(b2),
    K1 = first(K1),
    K2 = first(K2),
    .groups = "drop"
  )

write.csv(
  trace_summary,
  file.path(out_dir, "summary", "trace_summary_median_l2.csv"),
  row.names = FALSE
)

write.csv(
  budget_lines,
  file.path(out_dir, "summary", "budget_lines.csv"),
  row.names = FALSE
)

# Budget-matched diagnostic table: TS final point vs FD at the corresponding
# dotted vertical budget line. This is only a diagnostic summary; the main plot
# still uses the original iteration axis.
ts_final_by_rep <- raw_results %>%
  filter(arm_type == "TS", eligible == TRUE) %>%
  group_by(rep_id, Config) %>%
  slice_max(order_by = t, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  transmute(
    rep_id,
    Config,
    method = "TS final",
    eval_t = t,
    l2_error
  )

fd_at_budget_by_rep <- raw_results %>%
  filter(arm_type == "TS", eligible == TRUE) %>%
  group_by(rep_id, Config) %>%
  summarise(T_budget = first(T_fd_budget), .groups = "drop") %>%
  left_join(
    raw_results %>%
      filter(arm_type == "FD") %>%
      select(rep_id, fd_t = t, fd_l2_error = l2_error),
    by = "rep_id"
  ) %>%
  filter(fd_t == T_budget) %>%
  transmute(
    rep_id,
    Config,
    method = "FD at matched budget",
    eval_t = fd_t,
    l2_error = fd_l2_error
  )

budget_compare_summary <- bind_rows(ts_final_by_rep, fd_at_budget_by_rep) %>%
  group_by(Config, method) %>%
  summarise(
    eval_t = first(eval_t),
    median_l2 = median(l2_error, na.rm = TRUE),
    q25_l2 = quantile(l2_error, 0.25, na.rm = TRUE),
    q75_l2 = quantile(l2_error, 0.75, na.rm = TRUE),
    n_rep = dplyr::n_distinct(rep_id),
    .groups = "drop"
  )

write.csv(
  budget_compare_summary,
  file.path(out_dir, "summary", "budget_matched_TS_final_vs_FD.csv"),
  row.names = FALSE
)

all_configs <- unique(trace_summary$Config)
ts_configs <- sort(all_configs[all_configs != "FD (Baseline)"])

my_colors <- c("FD (Baseline)" = "grey50")
my_linetypes <- c("FD (Baseline)" = "dashed")
my_linewidths <- c("FD (Baseline)" = 1.2)

palette_colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
                    "#A65628", "#F781BF", "#999999", "#66C2A5", "#FC8D62")

for (i in seq_along(ts_configs)) {
  color_idx <- ((i - 1) %% length(palette_colors)) + 1
  my_colors[ts_configs[i]] <- palette_colors[color_idx]
  my_linetypes[ts_configs[i]] <- "solid"
  my_linewidths[ts_configs[i]] <- 0.8
}

y_max_label <- max(trace_summary$median_l2, na.rm = TRUE)

p_compare_without_title <- ggplot(
  trace_summary,
  aes(x = t, y = median_l2, color = Config, linetype = Config, linewidth = Config)
) +
  geom_line() +
  geom_vline(
    data = budget_lines,
    aes(xintercept = T_budget, color = Config),
    linetype = "dotted",
    linewidth = 1
  ) +
  geom_text(
    data = budget_lines,
    aes(x = T_budget, y = y_max_label, label = paste0("Budget=", T_budget)),
    angle = 90,
    vjust = -0.5,
    hjust = 1,
    size = 3.5,
    show.legend = FALSE
  ) +
  theme_bw() +
  scale_y_log10() +
  scale_color_manual(values = my_colors) +
  scale_linetype_manual(values = my_linetypes) +
  scale_linewidth_manual(values = my_linewidths) +
  labs(
    #title = "Trajectory Comparison: Two-Stage RIGHT vs Full-Data RIGHT",
    #subtitle = "Dotted vertical lines indicate the FD-equivalent sample-access budget for each TS configuration.",
    x = expression(Iteration~Step~(t)),
    y = expression(Median~L[2]~Error),
    color = "Algorithm Configuration",
    linetype = "Algorithm Configuration",
    linewidth = "Algorithm Configuration"
  ) +
  theme(
    legend.position = "bottom",
    legend.key.width = unit(3, "lines")
  )

print(p_compare_without_title)

ggsave(
  file.path(out_dir, "figures", "TS_vs_FD_trajectory_plot.pdf"),
  plot = p_compare,
  width = 10.5,
  height = 6
)


ggsave(
  file.path(out_dir, "figures", "TS_vs_FD_trajectory_plot(no_title).pdf"),
  plot = p_compare_without_title,
  width = 10.5,
  height = 6
)
