rm(list = ls())

source("./simulations/07_schedule_ablation/utils_schedule_linear.R")
 
out_dir <- "./results/07_schedule_ablation/E1_linear_TS_plateau_pilot"
dir.create(file.path(out_dir, "raw"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "summary"), recursive = TRUE, showWarnings = FALSE)

cfg <- list(
  n = 2000,
  p = 600,
  s_star = 10,
  s = 20,
  eta = 0.01,
  c_r = 0.5,
  T1_grid = c(128, 192, 256, 384),
  T2_grid = c(8, 12),
  m_min = 10,
  reps = 40,
  seed_base = 20260417L,
  df_X = 2.5,
  scale_X = 1,
  df_eps = 1.5
)

theta_star <- make_theta_star(cfg$p, cfg$s_star, magnitude = 5)
grid <- expand.grid(T1 = cfg$T1_grid, T2 = cfg$T2_grid, KEEP.OUT.ATTRS = FALSE)

num_cores <- max(1, parallel::detectCores() - 1)
cl <- parallel::makeCluster(num_cores)
doParallel::registerDoParallel(cl)

raw_results <- foreach(
  rep_id = 1:cfg$reps,
  .combine = dplyr::bind_rows,
  .packages = c("dplyr", "MASS", "mvtnorm")
) %dopar% {
  source("simulations/07_schedule_ablation/utils_schedule_linear.R")
  
  dat <- gen_one_linear_dataset(
    seed = cfg$seed_base + rep_id,
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
  
  ans <- lapply(seq_len(nrow(grid)), function(j) {
    pars <- grid[j, ]
    out <- run_ts_right_linear(
      X = X, y = y, theta_star = theta_star,
      s = cfg$s, eta = cfg$eta,
      c_r = cfg$c_r, T1 = pars$T1, T2 = pars$T2,
      m_min = cfg$m_min
    )
    cbind(data.frame(rep_id = rep_id), out)
  })
  
  bind_rows(ans)
}

parallel::stopCluster(cl)

saveRDS(raw_results, file.path(out_dir, "raw", "raw_TS_plateau_pilot.rds"))

summary_tbl <- raw_results %>%
  group_by(T1, T2) %>%
  summarise(
    n1 = first(n1),
    n2 = first(n2),
    b2 = first(b2),
    K1 = first(K1),
    K2 = first(K2),
    cost_proxy = first(cost_proxy),
    median_l2 = median(l2_error, na.rm = TRUE),
    mean_l2 = mean(l2_error, na.rm = TRUE),
    iqr_l2 = IQR(l2_error, na.rm = TRUE),
    median_runtime = median(runtime_sec, na.rm = TRUE),
    eligible = all(eligible),
    .groups = "drop"
  ) %>%
  arrange(median_l2, cost_proxy)

best_med <- min(summary_tbl$median_l2[summary_tbl$eligible], na.rm = TRUE)

plateau_pick <- summary_tbl %>%
  filter(eligible, median_l2 <= 1.05 * best_med) %>%
  arrange(cost_proxy, T1, T2) %>%
  slice(1)

write.csv(summary_tbl, file.path(out_dir, "summary", "summary_TS_plateau_pilot.csv"), row.names = FALSE)
write.csv(plateau_pick, file.path(out_dir, "summary", "default_TS_pick.csv"), row.names = FALSE)

print(summary_tbl)
print(plateau_pick)
