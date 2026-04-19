rm(list = ls())

library(here)
source(here("simulations","07_schedule_ablation","utils_schedule_linear.R"))

out_dir <- "results/07_schedule_ablation/E1_linear_FD_budget_pilot"
dir.create(file.path(out_dir, "raw"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "summary"), recursive = TRUE, showWarnings = FALSE)

cfg <- list(
  n_grid = c(1200, 2000, 2800),
  p = 600,
  s_star = 10,
  s = 20,
  eta = 0.01,
  m_min = 10,
  reps = 40,
  seed_base = 20260418L,
  df_X = 2.5,
  scale_X = 1,
  df_eps = 1.5
)

# -------- manually edit this after TS plateau pilot --------
TS_cfg <- list(
  c_r = 0.5,
  T1 = 384,
  T2 = 8
)

theta_star <- make_theta_star(cfg$p, cfg$s_star, magnitude = 5)

num_cores <- max(1, parallel::detectCores() - 1)
cl <- parallel::makeCluster(num_cores)
doParallel::registerDoParallel(cl)

raw_results <- foreach(
  n = cfg$n_grid,
  .combine = dplyr::bind_rows,
  .packages = c("dplyr", "MASS", "mvtnorm","here")
) %:%
  foreach(
    rep_id = 1:cfg$reps,
    .combine = dplyr::bind_rows,
    .packages = c("dplyr", "MASS", "mvtnorm")
  ) %dopar% {
    source(here("simulations","07_schedule_ablation","utils_schedule_linear.R"))
    
    dat <- gen_one_linear_dataset(
      seed = cfg$seed_base + 100000L * match(n, cfg$n_grid) + rep_id,
      n = n,
      p = cfg$p,
      theta_star = theta_star,
      df_X = cfg$df_X,
      Sigma_X = diag(cfg$p),
      scale_X = cfg$scale_X,
      df_eps = cfg$df_eps
    )
    
    X <- dat$X
    y <- dat$y
    
    budget <- compute_ts_budget(
      n = n, s_star = cfg$s_star, p = cfg$p,
      c_r = TS_cfg$c_r, T1 = TS_cfg$T1, T2 = TS_cfg$T2,
      m_min = cfg$m_min
    )
    
    T_budget <- budget$T_fd_budget
    T_fd_grid <- sort(unique(pmax(8, round(c(0.75, 1.00, 1.25, 1.50) * T_budget))))
    
    ans <- lapply(T_fd_grid, function(T_fd) {
      out <- run_fd_right_linear(
        X = X, y = y, theta_star = theta_star,
        s = cfg$s, eta = cfg$eta,
        T_fd = T_fd, m_min = cfg$m_min
      )
      cbind(
        data.frame(
          n = n, rep_id = rep_id,
          T_fd_budget = T_budget,
          c_r = TS_cfg$c_r, T1 = TS_cfg$T1, T2 = TS_cfg$T2
        ),
        out
      )
    })
    
    bind_rows(ans)
  }

parallel::stopCluster(cl)

saveRDS(raw_results, file.path(out_dir, "raw", "raw_FD_budget_pilot.rds"))

summary_tbl <- raw_results %>%
  group_by(n, T_fd_budget, T_fd) %>%
  summarise(
    K_fd = first(K_fd),
    median_l2 = median(l2_error, na.rm = TRUE),
    mean_l2 = mean(l2_error, na.rm = TRUE),
    iqr_l2 = IQR(l2_error, na.rm = TRUE),
    median_runtime = median(runtime_sec, na.rm = TRUE),
    ratio_to_budget = first(T_fd / T_fd_budget),
    .groups = "drop"
  ) %>%
  arrange(n, median_l2, T_fd)

write.csv(summary_tbl, file.path(out_dir, "summary", "summary_FD_budget_pilot.csv"), row.names = FALSE)
print(summary_tbl)
