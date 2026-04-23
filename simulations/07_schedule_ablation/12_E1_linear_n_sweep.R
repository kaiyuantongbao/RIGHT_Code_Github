rm(list = ls())
library(here)
source(here("simulations","07_schedule_ablation","utils_schedule_linear.R"))

out_dir <- "results/07_schedule_ablation/E1_linear_n_sweep"
dir.create(file.path(out_dir, "raw"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "summary"), recursive = TRUE, showWarnings = FALSE)

cfg <- list(
  n_grid = c(1200, 2000, 2800,3600,4400),
  p = 600,
  s_star = 10,
  s = 20,
  eta = 0.01,
  m_min = 10,
  reps = 200,
  seed_base = 20260419L,
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
    .packages = c("dplyr", "MASS", "mvtnorm","here")
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
    
    # TS
    ts_out <- run_ts_right_linear(
      X = X, y = y, theta_star = theta_star,
      s = cfg$s, eta = cfg$eta,
      c_r = TS_cfg$c_r, T1 = TS_cfg$T1, T2 = TS_cfg$T2,
      m_min = cfg$m_min
    )
    
    # FD budget-matched
    T_fd_budget <- unique(ts_out$T_fd_budget)
    fd_out <- run_fd_right_linear(
      X = X, y = y, theta_star = theta_star,
      s = cfg$s, eta = cfg$eta,
      T_fd = T_fd_budget, m_min = cfg$m_min
    )
    
    # FS default m_min*lop=70
    T_fs <- floor(n / 70)
    fs_out <- run_fs_right_linear(
      X = X, y = y, theta_star = theta_star,
      s = cfg$s, eta = cfg$eta,
      T_fs = T_fs, m_min = cfg$m_min
    )
    
    bind_rows(
      cbind(data.frame(n = n, rep_id = rep_id, method = "TS"), ts_out),
      cbind(data.frame(n = n, rep_id = rep_id, method = "FD_budget"), fd_out),
      cbind(data.frame(n = n, rep_id = rep_id, method = "FS"), fs_out)
    )
  }

parallel::stopCluster(cl)

saveRDS(raw_results, file.path(out_dir, "raw", "raw_E1_linear_n_sweep.rds"))

summary_tbl <- raw_results %>%
  group_by(n, method) %>%
  summarise(
    median_l2 = median(l2_error, na.rm = TRUE),
    mean_l2 = mean(l2_error, na.rm = TRUE),
    iqr_l2 = IQR(l2_error, na.rm = TRUE),
    median_runtime = median(runtime_sec, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(n, method)

write.csv(summary_tbl, file.path(out_dir, "summary", "summary_E1_linear_n_sweep.csv"), row.names = FALSE)
print(summary_tbl)

library(ggplot2)
library(dplyr)
p<-ggplot(data=summary_tbl,mapping=aes(x=log(n),y=log(median_l2),color=method))+
  geom_line()+
  labs(
    x="Log Sample Size",
    y="Log Median L2 Error"
  )
p

ggsave(file.path(out_dir,"method_vs_n.pdf"),plot=p,height=7,width=6.5)
