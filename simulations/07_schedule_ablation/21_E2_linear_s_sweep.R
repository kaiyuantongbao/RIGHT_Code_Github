rm(list = ls())
library(here)
source(here("simulations","07_schedule_ablation","utils_schedule_linear.R"))

out_dir <- "results/07_schedule_ablation/E2_linear_s_sweep"
#dir.create(file.path(out_dir, "raw"), recursive = TRUE, showWarnings = FALSE)
#dir.create(file.path(out_dir, "summary"), recursive = TRUE, showWarnings = FALSE)

cfg <- list(
  n = 3200,
  p = 600,
  s_star_grid = c(5, 15, 20,30),
  eta = 0.02,
  c_r = 0.5,
  T1 = 128,
  T2 = 8,
  m_min = 10,
  reps = 100,
  seed_base = 20260420L,
  df_X = 2.5,
  scale_X = 1,
  df_eps = 1.5
)

num_cores <- max(1, parallel::detectCores() - 1)
cl <- parallel::makeCluster(num_cores)
doParallel::registerDoParallel(cl)

raw_results <- foreach(
  s_star = cfg$s_star_grid,
  .combine = dplyr::bind_rows,
  .packages = c("dplyr", "MASS", "mvtnorm","here")
) %:%
  foreach(
    rep_id = 1:cfg$reps,
    .combine = dplyr::bind_rows,
    .packages = c("dplyr", "MASS", "mvtnorm","here")
  ) %dopar% {
    source(here("simulations","07_schedule_ablation","utils_schedule_linear.R"))
    
    s <- 2 * s_star
    theta_star <- make_theta_star(cfg$p, s_star, magnitude = 5)
    
    dat <- gen_one_linear_dataset(
      seed = cfg$seed_base + 100000L * match(s_star, cfg$s_star_grid) + rep_id,
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
    
    ts_budget <- compute_ts_budget(
      n = cfg$n,
      s_star = s_star,
      p = cfg$p,
      c_r = cfg$c_r,
      T1 = cfg$T1,
      T2 = cfg$T2,
      m_min = cfg$m_min
    )
    
    ts_out <- run_ts_right_linear(
      X = X, y = y, theta_star = theta_star,
      s = s, eta = cfg$eta,
      c_r = cfg$c_r, T1 = cfg$T1, T2 = cfg$T2,
      m_min = cfg$m_min
    )
    
    fd_out <- run_fd_right_linear(
      X = X, y = y, theta_star = theta_star,
      s = s, eta = cfg$eta,
      T_fd = ts_budget$T_fd_budget,
      m_min = cfg$m_min
    )
    
    T_fs <- floor(cfg$n / 70)
    fs_out <- run_fs_right_linear(
      X = X, y = y, theta_star = theta_star,
      s = s, eta = cfg$eta,
      T_fs = T_fs, m_min = cfg$m_min
    )
    
    
    ts_res <- ts_out %>% dplyr::mutate(method = "TS")
    fd_res <- fd_out %>% dplyr::mutate(method = "FD_budget")
    fs_res <- fs_out %>% dplyr::mutate(method = "FS")
    
    
    res_combined <- dplyr::bind_rows(ts_res, fd_res, fs_res)
    
   
    res_combined <- res_combined %>%
      dplyr::mutate(
        n = cfg$n, 
        s_star = s_star, 
        s = s, 
        rep_id = rep_id,
        T1 = cfg$T1, 
        T2 = cfg$T2, 
        c_r = cfg$c_r,
        T_fd_budget = ts_budget$T_fd_budget
      )
    

    
    return(res_combined)
  }

parallel::stopCluster(cl)

saveRDS(raw_results, file.path(out_dir, "raw", "raw_E2_linear_s_sweep.rds"))

summary_tbl <- raw_results %>%
  group_by(s_star, s, T_fd_budget, method) %>%
  summarise(
    median_l2 = median(l2_error, na.rm = TRUE),
    mean_l2 = mean(l2_error, na.rm = TRUE),
    iqr_l2 = IQR(l2_error, na.rm = TRUE),
    median_runtime = median(runtime_sec, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(s_star, method)

write.csv(summary_tbl, file.path(out_dir, "summary", "summary_E2_linear_s_sweep.csv"), row.names = FALSE)
print(summary_tbl)
summary_tbl<-read.csv(file.path(out_dir, "summary", "summary_E2_linear_s_sweep.csv"))
library(ggplot2)
p<-ggplot(data=summary_tbl,mapping = aes(x=s_star,y=log(median_l2),color=method))+
  geom_line()+
  labs(x="True Sparsity s star",
       y="Log Median L2 Error")
p

ggsave(file.path(out_dir, "summary", "linear_s_vs_methods.pdf"),p,height=7,width=6.5)
