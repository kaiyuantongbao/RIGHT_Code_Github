rm(list = ls())

library(dplyr)
library(ggplot2)  # 新增：用于画图
library(doParallel)
library(foreach)

source("./simulations/07_schedule_ablation/utils_schedule_linear.R")

 out_dir <- "./results/07_schedule_ablation/test_proper_T_and_eta_for_TS"
 dir.create(file.path(out_dir, "raw"), recursive = TRUE, showWarnings = FALSE)
 dir.create(file.path(out_dir, "summary"), recursive = TRUE, showWarnings = FALSE)
 dir.create(file.path(out_dir, "figure"), recursive = TRUE, showWarnings = FALSE)

cfg <- list(
  n = 2000,
  p = 600,
  s_star = 10,
  s = 20,
  eta = 0.02,
  c_r = 0.5,
  T1_grid = c(128, 192, 240),
  T2_grid = c(8, 12),
  m_min = 10,
  reps = 100,
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
  # 建议加上 .export，防止并行计算时找不到外层环境定义的函数
  .export = c("run_ts_right_linear"), 
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
      m_min = cfg$m_min,
      # 【修改 1】：打开追踪开关，获取带 t 维度的长格式数据
      record_trace = TRUE,      
      record_l2 = TRUE,         
      record_initial = TRUE     
    )
    cbind(data.frame(rep_id = rep_id), out)
  })
  
  bind_rows(ans)
}

parallel::stopCluster(cl)

saveRDS(raw_results, file.path(out_dir, "raw", "raw_TS_plateau_pilot.rds"))

# ==============================================================================
# 结果汇总与画图
# ==============================================================================

# 【修改 2】：只取最终收敛的那个点来计算表格（否则会把前面迭代步骤的大误差也平均进去）
final_results <- raw_results %>%
  filter(t == T1 + T2 | !eligible)

summary_tbl <- final_results %>%
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


# 过滤掉不符合预算的数据，然后计算每个 (T1, T2) 组合下，每一步 t 的中位数误差
trace_summary <- raw_results %>%
  filter(eligible) %>%
  group_by(T1, T2, t) %>%
  summarise(median_l2 = median(l2_error, na.rm = TRUE), .groups = "drop") %>%
  mutate(Config = factor(paste0("T1=", T1, ", T2=", T2)))

# 绘制误差下降轨迹图
p_trace <- ggplot(trace_summary, aes(x = t, y = median_l2, color = Config)) +
  geom_line(linewidth = 1) +
  theme_bw() +
  labs(
    title = "Median L2 Error Path over Iterations (100 Replications)",
    x = "Iteration Step (t)",
    y = expression(Median~L[2]~Error)
  ) +
  scale_y_log10() # 对数 Y 轴通常能更清晰地展示收敛率

print(p_trace)

ggsave(file.path(out_dir, "figure", "TS_error_trajectory.pdf"), plot = p_trace, width = 8, height = 6)
