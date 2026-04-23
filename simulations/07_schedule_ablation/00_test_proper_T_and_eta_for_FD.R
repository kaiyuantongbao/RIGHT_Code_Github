rm(list = ls())

library(here)
library(dplyr)
library(ggplot2) # 新增画图包
library(doParallel)
library(foreach)

source(here("simulations", "07_schedule_ablation", "utils_schedule_linear.R"))

out_dir <- "results/07_schedule_ablation/E1_linear_FD_budget_pilot"
dir.create(file.path(out_dir, "raw"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "summary"), recursive = TRUE, showWarnings = FALSE)


cfg <- list(
  n_grid = c(2000),        
  p = 600,
  s_star = 10,
  s = 20,
  eta = 0.02,
  m_min = 10,
  reps = 100,               # 修改为 100 次重复
  seed_base = 20260417L,
  df_X = 2.5,
  scale_X = 1,
  df_eps = 1.5
)

# -------- manually edit this after TS plateau pilot --------
TS_cfg <- list(
  c_r = 0.5,
  T1 = 192,
  T2 = 8
)

theta_star <- make_theta_star(cfg$p, cfg$s_star, magnitude = 5)

num_cores <- max(1, parallel::detectCores() - 1)
cl <- parallel::makeCluster(num_cores)
doParallel::registerDoParallel(cl)

# 确保导出带追踪功能的函数
raw_results <- foreach(
  n = cfg$n_grid,
  .combine = dplyr::bind_rows,
  .packages = c("dplyr", "MASS", "mvtnorm", "here")
) %:%
  foreach(
    rep_id = 1:cfg$reps,
    .combine = dplyr::bind_rows,
    .export = c("run_fd_right_linear"), 
    .packages = c("dplyr", "MASS", "mvtnorm")
  ) %dopar% {
    source(here("simulations", "07_schedule_ablation", "utils_schedule_linear.R"))
    
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
    
    # 计算理论 Budget
    budget <- compute_ts_budget(
      n = n, s_star = cfg$s_star, p = cfg$p,
      c_r = TS_cfg$c_r, T1 = TS_cfg$T1, T2 = TS_cfg$T2,
      m_min = cfg$m_min
    )
    
    T_budget <- budget$T_fd_budget
    
    # 【修改2】计算我们要跑的最大步数，即 1.5 倍的 budget
    T_max_run <- max(8, round(1.50 * T_budget))
    
    # 【修改3】抛弃循环，直接利用新版函数的 trace 机制，一口气跑到 1.5 倍
    out <- run_fd_right_linear(
      X = X, y = y, theta_star = theta_star,
      s = cfg$s, eta = cfg$eta,
      T_fd = T_max_run, m_min = cfg$m_min,
      record_trace = TRUE,      # 开启追踪
      record_l2 = TRUE,         # 记录误差
      record_initial = TRUE     # 从 t=0 开始记录
    )
    
    # 拼合一些元数据并返回
    cbind(
      data.frame(
        n = n, rep_id = rep_id,
        T_fd_budget = T_budget,
        c_r = TS_cfg$c_r, T1 = TS_cfg$T1, T2 = TS_cfg$T2
      ),
      out
    )
  }

parallel::stopCluster(cl)

saveRDS(raw_results, file.path(out_dir, "raw", "raw_FD_budget_trajectory.rds"))

# ==============================================================================
# 数据汇总与绘图 (画 100 次的中位数下降曲线和 Budget 竖线)
# ==============================================================================

# 提取这一组固定的 T_budget（因为 n 只有 2000，预算对于所有 rep 都是恒定的）
fixed_T_budget <- raw_results$T_fd_budget[1]

# 计算各步数下的中位数误差
trace_summary <- raw_results %>%
  group_by(t) %>%
  summarise(
    median_l2 = median(l2_error, na.rm = TRUE),
    .groups = "drop"
  )

# 计算我们要画竖线的位置
vlines_data <- data.frame(
  multiplier = c(0.75, 1.00, 1.25, 1.50),
  x_pos = round(c(0.75, 1.00, 1.25, 1.50) * fixed_T_budget)
)
vlines_data$label <- paste0(vlines_data$multiplier, "x Budget (T=", vlines_data$x_pos, ")")

# 画图
p_traj <- ggplot(trace_summary, aes(x = t, y = median_l2)) +
  geom_line(color = "steelblue", linewidth = 1.2) +
  # 添加竖线
  geom_vline(data = vlines_data, aes(xintercept = x_pos), 
             color = c("grey50", "red", "grey50", "grey50"), 
             linetype = c("dashed", "solid", "dashed", "dashed"),
             linewidth = c(0.8, 1.2, 0.8, 0.8)) +
  # 添加文本注释（稍微偏右上方一点避免挡住线）
  geom_text(data = vlines_data, aes(x = x_pos, y = max(trace_summary$median_l2), label = label),
            angle = 90, vjust = -0.5, hjust = 1, size = 3.5, color = "black") +
  theme_bw() +
  labs(
    title = paste0("FD Algorithm Convergence Trajectory (n = ", cfg$n_grid[1], ", Reps = ", cfg$reps, ")"),
    subtitle = paste0("Red line marks 1.00x Budget (T_budget = ", fixed_T_budget, ")"),
    x = expression(Iteration~Step~T),
    y = expression(Median~L[2]~Error)
  ) +
  scale_y_log10() # 使用对数刻度更好观察收敛率

print(p_traj)
 ggsave(file.path(out_dir, "summary", "FD_budget_trajectory.pdf"), plot = p_traj, width = 8, height = 6)