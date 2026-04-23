rm(list = ls())

library(here)
library(dplyr)
library(ggplot2)
library(doParallel)
library(foreach)

# 确保加载了你之前修改好的、带有追踪开关的两个核心函数
source(here("simulations", "07_schedule_ablation", "utils_schedule_linear.R"))

out_dir <- "results/07_schedule_ablation/TS_vs_FD_trajectory"
dir.create(file.path(out_dir, "raw"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "summary"), recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# 1. 实验参数配置
# ==============================================================================
cfg <- list(
  n = 2000,
  p = 600,
  s_star = 10,
  s = 20,
  eta = 0.02,       
  c_r = 0.5,
  T1_grid = c(128, 150, 200),
  T2 = 8,
  m_min = 10,
  reps = 100,        
  seed_base =20260417L,
  df_X = 2.5,
  scale_X = 1,
  df_eps = 1.5
)

theta_star <- make_theta_star(cfg$p, cfg$s_star, magnitude = 5)

num_cores <- max(1, parallel::detectCores() - 1)
cl <- parallel::makeCluster(num_cores)
doParallel::registerDoParallel(cl)

# ==============================================================================
# 2. 并行主循环 (公平对抗设计)
# ==============================================================================
raw_results <- foreach(
  rep_id = 1:cfg$reps,
  .combine = dplyr::bind_rows,
  .export = c("run_ts_right_linear", "run_fd_right_linear"),
  .packages = c("dplyr", "MASS", "mvtnorm", "here")
) %dopar% {
  source(here("simulations", "07_schedule_ablation", "utils_schedule_linear.R"))
  
  # 【公平性保证 1】: 每一次 rep_id 生成唯一一份数据
  dat <- gen_one_linear_dataset(
    seed = cfg$seed_base + rep_id,
    n = cfg$n, p = cfg$p, theta_star = theta_star,
    df_X = cfg$df_X, Sigma_X = diag(cfg$p), scale_X = cfg$scale_X, df_eps = cfg$df_eps
  )
  X <- dat$X
  y <- dat$y
  
  # -------- 跑 3 种 TS 算法 --------
  ans_ts <- lapply(cfg$T1_grid, function(t1) {
    out <- run_ts_right_linear(
      X = X, y = y, theta_star = theta_star,
      s = cfg$s, eta = cfg$eta, c_r = cfg$c_r, 
      T1 = t1, T2 = cfg$T2, m_min = cfg$m_min,
      record_trace = TRUE, record_l2 = TRUE, record_initial = TRUE
    )
    # 打上清晰的配置标签，方便画图
    out$Config <- paste0("TS (T1=", t1, ", T2=", cfg$T2, ")")
    out$arm_type <- "TS"
    return(out)
  })
  df_ts <- bind_rows(ans_ts)
  
  # -------- 跑 1 种 FD 算法 (基线) --------
  # 计算出这 3 种 TS 中，最大的等效 FD 预算是多少
  max_budget <- max(df_ts$T_fd_budget, na.rm = TRUE)
  
  # 为了让 FD 的曲线能在图上完全覆盖 TS 的横坐标跨度，我们让它多跑一点
  # 比如跑到 TS 最大步数 (240+8) 和 最大预算 中的较大值，再乘个 1.1
  T_fd_run <- ceiling(max(max(cfg$T1_grid) + cfg$T2, max_budget) * 1.1)
  
  # 【公平性保证 2】: 使用绝对相同的 X 和 y 跑 FD
  out_fd <- run_fd_right_linear(
    X = X, y = y, theta_star = theta_star,
    s = cfg$s, eta = cfg$eta, T_fd = T_fd_run, m_min = cfg$m_min,
    record_trace = TRUE, record_l2 = TRUE, record_initial = TRUE
  )
  
  out_fd$Config <- "FD (Baseline)"
  out_fd$arm_type <- "FD"
  # FD 没有预算字段，为了 bind_rows 不报错，我们给它补上 NA
  out_fd$T_fd_budget <- NA_real_ 
  
  # 合并当前 rep 的所有结果
  bind_rows(df_ts, out_fd) %>% mutate(rep_id = rep_id)
}

parallel::stopCluster(cl)
 saveRDS(raw_results, file.path(out_dir, "raw", "raw_TS_vs_FD_trajectory.rds"))

# ==============================================================================
# 3. 数据处理与可视化
# ==============================================================================

# 计算每一步的中位数误差
trace_summary <- raw_results %>%
  filter(eligible == TRUE | arm_type == "FD") %>%
  group_by(Config, arm_type, t) %>%
  summarise(median_l2 = median(l2_error, na.rm = TRUE), .groups = "drop")

# 提取每种 TS 配置的等效 Budget 竖线位置
budget_lines <- raw_results %>%
  filter(arm_type == "TS", eligible == TRUE) %>%
  group_by(Config) %>%
  summarise(
    T_budget = first(T_fd_budget),
    .groups = "drop"
  )

# 【核心修改区】：为 Color, Linetype, Linewidth 都建立严格映射
all_configs <- unique(trace_summary$Config)
ts_configs <- sort(all_configs[all_configs != "FD (Baseline)"])

# 建立字典：指定 FD 为灰色的加粗虚线
my_colors <- c("FD (Baseline)" = "grey50")
my_linetypes <- c("FD (Baseline)" = "dashed")
my_linewidths <- c("FD (Baseline)" = 1.2)

palette_colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00")

# 给每种 TS 分配颜色，并强制为常规的实线
for(i in seq_along(ts_configs)) {
  color_idx <- ifelse(i <= length(palette_colors), i, length(palette_colors))
  my_colors[ts_configs[i]] <- palette_colors[color_idx]
  my_linetypes[ts_configs[i]] <- "solid"
  my_linewidths[ts_configs[i]] <- 0.8
}

y_max_label <- max(trace_summary$median_l2, na.rm = TRUE)

# 画图：把 color, linetype, linewidth 全部映射给 Config！
p_compare <- ggplot(trace_summary, aes(x = t, y = median_l2, 
                                       color = Config, 
                                       linetype = Config, 
                                       linewidth = Config)) +
  geom_line() +
  geom_vline(data = budget_lines, aes(xintercept = T_budget, color = Config), 
             linetype = "dotted", linewidth = 1) +
  geom_text(data = budget_lines, aes(x = T_budget, y = y_max_label, 
                                     label = paste0("Budget=", T_budget)),
            angle = 90, vjust = -0.5, hjust = 1, size = 3.5, show.legend = FALSE) +
  theme_bw() +
  scale_y_log10() + 
  scale_color_manual(values = my_colors) +
  scale_linetype_manual(values = my_linetypes) +
  scale_linewidth_manual(values = my_linewidths) +
  labs(
    title = "Trajectory Comparison: Two-Stage (TS) vs Full-Dataset (FD)",
    subtitle = "Dotted vertical lines indicate the equivalent computational budget (T_budget) for each TS config.",
    x = expression(Iteration~Step~(t)),
    y = expression(Median~L[2]~Error),
    # 强制让三个图例合并为一个，只要标题一致，ggplot 就会自动拼合
    color = "Algorithm Configuration",
    linetype = "Algorithm Configuration",
    linewidth = "Algorithm Configuration"
  ) +
  theme(
    legend.position = "bottom",
    legend.key.width = unit(3, "lines") # 【小巧思】把图例画长一点，让虚线的间隔更明显
  )


print(p_compare)
 ggsave(file.path(out_dir, "summary", "TS_vs_FD_trajectory_plot.pdf"), plot = p_compare, width = 10.5, height = 6)
 