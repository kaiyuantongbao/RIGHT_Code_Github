rm(list = ls())

library(dplyr)
library(MASS)
library(mvtnorm)
library(foreach)
library(doParallel)

library(here)
source(here("R", "data_generator.R"))
source(here("R", "solvers.R"))

# --------------------------------------------------
# Common config helpers
# --------------------------------------------------
make_theta_star <- function(p, s_star, magnitude = 5) {
  theta_star <- rep(0, p)
  theta_star[1:s_star] <- magnitude
  theta_star
}

split_ratio_to_n1n2 <- function(n, s_star, c_r) {
  ratio <- c_r * s_star / log(s_star)
  rho <- ratio / (1 + ratio)
  n1 <- floor(rho * n)
  n2 <- n - n1
  list(n1 = n1, n2 = n2, rho = rho)
}

compute_ts_budget <- function(n, s_star, p, c_r, T1, T2, m_min = 10) {
  tmp <- split_ratio_to_n1n2(n, s_star, c_r)
  n1 <- tmp$n1
  n2 <- tmp$n2
  
  K1_tar <- ceiling(s_star * log(p))
  K1 <- min(K1_tar, floor(n1 / m_min))
  
  b2 <- floor(n2 / T2)
  K2_tar <- ceiling(log(p))
  K2 <- min(K2_tar, floor(b2 / m_min))
  
  cost_proxy <- T1 * n1 + T2 * b2
  T_fd_budget <- ceiling(cost_proxy / n)
  
  list(
    n1 = n1, n2 = n2, rho = tmp$rho,
    K1_tar = K1_tar, K1 = K1,
    b2 = b2, K2_tar = K2_tar, K2 = K2,
    cost_proxy = cost_proxy,
    T_fd_budget = T_fd_budget,
    eligible = (K1 >= 1) && (K2 >= 4) && (b2 >= m_min)
  )
}

gen_one_linear_dataset <- function(seed, n, p, theta_star,
                                   df_X = 2.5,
                                   Sigma_X = NULL,
                                   scale_X = 1,
                                   df_eps = 1.5) {
  set.seed(seed)
  if (is.null(Sigma_X)) Sigma_X <- diag(p)
  
  dat <- generate_linear_data(
    n = n,
    p = p,
    theta_star = theta_star,
    generator_X = generate_X_t,
    generator_epsilon = generate_epsilon_t,
    df_X = df_X,
    Sigma_X = Sigma_X,
    scale_X = scale_X,
    df_eps = df_eps
  )
  
  list(X = dat$X, y = as.numeric(dat$y))
}

run_fd_right_linear <- function(X, y, theta_star, s, eta, T_fd, m_min = 10,
                                record_trace = FALSE,      # 新增：主开关
                                record_l2 = FALSE,         # 新增：L2 开关
                                record_initial = FALSE) {  # 新增：初始值开关
  n <- nrow(X)
  p <- ncol(X)
  s_star <- sum(theta_star != 0)
  
  K_fd_tar <- ceiling(s_star * log(p))
  K_fd <- min(K_fd_tar, floor(n / m_min))
  
  t0 <- proc.time()[3]
  
  # 调用 solver_right，传入追踪开关和真实的 theta_star
  out_fd <- solver_right(
    X = X,
    y = y,
    s = s,
    eta = eta,
    T_max = T_fd,
    K = K_fd,
    theta_init = rep(0, p),
    grad_func_samplewise = grad_linear_regression_samplewise,
    theta_star = theta_star,         # 【关键新增】：底层算 L2 需要真实的参数
    record_trace = record_trace,     # 传入开关
    record_l2 = record_l2,           # 传入开关
    record_initial = record_initial  # 传入开关
  )
  
  runtime_sec <- proc.time()[3] - t0
  
  # 构建基础信息数据框（不包含误差和时间步 t）
  base_df <- data.frame(
    arm = "FD",
    T_fd = T_fd,
    K_fd_tar = K_fd_tar,
    K_fd = K_fd,
    runtime_sec = runtime_sec
  )
  
  # 根据追踪开关决定返回格式
  if (record_trace && record_l2) {
    # 追踪模式：解析底层返回的 list
    trace_t <- out_fd$trace$iteration
    trace_l2 <- out_fd$trace$l2_error
    
    # 扩展为长格式
    out_df <- base_df[rep(1, length(trace_t)), ] 
    out_df$t <- trace_t
    out_df$l2_error <- trace_l2
    return(out_df)
    
  } else {
    # 极速模式：解析底层返回的 theta_hat，手动算一次最终误差
    theta_hat <- if (record_trace) out_fd$theta else out_fd
    l2_error <- sqrt(sum((theta_hat - theta_star)^2))
    
    base_df$l2_error <- l2_error
    return(base_df)
  }
}

run_ts_right_linear <- function(X, y, theta_star, s, eta, c_r, T1, T2, m_min = 10,
                                record_trace = FALSE,      # 新增：主开关，默认为 FALSE
                                record_l2 = FALSE,         # 新增：L2 开关，默认为 FALSE
                                record_initial = FALSE) {  # 新增：初始值开关，默认为 FALSE
  n <- nrow(X)
  p <- ncol(X)
  s_star <- sum(theta_star != 0)
  
  budget <- compute_ts_budget(n, s_star, p, c_r, T1, T2, m_min)
  if (!budget$eligible) {
    # 如果不可行，无论是否开启追踪，都返回统一格式的 NA
    res <- data.frame(
      arm = "TS", c_r = c_r, T1 = T1, T2 = T2,
      n1 = budget$n1, n2 = budget$n2, b2 = budget$b2,
      K1_tar = budget$K1_tar, K1 = budget$K1,
      K2_tar = budget$K2_tar, K2 = budget$K2,
      cost_proxy = budget$cost_proxy, T_fd_budget = budget$T_fd_budget,
      eligible = FALSE, runtime_sec = NA_real_,
      l2_error = NA_real_
    )
    if (record_trace) res$t <- NA_integer_ # 追踪模式特有的 t 字段
    return(res)
  }
  
  idx <- sample.int(n)
  id1 <- idx[1:budget$n1]
  id2 <- idx[(budget$n1 + 1):n]
  
  X1 <- X[id1, , drop = FALSE]
  y1 <- y[id1]
  
  t0 <- proc.time()[3]
  
  # --- Stage 1 ---
  out_s1 <- solver_right(
    X = X1, y = y1, s = s, eta = eta, T_max = T1, K = budget$K1,
    theta_init = rep(0, p), grad_func_samplewise = grad_linear_regression_samplewise,
    theta_star = theta_star,
    record_trace = record_trace,     # 传入外部参数
    record_l2 = record_l2,           # 传入外部参数
    record_initial = record_initial  # 传入外部参数
  )
  
  # 如果开启了追踪，out_s1 是一个 list；如果没有，它直接就是 theta
  if (record_trace) {
    theta_cur <- out_s1$theta
    trace_t <- out_s1$trace$iteration
    trace_l2 <- out_s1$trace$l2_error
  } else {
    theta_cur <- out_s1
  }
  
  # --- Stage 2 ---
  use_n2 <- budget$b2 * T2
  id2 <- sample(id2)[seq_len(use_n2)]
  
  for (tt in seq_len(T2)) {
    batch_ids <- id2[((tt - 1) * budget$b2 + 1):(tt * budget$b2)]
    
    # Stage 2 无论如何每次只走1步，不在这里面调用内层追踪
    theta_cur <- solver_right(
      X = X[batch_ids, , drop = FALSE], y = y[batch_ids], s = s, eta = eta,
      T_max = 1, K = budget$K2, theta_init = theta_cur,
      grad_func_samplewise = grad_linear_regression_samplewise,
      record_trace = FALSE # 强制关掉内层，我们自己在外面拼接
    )
    
    # 如果外层开启了追踪，我们手动把 Stage 2 的每一步误差拼进去
    if (record_trace && record_l2) {
      trace_t <- c(trace_t, T1 + tt)
      trace_l2 <- c(trace_l2, sqrt(sum((theta_cur - theta_star)^2)))
    }
  }
  
  runtime_sec <- proc.time()[3] - t0
  final_l2_error <- sqrt(sum((theta_cur - theta_star)^2))
  
  # --- 返回结果 ---
  # 构建基础信息数据框
  base_df <- data.frame(
    arm = "TS", c_r = c_r, T1 = T1, T2 = T2,
    n1 = budget$n1, n2 = budget$n2, b2 = budget$b2,
    K1_tar = budget$K1_tar, K1 = budget$K1,
    K2_tar = budget$K2_tar, K2 = budget$K2,
    cost_proxy = budget$cost_proxy, T_fd_budget = budget$T_fd_budget,
    eligible = TRUE, runtime_sec = runtime_sec
  )
  
  if (record_trace && record_l2) {
    # 追踪模式：返回长格式（行数为 T1 + T2 + 1）
    out_df <- base_df[rep(1, length(trace_t)), ] # 复制基础信息行
    out_df$t <- trace_t
    out_df$l2_error <- trace_l2
    return(out_df)
  } else {
    # 极速模式（不追踪）：只返回 1 行最终结果，格式和原来的一模一样
    base_df$l2_error <- final_l2_error
    return(base_df)
  }
}

run_fs_right_linear <- function(X, y, theta_star, s, eta, T_fs, m_min = 10) {
  n <- nrow(X)
  p <- ncol(X)
  
  b_fs <- floor(n / T_fs)
  K_fs_tar <- ceiling(log(p))
  K_fs <- min(K_fs_tar, floor(b_fs / m_min))
  
  eligible <- (K_fs >= 4) && (b_fs >= m_min)
  if (!eligible) {
    return(data.frame(
      arm = "FS",
      T_fs = T_fs,
      b_fs = b_fs,
      K_fs_tar = K_fs_tar,
      K_fs = K_fs,
      eligible = FALSE,
      l2_error = NA_real_,
      runtime_sec = NA_real_
    ))
  }
  
  idx <- sample.int(n)
  use_n <- b_fs * T_fs
  idx <- idx[seq_len(use_n)]
  
  t0 <- proc.time()[3]
  theta_cur <- rep(0, p)
  
  for (tt in seq_len(T_fs)) {
    batch_ids <- idx[((tt - 1) * b_fs + 1):(tt * b_fs)]
    theta_cur <- solver_right(
      X = X[batch_ids, , drop = FALSE],
      y = y[batch_ids],
      s = s,
      eta = eta,
      T_max = 1,
      K = K_fs,
      theta_init = theta_cur,
      grad_func_samplewise = grad_linear_regression_samplewise
    )
  }
  
  runtime_sec <- proc.time()[3] - t0
  l2_error <- sqrt(sum((theta_cur - theta_star)^2))
  
  data.frame(
    arm = "FS",
    T_fs = T_fs,
    b_fs = b_fs,
    K_fs_tar = K_fs_tar,
    K_fs = K_fs,
    eligible = TRUE,
    l2_error = l2_error,
    runtime_sec = runtime_sec
  )
}