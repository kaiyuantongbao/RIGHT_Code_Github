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

run_fd_right_linear <- function(X, y, theta_star, s, eta, T_fd, m_min = 10) {
  n <- nrow(X)
  p <- ncol(X)
  s_star <- sum(theta_star != 0)
  
  K_fd_tar <- ceiling(s_star * log(p))
  K_fd <- min(K_fd_tar, floor(n / m_min))
  
  t0 <- proc.time()[3]
  theta_hat <- solver_right(
    X = X,
    y = y,
    s = s,
    eta = eta,
    T_max = T_fd,
    K = K_fd,
    theta_init = rep(0, p),
    grad_func_samplewise = grad_linear_regression_samplewise
  )
  runtime_sec <- proc.time()[3] - t0
  l2_error <- sqrt(sum((theta_hat - theta_star)^2))
  
  data.frame(
    arm = "FD",
    T_fd = T_fd,
    K_fd_tar = K_fd_tar,
    K_fd = K_fd,
    l2_error = l2_error,
    runtime_sec = runtime_sec
  )
}

run_ts_right_linear <- function(X, y, theta_star, s, eta, c_r, T1, T2, m_min = 10) {
  n <- nrow(X)
  p <- ncol(X)
  s_star <- sum(theta_star != 0)
  
  budget <- compute_ts_budget(n, s_star, p, c_r, T1, T2, m_min)
  if (!budget$eligible) {
    return(data.frame(
      arm = "TS",
      c_r = c_r, T1 = T1, T2 = T2,
      n1 = budget$n1, n2 = budget$n2,
      b2 = budget$b2,
      K1_tar = budget$K1_tar, K1 = budget$K1,
      K2_tar = budget$K2_tar, K2 = budget$K2,
      cost_proxy = budget$cost_proxy,
      T_fd_budget = budget$T_fd_budget,
      eligible = FALSE,
      l2_error = NA_real_,
      runtime_sec = NA_real_
    ))
  }
  
  idx <- sample.int(n)
  id1 <- idx[1:budget$n1]
  id2 <- idx[(budget$n1 + 1):n]
  
  X1 <- X[id1, , drop = FALSE]
  y1 <- y[id1]
  
  t0 <- proc.time()[3]
  
  theta_cur <- solver_right(
    X = X1,
    y = y1,
    s = s,
    eta = eta,
    T_max = T1,
    K = budget$K1,
    theta_init = rep(0, p),
    grad_func_samplewise = grad_linear_regression_samplewise
  )
  
  use_n2 <- budget$b2 * T2
  id2 <- sample(id2)[seq_len(use_n2)]
  
  for (tt in seq_len(T2)) {
    batch_ids <- id2[((tt - 1) * budget$b2 + 1):(tt * budget$b2)]
    theta_cur <- solver_right(
      X = X[batch_ids, , drop = FALSE],
      y = y[batch_ids],
      s = s,
      eta = eta,
      T_max = 1,
      K = budget$K2,
      theta_init = theta_cur,
      grad_func_samplewise = grad_linear_regression_samplewise
    )
  }
  
  runtime_sec <- proc.time()[3] - t0
  l2_error <- sqrt(sum((theta_cur - theta_star)^2))
  
  data.frame(
    arm = "TS",
    c_r = c_r, T1 = T1, T2 = T2,
    n1 = budget$n1, n2 = budget$n2,
    b2 = budget$b2,
    K1_tar = budget$K1_tar, K1 = budget$K1,
    K2_tar = budget$K2_tar, K2 = budget$K2,
    cost_proxy = budget$cost_proxy,
    T_fd_budget = budget$T_fd_budget,
    eligible = TRUE,
    l2_error = l2_error,
    runtime_sec = runtime_sec
  )
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