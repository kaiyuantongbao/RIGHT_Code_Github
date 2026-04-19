rm(list = ls())

library(dplyr)
library(foreach)
library(doParallel)

source("R/data_generator.R")
source("R/solvers.R")

# ============================================================
# 0. Config
# ============================================================
out_dir <- "results/07_schedule_ablation/E0_linear"
dir.create(file.path(out_dir, "raw"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "summary"), recursive = TRUE, showWarnings = FALSE)

cfg <- list(
  n = 1800,
  p = 600,
  s_star = 10,
  s = 20,
  eta = 0.01,
  df_X = 2.5,
  df_eps = 1.5,
  Sigma_X = diag(600),
  scale_X = 1,
  m_min = 10,
  reps = 80,
  seed_base = 20260416L,
  c_r_grid = c(1/4, 1/2, 1),
  T1_grid = c(64, 96, 128),
  T2_grid = c(4, 8, 12)
)

theta_star <- c(rep(5, cfg$s_star), rep(0, cfg$p - cfg$s_star))

# ============================================================
# 1. Helper: one dataset
# ============================================================
gen_one_dataset <- function(seed, cfg, theta_star) {
  set.seed(seed)
  dat <- generate_linear_data(
    n = cfg$n,
    p = cfg$p,
    theta_star = theta_star,
    generator_X = generate_X_t,
    generator_epsilon = generate_epsilon_t,
    df_X = cfg$df_X,
    Sigma_X = cfg$Sigma_X,
    scale_X = cfg$scale_X,
    df_eps = cfg$df_eps
  )
  list(
    X = dat$X,
    y = as.numeric(dat$y)
  )
}

# ============================================================
# 2. Helper: one two-stage run
# ============================================================
run_ts_right_linear <- function(X, y, theta_star, s, eta, c_r, T1, T2, m_min = 10) {
  n <- nrow(X)
  p <- ncol(X)
  s_star <- sum(theta_star != 0)
  
  # split ratio: n1 / n2 = c_r * s_star / log(s_star)
  ratio <- c_r * s_star / log(s_star)
  rho <- ratio / (1 + ratio)
  
  n1 <- floor(rho * n)
  n2 <- n - n1
  
  # sample split
  idx <- sample.int(n)
  id1 <- idx[1:n1]
  id2 <- idx[(n1 + 1):n]
  
  # phase 1
  K1_tar <- ceiling(s_star * log(p))
  K1 <- min(K1_tar, floor(n1 / m_min))
  
  # phase 2
  b2 <- floor(n2 / T2)
  K2_tar <- ceiling(log(p))
  K2 <- min(K2_tar, floor(b2 / m_min))
  
  feasible <- (n1 >= m_min) && (b2 >= m_min) && (K1 >= 1) && (K2 >= 1)
  
  if (!feasible) {
    return(data.frame(
      n1 = n1, n2 = n2, b2 = b2,
      K1_tar = K1_tar, K1 = K1,
      K2_tar = K2_tar, K2 = K2,
      feasible = FALSE,
      l2_error = NA_real_,
      runtime_sec = NA_real_
    ))
  }
  
  t0 <- proc.time()[3]
  
  # shuffle phase-1 subset once so MoM blocks are not tied to original row order
  X1 <- X[id1, , drop = FALSE]
  y1 <- y[id1]
  ord1 <- sample.int(n1)
  X1 <- X1[ord1, , drop = FALSE]
  y1 <- y1[ord1]
  
  theta_cur <- solver_right(
    X = X1,
    y = y1,
    s = s,
    eta = eta,
    T_max = T1,
    K = K1,
    theta_init = rep(0, p),
    grad_func_samplewise = grad_linear_regression_samplewise
  )
  
  # phase 2: split into T2 batches, one robust step per batch
  id2 <- sample(id2)
  use_n2 <- b2 * T2
  id2 <- id2[seq_len(use_n2)]
  
  for (tt in seq_len(T2)) {
    batch_ids <- id2[((tt - 1) * b2 + 1):(tt * b2)]
    
    theta_cur <- solver_right(
      X = X[batch_ids, , drop = FALSE],
      y = y[batch_ids],
      s = s,
      eta = eta,
      T_max = 1,
      K = K2,
      theta_init = theta_cur,
      grad_func_samplewise = grad_linear_regression_samplewise
    )
  }
  
  runtime_sec <- proc.time()[3] - t0
  l2_error <- sqrt(sum((theta_cur - theta_star)^2))
  
  data.frame(
    n1 = n1, n2 = n2, b2 = b2,
    K1_tar = K1_tar, K1 = K1,
    K2_tar = K2_tar, K2 = K2,
    feasible = TRUE,
    l2_error = l2_error,
    runtime_sec = runtime_sec
  )
}

# ============================================================
# 3. Grid
# ============================================================
grid <- expand.grid(
  c_r = cfg$c_r_grid,
  T1 = cfg$T1_grid,
  T2 = cfg$T2_grid,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

# ============================================================
# 4. Parallel over rep_id
# ============================================================
n_cores <- min(24, parallel::detectCores() - 1)
cl <- parallel::makeCluster(n_cores)
doParallel::registerDoParallel(cl)

raw_results <- foreach::foreach(
  rep_id = 1:cfg$reps,
  .combine = dplyr::bind_rows,
  .packages = c("dplyr", "mvtnorm", "MASS")
) %dopar% {
  
  source("R/data_generator.R")
  source("R/solvers.R")
  
  dat <- gen_one_dataset(cfg$seed_base + rep_id, cfg, theta_star)
  X <- dat$X
  y <- dat$y
  
  one_rep <- lapply(seq_len(nrow(grid)), function(j) {
    pars <- grid[j, ]
    
    out <- run_ts_right_linear(
      X = X,
      y = y,
      theta_star = theta_star,
      s = cfg$s,
      eta = cfg$eta,
      c_r = pars$c_r,
      T1 = pars$T1,
      T2 = pars$T2,
      m_min = cfg$m_min
    )
    
    cbind(
      data.frame(
        rep_id = rep_id,
        c_r = pars$c_r,
        T1 = pars$T1,
        T2 = pars$T2
      ),
      out
    )
  })
  
  dplyr::bind_rows(one_rep)
}

parallel::stopCluster(cl)

saveRDS(raw_results, file.path(out_dir, "raw", "E0_linear_raw_results.rds"))

# ============================================================
# 5. Summary + pick default config
# ============================================================
summary_tbl <- raw_results %>%
  group_by(c_r, T1, T2) %>%
  summarise(
    n1 = first(n1),
    n2 = first(n2),
    b2 = first(b2),
    K1_tar = first(K1_tar),
    K1 = first(K1),
    K2_tar = first(K2_tar),
    K2 = first(K2),
    feasible_all = all(feasible),
    median_l2 = median(l2_error, na.rm = TRUE),
    mean_l2 = mean(l2_error, na.rm = TRUE),
    iqr_l2 = IQR(l2_error, na.rm = TRUE),
    median_runtime = median(runtime_sec, na.rm = TRUE),
    cost_proxy = first(T1 * n1 + T2 * b2),
    .groups = "drop"
  ) %>%
  mutate(
    eligible = feasible_all & (K2 >= 4)
  ) %>%
  arrange(median_l2, T1, T2, c_r)

best_med <- min(summary_tbl$median_l2[summary_tbl$eligible], na.rm = TRUE)

default_pick <- summary_tbl %>%
  filter(eligible, median_l2 <= 1.05 * best_med) %>%
  arrange(cost_proxy, T1, T2, c_r) %>%
  slice(1)

saveRDS(summary_tbl, file.path(out_dir, "summary", "E0_linear_summary.rds"))
write.csv(summary_tbl, file.path(out_dir, "summary", "E0_linear_summary.csv"), row.names = FALSE)
write.csv(default_pick, file.path(out_dir, "summary", "E0_linear_default_pick.csv"), row.names = FALSE)

print(summary_tbl)
print(default_pick)
