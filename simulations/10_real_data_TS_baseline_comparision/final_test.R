# GDSC:277 / Linifanib: parallel repeated 5-fold TS tuned + baselines.
#
# Baselines:
#   LASSO
#   IHT
#   Adaptive Huber with tau quantile = 0.99
#   Shrinkage with tau_x quantile = 0.99 and tau_y quantile = 0.98
#   TS-RIGHT-TUNED over K1 cap in c(1, 5, 15, 25)
#
# Run from the project root:
#   "C:/Program Files/R/R-4.4.3/bin/Rscript.exe" simulations/14_GDSC_CCLE_GDSC277/gdsc277_eta02_repeated5fold_trajectory_compare.R

resolve_project_root <- function() {
  candidates <- c(
    getwd(),
    normalizePath(file.path(getwd(), ".."), mustWork = FALSE),
    normalizePath(file.path(getwd(), "..", ".."), mustWork = FALSE)
  )
  
  for (candidate in unique(candidates)) {
    if (file.exists(file.path(candidate, "data", "GDSC_CCLE_GDSC277", "GDSC277_Linifanib_processed.rds"))) {
      return(normalizePath(candidate, mustWork = TRUE))
    }
  }
  
  stop("Could not find data/GDSC_CCLE_GDSC277/GDSC277_Linifanib_processed.rds.")
}

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(glmnet)
  library(parallel)
})

n_cores <- parallel::detectCores(logical = TRUE)
if (!is.finite(n_cores)) n_cores <- 2L

cfg <- list(
  seed = 20260517L,
  n_repeats = 20L,
  n_folds = 5L,
  n_workers = max(1L, min(20L, n_cores - 1L)),
  lasso_nfolds = 5L,
  eta_scale = 0.5,
  q = 0.40,
  k1_cap_grid = c(2L, 5L, 15L, 25L),
  T2 = 4L,
  ts_seed_offset = 21L,
  m_min = 8L,
  t1_grid = seq(0L, 100L, by = 10L),
  compare_ts_T1 = 50L,
  compare_iht_T = 50L,
  huber_tau_q = 0.99,
  shrink_tau_x_q = 0.99,
  shrink_tau_y_q = 0.98,
  huber_T = 500L
)

project_root <- resolve_project_root()
source(file.path(project_root, "R", "solvers.R"))

data_path <- file.path(project_root, "data", "GDSC_CCLE_GDSC277", "GDSC277_Linifanib_processed.rds")
out_dir <- file.path(project_root, "results", "10_realdata_GDSC_CCLE_GDSC277", "eta05_repeated5fold_parallel_tuned_ts_huber_shrinkage")
dir.create(file.path(out_dir, "raw"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "summary"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "figures"), recursive = TRUE, showWarnings = FALSE)

task <- readRDS(data_path)
X_raw <- task$X
y_raw <- task$y

with_seed_preserved <- function(seed, expr) {
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old_seed <- if (had_seed) get(".Random.seed", envir = .GlobalEnv) else NULL
  
  on.exit({
    if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  
  set.seed(seed)
  force(expr)
}

make_folds <- function(n, k, seed) {
  with_seed_preserved(seed, {
    sample(rep(seq_len(k), length.out = n), size = n, replace = FALSE)
  })
}

fit_robust_x_scaler <- function(X_train) {
  centers <- apply(X_train, 2, median, na.rm = TRUE)
  centered_train <- sweep(X_train, 2, centers, FUN = "-")
  scales <- apply(centered_train, 2, mad, na.rm = TRUE)
  scales[!is.finite(scales) | scales == 0] <- 1
  list(center = centers, scale = scales)
}

apply_robust_x_scaler <- function(X, scaler) {
  sweep(sweep(X, 2, scaler$center, FUN = "-"), 2, scaler$scale, FUN = "/")
}

estimate_gradient_step <- function(X_train) {
  gram <- crossprod(X_train) / nrow(X_train)
  lambda_max <- suppressWarnings(eigen(gram, symmetric = TRUE, only.values = TRUE)$values[1])
  if (!is.finite(lambda_max) || lambda_max <= 0) return(1e-3)
  1 / lambda_max
}

hard_threshold <- function(v, s) {
  p <- length(v)
  if (s <= 0L) return(rep(0, p))
  if (s >= p) return(v)
  keep <- order(abs(v), decreasing = TRUE)[seq_len(s)]
  out <- rep(0, p)
  out[keep] <- v[keep]
  out
}

grad_linear_regression <- function(theta, X, y) {
  drop(t(X) %*% (X %*% theta - y)) / nrow(X)
}

grad_linear_regression_samplewise <- function(theta, X, y) {
  X * as.vector(X %*% theta - y)
}

robust_grad_MoM <- function(theta, X, y, K) {
  n <- nrow(X)
  K <- max(1L, min(as.integer(K), n))
  sample_grads <- grad_linear_regression_samplewise(theta, X, y)
  block_id <- cut(seq_len(n), breaks = K, labels = FALSE)
  
  block_grads <- vapply(seq_len(K), function(k) {
    colMeans(sample_grads[block_id == k, , drop = FALSE])
  }, numeric(ncol(X)))
  
  apply(block_grads, 1, median)
}

safe_rmse <- function(theta, X_val, y_val) {
  if (length(theta) != ncol(X_val) || any(!is.finite(theta))) return(NA_real_)
  sqrt(mean((as.numeric(X_val %*% theta) - y_val)^2))
}

safe_mae <- function(theta, X_val, y_val) {
  if (length(theta) != ncol(X_val) || any(!is.finite(theta))) return(NA_real_)
  mean(abs(as.numeric(X_val %*% theta) - y_val))
}

safe_r2 <- function(theta, X_val, y_val) {
  if (length(theta) != ncol(X_val) || any(!is.finite(theta))) return(NA_real_)
  pred <- as.numeric(X_val %*% theta)
  denom <- sum((y_val - mean(y_val))^2)
  if (!is.finite(denom) || denom == 0) return(NA_real_)
  1 - sum((y_val - pred)^2) / denom
}

fit_lasso_tuned <- function(X_train, y_train, seed) {
  fit <- with_seed_preserved(seed, {
    cv.glmnet(
      X_train,
      y_train,
      standardize = FALSE,
      intercept = FALSE,
      nfolds = cfg$lasso_nfolds
    )
  })
  
  theta <- as.vector(coef(fit, s = fit$lambda.min))[-1]
  list(theta = theta, lambda = fit$lambda.min, lambda_grid = fit$lambda, nonzero = sum(theta != 0))
}

fit_iht <- function(X_train, y_train, s, eta, T_max) {
  theta <- rep(0, ncol(X_train))
  for (tt in seq_len(T_max)) {
    theta <- hard_threshold(theta - eta * grad_linear_regression(theta, X_train, y_train), s)
  }
  theta
}

fit_huber_fixed_q <- function(X_train, y_train, lasso_theta, lambda_grid) {
  resid <- as.numeric(y_train - X_train %*% lasso_theta)
  tau <- as.numeric(quantile(abs(resid), probs = cfg$huber_tau_q, na.rm = TRUE))
  tau <- max(tau, 1e-8)
  
  lambda_grid <- sort(unique(as.numeric(lambda_grid)), decreasing = TRUE)
  lambda_grid <- lambda_grid[is.finite(lambda_grid) & lambda_grid > 0]
  if (length(lambda_grid) > 20L) {
    lambda_grid <- lambda_grid[round(seq(1, length(lambda_grid), length.out = 20L))]
  }
  
  candidates <- lapply(lambda_grid, function(lambda_val) {
    theta <- tryCatch(
      solver_huber_lasso(
        X = X_train,
        y = y_train,
        lambda = lambda_val,
        tau = tau,
        T_max = cfg$huber_T,
        beta_init = lasso_theta
      ),
      error = function(e) rep(NA_real_, ncol(X_train))
    )
    
    data.frame(
      lambda = lambda_val,
      objective = compute_huber_loss(y_train, X_train, theta, tau) + lambda_val * sum(abs(theta)),
      stringsAsFactors = FALSE
    )
  }) %>%
    bind_rows() %>%
    filter(is.finite(objective))
  
  if (nrow(candidates) == 0L) {
    return(list(theta = rep(NA_real_, ncol(X_train)), lambda = NA_real_, tau = tau))
  }
  
  chosen <- candidates %>%
    arrange(objective, desc(lambda)) %>%
    slice(1L)
  
  theta_final <- tryCatch(
    solver_huber_lasso(
      X = X_train,
      y = y_train,
      lambda = chosen$lambda[1],
      tau = tau,
      T_max = cfg$huber_T,
      beta_init = lasso_theta
    ),
    error = function(e) rep(NA_real_, ncol(X_train))
  )
  
  list(theta = as.numeric(theta_final), lambda = chosen$lambda[1], tau = tau)
}

fit_shrinkage_fixed_q <- function(X_train, y_train, seed) {
  tau_x <- as.numeric(quantile(abs(X_train), probs = cfg$shrink_tau_x_q, na.rm = TRUE))
  tau_y <- as.numeric(quantile(abs(y_train), probs = cfg$shrink_tau_y_q, na.rm = TRUE))
  tau_x <- max(tau_x, 1e-8)
  tau_y <- max(tau_y, 1e-8)
  
  X_shrunk <- truncate_operator(X_train, tau_x)
  y_shrunk <- truncate_operator(y_train, tau_y)
  
  fit <- tryCatch(
    with_seed_preserved(seed, {
      cv.glmnet(
        X_shrunk,
        y_shrunk,
        standardize = FALSE,
        intercept = FALSE,
        nfolds = cfg$lasso_nfolds
      )
    }),
    error = function(e) NULL
  )
  
  if (is.null(fit)) {
    return(list(theta = rep(NA_real_, ncol(X_train)), lambda = NA_real_, tau_x = tau_x, tau_y = tau_y))
  }
  
  theta <- tryCatch(
    solver_shrinkage_method(X_train, y_train, tau_x, tau_y, fit$lambda.min),
    error = function(e) rep(NA_real_, ncol(X_train))
  )
  
  list(theta = as.numeric(theta), lambda = fit$lambda.min, tau_x = tau_x, tau_y = tau_y)
}

compute_ts_budget <- function(n, p, s, q, T2, k1_cap, m_min = 8L, s_ref = s) {
  n2 <- floor(q * n)
  n1 <- n - n2
  
  K1_tar <- ceiling(s_ref * log(p))
  K1_default <- min(K1_tar, floor(n1 / m_min))
  K1 <- min(K1_default, as.integer(k1_cap))
  
  b2 <- floor(n2 / T2)
  K2_tar <- ceiling(log(p))
  K2 <- min(K2_tar, floor(b2 / m_min), 2L)
  
  eligible <- (n1 >= m_min) && (b2 >= m_min) && (K1 >= 1L) && (K2 >= 1L)
  
  list(
    n1 = n1, n2 = n2, b2 = b2,
    K1 = K1, K1_default = K1_default, K1_tar = K1_tar,
    K2 = K2, K2_tar = K2_tar,
    stage1_block_size = if (K1 > 0L) n1 / K1 else NA_real_,
    stage2_block_size = if (K2 > 0L) b2 / K2 else NA_real_,
    eligible = eligible
  )
}

make_ts_split <- function(n, budget, T2, seed) {
  with_seed_preserved(seed, {
    idx <- sample.int(n)
    id1 <- idx[seq_len(budget$n1)]
    id2_pool <- idx[(budget$n1 + 1L):n]
    use_n2 <- budget$b2 * T2
    id2 <- sample(id2_pool, size = length(id2_pool), replace = FALSE)[seq_len(use_n2)]
    list(id1 = id1, id2 = id2)
  })
}

run_stage2 <- function(theta, X_train, y_train, split, budget, eta, s, T2) {
  theta_out <- theta
  for (tt in seq_len(T2)) {
    batch_ids <- split$id2[((tt - 1L) * budget$b2 + 1L):(tt * budget$b2)]
    theta_out <- hard_threshold(
      theta_out - eta * robust_grad_MoM(theta_out, X_train[batch_ids, , drop = FALSE], y_train[batch_ids], K = budget$K2),
      s
    )
  }
  theta_out
}

make_metric_row <- function(method, config, repeat_id, outer_seed, fold, theta, X_val, y_val,
                            eta = NA_real_, eta_scale = NA_real_, q = NA_real_,
                            k1_cap = NA_integer_, T = NA_integer_, T1 = NA_integer_,
                            T2 = NA_integer_, lambda = NA_real_, s = NA_integer_,
                            tau_q = NA_real_, tau_x_q = NA_real_, tau_y_q = NA_real_,
                            tau = NA_real_, tau_x = NA_real_, tau_y = NA_real_,
                            budget = NULL, eta_base = NA_real_, y_center = NA_real_,
                            fold_seed = NA_integer_, lasso_seed = NA_integer_,
                            ts_seed = NA_integer_) {
  if (is.null(budget)) {
    budget <- list(
      n1 = NA_integer_, n2 = NA_integer_, b2 = NA_integer_,
      K1 = NA_integer_, K1_default = NA_integer_, K1_tar = NA_integer_,
      K2 = NA_integer_, K2_tar = NA_integer_,
      stage1_block_size = NA_real_, stage2_block_size = NA_real_,
      eligible = TRUE
    )
  }
  
  data.frame(
    method = method, config = config, repeat_id = repeat_id, outer_seed = outer_seed, fold = fold,
    rmse = safe_rmse(theta, X_val, y_val),
    mae = safe_mae(theta, X_val, y_val),
    r2 = safe_r2(theta, X_val, y_val),
    eta = eta, eta_scale = eta_scale, q = q,
    k1_cap = k1_cap, T = T, T1 = T1, T2 = T2,
    lambda = lambda, s = s,
    tau_q = tau_q, tau_x_q = tau_x_q, tau_y_q = tau_y_q,
    tau = tau, tau_x = tau_x, tau_y = tau_y,
    nonzero = if (any(!is.finite(theta))) NA_integer_ else sum(theta != 0),
    n1 = budget$n1, n2 = budget$n2, b2 = budget$b2,
    K1 = budget$K1, K1_default = budget$K1_default, K1_tar = budget$K1_tar,
    K2 = budget$K2, K2_tar = budget$K2_tar,
    stage1_block_size = budget$stage1_block_size,
    stage2_block_size = budget$stage2_block_size,
    eligible = budget$eligible,
    eta_base = eta_base, y_center = y_center,
    fold_seed = fold_seed, lasso_seed = lasso_seed, ts_seed = ts_seed,
    stringsAsFactors = FALSE
  )
}

run_ts_candidate <- function(k1_cap, repeat_id, outer_seed, fold, X_train, y_train, X_val, y_val,
                             s_fold, eta, eta_base, y_center, fold_seed, lasso_seed, ts_seed) {
  budget <- compute_ts_budget(nrow(X_train), ncol(X_train), s_fold, cfg$q, cfg$T2, k1_cap, cfg$m_min, s_fold)
  if (!budget$eligible) stop(sprintf("TS budget is not eligible in repeat %d fold %d K1 cap %d.", repeat_id, fold, k1_cap))
  
  split <- make_ts_split(nrow(X_train), budget, cfg$T2, seed = ts_seed)
  theta <- rep(0, ncol(X_train))
  trajectory_rows <- list()
  comparison_row <- NULL
  
  add_rows <- function(T1, theta_stage1) {
    theta_stage2 <- run_stage2(theta_stage1, X_train, y_train, split, budget, eta, s_fold, cfg$T2)
    
    trajectory_rows[[length(trajectory_rows) + 1L]] <<- cbind(
      stage = "stage1_only",
      make_metric_row(
        "TS-RIGHT-TUNED", "TS tuned K1cap", repeat_id, outer_seed, fold,
        theta_stage1, X_val, y_val,
        eta = eta, eta_scale = cfg$eta_scale, q = cfg$q,
        k1_cap = k1_cap, T1 = T1, T2 = cfg$T2, s = s_fold,
        budget = budget, eta_base = eta_base, y_center = y_center,
        fold_seed = fold_seed, lasso_seed = lasso_seed, ts_seed = ts_seed
      )
    )
    
    trajectory_rows[[length(trajectory_rows) + 1L]] <<- cbind(
      stage = "after_stage2",
      make_metric_row(
        "TS-RIGHT-TUNED", "TS tuned K1cap", repeat_id, outer_seed, fold,
        theta_stage2, X_val, y_val,
        eta = eta, eta_scale = cfg$eta_scale, q = cfg$q,
        k1_cap = k1_cap, T1 = T1, T2 = cfg$T2, s = s_fold,
        budget = budget, eta_base = eta_base, y_center = y_center,
        fold_seed = fold_seed, lasso_seed = lasso_seed, ts_seed = ts_seed
      )
    )
    
    if (T1 == cfg$compare_ts_T1) {
      comparison_row <<- make_metric_row(
        "TS-RIGHT-TUNED",
        sprintf("TS tuned K1cap in {%s}, T1=%d, T2=%d, eta=%.2g/L",
                paste(cfg$k1_cap_grid, collapse = ","), cfg$compare_ts_T1, cfg$T2, cfg$eta_scale),
        repeat_id, outer_seed, fold,
        theta_stage2, X_val, y_val,
        eta = eta, eta_scale = cfg$eta_scale, q = cfg$q,
        k1_cap = k1_cap, T1 = T1, T2 = cfg$T2, s = s_fold,
        budget = budget, eta_base = eta_base, y_center = y_center,
        fold_seed = fold_seed, lasso_seed = lasso_seed, ts_seed = ts_seed
      )
    }
  }
  
  add_rows(0L, theta)
  
  for (tt in seq_len(max(cfg$t1_grid))) {
    theta <- hard_threshold(
      theta - eta * robust_grad_MoM(theta, X_train[split$id1, , drop = FALSE], y_train[split$id1], K = budget$K1),
      s_fold
    )
    if (tt %in% cfg$t1_grid) add_rows(tt, theta)
  }
  
  list(trajectory = bind_rows(trajectory_rows), comparison = comparison_row)
}

run_repeat_fold <- function(job) {
  repeat_id <- job$repeat_id
  fold <- job$fold
  outer_seed <- cfg$seed + 100000L * (repeat_id - 1L)
  fold_id <- make_folds(nrow(X_raw), cfg$n_folds, seed = outer_seed)
  
  val_idx <- which(fold_id == fold)
  train_idx <- setdiff(seq_len(nrow(X_raw)), val_idx)
  
  scaler <- fit_robust_x_scaler(X_raw[train_idx, , drop = FALSE])
  X_train <- apply_robust_x_scaler(X_raw[train_idx, , drop = FALSE], scaler)
  X_val <- apply_robust_x_scaler(X_raw[val_idx, , drop = FALSE], scaler)
  
  y_center <- median(y_raw[train_idx], na.rm = TRUE)
  y_train <- y_raw[train_idx] - y_center
  y_val <- y_raw[val_idx] - y_center
  
  eta_base <- estimate_gradient_step(X_train)
  eta <- cfg$eta_scale * eta_base
  fold_seed <- outer_seed + 1000L * fold
  lasso_seed <- fold_seed + 1L
  huber_seed <- fold_seed + 31L
  shrink_seed <- fold_seed + 41L
  ts_seed <- fold_seed + 200L + cfg$ts_seed_offset
  
  lasso_fit <- fit_lasso_tuned(X_train, y_train, seed = lasso_seed)
  s_fold <- max(1L, min(as.integer(lasso_fit$nonzero), floor(nrow(X_train) / 2), ncol(X_train)))
  
  lasso_row <- make_metric_row(
    "LASSO", "LASSO cv.glmnet lambda.min", repeat_id, outer_seed, fold,
    lasso_fit$theta, X_val, y_val,
    lambda = lasso_fit$lambda, s = s_fold,
    eta_base = eta_base, y_center = y_center,
    fold_seed = fold_seed, lasso_seed = lasso_seed
  )
  
  iht_theta <- fit_iht(X_train, y_train, s_fold, eta, cfg$compare_iht_T)
  iht_row <- make_metric_row(
    "IHT", sprintf("IHT eta=%.2g/L, T=%d", cfg$eta_scale, cfg$compare_iht_T),
    repeat_id, outer_seed, fold,
    iht_theta, X_val, y_val,
    eta = eta, eta_scale = cfg$eta_scale, T = cfg$compare_iht_T, s = s_fold,
    eta_base = eta_base, y_center = y_center,
    fold_seed = fold_seed, lasso_seed = lasso_seed
  )
  
  huber_fit <- fit_huber_fixed_q(X_train, y_train, lasso_fit$theta, lasso_fit$lambda_grid)
  huber_row <- make_metric_row(
    "Adaptive Huber",
    sprintf("Adaptive Huber tau_q=%.2f", cfg$huber_tau_q),
    repeat_id, outer_seed, fold,
    huber_fit$theta, X_val, y_val,
    lambda = huber_fit$lambda, s = s_fold,
    tau_q = cfg$huber_tau_q, tau = huber_fit$tau,
    eta_base = eta_base, y_center = y_center,
    fold_seed = fold_seed, lasso_seed = huber_seed
  )
  
  shrink_fit <- fit_shrinkage_fixed_q(X_train, y_train, seed = shrink_seed)
  shrink_row <- make_metric_row(
    "Shrinkage",
    sprintf("Shrinkage tau_x_q=%.2f, tau_y_q=%.2f", cfg$shrink_tau_x_q, cfg$shrink_tau_y_q),
    repeat_id, outer_seed, fold,
    shrink_fit$theta, X_val, y_val,
    lambda = shrink_fit$lambda, s = s_fold,
    tau_x_q = cfg$shrink_tau_x_q, tau_y_q = cfg$shrink_tau_y_q,
    tau_x = shrink_fit$tau_x, tau_y = shrink_fit$tau_y,
    eta_base = eta_base, y_center = y_center,
    fold_seed = fold_seed, lasso_seed = shrink_seed
  )
  
  ts_outputs <- lapply(cfg$k1_cap_grid, function(k1_cap) {
    run_ts_candidate(
      k1_cap, repeat_id, outer_seed, fold,
      X_train, y_train, X_val, y_val,
      s_fold, eta, eta_base, y_center,
      fold_seed, lasso_seed, ts_seed
    )
  })
  
  ts_trajectory_candidates <- bind_rows(lapply(ts_outputs, `[[`, "trajectory"))
  ts_comparison_candidates <- bind_rows(lapply(ts_outputs, `[[`, "comparison"))
  
  tuned_comparison <- ts_comparison_candidates %>%
    arrange(rmse, k1_cap) %>%
    slice(1L)
  
  tuned_trajectory <- ts_trajectory_candidates %>%
    group_by(repeat_id, fold, stage, T1) %>%
    arrange(rmse, k1_cap, .by_group = TRUE) %>%
    slice(1L) %>%
    ungroup() %>%
    mutate(
      config = sprintf("TS tuned K1cap in {%s}, T1=%d, T2=%d, eta=%.2g/L",
                       paste(cfg$k1_cap_grid, collapse = ","), T1, T2, cfg$eta_scale)
    )
  
  list(
    trajectory = tuned_trajectory,
    comparison = bind_rows(lasso_row, iht_row, huber_row, shrink_row, tuned_comparison)
  )
}

jobs <- expand.grid(
  repeat_id = seq_len(cfg$n_repeats),
  fold = seq_len(cfg$n_folds),
  KEEP.OUT.ATTRS = FALSE
)
jobs <- split(jobs, seq_len(nrow(jobs)))

cat(sprintf("Loaded %s / %s: n = %d, p = %d.\n", task$drug_key, task$drug_name, nrow(X_raw), ncol(X_raw)))
cat(sprintf("Running %d repeat-fold jobs on %d worker(s).\n", length(jobs), cfg$n_workers))
cat("Outer split seeds:\n")
print(cfg$seed + 100000L * (seq_len(cfg$n_repeats) - 1L))

if (cfg$n_workers <= 1L) {
  all_results <- lapply(jobs, run_repeat_fold)
} else {
  cl <- parallel::makeCluster(cfg$n_workers)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  
  parallel::clusterEvalQ(cl, {
    suppressPackageStartupMessages({
      library(dplyr)
      library(glmnet)
    })
    NULL
  })
  
  parallel::clusterExport(
    cl,
    varlist = c(
      "cfg", "X_raw", "y_raw",
      "with_seed_preserved", "make_folds",
      "fit_robust_x_scaler", "apply_robust_x_scaler",
      "estimate_gradient_step", "hard_threshold",
      "grad_linear_regression", "grad_linear_regression_samplewise",
      "robust_grad_MoM", "safe_rmse", "safe_mae", "safe_r2",
      "compute_huber_loss", "grad_huber_loss", "soft_threshold",
      "truncate_operator", "solver_huber_lasso", "solver_shrinkage_method",
      "fit_lasso_tuned", "fit_iht", "fit_huber_fixed_q", "fit_shrinkage_fixed_q",
      "compute_ts_budget", "make_ts_split", "run_stage2",
      "make_metric_row", "run_ts_candidate", "run_repeat_fold"
    ),
    envir = environment()
  )
  
  all_results <- parallel::parLapply(cl, jobs, run_repeat_fold)
}

raw_trajectory <- bind_rows(lapply(all_results, `[[`, "trajectory")) %>%
  arrange(repeat_id, fold, stage, T1)

raw_comparison <- bind_rows(lapply(all_results, `[[`, "comparison")) %>%
  arrange(repeat_id, fold, method)

summary_trajectory <- raw_trajectory %>%
  group_by(method, config, stage, T1, T2, eta_scale, q) %>%
  summarise(
    median_rmse = median(rmse, na.rm = TRUE),
    mean_rmse = mean(rmse, na.rm = TRUE),
    sd_rmse = sd(rmse, na.rm = TRUE),
    iqr_rmse = IQR(rmse, na.rm = TRUE),
    min_rmse = min(rmse, na.rm = TRUE),
    max_rmse = max(rmse, na.rm = TRUE),
    median_mae = median(mae, na.rm = TRUE),
    mean_mae = mean(mae, na.rm = TRUE),
    median_r2 = median(r2, na.rm = TRUE),
    mean_r2 = mean(r2, na.rm = TRUE),
    fail_rate = mean(!is.finite(rmse) | is.na(rmse)),
    n_repeats = dplyr::n_distinct(repeat_id),
    n_folds_total = dplyr::n(),
    median_s = median(s, na.rm = TRUE),
    median_nonzero = median(nonzero, na.rm = TRUE),
    tuned_k1_counts = paste(paste0("K1=", names(table(k1_cap)), "(", as.numeric(table(k1_cap)), " times)"), collapse = " | "),
    K1_values = paste(sort(unique(K1[is.finite(K1)])), collapse = ","),
    K2_values = paste(sort(unique(K2[is.finite(K2)])), collapse = ","),
    .groups = "drop"
  ) %>%
  arrange(method, stage, median_rmse, T1)

summary_comparison <- raw_comparison %>%
  group_by(method, config, eta_scale, q, T, T1, T2, tau_q, tau_x_q, tau_y_q) %>%
  summarise(
    median_rmse = median(rmse, na.rm = TRUE),
    mean_rmse = mean(rmse, na.rm = TRUE),
    sd_rmse = sd(rmse, na.rm = TRUE),
    iqr_rmse = IQR(rmse, na.rm = TRUE),
    min_rmse = min(rmse, na.rm = TRUE),
    max_rmse = max(rmse, na.rm = TRUE),
    median_mae = median(mae, na.rm = TRUE),
    mean_mae = mean(mae, na.rm = TRUE),
    median_r2 = median(r2, na.rm = TRUE),
    mean_r2 = mean(r2, na.rm = TRUE),
    fail_rate = mean(!is.finite(rmse) | is.na(rmse)),
    n_repeats = dplyr::n_distinct(repeat_id),
    n_folds_total = dplyr::n(),
    median_s = median(s, na.rm = TRUE),
    median_nonzero = median(nonzero, na.rm = TRUE),
    tuned_k1_counts = paste(paste0("K1=", names(table(k1_cap)), "(", as.numeric(table(k1_cap)), " times)"), collapse = " | "),
    lambda_values = paste(signif(sort(unique(lambda[is.finite(lambda)])), 6), collapse = ","),
    K1_values = paste(sort(unique(K1[is.finite(K1)])), collapse = ","),
    K2_values = paste(sort(unique(K2[is.finite(K2)])), collapse = ","),
    tau_values = paste(signif(sort(unique(tau[is.finite(tau)])), 6), collapse = ","),
    tau_x_values = paste(signif(sort(unique(tau_x[is.finite(tau_x)])), 6), collapse = ","),
    tau_y_values = paste(signif(sort(unique(tau_y[is.finite(tau_y)])), 6), collapse = ","),
    .groups = "drop"
  ) %>%
  arrange(median_rmse, mean_rmse)

tuned_k1_counts <- raw_comparison %>%
  filter(method == "TS-RIGHT-TUNED") %>%
  count(k1_cap, name = "n_selected") %>%
  arrange(k1_cap)

best_trajectory_by_stage <- summary_trajectory %>%
  filter(method == "TS-RIGHT-TUNED") %>%
  group_by(stage) %>%
  slice_min(order_by = median_rmse, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(median_rmse)

write.csv(raw_trajectory, file.path(out_dir, "raw", "raw_tuned_ts_trajectory.csv"), row.names = FALSE)
write.csv(raw_comparison, file.path(out_dir, "raw", "raw_lasso_iht_huber_shrinkage_tuned_ts.csv"), row.names = FALSE)
write.csv(summary_trajectory, file.path(out_dir, "summary", "summary_tuned_ts_trajectory.csv"), row.names = FALSE)
write.csv(best_trajectory_by_stage, file.path(out_dir, "summary", "best_tuned_ts_trajectory_by_stage.csv"), row.names = FALSE)
write.csv(summary_comparison, file.path(out_dir, "summary", "summary_lasso_iht_huber_shrinkage_tuned_ts.csv"), row.names = FALSE)
write.csv(tuned_k1_counts, file.path(out_dir, "summary", "tuned_k1_selection_counts.csv"), row.names = FALSE)

plot_trajectory <- summary_trajectory %>%
  filter(method == "TS-RIGHT-TUNED") %>%
  mutate(stage = factor(stage, levels = c("stage1_only", "after_stage2"))) %>%
  ggplot(aes(x = T1, y = median_rmse, color = stage)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 1.9) +
  theme_bw(base_size = 12) +
  labs(
    x = "Stage-1 iterations T1",
    y = "Median held-out RMSE over repeated 5-fold results",
    color = "Evaluation",
    title = "Trajectory of TS-RIGHT-TUNED"
  ) +
  theme(legend.position = "bottom")
plot_trajectory

# 定义自定义调色板
method_palette <- c(
  "TS-RIGHT" = "#332288",
  "IHT" = "#CC6677",
  "Adaptive Huber" = "#3B82A1",
  "Shrinkage" = "#44AA99",
  "Lasso" = "#CAB54B"
)

plot_comparison <- raw_comparison %>%
  # 统一方法名称，确保能正确匹配你提供的调色板
  mutate(
    method = case_when(
      method == "LASSO" ~ "Lasso",
      method == "TS-RIGHT-TUNED" ~ "TS-RIGHT",
      TRUE ~ method
    )
  ) %>%
  # 筛选目标方法
  filter(method %in% c("Lasso", "IHT", "Adaptive Huber", "Shrinkage", "TS-RIGHT")) %>%
  # 核心调整：根据 100 次重复的 Median RMSE 对 method 进行从小到大的排序
  mutate(method = reorder(method, rmse, FUN = median, na.rm = TRUE)) %>%
  # x 映射为 method，配合 coord_flip 变成纵坐标，保证和图例文字一致
  ggplot(aes(x = method, y = rmse, fill = method)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.8) +
  geom_point(position = position_jitter(width = 0.1, height = 0), size = 1.2, alpha = 0.4, color = "black") +
  # 载入你指定的调色板
  scale_fill_manual(values = method_palette) +
  coord_flip() +
  theme_bw(base_size = 12) +
  labs(
    x = NULL,
    y = "Median RMSE of 100 repetitions", # 更新为你要求的横坐标文字
    fill = "Method"
  ) +
  theme(legend.position = "bottom")

plot_comparison

ggsave(file.path(out_dir, "figures", "tuned_ts_median_rmse_trajectory.png"), plot_trajectory, width = 8.8, height = 5.2, dpi = 300)
ggsave(file.path(out_dir, "figures", "tuned_ts_median_rmse_trajectory.pdf"), plot_trajectory, width = 8.8, height = 5.2)
ggsave(file.path(out_dir, "figures", "lasso_iht_huber_shrinkage_tuned_ts_rmse.png"), plot_comparison, width = 10.5, height = 6.2, dpi = 300)
ggsave(file.path(out_dir, "figures", "lasso_iht_huber_shrinkage_tuned_ts_rmse.pdf"), plot_comparison, width = 10.5, height = 6.2)

print(best_trajectory_by_stage)
print(summary_comparison %>% select(method, config, median_rmse, mean_rmse, tuned_k1_counts))
print(tuned_k1_counts)

cat("\nParallel repeated 5-fold comparison with Huber, shrinkage, and tuned TS complete.\n")
cat("Trajectory summary: ", file.path(out_dir, "summary", "summary_tuned_ts_trajectory.csv"), "\n", sep = "")
cat("Comparison summary: ", file.path(out_dir, "summary", "summary_lasso_iht_huber_shrinkage_tuned_ts.csv"), "\n", sep = "")
cat("K1 selection counts: ", file.path(out_dir, "summary", "tuned_k1_selection_counts.csv"), "\n", sep = "")
