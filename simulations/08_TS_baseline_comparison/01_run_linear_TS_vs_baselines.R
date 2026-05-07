# ======================================================================
# Linear simulation: TS-RIGHT vs IHT / Lasso / Huber-Lasso / Shrinkage
#
# This script is intended to live at:
#   simulations/08_TS_baseline_comparison/01_run_linear_TS_vs_baselines.R
#
# It assumes that the updated schedule utility exists at:
#   simulations/07_schedule_ablation/utils_schedule_right.R
# and that this utility has the seed-enabled TS/FS wrappers.
#
# Main tuning protocol:
#   - TS-RIGHT: per-n pilot tuning over q,T1,T2, then frozen for final reps.
#   - IHT: per-n pilot tuning over eta,T, then frozen for final reps.
#   - Lasso: cv.glmnet.
#   - Huber-Lasso: train/validation tuning over lambda and residual quantile tau.
#   - Shrinkage: train/validation tuning over truncation quantile; lambda by cv.glmnet.
#
# Validation metric throughout tuning:
#   median_i | y_i - x_i^T theta |.
# ======================================================================

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(ggplot2)
  library(glmnet)
  library(foreach)
  library(doParallel)
})

source(here("simulations", "07_schedule_ablation", "utils_schedule_right.R"))

# ----------------------------------------------------------------------
# Configuration
# ----------------------------------------------------------------------

cfg <- list(
  experiment_id = "linear_TS_vs_baselines_n_sweep",

  p = 800,
  n_grid = c(800, 1200, 1600, 2000, 2400),

  s_star = 5,
  s = 10,
  m_min = 10,

  reps = 100,
  n_pilot = 15,
  seed_base = 20260507L,

  df_X = 2.5,
  scale_X = 1,
  df_eps = 1.5,
  scale_eps = 1,
  theta_magnitude = 5,

  # TS-RIGHT schedule grid requested by the user.
  q_grid = c(0.50),
  T1_grid = c(150, 200),
  T2_grid = c(8),
  eta_ts = 0.02,

  # IHT grid requested by the user.
  iht_eta_grid = c(0.005, 0.01, 0.02),
  iht_T_grid = c(200, 400, 800),

  # Baseline tuning grids.
  huber_lambda_grid = 10^seq(-3, 0, length.out = 8),
  huber_tau_quantile_grid = c(0.90, 0.95, 0.98),
  shrink_tau_quantile_grid = c(0.90, 0.95, 0.98),

  # Solver iteration lengths for Huber-Lasso. These match the spirit of the old
  # comparison script but can be lowered for smoke tests.
  T_huber_tune = 250,
  T_huber_final = 500,

  val_frac = 0.20,
  cv_nfolds = 5,
  selection_alpha = 1.05,

  # Parallel control.
  use_parallel = TRUE,
  n_cores = max(1L, parallel::detectCores() - 1L)
)

# For smoke tests, temporarily override, e.g.:
# cfg$reps <- 2
# cfg$n_pilot <- 2
# cfg$n_grid <- c(800)

out_dir <- here("results", "08_TS_baseline_comparison", "linear")
raw_dir <- file.path(out_dir, "raw")
summary_dir <- file.path(out_dir, "summary")
fig_dir <- file.path(out_dir, "figures")
registry_dir <- here("results", "00_tuning_registry")
candidate_dir <- file.path(registry_dir, "candidates_linear_baseline")

dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(registry_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(candidate_dir, recursive = TRUE, showWarnings = FALSE)

registry_path <- file.path(registry_dir, "linear_TS_vs_baselines_tuning_registry.csv")

# ----------------------------------------------------------------------
# Small utilities
# ----------------------------------------------------------------------

stable_hash <- function(x) {
  bytes <- utf8ToInt(paste(x, collapse = "|"))
  h <- 0L
  for (b in bytes) {
    h <- as.integer((as.numeric(h) * 131 + b) %% 1000000007)
  }
  sprintf("h%08x", h)
}

make_short_config_id <- function(prefix, signature, n) {
  paste0(prefix, "_n", n, "_", stable_hash(signature))
}

read_registry <- function(path) {
  if (!file.exists(path)) {
    return(data.frame())
  }
  read.csv(path, stringsAsFactors = FALSE)
}

append_registry_row <- function(row, path) {
  old <- read_registry(path)
  out <- dplyr::bind_rows(old, row)
  write.csv(out, path, row.names = FALSE)
  invisible(out)
}

seed_for <- function(base, n_index, rep_id, offset = 0L) {
  as.integer(base + offset + 10000L * n_index + rep_id)
}

with_temp_seed <- function(seed) {
  if (is.null(seed)) return(NULL)
  old_exists <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old_state <- if (old_exists) get(".Random.seed", envir = .GlobalEnv, inherits = FALSE) else NULL
  set.seed(seed)
  list(old_exists = old_exists, old_state = old_state)
}

restore_temp_seed <- function(state) {
  if (is.null(state)) return(invisible(NULL))
  if (state$old_exists) {
    assign(".Random.seed", state$old_state, envir = .GlobalEnv)
  } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    rm(".Random.seed", envir = .GlobalEnv)
  }
  invisible(NULL)
}

make_train_val_split <- function(n, val_frac = 0.2, seed) {
  rng <- with_temp_seed(seed)
  on.exit(restore_temp_seed(rng), add = TRUE)

  n_val <- max(1L, floor(val_frac * n))
  val_idx <- sample.int(n, size = n_val, replace = FALSE)
  train_idx <- setdiff(seq_len(n), val_idx)
  list(train = train_idx, val = val_idx)
}

make_foldid <- function(n, nfolds = 5, seed = 1L) {
  nfolds <- min(as.integer(nfolds), n)
  rng <- with_temp_seed(seed)
  on.exit(restore_temp_seed(rng), add = TRUE)

  foldid <- rep(seq_len(nfolds), length.out = n)
  sample(foldid, size = n, replace = FALSE)
}

median_abs_prediction_error <- function(theta, X_val, y_val) {
  if (any(!is.finite(theta))) return(Inf)
  pred <- as.numeric(X_val %*% theta)
  median(abs(y_val - pred), na.rm = TRUE)
}

safe_l2 <- function(theta, theta_star) {
  if (any(!is.finite(theta))) return(NA_real_)
  sqrt(sum((theta - theta_star)^2))
}

# ----------------------------------------------------------------------
# Fitters that return theta_hat
# ----------------------------------------------------------------------

fit_ts_right_linear_theta <- function(
  X,
  y,
  s,
  eta,
  q,
  T1,
  T2,
  m_min,
  s_ref,
  seed = NULL
) {
  rng <- with_temp_seed(seed)
  on.exit(restore_temp_seed(rng), add = TRUE)

  n <- nrow(X)
  p <- ncol(X)

  budget <- compute_ts_budget(
    n = n,
    p = p,
    s = s,
    q = q,
    T1 = T1,
    T2 = T2,
    s_ref = s_ref,
    m_min = m_min
  )

  if (!budget$eligible) {
    return(list(theta = rep(NA_real_, p), budget = budget, eligible = FALSE))
  }

  idx <- sample.int(n)
  id1 <- idx[seq_len(budget$n1)]
  id2 <- idx[(budget$n1 + 1L):n]

  theta_cur <- solver_right(
    X = X[id1, , drop = FALSE],
    y = y[id1],
    s = s,
    eta = eta,
    T_max = T1,
    K = budget$K1,
    theta_init = rep(0, p),
    grad_func_samplewise = grad_linear_regression_samplewise,
    record_trace = FALSE
  )

  use_n2 <- budget$b2 * T2
  id2 <- sample(id2, size = length(id2), replace = FALSE)[seq_len(use_n2)]

  for (tt in seq_len(T2)) {
    batch_ids <- id2[((tt - 1L) * budget$b2 + 1L):(tt * budget$b2)]
    theta_cur <- solver_right(
      X = X[batch_ids, , drop = FALSE],
      y = y[batch_ids],
      s = s,
      eta = eta,
      T_max = 1L,
      K = budget$K2,
      theta_init = theta_cur,
      grad_func_samplewise = grad_linear_regression_samplewise,
      record_trace = FALSE
    )
  }

  list(theta = as.numeric(theta_cur), budget = budget, eligible = TRUE)
}

fit_iht_linear_theta <- function(X, y, s, eta, T) {
  as.numeric(solver_iht(
    X = X,
    y = y,
    s = s,
    eta = eta,
    T_max = T,
    grad_func = grad_linear_regression
  ))
}

fit_lasso_cv <- function(X, y, seed, nfolds = 5) {
  foldid <- make_foldid(nrow(X), nfolds = nfolds, seed = seed)
  cv_fit <- cv.glmnet(
    X,
    y,
    standardize = FALSE,
    intercept = FALSE,
    foldid = foldid
  )
  lambda_best <- cv_fit$lambda.min
  final_fit <- glmnet(
    X,
    y,
    lambda = lambda_best,
    standardize = FALSE,
    intercept = FALSE
  )
  list(
    theta = as.numeric(final_fit$beta),
    lambda = as.numeric(lambda_best)
  )
}

fit_huber_lasso_tuned <- function(
  X_full,
  y_full,
  train_idx,
  lambda_grid,
  tau_quantile_grid,
  val_alpha = 1.05,
  T_tune = 250,
  T_final = 500,
  seed = 1L,
  nfolds = 5
) {
  X_train <- X_full[train_idx, , drop = FALSE]
  y_train <- y_full[train_idx]
  val_idx <- setdiff(seq_len(nrow(X_full)), train_idx)
  X_val <- X_full[val_idx, , drop = FALSE]
  y_val <- y_full[val_idx]

  # Lasso initialization for residual-adaptive Huber thresholds.
  init_train <- fit_lasso_cv(X_train, y_train, seed = seed + 11L, nfolds = nfolds)
  beta_init_train <- init_train$theta
  resid_train <- as.numeric(y_train - X_train %*% beta_init_train)

  cand <- list()
  cc <- 0L

  for (tau_q in tau_quantile_grid) {
    tau_val <- as.numeric(quantile(abs(resid_train), probs = tau_q, na.rm = TRUE))
    tau_val <- max(tau_val, 1e-8)

    for (lambda_val in lambda_grid) {
      cc <- cc + 1L
      theta_tmp <- tryCatch(
        solver_huber_lasso(
          X = X_train,
          y = y_train,
          lambda = lambda_val,
          tau = tau_val,
          T_max = T_tune,
          beta_init = beta_init_train
        ),
        error = function(e) rep(NA_real_, ncol(X_full))
      )

      cand[[cc]] <- data.frame(
        tau_q = tau_q,
        tau = tau_val,
        lambda = lambda_val,
        val_loss = median_abs_prediction_error(theta_tmp, X_val, y_val)
      )
    }
  }

  cand_df <- dplyr::bind_rows(cand)
  cand_df <- cand_df[is.finite(cand_df$val_loss), , drop = FALSE]

  if (nrow(cand_df) == 0L) {
    return(list(theta = rep(NA_real_, ncol(X_full)), selected = data.frame()))
  }

  best <- min(cand_df$val_loss, na.rm = TRUE)
  near <- cand_df[cand_df$val_loss <= val_alpha * best, , drop = FALSE]
  # Simpler model preference: choose the largest lambda among near-best candidates;
  # then choose tau_q closest to 0.95.
  near$dist_tau_center <- abs(near$tau_q - 0.95)
  near <- near[order(-near$lambda, near$dist_tau_center, near$val_loss), , drop = FALSE]
  chosen <- near[1, , drop = FALSE]

  init_full <- fit_lasso_cv(X_full, y_full, seed = seed + 17L, nfolds = nfolds)
  beta_init_full <- init_full$theta
  resid_full <- as.numeric(y_full - X_full %*% beta_init_full)
  tau_full <- as.numeric(quantile(abs(resid_full), probs = chosen$tau_q, na.rm = TRUE))
  tau_full <- max(tau_full, 1e-8)

  theta_final <- tryCatch(
    solver_huber_lasso(
      X = X_full,
      y = y_full,
      lambda = chosen$lambda,
      tau = tau_full,
      T_max = T_final,
      beta_init = beta_init_full
    ),
    error = function(e) rep(NA_real_, ncol(X_full))
  )

  chosen$tau_final <- tau_full
  list(theta = as.numeric(theta_final), selected = chosen)
}

fit_shrinkage_tuned <- function(
  X_full,
  y_full,
  train_idx,
  tau_quantile_grid,
  val_alpha = 1.05,
  seed = 1L,
  nfolds = 5
) {
  X_train <- X_full[train_idx, , drop = FALSE]
  y_train <- y_full[train_idx]
  val_idx <- setdiff(seq_len(nrow(X_full)), train_idx)
  X_val <- X_full[val_idx, , drop = FALSE]
  y_val <- y_full[val_idx]

  cand <- list()
  cc <- 0L

  for (tau_q in tau_quantile_grid) {
    tau_x <- as.numeric(quantile(abs(X_train), probs = tau_q, na.rm = TRUE))
    tau_y <- as.numeric(quantile(abs(y_train), probs = tau_q, na.rm = TRUE))
    tau_x <- max(tau_x, 1e-8)
    tau_y <- max(tau_y, 1e-8)

    X_train_shrunk <- truncate_operator(X_train, tau_x)
    y_train_shrunk <- truncate_operator(y_train, tau_y)

    theta_tmp <- tryCatch({
      foldid <- make_foldid(nrow(X_train), nfolds = nfolds, seed = seed + as.integer(1000 * tau_q))
      cv_fit <- cv.glmnet(
        X_train_shrunk,
        y_train_shrunk,
        standardize = FALSE,
        intercept = FALSE,
        foldid = foldid
      )
      lambda_best <- cv_fit$lambda.min
      theta_fit <- solver_shrinkage_method(X_train, y_train, tau_x, tau_y, lambda_best)
      list(theta = as.numeric(theta_fit), lambda = as.numeric(lambda_best))
    }, error = function(e) {
      list(theta = rep(NA_real_, ncol(X_full)), lambda = NA_real_)
    })

    cc <- cc + 1L
    cand[[cc]] <- data.frame(
      tau_q = tau_q,
      tau_x = tau_x,
      tau_y = tau_y,
      lambda = theta_tmp$lambda,
      val_loss = median_abs_prediction_error(theta_tmp$theta, X_val, y_val)
    )
  }

  cand_df <- dplyr::bind_rows(cand)
  cand_df <- cand_df[is.finite(cand_df$val_loss), , drop = FALSE]

  if (nrow(cand_df) == 0L) {
    return(list(theta = rep(NA_real_, ncol(X_full)), selected = data.frame()))
  }

  best <- min(cand_df$val_loss, na.rm = TRUE)
  near <- cand_df[cand_df$val_loss <= val_alpha * best, , drop = FALSE]
  near$dist_tau_center <- abs(near$tau_q - 0.95)
  near <- near[order(near$dist_tau_center, near$val_loss), , drop = FALSE]
  chosen <- near[1, , drop = FALSE]

  tau_x_full <- as.numeric(quantile(abs(X_full), probs = chosen$tau_q, na.rm = TRUE))
  tau_y_full <- as.numeric(quantile(abs(y_full), probs = chosen$tau_q, na.rm = TRUE))
  tau_x_full <- max(tau_x_full, 1e-8)
  tau_y_full <- max(tau_y_full, 1e-8)

  theta_final <- tryCatch(
    solver_shrinkage_method(X_full, y_full, tau_x_full, tau_y_full, chosen$lambda),
    error = function(e) rep(NA_real_, ncol(X_full))
  )

  chosen$tau_x_final <- tau_x_full
  chosen$tau_y_final <- tau_y_full
  list(theta = as.numeric(theta_final), selected = chosen)
}

# ----------------------------------------------------------------------
# Per-n pilot tuning for TS-RIGHT and IHT
# ----------------------------------------------------------------------

make_ts_signature <- function(n, cfg) {
  paste(
    cfg$experiment_id,
    "method=TS-RIGHT",
    paste0("n=", n),
    paste0("p=", cfg$p),
    paste0("s_star=", cfg$s_star),
    paste0("s=", cfg$s),
    paste0("m_min=", cfg$m_min),
    paste0("eta_ts=", cfg$eta_ts),
    paste0("df_X=", cfg$df_X),
    paste0("scale_X=", cfg$scale_X),
    paste0("df_eps=", cfg$df_eps),
    paste0("theta_magnitude=", cfg$theta_magnitude),
    paste0("n_pilot=", cfg$n_pilot),
    paste0("q_grid=", paste(cfg$q_grid, collapse = ",")),
    paste0("T1_grid=", paste(cfg$T1_grid, collapse = ",")),
    paste0("T2_grid=", paste(cfg$T2_grid, collapse = ",")),
    paste0("metric=median_abs_prediction_error"),
    paste0("alpha=", cfg$selection_alpha),
    sep = "__"
  )
}

make_iht_signature <- function(n, cfg) {
  paste(
    cfg$experiment_id,
    "method=IHT",
    paste0("n=", n),
    paste0("p=", cfg$p),
    paste0("s_star=", cfg$s_star),
    paste0("s=", cfg$s),
    paste0("df_X=", cfg$df_X),
    paste0("scale_X=", cfg$scale_X),
    paste0("df_eps=", cfg$df_eps),
    paste0("theta_magnitude=", cfg$theta_magnitude),
    paste0("n_pilot=", cfg$n_pilot),
    paste0("eta_grid=", paste(cfg$iht_eta_grid, collapse = ",")),
    paste0("T_grid=", paste(cfg$iht_T_grid, collapse = ",")),
    paste0("metric=median_abs_prediction_error"),
    paste0("alpha=", cfg$selection_alpha),
    sep = "__"
  )
}

tune_ts_for_n <- function(n, n_index, theta_star, cfg, registry_path, candidate_dir) {
  signature <- make_ts_signature(n, cfg)
  config_id <- make_short_config_id("TS", signature, n)

  registry <- read_registry(registry_path)
  if (nrow(registry) > 0L) {
    hit <- registry[registry$config_id == config_id & registry$method == "TS-RIGHT", , drop = FALSE]
    if (nrow(hit) > 0L) {
      message("Using frozen TS schedule for n = ", n, ": ", config_id)
      return(hit[1, , drop = FALSE])
    }
  }

  message("Running pilot TS tuning for n = ", n, ": ", config_id)

  grid <- make_ts_grid(
    q_grid = cfg$q_grid,
    T1_grid = cfg$T1_grid,
    T2_grid = cfg$T2_grid
  )

  # Budget diagnostics on the full n. Actual pilot fitting uses the training split.
  budget_diag <- lapply(seq_len(nrow(grid)), function(ii) {
    b <- compute_ts_budget(
      n = n,
      p = cfg$p,
      s = cfg$s,
      q = grid$q[ii],
      T1 = grid$T1[ii],
      T2 = grid$T2[ii],
      s_ref = cfg$s_star,
      m_min = cfg$m_min
    )
    data.frame(
      q = grid$q[ii],
      T1 = grid$T1[ii],
      T2 = grid$T2[ii],
      n1 = b$n1,
      n2 = b$n2,
      b2 = b$b2,
      K1 = b$K1,
      K2 = b$K2,
      T_fd_budget = b$T_fd_budget,
      eligible = b$eligible
    )
  }) %>% dplyr::bind_rows()

  write.csv(
    budget_diag,
    file.path(candidate_dir, paste0("ts_budget_diag_", config_id, ".csv")),
    row.names = FALSE
  )

  pilot_rows <- list()
  rr <- 0L

  for (ii in seq_len(nrow(grid))) {
    q_i <- grid$q[ii]
    T1_i <- grid$T1[ii]
    T2_i <- grid$T2[ii]

    for (pilot_id in seq_len(cfg$n_pilot)) {
      seed_data <- seed_for(cfg$seed_base, n_index, pilot_id, offset = 100000L)
      seed_split <- seed_for(cfg$seed_base, n_index, pilot_id, offset = 200000L)
      seed_algo <- seed_for(cfg$seed_base, n_index, pilot_id, offset = 300000L)

      dat <- gen_one_linear_dataset(
        seed = seed_data,
        n = n,
        p = cfg$p,
        theta_star = theta_star,
        df_X = cfg$df_X,
        Sigma_X = diag(cfg$p),
        scale_X = cfg$scale_X,
        df_eps = cfg$df_eps,
        scale_eps = cfg$scale_eps
      )

      split <- make_train_val_split(n = n, val_frac = cfg$val_frac, seed = seed_split)
      X_train <- dat$X[split$train, , drop = FALSE]
      y_train <- dat$y[split$train]
      X_val <- dat$X[split$val, , drop = FALSE]
      y_val <- dat$y[split$val]

      fit <- fit_ts_right_linear_theta(
        X = X_train,
        y = y_train,
        s = cfg$s,
        eta = cfg$eta_ts,
        q = q_i,
        T1 = T1_i,
        T2 = T2_i,
        m_min = cfg$m_min,
        s_ref = cfg$s_star,
        seed = seed_algo
      )

      rr <- rr + 1L
      pilot_rows[[rr]] <- data.frame(
        config_id = config_id,
        n = n,
        pilot_id = pilot_id,
        q = q_i,
        T1 = T1_i,
        T2 = T2_i,
        eligible = fit$eligible,
        val_loss = median_abs_prediction_error(fit$theta, X_val, y_val),
        seed_data = seed_data,
        seed_split = seed_split,
        seed_algo = seed_algo
      )
    }
  }

  pilot_raw <- dplyr::bind_rows(pilot_rows)
  write.csv(
    pilot_raw,
    file.path(candidate_dir, paste0("ts_candidate_raw_", config_id, ".csv")),
    row.names = FALSE
  )

  pilot_summary <- pilot_raw %>%
    group_by(q, T1, T2) %>%
    summarise(
      pilot_median = median(val_loss, na.rm = TRUE),
      pilot_iqr = IQR(val_loss, na.rm = TRUE),
      fail_rate = mean(!is.finite(val_loss) | !eligible),
      .groups = "drop"
    ) %>%
    left_join(budget_diag, by = c("q", "T1", "T2"))

  write.csv(
    pilot_summary,
    file.path(candidate_dir, paste0("ts_candidate_summary_", config_id, ".csv")),
    row.names = FALSE
  )

  valid <- pilot_summary[is.finite(pilot_summary$pilot_median) & pilot_summary$eligible, , drop = FALSE]
  if (nrow(valid) == 0L) stop("No valid TS candidate for n = ", n)

  best <- min(valid$pilot_median, na.rm = TRUE)
  near <- valid[valid$pilot_median <= cfg$selection_alpha * best, , drop = FALSE]
  near$dist_center <- abs(near$q - 0.5) + abs(near$T1 - 150) / 150 + abs(near$T2 - 8) / 8
  near <- near[order(near$T_fd_budget, near$dist_center, near$pilot_median), , drop = FALSE]
  chosen <- near[1, , drop = FALSE]

  row <- data.frame(
    config_id = config_id,
    config_signature = signature,
    experiment_id = cfg$experiment_id,
    method = "TS-RIGHT",
    model = "linear",
    n = n,
    p = cfg$p,
    s_star = cfg$s_star,
    s = cfg$s,
    s_ref = cfg$s_star,
    df_X = cfg$df_X,
    scale_X = cfg$scale_X,
    df_eps = cfg$df_eps,
    theta_magnitude = cfg$theta_magnitude,
    validation_metric = "median_abs_prediction_error",
    n_pilot = cfg$n_pilot,
    selection_alpha = cfg$selection_alpha,
    eta_ts = cfg$eta_ts,
    q_selected = chosen$q,
    T1_selected = chosen$T1,
    T2_selected = chosen$T2,
    T_fd_budget = chosen$T_fd_budget,
    n1 = chosen$n1,
    n2 = chosen$n2,
    b2 = chosen$b2,
    K1 = chosen$K1,
    K2 = chosen$K2,
    pilot_median = chosen$pilot_median,
    pilot_iqr = chosen$pilot_iqr,
    fail_rate = chosen$fail_rate,
    selection_rule = "within_1.05_best_then_min_budget",
    seed_base = cfg$seed_base
  )

  append_registry_row(row, registry_path)
  row
}

tune_iht_for_n <- function(n, n_index, theta_star, cfg, registry_path, candidate_dir) {
  signature <- make_iht_signature(n, cfg)
  config_id <- make_short_config_id("IHT", signature, n)

  registry <- read_registry(registry_path)
  if (nrow(registry) > 0L) {
    hit <- registry[registry$config_id == config_id & registry$method == "IHT", , drop = FALSE]
    if (nrow(hit) > 0L) {
      message("Using frozen IHT schedule for n = ", n, ": ", config_id)
      return(hit[1, , drop = FALSE])
    }
  }

  message("Running pilot IHT tuning for n = ", n, ": ", config_id)

  grid <- expand.grid(
    eta = cfg$iht_eta_grid,
    T = cfg$iht_T_grid,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )

  pilot_rows <- list()
  rr <- 0L

  for (ii in seq_len(nrow(grid))) {
    eta_i <- grid$eta[ii]
    T_i <- grid$T[ii]

    for (pilot_id in seq_len(cfg$n_pilot)) {
      seed_data <- seed_for(cfg$seed_base, n_index, pilot_id, offset = 400000L)
      seed_split <- seed_for(cfg$seed_base, n_index, pilot_id, offset = 500000L)

      dat <- gen_one_linear_dataset(
        seed = seed_data,
        n = n,
        p = cfg$p,
        theta_star = theta_star,
        df_X = cfg$df_X,
        Sigma_X = diag(cfg$p),
        scale_X = cfg$scale_X,
        df_eps = cfg$df_eps,
        scale_eps = cfg$scale_eps
      )

      split <- make_train_val_split(n = n, val_frac = cfg$val_frac, seed = seed_split)
      X_train <- dat$X[split$train, , drop = FALSE]
      y_train <- dat$y[split$train]
      X_val <- dat$X[split$val, , drop = FALSE]
      y_val <- dat$y[split$val]

      theta_hat <- tryCatch(
        fit_iht_linear_theta(X_train, y_train, s = cfg$s, eta = eta_i, T = T_i),
        error = function(e) rep(NA_real_, cfg$p)
      )

      rr <- rr + 1L
      pilot_rows[[rr]] <- data.frame(
        config_id = config_id,
        n = n,
        pilot_id = pilot_id,
        eta = eta_i,
        T = T_i,
        val_loss = median_abs_prediction_error(theta_hat, X_val, y_val),
        seed_data = seed_data,
        seed_split = seed_split
      )
    }
  }

  pilot_raw <- dplyr::bind_rows(pilot_rows)
  write.csv(
    pilot_raw,
    file.path(candidate_dir, paste0("iht_candidate_raw_", config_id, ".csv")),
    row.names = FALSE
  )

  pilot_summary <- pilot_raw %>%
    group_by(eta, T) %>%
    summarise(
      pilot_median = median(val_loss, na.rm = TRUE),
      pilot_iqr = IQR(val_loss, na.rm = TRUE),
      fail_rate = mean(!is.finite(val_loss)),
      .groups = "drop"
    )

  write.csv(
    pilot_summary,
    file.path(candidate_dir, paste0("iht_candidate_summary_", config_id, ".csv")),
    row.names = FALSE
  )

  valid <- pilot_summary[is.finite(pilot_summary$pilot_median), , drop = FALSE]
  if (nrow(valid) == 0L) stop("No valid IHT candidate for n = ", n)

  best <- min(valid$pilot_median, na.rm = TRUE)
  near <- valid[valid$pilot_median <= cfg$selection_alpha * best, , drop = FALSE]
  near$dist_eta_center <- abs(log(near$eta) - log(0.01))
  near <- near[order(near$T, near$dist_eta_center, near$pilot_median), , drop = FALSE]
  chosen <- near[1, , drop = FALSE]

  row <- data.frame(
    config_id = config_id,
    config_signature = signature,
    experiment_id = cfg$experiment_id,
    method = "IHT",
    model = "linear",
    n = n,
    p = cfg$p,
    s_star = cfg$s_star,
    s = cfg$s,
    df_X = cfg$df_X,
    scale_X = cfg$scale_X,
    df_eps = cfg$df_eps,
    theta_magnitude = cfg$theta_magnitude,
    validation_metric = "median_abs_prediction_error",
    n_pilot = cfg$n_pilot,
    selection_alpha = cfg$selection_alpha,
    eta_selected = chosen$eta,
    T_selected = chosen$T,
    pilot_median = chosen$pilot_median,
    pilot_iqr = chosen$pilot_iqr,
    fail_rate = chosen$fail_rate,
    selection_rule = "within_1.05_best_then_min_T",
    seed_base = cfg$seed_base
  )

  append_registry_row(row, registry_path)
  row
}

# ----------------------------------------------------------------------
# Stage A: pilot tuning registry
# ----------------------------------------------------------------------

theta_star <- make_theta_star(
  p = cfg$p,
  s_star = cfg$s_star,
  magnitude = cfg$theta_magnitude
)

selected_ts <- list()
selected_iht <- list()

for (jj in seq_along(cfg$n_grid)) {
  n_val <- cfg$n_grid[jj]
  selected_ts[[as.character(n_val)]] <- tune_ts_for_n(
    n = n_val,
    n_index = jj,
    theta_star = theta_star,
    cfg = cfg,
    registry_path = registry_path,
    candidate_dir = candidate_dir
  )

  selected_iht[[as.character(n_val)]] <- tune_iht_for_n(
    n = n_val,
    n_index = jj,
    theta_star = theta_star,
    cfg = cfg,
    registry_path = registry_path,
    candidate_dir = candidate_dir
  )
}

selected_ts_df <- dplyr::bind_rows(selected_ts)
selected_iht_df <- dplyr::bind_rows(selected_iht)
selected_tuning_used <- dplyr::bind_rows(selected_ts_df, selected_iht_df)
write.csv(
  selected_tuning_used,
  file.path(summary_dir, "selected_tuning_used_linear_TS_vs_baselines.csv"),
  row.names = FALSE
)

# ----------------------------------------------------------------------
# Stage B: final Monte Carlo evaluation
# ----------------------------------------------------------------------

fit_one_final_rep <- function(n, n_index, rep_id, theta_star, cfg, selected_ts_df, selected_iht_df) {
  seed_data <- seed_for(cfg$seed_base, n_index, rep_id, offset = 600000L)
  seed_split <- seed_for(cfg$seed_base, n_index, rep_id, offset = 700000L)
  seed_algo_ts <- seed_for(cfg$seed_base, n_index, rep_id, offset = 800000L)
  seed_lasso <- seed_for(cfg$seed_base, n_index, rep_id, offset = 900000L)
  seed_huber <- seed_for(cfg$seed_base, n_index, rep_id, offset = 1000000L)
  seed_shrink <- seed_for(cfg$seed_base, n_index, rep_id, offset = 1100000L)

  dat <- gen_one_linear_dataset(
    seed = seed_data,
    n = n,
    p = cfg$p,
    theta_star = theta_star,
    df_X = cfg$df_X,
    Sigma_X = diag(cfg$p),
    scale_X = cfg$scale_X,
    df_eps = cfg$df_eps,
    scale_eps = cfg$scale_eps
  )

  X <- dat$X
  y <- dat$y
  split <- make_train_val_split(n = n, val_frac = cfg$val_frac, seed = seed_split)
  train_idx <- split$train

  rows <- list()

  # TS-RIGHT, frozen per n.
  ts_par <- selected_ts_df[selected_ts_df$n == n, , drop = FALSE][1, ]
  t0 <- proc.time()[3]
  fit_ts <- fit_ts_right_linear_theta(
    X = X,
    y = y,
    s = cfg$s,
    eta = cfg$eta_ts,
    q = ts_par$q_selected,
    T1 = ts_par$T1_selected,
    T2 = ts_par$T2_selected,
    m_min = cfg$m_min,
    s_ref = cfg$s_star,
    seed = seed_algo_ts
  )
  runtime_ts <- proc.time()[3] - t0
  rows[[length(rows) + 1L]] <- data.frame(
    n = n,
    p = cfg$p,
    rep_id = rep_id,
    method = "TS-RIGHT",
    l2_error = safe_l2(fit_ts$theta, theta_star),
    runtime_sec = runtime_ts,
    seed_data = seed_data,
    seed_split = seed_split,
    seed_algo = seed_algo_ts,
    tuning_config_id = ts_par$config_id,
    q = ts_par$q_selected,
    T1 = ts_par$T1_selected,
    T2 = ts_par$T2_selected,
    T_fd_budget = ts_par$T_fd_budget,
    eta = cfg$eta_ts,
    T = NA_integer_,
    lambda = NA_real_,
    tau = NA_real_,
    tau_q = NA_real_
  )

  # IHT, frozen per n.
  iht_par <- selected_iht_df[selected_iht_df$n == n, , drop = FALSE][1, ]
  t0 <- proc.time()[3]
  theta_iht <- tryCatch(
    fit_iht_linear_theta(
      X = X,
      y = y,
      s = cfg$s,
      eta = iht_par$eta_selected,
      T = iht_par$T_selected
    ),
    error = function(e) rep(NA_real_, cfg$p)
  )
  runtime_iht <- proc.time()[3] - t0
  rows[[length(rows) + 1L]] <- data.frame(
    n = n,
    p = cfg$p,
    rep_id = rep_id,
    method = "IHT",
    l2_error = safe_l2(theta_iht, theta_star),
    runtime_sec = runtime_iht,
    seed_data = seed_data,
    seed_split = seed_split,
    seed_algo = NA_integer_,
    tuning_config_id = iht_par$config_id,
    q = NA_real_,
    T1 = NA_integer_,
    T2 = NA_integer_,
    T_fd_budget = NA_integer_,
    eta = iht_par$eta_selected,
    T = iht_par$T_selected,
    lambda = NA_real_,
    tau = NA_real_,
    tau_q = NA_real_
  )

  # Lasso via cv.glmnet on the full dataset.
  t0 <- proc.time()[3]
  fit_lasso <- tryCatch(
    fit_lasso_cv(X, y, seed = seed_lasso, nfolds = cfg$cv_nfolds),
    error = function(e) list(theta = rep(NA_real_, cfg$p), lambda = NA_real_)
  )
  runtime_lasso <- proc.time()[3] - t0
  rows[[length(rows) + 1L]] <- data.frame(
    n = n,
    p = cfg$p,
    rep_id = rep_id,
    method = "Lasso",
    l2_error = safe_l2(fit_lasso$theta, theta_star),
    runtime_sec = runtime_lasso,
    seed_data = seed_data,
    seed_split = seed_split,
    seed_algo = seed_lasso,
    tuning_config_id = NA_character_,
    q = NA_real_,
    T1 = NA_integer_,
    T2 = NA_integer_,
    T_fd_budget = NA_integer_,
    eta = NA_real_,
    T = NA_integer_,
    lambda = fit_lasso$lambda,
    tau = NA_real_,
    tau_q = NA_real_
  )

  # Adaptive Huber-Lasso via validation median absolute prediction error.
  t0 <- proc.time()[3]
  fit_huber <- fit_huber_lasso_tuned(
    X_full = X,
    y_full = y,
    train_idx = train_idx,
    lambda_grid = cfg$huber_lambda_grid,
    tau_quantile_grid = cfg$huber_tau_quantile_grid,
    val_alpha = cfg$selection_alpha,
    T_tune = cfg$T_huber_tune,
    T_final = cfg$T_huber_final,
    seed = seed_huber,
    nfolds = cfg$cv_nfolds
  )
  runtime_huber <- proc.time()[3] - t0
  huber_sel <- fit_huber$selected
  if (nrow(huber_sel) == 0L) {
    huber_lambda <- NA_real_; huber_tau <- NA_real_; huber_tau_q <- NA_real_
  } else {
    huber_lambda <- huber_sel$lambda[1]
    huber_tau <- huber_sel$tau_final[1]
    huber_tau_q <- huber_sel$tau_q[1]
  }
  rows[[length(rows) + 1L]] <- data.frame(
    n = n,
    p = cfg$p,
    rep_id = rep_id,
    method = "Adaptive Huber",
    l2_error = safe_l2(fit_huber$theta, theta_star),
    runtime_sec = runtime_huber,
    seed_data = seed_data,
    seed_split = seed_split,
    seed_algo = seed_huber,
    tuning_config_id = NA_character_,
    q = NA_real_,
    T1 = NA_integer_,
    T2 = NA_integer_,
    T_fd_budget = NA_integer_,
    eta = NA_real_,
    T = NA_integer_,
    lambda = huber_lambda,
    tau = huber_tau,
    tau_q = huber_tau_q
  )

  # Shrinkage/truncation baseline.
  t0 <- proc.time()[3]
  fit_shrink <- fit_shrinkage_tuned(
    X_full = X,
    y_full = y,
    train_idx = train_idx,
    tau_quantile_grid = cfg$shrink_tau_quantile_grid,
    val_alpha = cfg$selection_alpha,
    seed = seed_shrink,
    nfolds = cfg$cv_nfolds
  )
  runtime_shrink <- proc.time()[3] - t0
  shrink_sel <- fit_shrink$selected
  if (nrow(shrink_sel) == 0L) {
    shrink_lambda <- NA_real_; shrink_tau <- NA_real_; shrink_tau_q <- NA_real_
  } else {
    shrink_lambda <- shrink_sel$lambda[1]
    shrink_tau <- shrink_sel$tau_x_final[1]
    shrink_tau_q <- shrink_sel$tau_q[1]
  }
  rows[[length(rows) + 1L]] <- data.frame(
    n = n,
    p = cfg$p,
    rep_id = rep_id,
    method = "Shrinkage",
    l2_error = safe_l2(fit_shrink$theta, theta_star),
    runtime_sec = runtime_shrink,
    seed_data = seed_data,
    seed_split = seed_split,
    seed_algo = seed_shrink,
    tuning_config_id = NA_character_,
    q = NA_real_,
    T1 = NA_integer_,
    T2 = NA_integer_,
    T_fd_budget = NA_integer_,
    eta = NA_real_,
    T = NA_integer_,
    lambda = shrink_lambda,
    tau = shrink_tau,
    tau_q = shrink_tau_q
  )

  dplyr::bind_rows(rows)
}

param_grid <- expand.grid(
  n = cfg$n_grid,
  rep_id = seq_len(cfg$reps),
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)
param_grid$n_index <- match(param_grid$n, cfg$n_grid)

message("Running final linear TS-vs-baselines experiment: ", nrow(param_grid), " data replications.")

helper_exports <- c(
  "cfg", "param_grid", "theta_star", "selected_ts_df", "selected_iht_df",
  "seed_for", "with_temp_seed", "restore_temp_seed", "make_train_val_split", "make_foldid",
  "median_abs_prediction_error", "safe_l2",
  "fit_ts_right_linear_theta", "fit_iht_linear_theta", "fit_lasso_cv",
  "fit_huber_lasso_tuned", "fit_shrinkage_tuned", "fit_one_final_rep"
)

if (cfg$use_parallel && cfg$n_cores > 1L) {
  cl <- parallel::makeCluster(cfg$n_cores)
  doParallel::registerDoParallel(cl)
  on.exit({
    try(parallel::stopCluster(cl), silent = TRUE)
  }, add = TRUE)

  raw_results <- foreach(
    ii = seq_len(nrow(param_grid)),
    .combine = dplyr::bind_rows,
    .packages = c("dplyr", "glmnet", "here"),
    .export = helper_exports
  ) %dopar% {
    source(here("simulations", "07_schedule_ablation", "utils_schedule_right.R"))
    row <- param_grid[ii, ]
    fit_one_final_rep(
      n = row$n,
      n_index = row$n_index,
      rep_id = row$rep_id,
      theta_star = theta_star,
      cfg = cfg,
      selected_ts_df = selected_ts_df,
      selected_iht_df = selected_iht_df
    )
  }
} else {
  raw_results <- lapply(seq_len(nrow(param_grid)), function(ii) {
    row <- param_grid[ii, ]
    fit_one_final_rep(
      n = row$n,
      n_index = row$n_index,
      rep_id = row$rep_id,
      theta_star = theta_star,
      cfg = cfg,
      selected_ts_df = selected_ts_df,
      selected_iht_df = selected_iht_df
    )
  }) %>% dplyr::bind_rows()
}

# ----------------------------------------------------------------------
# Save raw results, summaries, and figures
# ----------------------------------------------------------------------

saveRDS(raw_results, file.path(raw_dir, "raw_linear_TS_vs_baselines.rds"))
write.csv(raw_results, file.path(raw_dir, "raw_linear_TS_vs_baselines.csv"), row.names = FALSE)

summary_results <- raw_results %>%
  group_by(n, method) %>%
  summarise(
    median_l2 = median(l2_error, na.rm = FALSE),
    q25_l2 = quantile(l2_error, 0.25, na.rm = FALSE),
    q75_l2 = quantile(l2_error, 0.75, na.rm = FALSE),
    iqr_l2 = IQR(l2_error, na.rm = FALSE),
    mean_l2 = mean(l2_error, na.rm = FALSE),
    sd_l2 = sd(l2_error, na.rm = FALSE),
    fail_rate = mean(!is.finite(l2_error) | is.na(l2_error)),
    median_runtime_sec = median(runtime_sec, na.rm = FALSE),
    .groups = "drop"
  )

write.csv(summary_results, file.path(summary_dir, "summary_linear_TS_vs_baselines.csv"), row.names = FALSE)
saveRDS(summary_results, file.path(summary_dir, "summary_linear_TS_vs_baselines.rds"))

method_levels <- c("TS-RIGHT", "Adaptive Huber", "Shrinkage", "Lasso", "IHT")
summary_results$method <- factor(summary_results$method, levels = method_levels)

plot_all <- ggplot(
  summary_results,
  aes(x = log(n), y = log(median_l2), color = method, shape = method, group = method)
) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.8) +
  labs(
    x = "Log(Sample Size, n)",
    y = "Log(Median L2 Error)",
    color = "Method",
    shape = "Method"
  ) +
  theme_bw(base_size = 15) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )

print(plot_all)
ggsave(file.path(fig_dir, "linear_TS_vs_baselines_all_methods.pdf"), plot_all, width = 9.5, height = 5)

plot_no_iht <- summary_results %>%
  filter(method != "IHT") %>%
  ggplot(aes(x = log(n), y = log(median_l2), color = method, shape = method, group = method)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.8) +
  labs(
    x = "Log(Sample Size, n)",
    y = "Log(Median L2 Error)",
    color = "Method",
    shape = "Method"
  ) +
  theme_bw(base_size = 15) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )

print(plot_no_iht)
ggsave(file.path(fig_dir, "linear_TS_vs_baselines_without_IHT.pdf"), plot_no_iht, width = 9.5, height = 5)

message("Done. Results saved under: ", out_dir)







############################################################################################################
#box plot 
library(dplyr)
library(ggplot2)
library(here)
plot_saved_out_dir <- here("results", "08_TS_baseline_comparison", "linear")
plot_saved_raw_dir <- file.path(plot_saved_out_dir, "raw")
plot_saved_summary_dir <- file.path(plot_saved_out_dir, "summary")
plot_saved_fig_dir <- file.path(plot_saved_out_dir, "figures")

dir.create(plot_saved_summary_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_saved_fig_dir, recursive = TRUE, showWarnings = FALSE)

raw_results_saved <- readRDS(file.path(plot_saved_raw_dir, "raw_linear_TS_vs_baselines.rds"))

method_palette <- c(
  "RIGHT" = "#332288",
  "Adaptive Huber" = "#3B82A1",
  "Shrinkage" = "#44AA99",
  "Lasso" = "#CAB54B",
  "IHT" = "#CC6677"
)
method_levels_display <- names(method_palette)
n_levels_display <- sort(unique(raw_results_saved$n))

quantile_plot_raw <- raw_results_saved %>%
  filter(is.finite(l2_error), !is.na(l2_error), l2_error > 0) %>%
  mutate(
    method_display = dplyr::recode(method, "TS-RIGHT" = "RIGHT"),
    method_display = factor(method_display, levels = method_levels_display),
    n_factor = factor(n, levels = n_levels_display)
  )

non_iht_l2 <- quantile_plot_raw %>%
  filter(method_display != "IHT") %>%
  pull(l2_error)

display_cap_l2 <- as.numeric(quantile(non_iht_l2, 0.99, na.rm = TRUE))
display_floor_l2 <- as.numeric(quantile(quantile_plot_raw$l2_error, 0.005, na.rm = TRUE))
display_floor_l2 <- max(display_floor_l2, 1e-8)

dodge_step <- 0.13
tick_half_width <- 0.045
method_mid <- (length(method_levels_display) + 1) / 2

grouped_quantiles <- quantile_plot_raw %>%
  group_by(n, n_factor, method_display) %>%
  summarise(
    q25_l2 = quantile(l2_error, 0.25, na.rm = TRUE),
    q50_l2 = quantile(l2_error, 0.50, na.rm = TRUE),
    q75_l2 = quantile(l2_error, 0.75, na.rm = TRUE),
    n_obs = dplyr::n(),
    .groups = "drop"
  ) %>%
  mutate(
    x_base = as.numeric(n_factor),
    method_index = as.numeric(method_display),
    x_pos = x_base + (method_index - method_mid) * dodge_step
  )

iht_extreme_marks <- quantile_plot_raw %>%
  filter(method_display == "IHT") %>%
  group_by(n, n_factor, method_display) %>%
  summarise(
    n_obs = dplyr::n(),
    n_above_cap = sum(l2_error > display_cap_l2, na.rm = TRUE),
    max_l2 = max(l2_error, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(n_above_cap > 0) %>%
  mutate(
    x_base = as.numeric(n_factor),
    method_index = as.numeric(method_display),
    x_pos = x_base + (method_index - method_mid) * dodge_step,
    y_pos = display_cap_l2,
    label = paste0(n_above_cap, "/", n_obs)
  )

write.csv(
  grouped_quantiles,
  file.path(plot_saved_summary_dir, "grouped_quantile_lines_linear_TS_vs_baselines.csv"),
  row.names = FALSE
)
write.csv(
  iht_extreme_marks,
  file.path(plot_saved_summary_dir, "iht_extreme_counts_grouped_quantile_plot.csv"),
  row.names = FALSE
)

plot_grouped_quantile_lines <- ggplot(grouped_quantiles, aes(color = method_display)) +
  geom_linerange(
    aes(x = x_pos, ymin = q25_l2, ymax = q75_l2),
    linewidth = 0.95
  ) +
  geom_segment(
    aes(
      x = x_pos - tick_half_width,
      xend = x_pos + tick_half_width,
      y = q25_l2,
      yend = q25_l2
    ),
    linewidth = 0.65
  ) +
  geom_segment(
    aes(
      x = x_pos - tick_half_width,
      xend = x_pos + tick_half_width,
      y = q50_l2,
      yend = q50_l2
    ),
    linewidth = 1.25
  ) +
  geom_segment(
    aes(
      x = x_pos - tick_half_width,
      xend = x_pos + tick_half_width,
      y = q75_l2,
      yend = q75_l2
    ),
    linewidth = 0.65
  ) +
  geom_point(
    data = iht_extreme_marks,
    aes(x = x_pos, y = y_pos),
    inherit.aes = FALSE,
    shape = 24,
    size = 3.1,
    stroke = 0.4,
    color = method_palette["IHT"],
    fill = method_palette["IHT"]
  ) +
  geom_text(
    data = iht_extreme_marks,
    aes(x = x_pos, y = y_pos, label = label),
    inherit.aes = FALSE,
    vjust = -0.7,
    size = 3.0,
    color = method_palette["IHT"]
  ) +
  scale_x_continuous(
    breaks = seq_along(n_levels_display),
    labels = n_levels_display,
    expand = expansion(mult = c(0.05, 0.07))
  ) +
  scale_y_log10(
    labels = function(x) formatC(x, format = "fg", digits = 3)
  ) +
  coord_cartesian(
    ylim = c(display_floor_l2, display_cap_l2),
    clip = "off"
  ) +
  scale_color_manual(
    values = method_palette,
    breaks = method_levels_display,
    drop = FALSE
  ) +
  labs(
    x = "Sample Size, n",
    y = "L2 Error (log scale)",
    color = "Method",
    #title = "Linear Simulation: Grouped L2 Quantiles",
    #subtitle = "Vertical lines show q25-q75; thick center ticks show q50. IHT triangles mark values above the display cap.",
  ) +
  theme_bw(base_size = 15) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    plot.margin = margin(t = 18, r = 18, b = 8, l = 8)
  )

print(plot_grouped_quantile_lines)
ggsave(
  file.path(plot_saved_fig_dir, "linear_TS_vs_baselines_grouped_quantile_lines.pdf"),
  plot_grouped_quantile_lines,
  width = 10.5,
  height = 6
)
ggsave(
  file.path(plot_saved_fig_dir, "linear_TS_vs_baselines_grouped_quantile_lines.png"),
  plot_grouped_quantile_lines,
  width = 10.5,
  height = 6,
  dpi = 300
)
saveRDS(plot_grouped_quantile_lines,file.path(plot_saved_fig_dir,"linear_TS_vs_baselines_grouped_quantile_lines.rds"))
