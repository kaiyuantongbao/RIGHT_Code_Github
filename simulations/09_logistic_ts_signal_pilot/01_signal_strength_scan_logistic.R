# ======================================================================
# Logistic pilot: signal-strength scan for TS-RIGHT / IHT / Lasso / Shrinkage
#
# Goal:
#   First-stage pre-experiment for choosing a useful signal-to-noise regime.
#   We deliberately freeze TS-RIGHT at q = 0.5 and T1 = 200 by default,
#   so the only scientific axis in this script is theta_magnitude.
#
# Main metric:
#   ||theta_hat - theta_star||_2
#
# Outputs:
#   results/09_logistic_ts_signal_pilot/signal_strength_scan/raw/
#   results/09_logistic_ts_signal_pilot/signal_strength_scan/summary/
#   results/09_logistic_ts_signal_pilot/signal_strength_scan/figures/
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
  experiment_id = "logistic_signal_strength_scan",
  
  p = 800,
  n_grid = c(1000, 1400),
  s_star = 10,
  s = 20,
  m_min = 15,
  c_K1 = 0.10,
  c_K2 = 0.50,
  
  # For a quick smoke test use 2-5 reps. For a real pilot use 100-200.
  reps = 25,
  seed_base = 20260508L,
  
  df_X = 2.1,
  scale_X = 5,
  theta_magnitude_grid = c( 0.50),
  
  # Frozen TS-RIGHT schedule for the first pilot.
  # q is the Stage-II fraction n2 / n. T2 is needed by the existing TS wrapper.
  q_grid = c(0.50),
  T1_grid = c(200),
  # With m_min = 15 and c_K2 = 0.50, T2 = 8 is feasible at n = 600.
  T2_grid = c(8),
  eta_ts = 0.02,
  
  # IHT is intentionally simple here. Tune eta/T later if IHT is unstable.
  iht_eta = 0.02,
  iht_T = 200,
  
  # Shrinkage tuning: lambda by CV after truncating X at these quantiles.
  shrink_tau_quantile_grid = c(0.90, 0.95, 0.98),
  
  cv_nfolds = 5,
  use_parallel = TRUE,
  n_cores = max(1L, parallel::detectCores() - 1L)
)

env_int <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) return(default)
  as.integer(value)
}

env_num_vec <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) return(default)
  as.numeric(strsplit(value, ",", fixed = TRUE)[[1]])
}

env_bool <- function(name, default) {
  value <- tolower(Sys.getenv(name, unset = ""))
  if (!nzchar(value)) return(default)
  value %in% c("1", "true", "t", "yes", "y")
}

# Optional smoke-test overrides, for example:
#   $env:PILOT_REPS='1'
#   $env:PILOT_N_GRID='600'
#   $env:PILOT_THETA_GRID='0.5'
#   $env:PILOT_USE_PARALLEL='false'
cfg$reps <- env_int("PILOT_REPS", cfg$reps)
cfg$n_grid <- as.integer(env_num_vec("PILOT_N_GRID", cfg$n_grid))
cfg$theta_magnitude_grid <- env_num_vec("PILOT_THETA_GRID", cfg$theta_magnitude_grid)
cfg$q_grid <- env_num_vec("PILOT_Q_GRID", cfg$q_grid)
cfg$T1_grid <- as.integer(env_num_vec("PILOT_T1_GRID", cfg$T1_grid))
cfg$T2_grid <- as.integer(env_num_vec("PILOT_T2_GRID", cfg$T2_grid))
cfg$m_min <- env_int("PILOT_M_MIN", cfg$m_min)
cfg$c_K1 <- env_num_vec("PILOT_C_K1", cfg$c_K1)[1]
cfg$c_K2 <- env_num_vec("PILOT_C_K2", cfg$c_K2)[1]
cfg$use_parallel <- env_bool("PILOT_USE_PARALLEL", cfg$use_parallel)
cfg$n_cores <- env_int("PILOT_N_CORES", cfg$n_cores)

out_dir <- here("results", "09_logistic_ts_signal_pilot", "signal_strength_scan")
raw_dir <- file.path(out_dir, "raw")
summary_dir <- file.path(out_dir, "summary")
fig_dir <- file.path(out_dir, "figures")

dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# ----------------------------------------------------------------------
# Utilities
# ----------------------------------------------------------------------

seed_for <- function(base, n_index, mag_index, rep_id, offset = 0L) {
  as.integer(base + offset + 100000L * n_index + 1000L * mag_index + rep_id)
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
  if (isTRUE(state$old_exists)) {
    assign(".Random.seed", state$old_state, envir = .GlobalEnv)
  } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    rm(".Random.seed", envir = .GlobalEnv)
  }
  invisible(NULL)
}

make_theta_star_signed <- function(p, s_star, magnitude) {
  theta_star <- rep(0, p)
  signs <- rep(c(1, -1), length.out = s_star)
  theta_star[seq_len(s_star)] <- signs * magnitude
  theta_star
}

make_foldid <- function(n, nfolds = 5, seed = 1L) {
  nfolds <- min(as.integer(nfolds), n)
  rng <- with_temp_seed(seed)
  on.exit(restore_temp_seed(rng), add = TRUE)
  
  foldid <- rep(seq_len(nfolds), length.out = n)
  sample(foldid, size = n, replace = FALSE)
}

safe_l2 <- function(theta, theta_star) {
  if (length(theta) != length(theta_star) || any(!is.finite(theta))) return(NA_real_)
  sqrt(sum((theta - theta_star)^2))
}

logistic_nll <- function(theta, X, y) {
  if (any(!is.finite(theta))) return(Inf)
  eta <- as.numeric(X %*% theta)
  mean(log1p(exp(-abs(eta))) + pmax(eta, 0) - y * eta)
}

eta_diagnostics <- function(X, theta_star) {
  eta <- as.numeric(X %*% theta_star)
  prob <- plogis(eta)
  data.frame(
    eta_sd = sd(eta),
    eta_median_abs = median(abs(eta)),
    eta_q90_abs = as.numeric(quantile(abs(eta), 0.90, na.rm = TRUE)),
    eta_gt3_rate = mean(abs(eta) > 3),
    mean_bernoulli_var = mean(prob * (1 - prob))
  )
}

generate_one_logistic_dataset <- function(seed, n, p, theta_star, df_X, scale_X) {
  set.seed(seed)
  dat <- generate_logistic_data(
    n = n,
    p = p,
    theta_star = theta_star,
    generator_X = generate_X_t,
    df_X = df_X,
    Sigma_X = diag(p),
    scale_X = scale_X
  )
  list(X = dat$X, y = as.numeric(dat$y))
}

fit_ts_right_logistic_theta <- function(
    X,
    y,
    theta_star,
    s,
    eta,
    q,
    T1,
    T2,
    m_min,
    c_K1,
    c_K2,
    s_ref,
    seed
) {
  budget <- compute_ts_budget(
    n = nrow(X),
    p = ncol(X),
    s = s,
    q = q,
    T1 = T1,
    T2 = T2,
    s_ref = s_ref,
    m_min = m_min,
    c_K1 = c_K1,
    c_K2 = c_K2
  )
  
  meta <- data.frame(
    eligible = budget$eligible,
    n1 = budget$n1,
    n2 = budget$n2,
    b2 = budget$b2,
    m_min = m_min,
    c_K1 = c_K1,
    c_K2 = c_K2,
    K1 = budget$K1,
    K2 = budget$K2,
    K1_block_size = budget$n1 / budget$K1,
    K2_block_size = budget$b2 / budget$K2,
    T_fd_budget = budget$T_fd_budget
  )
  
  if (!isTRUE(budget$eligible)) {
    return(list(theta = rep(NA_real_, ncol(X)), meta = meta))
  }
  
  rng <- with_temp_seed(seed)
  on.exit(restore_temp_seed(rng), add = TRUE)
  
  idx <- sample.int(nrow(X))
  id1 <- idx[seq_len(budget$n1)]
  id2 <- idx[(budget$n1 + 1L):nrow(X)]
  
  theta_cur <- solver_right(
    X = X[id1, , drop = FALSE],
    y = y[id1],
    s = s,
    eta = eta,
    T_max = T1,
    K = budget$K1,
    theta_init = rep(0, ncol(X)),
    grad_func_samplewise = grad_logistic_regression_samplewise,
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
      grad_func_samplewise = grad_logistic_regression_samplewise,
      record_trace = FALSE
    )
  }
  
  list(theta = as.numeric(theta_cur), meta = meta)
}

fit_iht_logistic_theta <- function(X, y, s, eta, T) {
  as.numeric(solver_iht(
    X = X,
    y = y,
    s = s,
    eta = eta,
    T_max = T,
    grad_func = grad_logistic_regression
  ))
}

fit_lasso_logistic_cv <- function(X, y, seed, nfolds = 5) {
  foldid <- make_foldid(nrow(X), nfolds = nfolds, seed = seed)
  cv_fit <- cv.glmnet(
    X,
    y,
    family = "binomial",
    standardize = FALSE,
    intercept = FALSE,
    foldid = foldid
  )
  lambda_best <- cv_fit$lambda.min
  final_fit <- glmnet(
    X,
    y,
    family = "binomial",
    lambda = lambda_best,
    standardize = FALSE,
    intercept = FALSE
  )
  list(theta = as.numeric(final_fit$beta), lambda = as.numeric(lambda_best))
}

fit_shrinkage_logistic_tuned <- function(
    X,
    y,
    tau_quantile_grid,
    seed,
    nfolds = 5
) {
  cand <- vector("list", length(tau_quantile_grid))
  
  for (jj in seq_along(tau_quantile_grid)) {
    tau_q <- tau_quantile_grid[jj]
    tau_x <- as.numeric(quantile(abs(X), probs = tau_q, na.rm = TRUE))
    tau_x <- max(tau_x, 1e-8)
    
    X_shrunk <- truncate_operator(X, tau_x)
    
    fit_j <- tryCatch({
      foldid <- make_foldid(nrow(X), nfolds = nfolds, seed = seed + as.integer(1000 * tau_q))
      cv_fit <- cv.glmnet(
        X_shrunk,
        y,
        family = "binomial",
        standardize = FALSE,
        intercept = FALSE,
        foldid = foldid
      )
      lambda_best <- cv_fit$lambda.min
      theta_hat <- solver_shrinkage_method(
        X = X,
        y = y,
        tau_x = tau_x,
        tau_y = 999,
        lambda = lambda_best,
        family = "binomial"
      )
      lambda_idx <- which.min(abs(cv_fit$lambda - lambda_best))
      list(theta = as.numeric(theta_hat), lambda = as.numeric(lambda_best), val = cv_fit$cvm[lambda_idx])
    }, error = function(e) {
      list(theta = rep(NA_real_, ncol(X)), lambda = NA_real_, val = Inf)
    })
    
    cand[[jj]] <- data.frame(
      tau_q = tau_q,
      tau_x = tau_x,
      lambda = fit_j$lambda,
      cv_nll = fit_j$val,
      theta_index = jj
    )
    attr(cand[[jj]], "theta") <- fit_j$theta
  }
  
  cand_df <- dplyr::bind_rows(cand)
  valid <- cand_df[is.finite(cand_df$cv_nll), , drop = FALSE]
  
  if (nrow(valid) == 0L) {
    return(list(theta = rep(NA_real_, ncol(X)), selected = data.frame()))
  }
  
  valid$dist_tau_center <- abs(valid$tau_q - 0.95)
  valid <- valid[order(valid$cv_nll, valid$dist_tau_center), , drop = FALSE]
  chosen <- valid[1, , drop = FALSE]
  theta_chosen <- attr(cand[[chosen$theta_index[1]]], "theta")
  
  list(theta = theta_chosen, selected = chosen)
}

fit_one_rep <- function(n, n_index, theta_magnitude, mag_index, rep_id, cfg) {
  theta_star <- make_theta_star_signed(
    p = cfg$p,
    s_star = cfg$s_star,
    magnitude = theta_magnitude
  )
  
  seed_data <- seed_for(cfg$seed_base, n_index, mag_index, rep_id, offset = 1000000L)
  seed_ts <- seed_for(cfg$seed_base, n_index, mag_index, rep_id, offset = 2000000L)
  seed_lasso <- seed_for(cfg$seed_base, n_index, mag_index, rep_id, offset = 3000000L)
  seed_shrink <- seed_for(cfg$seed_base, n_index, mag_index, rep_id, offset = 4000000L)
  
  dat <- generate_one_logistic_dataset(
    seed = seed_data,
    n = n,
    p = cfg$p,
    theta_star = theta_star,
    df_X = cfg$df_X,
    scale_X = cfg$scale_X
  )
  
  X <- dat$X
  y <- dat$y
  eta_diag <- eta_diagnostics(X, theta_star)
  
  rows <- list()
  
  for (q in cfg$q_grid) {
    for (T1 in cfg$T1_grid) {
      for (T2 in cfg$T2_grid) {
        t0 <- proc.time()[3]
        fit_ts <- tryCatch(
          fit_ts_right_logistic_theta(
            X = X,
            y = y,
            theta_star = theta_star,
            s = cfg$s,
            eta = cfg$eta_ts,
            q = q,
            T1 = T1,
            T2 = T2,
            m_min = cfg$m_min,
            c_K1 = cfg$c_K1,
            c_K2 = cfg$c_K2,
            s_ref = cfg$s_star,
            seed = seed_ts
          ),
          error = function(e) {
            list(theta = rep(NA_real_, cfg$p), meta = data.frame(eligible = FALSE))
          }
        )
        runtime_sec <- proc.time()[3] - t0
        meta <- fit_ts$meta
        
        rows[[length(rows) + 1L]] <- data.frame(
          n = n,
          p = cfg$p,
          rep_id = rep_id,
          theta_magnitude = theta_magnitude,
          method = "TS-RIGHT",
          l2_error = safe_l2(fit_ts$theta, theta_star),
          train_nll = logistic_nll(fit_ts$theta, X, y),
          runtime_sec = runtime_sec,
          q = q,
          T1 = T1,
          T2 = T2,
          eta = cfg$eta_ts,
          T = NA_integer_,
          lambda = NA_real_,
          tau_q = NA_real_,
          eligible = if ("eligible" %in% names(meta)) isTRUE(meta$eligible[1]) else FALSE,
          m_min = if ("m_min" %in% names(meta)) meta$m_min[1] else cfg$m_min,
          c_K1 = if ("c_K1" %in% names(meta)) meta$c_K1[1] else cfg$c_K1,
          c_K2 = if ("c_K2" %in% names(meta)) meta$c_K2[1] else cfg$c_K2,
          K1 = if ("K1" %in% names(meta)) meta$K1[1] else NA_integer_,
          K2 = if ("K2" %in% names(meta)) meta$K2[1] else NA_integer_,
          K1_block_size = if ("K1_block_size" %in% names(meta)) meta$K1_block_size[1] else NA_real_,
          K2_block_size = if ("K2_block_size" %in% names(meta)) meta$K2_block_size[1] else NA_real_,
          seed_data = seed_data,
          seed_algo = seed_ts
        )
      }
    }
  }
  
  t0 <- proc.time()[3]
  theta_iht <- tryCatch(
    fit_iht_logistic_theta(
      X = X,
      y = y,
      s = cfg$s,
      eta = cfg$iht_eta,
      T = cfg$iht_T
    ),
    error = function(e) rep(NA_real_, cfg$p)
  )
  runtime_sec <- proc.time()[3] - t0
  rows[[length(rows) + 1L]] <- data.frame(
    n = n,
    p = cfg$p,
    rep_id = rep_id,
    theta_magnitude = theta_magnitude,
    method = "IHT",
    l2_error = safe_l2(theta_iht, theta_star),
    train_nll = logistic_nll(theta_iht, X, y),
    runtime_sec = runtime_sec,
    q = NA_real_,
    T1 = NA_integer_,
    T2 = NA_integer_,
    eta = cfg$iht_eta,
    T = cfg$iht_T,
    lambda = NA_real_,
    tau_q = NA_real_,
    eligible = TRUE,
    m_min = NA_integer_,
    c_K1 = NA_real_,
    c_K2 = NA_real_,
    K1 = NA_integer_,
    K2 = NA_integer_,
    K1_block_size = NA_real_,
    K2_block_size = NA_real_,
    seed_data = seed_data,
    seed_algo = NA_integer_
  )
  
  t0 <- proc.time()[3]
  fit_lasso <- tryCatch(
    fit_lasso_logistic_cv(X, y, seed = seed_lasso, nfolds = cfg$cv_nfolds),
    error = function(e) list(theta = rep(NA_real_, cfg$p), lambda = NA_real_)
  )
  runtime_sec <- proc.time()[3] - t0
  rows[[length(rows) + 1L]] <- data.frame(
    n = n,
    p = cfg$p,
    rep_id = rep_id,
    theta_magnitude = theta_magnitude,
    method = "Lasso",
    l2_error = safe_l2(fit_lasso$theta, theta_star),
    train_nll = logistic_nll(fit_lasso$theta, X, y),
    runtime_sec = runtime_sec,
    q = NA_real_,
    T1 = NA_integer_,
    T2 = NA_integer_,
    eta = NA_real_,
    T = NA_integer_,
    lambda = fit_lasso$lambda,
    tau_q = NA_real_,
    eligible = TRUE,
    m_min = NA_integer_,
    c_K1 = NA_real_,
    c_K2 = NA_real_,
    K1 = NA_integer_,
    K2 = NA_integer_,
    K1_block_size = NA_real_,
    K2_block_size = NA_real_,
    seed_data = seed_data,
    seed_algo = seed_lasso
  )
  
  t0 <- proc.time()[3]
  fit_shrink <- tryCatch(
    fit_shrinkage_logistic_tuned(
      X = X,
      y = y,
      tau_quantile_grid = cfg$shrink_tau_quantile_grid,
      seed = seed_shrink,
      nfolds = cfg$cv_nfolds
    ),
    error = function(e) list(theta = rep(NA_real_, cfg$p), selected = data.frame())
  )
  runtime_sec <- proc.time()[3] - t0
  shrink_sel <- fit_shrink$selected
  rows[[length(rows) + 1L]] <- data.frame(
    n = n,
    p = cfg$p,
    rep_id = rep_id,
    theta_magnitude = theta_magnitude,
    method = "Shrinkage",
    l2_error = safe_l2(fit_shrink$theta, theta_star),
    train_nll = logistic_nll(fit_shrink$theta, X, y),
    runtime_sec = runtime_sec,
    q = NA_real_,
    T1 = NA_integer_,
    T2 = NA_integer_,
    eta = NA_real_,
    T = NA_integer_,
    lambda = if (nrow(shrink_sel) > 0L) shrink_sel$lambda[1] else NA_real_,
    tau_q = if (nrow(shrink_sel) > 0L) shrink_sel$tau_q[1] else NA_real_,
    eligible = TRUE,
    m_min = NA_integer_,
    c_K1 = NA_real_,
    c_K2 = NA_real_,
    K1 = NA_integer_,
    K2 = NA_integer_,
    K1_block_size = NA_real_,
    K2_block_size = NA_real_,
    seed_data = seed_data,
    seed_algo = seed_shrink
  )
  
  dplyr::bind_cols(dplyr::bind_rows(rows), eta_diag[rep(1, length(rows)), ])
}

# ----------------------------------------------------------------------
# Main loop
# ----------------------------------------------------------------------

param_grid <- expand.grid(
  n = cfg$n_grid,
  theta_magnitude = cfg$theta_magnitude_grid,
  rep_id = seq_len(cfg$reps),
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)
param_grid$n_index <- match(param_grid$n, cfg$n_grid)
param_grid$mag_index <- match(param_grid$theta_magnitude, cfg$theta_magnitude_grid)

message("Running logistic signal-strength pilot: ", nrow(param_grid), " data replications.")
message("p = ", cfg$p, ", s_star = ", cfg$s_star, ", s = ", cfg$s)
message("df_X = ", cfg$df_X, ", scale_X = ", signif(cfg$scale_X, 5))
message("TS grid: q = ", paste(cfg$q_grid, collapse = ","), 
        "; T1 = ", paste(cfg$T1_grid, collapse = ","),
        "; T2 = ", paste(cfg$T2_grid, collapse = ","))
message("TS blocks: m_min = ", cfg$m_min, ", c_K1 = ", cfg$c_K1, ", c_K2 = ", cfg$c_K2)

helper_exports <- c(
  "cfg", "param_grid", "seed_for", "with_temp_seed", "restore_temp_seed",
  "make_theta_star_signed", "make_foldid", "safe_l2", "logistic_nll",
  "eta_diagnostics", "generate_one_logistic_dataset",
  "fit_ts_right_logistic_theta", "fit_iht_logistic_theta",
  "fit_lasso_logistic_cv", "fit_shrinkage_logistic_tuned", "fit_one_rep"
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
    fit_one_rep(
      n = row$n,
      n_index = row$n_index,
      theta_magnitude = row$theta_magnitude,
      mag_index = row$mag_index,
      rep_id = row$rep_id,
      cfg = cfg
    )
  }
} else {
  raw_results <- lapply(seq_len(nrow(param_grid)), function(ii) {
    row <- param_grid[ii, ]
    fit_one_rep(
      n = row$n,
      n_index = row$n_index,
      theta_magnitude = row$theta_magnitude,
      mag_index = row$mag_index,
      rep_id = row$rep_id,
      cfg = cfg
    )
  }) %>% dplyr::bind_rows()
}

# ----------------------------------------------------------------------
# Save summaries and plots
# ----------------------------------------------------------------------

saveRDS(raw_results, file.path(raw_dir, "raw_signal_strength_scan_logistic.rds"))
write.csv(raw_results, file.path(raw_dir, "raw_signal_strength_scan_logistic.csv"), row.names = FALSE)

summary_results <- raw_results %>%
  group_by(n, theta_magnitude, method) %>%
  summarise(
    median_l2 = median(l2_error, na.rm = TRUE),
    mean_l2 = mean(l2_error, na.rm = TRUE),
    q25_l2 = quantile(l2_error, 0.25, na.rm = TRUE),
    q75_l2 = quantile(l2_error, 0.75, na.rm = TRUE),
    fail_rate = mean(!is.finite(l2_error) | is.na(l2_error)),
    median_eta_sd = median(eta_sd, na.rm = TRUE),
    median_eta_abs = median(eta_median_abs, na.rm = TRUE),
    median_eta_gt3_rate = median(eta_gt3_rate, na.rm = TRUE),
    median_bernoulli_var = median(mean_bernoulli_var, na.rm = TRUE),
    median_runtime_sec = median(runtime_sec, na.rm = TRUE),
    K1 = suppressWarnings(first(stats::na.omit(K1))),
    K2 = suppressWarnings(first(stats::na.omit(K2))),
    K1_block_size = suppressWarnings(first(stats::na.omit(K1_block_size))),
    K2_block_size = suppressWarnings(first(stats::na.omit(K2_block_size))),
    .groups = "drop"
  )

write.csv(summary_results, file.path(summary_dir, "summary_signal_strength_scan_logistic.csv"), row.names = FALSE)
saveRDS(summary_results, file.path(summary_dir, "summary_signal_strength_scan_logistic.rds"))

method_levels <- c("TS-RIGHT", "IHT", "Shrinkage", "Lasso")
summary_results$method <- factor(summary_results$method, levels = method_levels)

plot_median <- summary_results %>%
  filter(is.finite(median_l2), median_l2 > 0) %>%
  ggplot(aes(x = theta_magnitude, y = median_l2, color = method, shape = method, group = method)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.5) +
  facet_wrap(~ n, labeller = label_both) +
  scale_y_log10() +
  labs(
    x = expression(theta^"* magnitude"),
    y = "Median L2 Error (log scale)",
    color = "Method",
    shape = "Method"
  ) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )

print(plot_median)
ggsave(file.path(fig_dir, "logistic_signal_scan_median_l2.pdf"), plot_median, width = 8.5, height = 5)
ggsave(file.path(fig_dir, "logistic_signal_scan_median_l2.png"), plot_median, width = 8.5, height = 5, dpi = 300)

plot_eta <- summary_results %>%
  distinct(n, theta_magnitude, median_eta_sd, median_eta_abs, median_eta_gt3_rate, median_bernoulli_var) %>%
  ggplot(aes(x = theta_magnitude, y = median_eta_sd, group = factor(n), color = factor(n))) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.3) +
  labs(
    x = expression(theta^"* magnitude"),
    y = "Median sd(X theta*)",
    color = "n"
  ) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )

print(plot_eta)
ggsave(file.path(fig_dir, "logistic_signal_scan_eta_sd.pdf"), plot_eta, width = 7, height = 4.5)

message("Done. Results saved under: ", out_dir)
