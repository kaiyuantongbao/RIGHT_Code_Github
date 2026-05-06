rm(list = ls())

library(here)
library(dplyr)
library(foreach)
library(doParallel)
library(ggplot2)

source(here("simulations", "07_schedule_ablation", "utils_schedule_right.R"))

# ======================================================================
# E2: linear s-sweep, TS-RIGHT vs budget-matched FD-RIGHT vs FS-RIGHT
#
# Main changes relative to the old script:
#   1. TS split uses q = n2 / n, not c_r * s_star / log(s_star).
#   2. For each s_star configuration, q, T1, T2 are selected by an
#      independent pilot tuning procedure and then frozen.
#   3. The selected schedules are saved in a registry and reused by the
#      final Monte Carlo runs.
#   4. Linear pilot tuning uses robust validation loss:
#          median_i | y_i - x_i^T theta |.
# ======================================================================

# ----------------------------------------------------------------------
# Output locations
# ----------------------------------------------------------------------

out_dir <- "results/07_schedule_ablation/E2_linear_s_sweep"
registry_dir <- "results/00_tuning_registry"
registry_file <- file.path(registry_dir, "ts_schedule_registry.csv")
candidate_dir <- file.path(registry_dir, "candidates")

dir.create(file.path(out_dir, "raw"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "summary"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(registry_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(candidate_dir, recursive = TRUE, showWarnings = FALSE)

# ----------------------------------------------------------------------
# Experiment configuration
# ----------------------------------------------------------------------

cfg <- list(
  experiment_id = "E2_linear_s_sweep",
  model = "linear",
  
  n = 3200,
  p = 600,
  s_star_grid = c(5, 15, 20, 30),
  s_multiplier = 2,
  theta_magnitude = 5,
  
  eta = 0.02,
  m_min = 10,
  
  # TS schedule tuning grid.
  q_grid = c(0.1, 0.2, 0.25, 0.5),
  T1_grid = c(128, 150, 200),
  T2_grid = c(8, 12),
  
  # Tuning rule: select the minimum-budget schedule among candidates whose
  # pilot median validation loss is within 5% of the best pilot median.
  n_pilot = 15,
  selection_alpha = 1.05,
  selection_rule = "within_1.05_best_then_min_budget",
  validation_metric = "median_abs_prediction_error",
  enforce_tail_lower = TRUE,
  relax_tail_if_empty = TRUE,
  
  # Main Monte Carlo repetitions after schedule selection is frozen.
  reps = 100,
  
  # Heavy-tailed design/noise.
  df_X = 2.5,
  scale_X = 1,
  df_eps = 1.5,
  scale_eps = 1,
  
  # Seed blocks are deliberately separated so pilot tuning and final
  # evaluation do not reuse the same data or algorithmic randomness.
  pilot_data_seed_base = 20260420L,
  pilot_val_seed_base = 20270420L,
  pilot_algo_seed_base = 20280420L,
  final_data_seed_base = 20290420L,
  final_algo_seed_base = 20300420L,
  
  # Parallelism. Reduce this manually if memory is tight.
  num_cores = max(1L, parallel::detectCores() - 1L)
)

# ----------------------------------------------------------------------
# Tuning helper functions
# ----------------------------------------------------------------------

sanitize_id <- function(x) {
  x <- gsub("\\.", "p", as.character(x))
  x <- gsub("[^A-Za-z0-9_\\u002D]+", "_", x)
  x
}

# A small base-R hash to keep file names short and Windows-safe.
# We avoid using all tuning parameters in candidate file names because Windows
# may reject long absolute paths. The full metadata are stored in
# config_signature and in the registry CSV.
stable_hash <- function(x) {
  txt <- paste(as.character(x), collapse = "||")
  ints <- utf8ToInt(txt)
  h <- 5381
  mod <- 2147483647
  if (length(ints) > 0L) {
    for (v in ints) {
      h <- (h * 33 + v) %% mod
    }
  }
  sprintf("%08x", as.integer(h))
}

safe_write_csv <- function(x, file) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  write.csv(x, file, row.names = FALSE)
}

make_grid_signature <- function(cfg) {
  paste0(
    "q", paste(sanitize_id(cfg$q_grid), collapse = "-"),
    "__T1", paste(sanitize_id(cfg$T1_grid), collapse = "-"),
    "__T2", paste(sanitize_id(cfg$T2_grid), collapse = "-"),
    "__tail", as.integer(isTRUE(cfg$enforce_tail_lower)),
    "__alpha", sanitize_id(cfg$selection_alpha)
  )
}

make_config_signature <- function(cfg, s_star, s, s_ref) {
  paste(
    paste0("experiment_id=", cfg$experiment_id),
    paste0("model=", cfg$model),
    paste0("n=", cfg$n),
    paste0("p=", cfg$p),
    paste0("s_star=", s_star),
    paste0("s=", s),
    paste0("s_ref=", s_ref),
    paste0("df_X=", cfg$df_X),
    paste0("df_eps=", cfg$df_eps),
    paste0("scale_X=", cfg$scale_X),
    paste0("scale_eps=", cfg$scale_eps),
    paste0("theta_magnitude=", cfg$theta_magnitude),
    paste0("eta=", cfg$eta),
    paste0("m_min=", cfg$m_min),
    paste0("validation_metric=", cfg$validation_metric),
    paste0("n_pilot=", cfg$n_pilot),
    paste0("selection_alpha=", cfg$selection_alpha),
    paste0("selection_rule=", cfg$selection_rule),
    paste0("grid_signature=", make_grid_signature(cfg)),
    sep = "||"
  )
}

make_config_id <- function(cfg, s_star, s, s_ref) {
  signature <- make_config_signature(cfg = cfg, s_star = s_star, s = s, s_ref = s_ref)
  paste0(
    sanitize_id(cfg$experiment_id),
    "_n", cfg$n,
    "_p", cfg$p,
    "_sstar", s_star,
    "_s", s,
    "_h", stable_hash(signature)
  )
}

make_candidate_file <- function(candidate_dir, prefix, config_id) {
  file.path(candidate_dir, paste0(prefix, "_", config_id, ".csv"))
}

diagnose_ts_grid <- function(
    grid,
    n,
    p,
    s,
    s_ref = s,
    m_min = 10,
    c_K1 = 1,
    c_K2 = 1,
    enforce_tail_lower = TRUE
) {
  required_cols <- c("q", "T1", "T2")
  missing_cols <- setdiff(required_cols, names(grid))
  if (length(missing_cols) > 0L) {
    stop("grid is missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  s_ref_eff <- as.integer(s_ref)
  
  out <- lapply(seq_len(nrow(grid)), function(ii) {
    q_i <- grid$q[ii]
    T1_i <- grid$T1[ii]
    T2_i <- grid$T2[ii]
    
    budget <- compute_ts_budget(
      n = n,
      p = p,
      s = s,
      q = q_i,
      T1 = T1_i,
      T2 = T2_i,
      s_ref = s_ref_eff,
      m_min = m_min,
      c_K1 = c_K1,
      c_K2 = c_K2
    )
    
    q_min_tail <- log(max(s_ref_eff, 2L)) / max(s_ref_eff, 2L)
    q_min_K2 <- T2_i * m_min * ceiling(c_K2 * log(p)) / n
    q_max_K1 <- 1 - m_min * ceiling(c_K1 * s_ref_eff * log(p)) / n
    
    feasible_tail <- if (isTRUE(enforce_tail_lower)) q_i >= q_min_tail else TRUE
    feasible_K2 <- q_i >= q_min_K2
    feasible_K1 <- q_i <= q_max_K1
    feasible_blocks <- feasible_K1 && feasible_K2
    feasible <- feasible_tail && feasible_blocks && isTRUE(budget$eligible)
    
    data.frame(
      q = q_i,
      T1 = T1_i,
      T2 = T2_i,
      s = as.integer(s),
      s_ref = s_ref_eff,
      q_min_tail = q_min_tail,
      q_min_K2 = q_min_K2,
      q_max_K1 = q_max_K1,
      feasible_tail = feasible_tail,
      feasible_K1 = feasible_K1,
      feasible_K2 = feasible_K2,
      feasible_blocks = feasible_blocks,
      feasible = feasible,
      n1 = budget$n1,
      n2 = budget$n2,
      b2 = budget$b2,
      K1_tar = budget$K1_tar,
      K1 = budget$K1,
      K2_tar = budget$K2_tar,
      K2 = budget$K2,
      T_fd_budget = budget$T_fd_budget,
      cost_proxy = budget$cost_proxy,
      budget_eligible = budget$eligible
    )
  })
  
  do.call(rbind, out)
}

filter_ts_grid_for_tuning <- function(
    grid,
    n,
    p,
    s,
    s_ref = s,
    m_min = 10,
    c_K1 = 1,
    c_K2 = 1,
    enforce_tail_lower = TRUE,
    relax_tail_if_empty = TRUE
) {
  diag <- diagnose_ts_grid(
    grid = grid,
    n = n,
    p = p,
    s = s,
    s_ref = s_ref,
    m_min = m_min,
    c_K1 = c_K1,
    c_K2 = c_K2,
    enforce_tail_lower = enforce_tail_lower
  )
  
  out <- diag[diag$feasible, , drop = FALSE]
  
  if (nrow(out) == 0L && enforce_tail_lower && relax_tail_if_empty) {
    warning(
      "No q candidate passed the tail lower bound for this configuration. ",
      "Relaxing only the tail lower bound, while keeping K1/K2 block feasibility."
    )
    
    diag_relaxed <- diagnose_ts_grid(
      grid = grid,
      n = n,
      p = p,
      s = s,
      s_ref = s_ref,
      m_min = m_min,
      c_K1 = c_K1,
      c_K2 = c_K2,
      enforce_tail_lower = FALSE
    )
    
    out <- diag_relaxed[diag_relaxed$feasible_blocks & diag_relaxed$budget_eligible, , drop = FALSE]
  }
  
  if (nrow(out) == 0L) {
    stop(
      "No feasible TS schedule remains after filtering. ",
      "Consider increasing n, decreasing T2, lowering m_min, or enlarging q_grid."
    )
  }
  
  out
}

linear_val_median_abs <- function(theta, X_val, y_val) {
  pred <- as.vector(X_val %*% theta)
  median(abs(y_val - pred), na.rm = TRUE)
}

# The current runners return a data frame, not theta_hat. For validation-based
# pilot tuning, we need theta_hat, so this local fitter mirrors run_ts_right()
# for the linear model and returns the final vector.
fit_ts_right_linear_theta <- function(
    X,
    y,
    s,
    eta,
    q,
    T1,
    T2,
    m_min = 10,
    s_ref,
    seed = NULL
) {
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
  
  if (!isTRUE(budget$eligible)) {
    return(list(
      theta = rep(NA_real_, p),
      eligible = FALSE,
      budget = budget,
      runtime_sec = NA_real_
    ))
  }
  
  rng_state <- with_local_seed(seed)
  on.exit(restore_local_seed(rng_state), add = TRUE)
  
  t0 <- proc.time()[3]
  
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
      T_max = 1,
      K = budget$K2,
      theta_init = theta_cur,
      grad_func_samplewise = grad_linear_regression_samplewise,
      record_trace = FALSE
    )
  }
  
  runtime_sec <- proc.time()[3] - t0
  
  list(
    theta = theta_cur,
    eligible = TRUE,
    budget = budget,
    runtime_sec = runtime_sec
  )
}

read_schedule_registry <- function(registry_file) {
  if (!file.exists(registry_file)) {
    return(data.frame())
  }
  read.csv(registry_file, stringsAsFactors = FALSE)
}

append_schedule_registry <- function(row, registry_file) {
  if (!file.exists(registry_file)) {
    safe_write_csv(row, registry_file)
  } else {
    existing <- read.csv(registry_file, stringsAsFactors = FALSE)
    out <- dplyr::bind_rows(existing, row)
    safe_write_csv(out, registry_file)
  }
}

select_ts_schedule <- function(candidate_summary, cfg) {
  finite_summary <- candidate_summary %>%
    dplyr::filter(is.finite(pilot_median_val_loss))
  
  if (nrow(finite_summary) == 0L) {
    stop("All TS tuning candidates failed or produced non-finite validation losses.")
  }
  
  best_loss <- min(finite_summary$pilot_median_val_loss, na.rm = TRUE)
  
  finite_summary %>%
    dplyr::mutate(
      within_best = pilot_median_val_loss <= cfg$selection_alpha * best_loss,
      tie_distance = abs(q - 0.25) / 0.25 +
        abs(T1 - 150) / 150 +
        abs(T2 - 8) / 8
    ) %>%
    dplyr::filter(within_best) %>%
    dplyr::arrange(T_fd_budget, tie_distance, q, T1, T2) %>%
    dplyr::slice(1) %>%
    dplyr::mutate(best_pilot_median_val_loss = best_loss)
}

run_or_load_ts_tuning <- function(
    cfg,
    s_star,
    registry_file,
    candidate_dir
) {
  s <- as.integer(cfg$s_multiplier * s_star)
  theta_star <- make_theta_star(cfg$p, s_star, magnitude = cfg$theta_magnitude)
  s_ref <- choose_s_ref(s = s, theta_star = theta_star)
  config_signature <- make_config_signature(cfg = cfg, s_star = s_star, s = s, s_ref = s_ref)
  config_id <- make_config_id(cfg = cfg, s_star = s_star, s = s, s_ref = s_ref)
  
  registry <- read_schedule_registry(registry_file)
  if (nrow(registry) > 0L && "config_id" %in% names(registry)) {
    matches <- registry[registry$config_id == config_id, , drop = FALSE]
    if (nrow(matches) > 0L && "config_signature" %in% names(matches)) {
      matches <- matches[matches$config_signature == config_signature, , drop = FALSE]
    }
    if (nrow(matches) > 0L) {
      message("Using frozen TS schedule from registry for config_id = ", config_id)
      return(matches[1, , drop = FALSE])
    }
  }
  
  message("Running pilot TS tuning for config_id = ", config_id)
  
  full_grid <- make_ts_grid(
    q_grid = cfg$q_grid,
    T1_grid = cfg$T1_grid,
    T2_grid = cfg$T2_grid
  )
  
  grid_filt <- filter_ts_grid_for_tuning(
    grid = full_grid,
    n = cfg$n,
    p = cfg$p,
    s = s,
    s_ref = s_ref,
    m_min = cfg$m_min,
    enforce_tail_lower = cfg$enforce_tail_lower,
    relax_tail_if_empty = cfg$relax_tail_if_empty
  )
  
  grid_diag_file <- make_candidate_file(candidate_dir, "grid_diagnostics", config_id)
  safe_write_csv(grid_filt, grid_diag_file)
  
  pilot_results <- foreach(
    rep_id = seq_len(cfg$n_pilot),
    .combine = dplyr::bind_rows,
    .packages = c("dplyr", "MASS", "mvtnorm", "here"),
    .export = c(
      "fit_ts_right_linear_theta",
      "linear_val_median_abs",
      "cfg",
      "config_id",
      "config_signature",
      "grid_filt",
      "s_star",
      "s",
      "s_ref"
    )
  ) %dopar% {
    source(here("simulations", "07_schedule_ablation", "utils_schedule_right.R"))
    
    theta_star_local <- make_theta_star(cfg$p, s_star, magnitude = cfg$theta_magnitude)
    s_index <- match(s_star, cfg$s_star_grid)
    
    dat_train <- gen_one_linear_dataset(
      seed = cfg$pilot_data_seed_base + 100000L * s_index + rep_id,
      n = cfg$n,
      p = cfg$p,
      theta_star = theta_star_local,
      df_X = cfg$df_X,
      Sigma_X = diag(cfg$p),
      scale_X = cfg$scale_X,
      df_eps = cfg$df_eps,
      scale_eps = cfg$scale_eps
    )
    
    dat_val <- gen_one_linear_dataset(
      seed = cfg$pilot_val_seed_base + 100000L * s_index + rep_id,
      n = cfg$n,
      p = cfg$p,
      theta_star = theta_star_local,
      df_X = cfg$df_X,
      Sigma_X = diag(cfg$p),
      scale_X = cfg$scale_X,
      df_eps = cfg$df_eps,
      scale_eps = cfg$scale_eps
    )
    
    one_rep <- lapply(seq_len(nrow(grid_filt)), function(ii) {
      cand <- grid_filt[ii, , drop = FALSE]
      algo_seed <- cfg$pilot_algo_seed_base + 100000L * s_index + rep_id
      
      fit <- tryCatch(
        fit_ts_right_linear_theta(
          X = dat_train$X,
          y = dat_train$y,
          s = s,
          eta = cfg$eta,
          q = cand$q,
          T1 = cand$T1,
          T2 = cand$T2,
          m_min = cfg$m_min,
          s_ref = s_ref,
          seed = algo_seed
        ),
        error = function(e) {
          list(
            theta = rep(NA_real_, cfg$p),
            eligible = FALSE,
            budget = NULL,
            runtime_sec = NA_real_,
            error_message = conditionMessage(e)
          )
        }
      )
      
      val_loss <- if (isTRUE(fit$eligible) && all(is.finite(fit$theta))) {
        linear_val_median_abs(fit$theta, dat_val$X, dat_val$y)
      } else {
        NA_real_
      }
      
      data.frame(
        config_id = config_id,
        config_signature = config_signature,
        s_star = s_star,
        s = s,
        s_ref = s_ref,
        pilot_rep_id = rep_id,
        q = cand$q,
        T1 = cand$T1,
        T2 = cand$T2,
        n1 = cand$n1,
        n2 = cand$n2,
        b2 = cand$b2,
        K1 = cand$K1,
        K2 = cand$K2,
        T_fd_budget = cand$T_fd_budget,
        cost_proxy = cand$cost_proxy,
        val_loss = val_loss,
        eligible = isTRUE(fit$eligible),
        runtime_sec = fit$runtime_sec,
        seed_data = cfg$pilot_data_seed_base + 100000L * s_index + rep_id,
        seed_val = cfg$pilot_val_seed_base + 100000L * s_index + rep_id,
        seed_algo = algo_seed,
        error_message = if (!is.null(fit$error_message)) fit$error_message else NA_character_
      )
    })
    
    dplyr::bind_rows(one_rep)
  }
  
  candidate_raw_file <- make_candidate_file(candidate_dir, "candidate_raw", config_id)
  safe_write_csv(pilot_results, candidate_raw_file)
  
  candidate_summary <- pilot_results %>%
    dplyr::group_by(q, T1, T2, n1, n2, b2, K1, K2, T_fd_budget, cost_proxy) %>%
    dplyr::summarise(
      n_pilot = dplyr::n(),
      pilot_median_val_loss = median(val_loss, na.rm = TRUE),
      pilot_mean_val_loss = mean(val_loss, na.rm = TRUE),
      pilot_iqr_val_loss = IQR(val_loss, na.rm = TRUE),
      pilot_fail_rate = mean(!eligible | !is.finite(val_loss)),
      pilot_median_runtime = median(runtime_sec, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::arrange(pilot_median_val_loss, T_fd_budget, q, T1, T2)
  
  candidate_summary_file <- make_candidate_file(candidate_dir, "candidate_summary", config_id)
  safe_write_csv(candidate_summary, candidate_summary_file)
  
  selected <- select_ts_schedule(candidate_summary, cfg)
  
  registry_row <- data.frame(
    config_id = config_id,
    config_signature = config_signature,
    experiment_id = cfg$experiment_id,
    model = cfg$model,
    n = cfg$n,
    p = cfg$p,
    s_star = s_star,
    s = s,
    s_ref = s_ref,
    df_X = cfg$df_X,
    scale_X = cfg$scale_X,
    df_eps = cfg$df_eps,
    scale_eps = cfg$scale_eps,
    theta_magnitude = cfg$theta_magnitude,
    eta = cfg$eta,
    m_min = cfg$m_min,
    validation_metric = cfg$validation_metric,
    n_pilot = cfg$n_pilot,
    selection_alpha = cfg$selection_alpha,
    selection_rule = cfg$selection_rule,
    enforce_tail_lower = cfg$enforce_tail_lower,
    q_selected = selected$q,
    T1_selected = selected$T1,
    T2_selected = selected$T2,
    n1 = selected$n1,
    n2 = selected$n2,
    b2 = selected$b2,
    K1 = selected$K1,
    K2 = selected$K2,
    T_fd_budget = selected$T_fd_budget,
    cost_proxy = selected$cost_proxy,
    pilot_median_val_loss = selected$pilot_median_val_loss,
    best_pilot_median_val_loss = selected$best_pilot_median_val_loss,
    pilot_iqr_val_loss = selected$pilot_iqr_val_loss,
    pilot_fail_rate = selected$pilot_fail_rate,
    pilot_median_runtime = selected$pilot_median_runtime,
    grid_signature = make_grid_signature(cfg),
    candidate_raw_file = candidate_raw_file,
    candidate_summary_file = candidate_summary_file,
    grid_diag_file = grid_diag_file,
    stringsAsFactors = FALSE
  )
  
  append_schedule_registry(registry_row, registry_file)
  registry_row
}

# ----------------------------------------------------------------------
# Parallel backend
# ----------------------------------------------------------------------

num_cores <- as.integer(cfg$num_cores)
message("Using ", num_cores, " parallel workers.")
cl <- parallel::makeCluster(num_cores)
doParallel::registerDoParallel(cl)
on.exit(parallel::stopCluster(cl), add = TRUE)

# ----------------------------------------------------------------------
# Stage A: tune TS schedule for each s_star, or load frozen schedule
# ----------------------------------------------------------------------

schedule_tbl <- dplyr::bind_rows(lapply(cfg$s_star_grid, function(s_star) {
  run_or_load_ts_tuning(
    cfg = cfg,
    s_star = s_star,
    registry_file = registry_file,
    candidate_dir = candidate_dir
  )
}))

schedule_out_file <- file.path(out_dir, "summary", "selected_TS_schedules_E2_linear_s_sweep.csv")
safe_write_csv(schedule_tbl, schedule_out_file)
print(schedule_tbl[, c(
  "s_star", "s", "s_ref", "q_selected", "T1_selected", "T2_selected",
  "n1", "n2", "b2", "K1", "K2", "T_fd_budget",
  "pilot_median_val_loss", "best_pilot_median_val_loss", "pilot_fail_rate"
)])

# ----------------------------------------------------------------------
# Stage B: final Monte Carlo with frozen schedules
# ----------------------------------------------------------------------

raw_results <- foreach(
  s_star = cfg$s_star_grid,
  .combine = dplyr::bind_rows,
  .packages = c("dplyr", "MASS", "mvtnorm", "here")
) %:% foreach(
  rep_id = seq_len(cfg$reps),
  .combine = dplyr::bind_rows,
  .packages = c("dplyr", "MASS", "mvtnorm", "here")
) %dopar% {
  source(here("simulations", "07_schedule_ablation", "utils_schedule_right.R"))
  
  s_index <- match(s_star, cfg$s_star_grid)
  sched <- schedule_tbl[schedule_tbl$s_star == s_star, , drop = FALSE]
  if (nrow(sched) != 1L) stop("Expected exactly one frozen schedule for s_star = ", s_star)
  
  s <- as.integer(sched$s)
  s_ref <- as.integer(sched$s_ref)
  theta_star <- make_theta_star(cfg$p, s_star, magnitude = cfg$theta_magnitude)
  
  seed_data <- cfg$final_data_seed_base + 100000L * s_index + rep_id
  seed_algo_ts <- cfg$final_algo_seed_base + 100000L * s_index + rep_id
  seed_algo_fs <- cfg$final_algo_seed_base + 200000L * s_index + rep_id
  
  dat <- gen_one_linear_dataset(
    seed = seed_data,
    n = cfg$n,
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
  
  ts_out <- run_ts_right_linear(
    X = X,
    y = y,
    theta_star = theta_star,
    s = s,
    eta = cfg$eta,
    q = sched$q_selected,
    T1 = sched$T1_selected,
    T2 = sched$T2_selected,
    m_min = cfg$m_min,
    s_ref = s_ref,
    record_trace = FALSE,
    record_l2 = TRUE,
    record_initial = FALSE,
    seed = seed_algo_ts
  ) %>%
    dplyr::mutate(method = "TS")
  
  fd_out <- run_fd_right_linear(
    X = X,
    y = y,
    theta_star = theta_star,
    s = s,
    eta = cfg$eta,
    T_fd = sched$T_fd_budget,
    m_min = cfg$m_min,
    s_ref = s_ref,
    record_trace = FALSE,
    record_l2 = TRUE,
    record_initial = FALSE
  ) %>%
    dplyr::mutate(method = "FD_budget")
  
  # Fully split baseline: use one fresh batch per iteration. The default
  # batch size is approximately m_min * ceil(log p), so K_fs can reach log p.
  T_fs <- max(1L, floor(cfg$n / (cfg$m_min * ceiling(log(cfg$p)))))
  
  fs_out <- run_fs_right_linear(
    X = X,
    y = y,
    theta_star = theta_star,
    s = s,
    eta = cfg$eta,
    T_fs = T_fs,
    m_min = cfg$m_min,
    record_trace = FALSE,
    record_l2 = TRUE,
    record_initial = FALSE,
    seed = seed_algo_fs
  ) %>%
    dplyr::mutate(method = "FS")
  
  dplyr::bind_rows(ts_out, fd_out, fs_out) %>%
    dplyr::mutate(
      experiment_id = cfg$experiment_id,
      config_id = sched$config_id,
      n = cfg$n,
      p = cfg$p,
      s_star = s_star,
      s = s,
      s_ref = s_ref,
      rep_id = rep_id,
      q_selected = sched$q_selected,
      T1_selected = sched$T1_selected,
      T2_selected = sched$T2_selected,
      T_fd_budget = sched$T_fd_budget,
      T_fs = T_fs,
      seed_data = seed_data,
      seed_algo_ts = seed_algo_ts,
      seed_algo_fs = seed_algo_fs,
      df_X = cfg$df_X,
      df_eps = cfg$df_eps,
      theta_magnitude = cfg$theta_magnitude,
      validation_metric = cfg$validation_metric,
      tuning_selection_rule = cfg$selection_rule
    )
}

raw_file <- file.path(out_dir, "raw", "raw_E2_linear_s_sweep_TS_tuned.rds")
saveRDS(raw_results, raw_file)

raw_csv_file <- file.path(out_dir, "raw", "raw_E2_linear_s_sweep_TS_tuned.csv")
safe_write_csv(raw_results, raw_csv_file)

# ----------------------------------------------------------------------
# Summary and figure
# ----------------------------------------------------------------------

summary_tbl <- raw_results %>%
  dplyr::group_by(
    s_star, s, s_ref, method,
    q_selected, T1_selected, T2_selected, T_fd_budget
  ) %>%
  dplyr::summarise(
    reps = dplyr::n(),
    fail_rate = mean(!eligible | !is.finite(l2_error)),
    median_l2 = median(l2_error, na.rm = TRUE),
    mean_l2 = mean(l2_error, na.rm = TRUE),
    q25_l2 = as.numeric(stats::quantile(l2_error, 0.25, na.rm = TRUE)),
    q75_l2 = as.numeric(stats::quantile(l2_error, 0.75, na.rm = TRUE)),
    iqr_l2 = IQR(l2_error, na.rm = TRUE),
    median_runtime = median(runtime_sec, na.rm = TRUE),
    mean_runtime = mean(runtime_sec, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(s_star, method)

summary_file <- file.path(out_dir, "summary", "summary_E2_linear_s_sweep_TS_tuned.csv")
safe_write_csv(summary_tbl, summary_file)
print(summary_tbl)

# Optional: TS improvement over budget-matched FD. Avoid tidyr dependency.
improvement_tbl <- Reduce(
  function(left, right) merge(left, right, by = "s_star", all = TRUE),
  list(
    summary_tbl[summary_tbl$method == "TS", c("s_star", "median_l2")],
    summary_tbl[summary_tbl$method == "FD_budget", c("s_star", "median_l2")],
    summary_tbl[summary_tbl$method == "FS", c("s_star", "median_l2")]
  )
)
names(improvement_tbl) <- c("s_star", "TS", "FD_budget", "FS")
improvement_tbl$TS_vs_FD_budget_reduction <-
  (improvement_tbl$FD_budget - improvement_tbl$TS) / improvement_tbl$FD_budget
improvement_tbl$TS_vs_FS_reduction <-
  (improvement_tbl$FS - improvement_tbl$TS) / improvement_tbl$FS

improvement_file <- file.path(out_dir, "summary", "improvement_E2_linear_s_sweep_TS_tuned.csv")
safe_write_csv(improvement_tbl, improvement_file)
print(improvement_tbl)

p <- ggplot(summary_tbl, aes(x = s_star, y = median_l2, color = method, group = method)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.2) +
  scale_y_log10() +
  labs(
    #title = "Linear s-sweep: TS-RIGHT vs FD-budget vs FS-RIGHT",
    #subtitle = "TS schedules are selected by independent pilot validation and frozen before final evaluation.",
    x = "True sparsity s*",
    y = "Median L2 error (log scale)",
    color = "Method"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "bottom"
  )

print(p)

ggsave(
  filename = file.path(out_dir, "figures", "linear_s_sweep_TS_tuned_median_l2.pdf"),
  plot = p,
  width = 7,
  height = 5
)

ggsave(
  filename = file.path(out_dir, "figures", "linear_s_sweep_TS_tuned_median_l2.png"),
  plot = p,
  width = 7,
  height = 5,
  dpi = 300
)

message("Done. Raw results saved to: ", raw_file)
message("Summary saved to: ", summary_file)
message("Selected schedules saved to: ", schedule_out_file)
message("Registry file: ", registry_file)
