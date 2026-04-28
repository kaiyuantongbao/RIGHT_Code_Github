# ======================================================================
# Utilities for schedule ablation experiments: FD-RIGHT / TS-RIGHT / FS-RIGHT
# Supports both sparse linear regression and sparse logistic regression.
#
# Main changes relative to the old utils_schedule_linear.R:
#   1. TS sample split uses fixed-fraction q: n2 = floor(q * n), n1 = n - n2.
#   2. Linear/logistic model differences are handled through get_model_spec().
#   3. FD/TS/FS runners are generic, with thin linear/logistic aliases.
#   4. FS supports record_trace / record_l2 / record_initial.
#   5. Trace output includes stage, sample_budget_used, and fd_equiv_iter.
# ======================================================================

library(dplyr)
library(MASS)
library(mvtnorm)
library(foreach)
library(doParallel)
library(here)

source(here("R", "data_generator.R"))
source(here("R", "solvers.R"))

# ----------------------------------------------------------------------
# Common config helpers
# ----------------------------------------------------------------------

make_theta_star <- function(p, s_star, magnitude = 5) {
  if (!is.finite(p) || p < 1L) stop("p must be a positive integer.")
  if (!is.finite(s_star) || s_star < 0L || s_star > p) {
    stop("s_star must be an integer in [0, p].")
  }
  
  theta_star <- rep(0, p)
  if (s_star > 0L) {
    theta_star[seq_len(s_star)] <- magnitude
  }
  theta_star
}

choose_s_ref <- function(s, theta_star = NULL, s_ref = NULL, tol = 1e-12) {
  if (!is.null(s_ref)) {
    out <- as.integer(s_ref)
  } else if (!is.null(theta_star)) {
    out <- as.integer(sum(abs(theta_star) > tol))
  } else {
    out <- as.integer(s)
  }
  
  if (!is.finite(out) || out < 1L) {
    stop("s_ref must be a positive integer. In simulation, pass theta_star; in real data, use s_ref = s.")
  }
  
  out
}

safe_l2_error <- function(theta_hat, theta_star) {
  if (is.null(theta_star)) return(NA_real_)
  sqrt(sum((theta_hat - theta_star)^2))
}

with_local_seed <- function(seed) {
  if (is.null(seed)) return(NULL)
  
  old_exists <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old_state <- if (old_exists) get(".Random.seed", envir = .GlobalEnv, inherits = FALSE) else NULL
  
  set.seed(seed)
  
  list(
    old_exists = old_exists,
    old_state = old_state
  )
}

restore_local_seed <- function(state) {
  if (is.null(state)) return(invisible(NULL))
  
  if (isTRUE(state$old_exists)) {
    assign(".Random.seed", state$old_state, envir = .GlobalEnv)
  } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    rm(".Random.seed", envir = .GlobalEnv)
  }
  
  invisible(NULL)
}

# ----------------------------------------------------------------------
# Model-specific pieces
# ----------------------------------------------------------------------

# The current R/solvers.R already defines grad_linear_regression_samplewise.
# Logistic sample-wise gradient is defined here so that logistic TS/FD/FS
# wrappers can share the same solver_right() interface.
grad_logistic_regression_samplewise <- function(theta, X, y) {
  linear_predictor <- X %*% theta
  probabilities <- plogis(linear_predictor)
  sample_gradients <- X * as.vector(probabilities - y)
  sample_gradients
}

get_model_spec <- function(model = c("linear", "logistic")) {
  model <- match.arg(model)
  
  if (model == "linear") {
    return(list(
      model = "linear",
      grad_func_samplewise = grad_linear_regression_samplewise
    ))
  }
  
  if (model == "logistic") {
    return(list(
      model = "logistic",
      grad_func_samplewise = grad_logistic_regression_samplewise
    ))
  }
  
  stop("Unknown model.")
}

# ----------------------------------------------------------------------
# Data generators used by simulations
# ----------------------------------------------------------------------

gen_one_linear_dataset <- function(
    seed,
    n,
    p,
    theta_star,
    df_X = 2.5,
    Sigma_X = NULL,
    scale_X = 1,
    df_eps = 1.5,
    scale_eps = 1
) {
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
    df_eps = df_eps,
    scale_eps = scale_eps
  )
  
  list(X = dat$X, y = as.numeric(dat$y))
}

gen_one_logistic_dataset <- function(
    seed,
    n,
    p,
    theta_star,
    df_X = 2.5,
    Sigma_X = NULL,
    scale_X = 1
) {
  set.seed(seed)
  if (is.null(Sigma_X)) Sigma_X <- diag(p)
  
  dat <- generate_logistic_data(
    n = n,
    p = p,
    theta_star = theta_star,
    generator_X = generate_X_t,
    df_X = df_X,
    Sigma_X = Sigma_X,
    scale_X = scale_X
  )
  
  list(X = dat$X, y = as.numeric(dat$y))
}

# ----------------------------------------------------------------------
# Schedule helpers
# ----------------------------------------------------------------------

split_q_to_n1n2 <- function(n, q) {
  if (!is.finite(n) || n < 2L) {
    stop("n must be at least 2.")
  }
  if (!is.finite(q) || q <= 0 || q >= 1) {
    stop("q must be in (0, 1). Here q is the Stage II refinement fraction n2 / n.")
  }
  
  n2 <- floor(q * n)
  n1 <- n - n2
  
  list(
    n1 = n1,
    n2 = n2,
    q = q,
    frac_stage1 = n1 / n,
    frac_stage2 = n2 / n
  )
}

compute_ts_budget <- function(
    n,
    p,
    s,
    q,
    T1,
    T2,
    s_ref = s,
    m_min = 10,
    c_K1 = 1,
    c_K2 = 1,
    K1_formula = c("s_log_p", "s_log_ep_over_s")
) {
  K1_formula <- match.arg(K1_formula)
  
  if (!is.finite(p) || p < 2L) stop("p must be at least 2.")
  if (!is.finite(s) || s < 1L) stop("s must be a positive integer.")
  if (!is.finite(s_ref) || s_ref < 1L) stop("s_ref must be a positive integer.")
  if (!is.finite(T1) || T1 < 1L) stop("T1 must be a positive integer.")
  if (!is.finite(T2) || T2 < 1L) stop("T2 must be a positive integer.")
  if (!is.finite(m_min) || m_min < 1L) stop("m_min must be a positive integer.")
  
  tmp <- split_q_to_n1n2(n = n, q = q)
  n1 <- tmp$n1
  n2 <- tmp$n2
  
  K1_tar <- switch(
    K1_formula,
    s_log_p = ceiling(c_K1 * s_ref * log(p)),
    s_log_ep_over_s = ceiling(c_K1 * s_ref * log(exp(1) * p / max(s_ref, 1)))
  )
  K1 <- min(K1_tar, floor(n1 / m_min))
  
  b2 <- floor(n2 / T2)
  
  K2_tar <- ceiling(c_K2 * log(p))
  K2 <- min(K2_tar, floor(b2 / m_min))
  
  cost_proxy <- T1 * n1 + T2 * b2
  T_fd_budget <- ceiling(cost_proxy / n)
  
  eligible <- (n1 >= m_min) &&
    (b2 >= m_min) &&
    (K1 >= 1L) &&
    (K2 >= 4L)
  
  list(
    n1 = n1,
    n2 = n2,
    q = q,
    frac_stage1 = tmp$frac_stage1,
    frac_stage2 = tmp$frac_stage2,
    
    s = as.integer(s),
    s_ref = as.integer(s_ref),
    
    K1_tar = K1_tar,
    K1 = K1,
    
    b2 = b2,
    K2_tar = K2_tar,
    K2 = K2,
    
    cost_proxy = cost_proxy,
    T_fd_budget = T_fd_budget,
    eligible = eligible
  )
}

make_ts_grid <- function(
    q_grid = c(0.1, 0.2, 0.25, 0.5),
    T1_grid = c(128, 150, 200),
    T2_grid = c(8, 12)
) {
  expand.grid(
    q = q_grid,
    T1 = T1_grid,
    T2 = T2_grid,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
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
    q_i  <- grid$q[ii]
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
    
    feasible_tail <- if (enforce_tail_lower) {
      q_i >= q_min_tail
    } else {
      TRUE
    }
    
    feasible_K2 <- q_i >= q_min_K2
    feasible_K1 <- q_i <= q_max_K1
    
    feasible_blocks <- feasible_K1 && feasible_K2
    feasible <- feasible_tail && feasible_blocks
    
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
      "No q candidate passed the tail lower bound. ",
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
    
    out <- diag_relaxed[diag_relaxed$feasible_blocks, , drop = FALSE]
  }
  
  if (nrow(out) == 0L) {
    stop(
      "No feasible TS schedule remains after filtering. ",
      "Consider increasing n, decreasing T2, lowering m_min, or enlarging q_grid."
    )
  }
  
  out
}
# ----------------------------------------------------------------------
# Trace helper: adds sample-access budget axis
# ----------------------------------------------------------------------

add_budget_axis_to_trace <- function(
    trace_t,
    n,
    arm,
    n1 = NULL,
    T1 = NULL,
    b2 = NULL,
    b_fs = NULL
) {
  trace_t <- as.numeric(trace_t)
  
  if (arm == "FD") {
    stage <- rep("full_data", length(trace_t))
    sample_budget_used <- trace_t * n
  } else if (arm == "TS") {
    if (is.null(n1) || is.null(T1) || is.null(b2)) {
      stop("For TS trace, n1, T1, and b2 must be provided.")
    }
    
    stage <- ifelse(
      trace_t == 0,
      "init",
      ifelse(trace_t <= T1, "stage1", "stage2")
    )
    
    sample_budget_used <- ifelse(
      trace_t == 0,
      0,
      ifelse(
        trace_t <= T1,
        trace_t * n1,
        T1 * n1 + (trace_t - T1) * b2
      )
    )
  } else if (arm == "FS") {
    if (is.null(b_fs)) {
      stop("For FS trace, b_fs must be provided.")
    }
    
    stage <- ifelse(trace_t == 0, "init", "fresh_batch")
    sample_budget_used <- trace_t * b_fs
  } else {
    stop("Unknown arm. Expected one of: FD, TS, FS.")
  }
  
  data.frame(
    stage = stage,
    sample_budget_used = as.numeric(sample_budget_used),
    fd_equiv_iter = as.numeric(sample_budget_used) / n
  )
}

# ----------------------------------------------------------------------
# Generic runners
# ----------------------------------------------------------------------

run_fd_right <- function(
    X,
    y,
    theta_star = NULL,
    s,
    eta,
    T_fd,
    model = c("linear", "logistic"),
    m_min = 10,
    s_ref = NULL,
    record_trace = FALSE,
    record_l2 = !is.null(theta_star),
    record_initial = FALSE,
    record_every = 1L,
    seed = NULL
) {
  model_spec <- get_model_spec(model)
  n <- nrow(X)
  p <- ncol(X)
  
  if (!is.finite(T_fd) || T_fd < 1L) stop("T_fd must be a positive integer.")
  if (!is.finite(m_min) || m_min < 1L) stop("m_min must be a positive integer.")
  record_every <- as.integer(record_every)
  if (record_every <= 0L) stop("record_every must be a positive integer.")
  
  s_ref_use <- choose_s_ref(s = s, theta_star = theta_star, s_ref = s_ref)
  
  K_fd_tar <- ceiling(s_ref_use * log(p))
  K_fd <- min(K_fd_tar, floor(n / m_min))
  eligible <- (n >= m_min) && (K_fd >= 1L)
  
  base_df <- data.frame(
    arm = "FD",
    model = model_spec$model,
    T_fd = T_fd,
    K_fd_tar = K_fd_tar,
    K_fd = K_fd,
    s_ref = s_ref_use,
    eligible = eligible,
    T_fd_budget = NA_integer_
  )
  
  if (!eligible) {
    base_df$runtime_sec <- NA_real_
    base_df$l2_error <- NA_real_
    if (record_trace) {
      base_df$t <- NA_integer_
      base_df$stage <- NA_character_
      base_df$sample_budget_used <- NA_real_
      base_df$fd_equiv_iter <- NA_real_
    }
    return(base_df)
  }
  
  if (record_l2 && is.null(theta_star)) {
    stop("record_l2 = TRUE requires theta_star.")
  }
  
  t0 <- proc.time()[3]
  out_fd <- solver_right(
    X = X,
    y = y,
    s = s,
    eta = eta,
    T_max = T_fd,
    K = K_fd,
    theta_init = rep(0, p),
    grad_func_samplewise = model_spec$grad_func_samplewise,
    theta_star = if (record_l2) theta_star else NULL,
    record_trace = record_trace,
    record_l2 = record_l2,
    record_initial = record_initial,
    record_every = record_every
  )
  runtime_sec <- proc.time()[3] - t0
  base_df$runtime_sec <- runtime_sec
  
  if (record_trace) {
    theta_hat <- out_fd$theta
    trace_t <- out_fd$trace$iteration
    trace_l2 <- if (record_l2 && "l2_error" %in% names(out_fd$trace)) {
      out_fd$trace$l2_error
    } else {
      rep(NA_real_, length(trace_t))
    }
    
    budget_trace <- add_budget_axis_to_trace(
      trace_t = trace_t,
      n = n,
      arm = "FD"
    )
    
    out_df <- base_df[rep(1, length(trace_t)), ]
    out_df$t <- trace_t
    out_df$stage <- budget_trace$stage
    out_df$sample_budget_used <- budget_trace$sample_budget_used
    out_df$fd_equiv_iter <- budget_trace$fd_equiv_iter
    out_df$l2_error <- trace_l2
    return(out_df)
  }
  
  theta_hat <- out_fd
  base_df$l2_error <- safe_l2_error(theta_hat, theta_star)
  base_df
}

run_ts_right <- function(
    X,
    y,
    theta_star = NULL,
    s,
    eta,
    q,
    T1,
    T2,
    model = c("linear", "logistic"),
    m_min = 10,
    s_ref = NULL,
    record_trace = FALSE,
    record_l2 = !is.null(theta_star),
    record_initial = FALSE,
    record_every = 1L,
    seed = NULL
) {
  model_spec <- get_model_spec(model)
  n <- nrow(X)
  p <- ncol(X)
  
  record_every <- as.integer(record_every)
  if (record_every <= 0L) stop("record_every must be a positive integer.")
  s_ref_use <- choose_s_ref(s = s, theta_star = theta_star, s_ref = s_ref)
  
  budget <- compute_ts_budget(
    n = n,
    p = p,
    s = s,
    q = q,
    T1 = T1,
    T2 = T2,
    s_ref = s_ref_use,
    m_min = m_min
  )
  
  base_df <- data.frame(
    arm = "TS",
    model = model_spec$model,
    q = q,
    frac_stage1 = budget$frac_stage1,
    frac_stage2 = budget$frac_stage2,
    T1 = T1,
    T2 = T2,
    n1 = budget$n1,
    n2 = budget$n2,
    b2 = budget$b2,
    s_ref = s_ref_use,
    K1_tar = budget$K1_tar,
    K1 = budget$K1,
    K2_tar = budget$K2_tar,
    K2 = budget$K2,
    cost_proxy = budget$cost_proxy,
    T_fd_budget = budget$T_fd_budget,
    eligible = budget$eligible
  )
  
  if (!budget$eligible) {
    base_df$runtime_sec <- NA_real_
    base_df$l2_error <- NA_real_
    if (record_trace) {
      base_df$t <- NA_integer_
      base_df$stage <- NA_character_
      base_df$sample_budget_used <- NA_real_
      base_df$fd_equiv_iter <- NA_real_
    }
    return(base_df)
  }
  
  if (record_l2 && is.null(theta_star)) {
    stop("record_l2 = TRUE requires theta_star.")
  }
  
  rng_state <- with_local_seed(seed)
  on.exit(restore_local_seed(rng_state), add = TRUE)
  
  idx <- sample.int(n)
  id1 <- idx[seq_len(budget$n1)]
  id2 <- idx[(budget$n1 + 1L):n]
  
  X1 <- X[id1, , drop = FALSE]
  y1 <- y[id1]
  
  t0 <- proc.time()[3]
  
  # Stage I: full-data localization on D1.
  out_s1 <- solver_right(
    X = X1,
    y = y1,
    s = s,
    eta = eta,
    T_max = T1,
    K = budget$K1,
    theta_init = rep(0, p),
    grad_func_samplewise = model_spec$grad_func_samplewise,
    theta_star = if (record_l2) theta_star else NULL,
    record_trace = record_trace,
    record_l2 = record_l2,
    record_initial = record_initial,
    record_every = record_every
  )
  
  if (record_trace) {
    theta_cur <- out_s1$theta
    trace_t <- out_s1$trace$iteration
    trace_l2 <- if (record_l2 && "l2_error" %in% names(out_s1$trace)) {
      out_s1$trace$l2_error
    } else {
      rep(NA_real_, length(trace_t))
    }
  } else {
    theta_cur <- out_s1
    trace_t <- integer(0)
    trace_l2 <- numeric(0)
  }
  
  # Stage II: fresh-batch refinement on D2.
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
      grad_func_samplewise = model_spec$grad_func_samplewise,
      record_trace = FALSE
    )
    
    if (record_trace && (tt %% record_every == 0L || tt == T2)) {
      trace_t <- c(trace_t, T1 + tt)
      trace_l2 <- c(trace_l2, if (record_l2) safe_l2_error(theta_cur, theta_star) else NA_real_)
    }
  }
  
  runtime_sec <- proc.time()[3] - t0
  final_l2_error <- safe_l2_error(theta_cur, theta_star)
  
  base_df$runtime_sec <- runtime_sec
  
  if (record_trace) {
    budget_trace <- add_budget_axis_to_trace(
      trace_t = trace_t,
      n = n,
      arm = "TS",
      n1 = budget$n1,
      T1 = T1,
      b2 = budget$b2
    )
    
    out_df <- base_df[rep(1, length(trace_t)), ]
    out_df$t <- trace_t
    out_df$stage <- budget_trace$stage
    out_df$sample_budget_used <- budget_trace$sample_budget_used
    out_df$fd_equiv_iter <- budget_trace$fd_equiv_iter
    out_df$l2_error <- trace_l2
    return(out_df)
  }
  
  base_df$l2_error <- final_l2_error
  base_df
}

run_fs_right <- function(
    X,
    y,
    theta_star = NULL,
    s,
    eta,
    T_fs,
    model = c("linear", "logistic"),
    m_min = 10,
    record_trace = FALSE,
    record_l2 = !is.null(theta_star),
    record_initial = FALSE,
    record_every = 1L,
    seed = NULL
) {
  model_spec <- get_model_spec(model)
  n <- nrow(X)
  p <- ncol(X)
  
  if (!is.finite(T_fs) || T_fs < 1L) stop("T_fs must be a positive integer.")
  if (!is.finite(m_min) || m_min < 1L) stop("m_min must be a positive integer.")
  record_every <- as.integer(record_every)
  if (record_every <= 0L) stop("record_every must be a positive integer.")
  
  
  b_fs <- floor(n / T_fs)
  K_fs_tar <- ceiling(log(p))
  K_fs <- min(K_fs_tar, floor(b_fs / m_min))
  eligible <- (b_fs >= m_min) && (K_fs >= 4L)
  
  base_df <- data.frame(
    arm = "FS",
    model = model_spec$model,
    T_fs = T_fs,
    b_fs = b_fs,
    K_fs_tar = K_fs_tar,
    K_fs = K_fs,
    eligible = eligible,
    T_fd_budget = NA_integer_
  )
  
  if (!eligible) {
    base_df$runtime_sec <- NA_real_
    base_df$l2_error <- NA_real_
    if (record_trace) {
      base_df$t <- NA_integer_
      base_df$stage <- NA_character_
      base_df$sample_budget_used <- NA_real_
      base_df$fd_equiv_iter <- NA_real_
    }
    return(base_df)
  }
  
  if (record_l2 && is.null(theta_star)) {
    stop("record_l2 = TRUE requires theta_star.")
  }
  
  rng_state <- with_local_seed(seed)
  on.exit(restore_local_seed(rng_state), add = TRUE)
  
  idx <- sample.int(n)
  use_n <- b_fs * T_fs
  idx <- idx[seq_len(use_n)]
  
  theta_cur <- rep(0, p)
  trace_t <- integer(0)
  trace_l2 <- numeric(0)
  
  if (record_trace && record_initial) {
    trace_t <- c(trace_t, 0L)
    trace_l2 <- c(trace_l2, if (record_l2) safe_l2_error(theta_cur, theta_star) else NA_real_)
  }
  
  t0 <- proc.time()[3]
  
  for (tt in seq_len(T_fs)) {
    batch_ids <- idx[((tt - 1L) * b_fs + 1L):(tt * b_fs)]
    
    theta_cur <- solver_right(
      X = X[batch_ids, , drop = FALSE],
      y = y[batch_ids],
      s = s,
      eta = eta,
      T_max = 1,
      K = K_fs,
      theta_init = theta_cur,
      grad_func_samplewise = model_spec$grad_func_samplewise,
      record_trace = FALSE
    )
    
    if (record_trace && (tt %% record_every == 0L || tt == T_fs)) {
      trace_t <- c(trace_t, tt)
      trace_l2 <- c(trace_l2, if (record_l2) safe_l2_error(theta_cur, theta_star) else NA_real_)
    }
  }
  
  runtime_sec <- proc.time()[3] - t0
  final_l2_error <- safe_l2_error(theta_cur, theta_star)
  base_df$runtime_sec <- runtime_sec
  
  if (record_trace) {
    budget_trace <- add_budget_axis_to_trace(
      trace_t = trace_t,
      n = n,
      arm = "FS",
      b_fs = b_fs
    )
    
    out_df <- base_df[rep(1, length(trace_t)), ]
    out_df$t <- trace_t
    out_df$stage <- budget_trace$stage
    out_df$sample_budget_used <- budget_trace$sample_budget_used
    out_df$fd_equiv_iter <- budget_trace$fd_equiv_iter
    out_df$l2_error <- trace_l2
    return(out_df)
  }
  
  base_df$l2_error <- final_l2_error
  base_df
}

# ----------------------------------------------------------------------
# Thin aliases: linear regression
# ----------------------------------------------------------------------

run_fd_right_linear <- function(
    X,
    y,
    theta_star,
    s,
    eta,
    T_fd,
    m_min = 10,
    s_ref = NULL,
    record_trace = FALSE,
    record_l2 = !is.null(theta_star),
    record_initial = FALSE,
    record_every = 1L,
    seed = NULL
) {
  run_fd_right(
    X = X,
    y = y,
    theta_star = theta_star,
    s = s,
    eta = eta,
    T_fd = T_fd,
    model = "linear",
    m_min = m_min,
    s_ref = s_ref,
    record_trace = record_trace,
    record_l2 = record_l2,
    record_initial = record_initial,
    record_every = record_every,
    seed = seed
  )
}

run_ts_right_linear <- function(
    X,
    y,
    theta_star,
    s,
    eta,
    q,
    T1,
    T2,
    m_min = 10,
    s_ref = NULL,
    record_trace = FALSE,
    record_l2 = !is.null(theta_star),
    record_initial = FALSE,
    record_every = 1L,
    seed = NULL
) {
  run_ts_right(
    X = X,
    y = y,
    theta_star = theta_star,
    s = s,
    eta = eta,
    q = q,
    T1 = T1,
    T2 = T2,
    model = "linear",
    m_min = m_min,
    s_ref = s_ref,
    record_trace = record_trace,
    record_l2 = record_l2,
    record_initial = record_initial,
    record_every = record_every,
    seed = seed
  )
}

run_fs_right_linear <- function(
    X,
    y,
    theta_star,
    s,
    eta,
    T_fs,
    m_min = 10,
    record_trace = FALSE,
    record_l2 = !is.null(theta_star),
    record_initial = FALSE,
    record_every = 1L,
    seed = NULL
) {
  run_fs_right(
    X = X,
    y = y,
    theta_star = theta_star,
    s = s,
    eta = eta,
    T_fs = T_fs,
    model = "linear",
    m_min = m_min,
    record_trace = record_trace,
    record_l2 = record_l2,
    record_initial = record_initial,
    record_every = record_every,
    seed = seed
  )
}

# ----------------------------------------------------------------------
# Thin aliases: logistic regression
# ----------------------------------------------------------------------

run_fd_right_logistic <- function(
    X,
    y,
    theta_star,
    s,
    eta,
    T_fd,
    m_min = 10,
    s_ref = NULL,
    record_trace = FALSE,
    record_l2 = !is.null(theta_star),
    record_initial = FALSE,
    record_every = 1L,
    seed = NULL
) {
  run_fd_right(
    X = X,
    y = y,
    theta_star = theta_star,
    s = s,
    eta = eta,
    T_fd = T_fd,
    model = "logistic",
    m_min = m_min,
    s_ref = s_ref,
    record_trace = record_trace,
    record_l2 = record_l2,
    record_initial = record_initial,
    record_every = record_every,
    seed = seed
  )
}

run_ts_right_logistic <- function(
    X,
    y,
    theta_star,
    s,
    eta,
    q,
    T1,
    T2,
    m_min = 10,
    s_ref = NULL,
    record_trace = FALSE,
    record_l2 = !is.null(theta_star),
    record_initial = FALSE,
    record_every = 1L,
    seed = NULL
) {
  run_ts_right(
    X = X,
    y = y,
    theta_star = theta_star,
    s = s,
    eta = eta,
    q = q,
    T1 = T1,
    T2 = T2,
    model = "logistic",
    m_min = m_min,
    s_ref = s_ref,
    record_trace = record_trace,
    record_l2 = record_l2,
    record_initial = record_initial,
    record_every = record_every,
    seed = seed
  )
}

run_fs_right_logistic <- function(
    X,
    y,
    theta_star,
    s,
    eta,
    T_fs,
    m_min = 10,
    record_trace = FALSE,
    record_l2 = !is.null(theta_star),
    record_initial = FALSE,
    record_every = 1L,
    seed = NULL
) {
  run_fs_right(
    X = X,
    y = y,
    theta_star = theta_star,
    s = s,
    eta = eta,
    T_fs = T_fs,
    model = "logistic",
    m_min = m_min,
    record_trace = record_trace,
    record_l2 = record_l2,
    record_initial = record_initial,
    record_every = record_every,
    seed = seed
  )
}
