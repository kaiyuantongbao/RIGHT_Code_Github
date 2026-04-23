# Parallel diagnostic helpers for RIGHT trajectory-localization and gradient-complexity experiments.
#
# This file intentionally reuses RIGHT_Code_Github/R/data_generator.R and R/solvers.R.
# It does not redefine hard_threshold(), robust_grad_MoM(), solver_right(),
# generate_linear_data(), generate_logistic_data(), generate_X_t(), generate_epsilon_t(),
# grad_linear_regression_samplewise(), or grad_logistic_regression_samplewise().
#
# New code here is limited to diagnostic wrappers, parallel execution helpers,
# population/reference gradients, summaries, and plotting.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

find_right_root <- function(start = getwd()) {
  cur <- normalizePath(start, mustWork = TRUE)
  for (i in seq_len(30)) {
    if (file.exists(file.path(cur, "R", "solvers.R")) &&
        file.exists(file.path(cur, "R", "data_generator.R"))) {
      return(cur)
    }
    parent <- dirname(cur)
    if (identical(parent, cur)) break
    cur <- parent
  }
  stop(
    "Cannot locate RIGHT_Code_Github root. Run from inside the repo, ",
    "or place these scripts under simulations/07_schedule_ablation/."
  )
}

source_right_core <- function(root = find_right_root()) {
  source(file.path(root, "R", "data_generator.R"))
  source(file.path(root, "R", "solvers.R"))
  invisible(root)
}

RIGHT_ROOT <- source_right_core()

l2_norm <- function(x) sqrt(sum(as.numeric(x)^2))

safe_plogis <- function(z) stats::plogis(pmax(pmin(z, 35), -35))

make_sparse_theta <- function(p, s_star, seed = 1L, lower = 1.0, upper = 2.0,
                              magnitude = NULL, random_signs = FALSE) {
  stopifnot(p >= s_star, s_star >= 1)
  old_seed <- if (exists(".Random.seed", envir = .GlobalEnv)) {
    get(".Random.seed", envir = .GlobalEnv)
  } else {
    NULL
  }
  set.seed(seed)
  theta <- rep(0, p)
  vals <- if (is.null(magnitude)) stats::runif(s_star, lower, upper) else rep(magnitude, s_star)
  if (random_signs) vals <- vals * sample(c(-1, 1), s_star, replace = TRUE)
  theta[seq_len(s_star)] <- vals
  if (!is.null(old_seed)) assign(".Random.seed", old_seed, envir = .GlobalEnv)
  theta
}

write_config_csv <- function(cfg, path) {
  rows <- lapply(names(cfg), function(k) {
    v <- cfg[[k]]
    data.frame(parameter = k, value = paste(v, collapse = ","), stringsAsFactors = FALSE)
  })
  utils::write.csv(dplyr::bind_rows(rows), path, row.names = FALSE)
}

get_samplewise_gradient <- function(model) {
  model <- match.arg(model, c("linear", "logistic"))
  if (model == "linear") return(grad_linear_regression_samplewise)
  grad_logistic_regression_samplewise
}

generate_diag_data <- function(seed, model, n, p, theta_star,
                               df_X = 2.5, scale_X = 1.0, Sigma_X = NULL,
                               df_eps = 1.5, scale_eps = 1.0) {
  model <- match.arg(model, c("linear", "logistic"))
  set.seed(seed)
  if (is.null(Sigma_X)) Sigma_X <- diag(p)

  if (model == "linear") {
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
    return(list(X = dat$X, y = as.numeric(dat$y)))
  }

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

make_record_iters <- function(T_max, record_every = 5L, dense_until = 50L) {
  dense_until <- min(T_max, max(0L, as.integer(dense_until)))
  dense <- seq.int(0L, dense_until, by = 1L)
  coarse <- seq.int(0L, T_max, by = max(1L, as.integer(record_every)))
  sort(unique(c(dense, coarse, T_max)))
}

compute_diag_K <- function(n, p, s_star, m_min = 10L,
                           K1_multiplier = 1.0, K2_multiplier = 1.0,
                           K1_formula = c("s_log_p", "s_log_ep_over_s")) {
  K1_formula <- match.arg(K1_formula)
  K1_base <- if (K1_formula == "s_log_p") {
    s_star * log(p)
  } else {
    s_star * log(exp(1) * p / s_star)
  }
  K1_tar <- ceiling(K1_multiplier * K1_base)
  K2_tar <- ceiling(K2_multiplier * log(p))
  max_K <- floor(n / m_min)
  list(
    K1_tar = K1_tar,
    K1 = min(K1_tar, max_K),
    K2_tar = K2_tar,
    K2 = min(K2_tar, max_K),
    max_K = max_K
  )
}

make_K_schemes <- function(n, p, s_star, m_min = 10L,
                           K1_multiplier = 1.0,
                           K2_multipliers = c(1, 2, 4, 8),
                           K1_formula = c("s_log_p", "s_log_ep_over_s")) {
  K1_formula <- match.arg(K1_formula)
  max_K <- floor(n / m_min)
  if (max_K < 1L) stop("floor(n / m_min) < 1; increase n or decrease m_min.")

  K1_base <- if (K1_formula == "s_log_p") {
    s_star * log(p)
  } else {
    s_star * log(exp(1) * p / s_star)
  }
  K1_tar <- ceiling(K1_multiplier * K1_base)
  out <- data.frame(
    K_name = "K1_global",
    regime = "global",
    multiplier = K1_multiplier,
    K_tar = K1_tar,
    K = min(K1_tar, max_K),
    block_size = floor(n / min(K1_tar, max_K)),
    stringsAsFactors = FALSE
  )
  for (cK in K2_multipliers) {
    K2_tar <- ceiling(cK * log(p))
    K2 <- min(K2_tar, max_K)
    out <- rbind(out, data.frame(
      K_name = paste0("K2_c", gsub("\\.", "p", as.character(cK))),
      regime = "local",
      multiplier = cK,
      K_tar = K2_tar,
      K = K2,
      block_size = floor(n / K2),
      stringsAsFactors = FALSE
    ))
  }
  out <- out[!duplicated(out$K_name), ]
  rownames(out) <- NULL
  out
}

parallel_lapply <- function(X, FUN, n_cores = 1L, label = "tasks") {
  n_cores <- as.integer(n_cores)
  if (!is.finite(n_cores) || n_cores < 1L) n_cores <- 1L
  if (n_cores == 1L || length(X) <= 1L) {
    return(lapply(X, FUN))
  }
  if (.Platform$OS.type == "windows") {
    message("Parallel execution via fork is unavailable on Windows; running sequentially.")
    return(lapply(X, FUN))
  }
  message("Running ", length(X), " ", label, " with ", n_cores, " forked workers.")
  parallel::mclapply(X, FUN, mc.cores = n_cores, mc.preschedule = FALSE)
}

support_overlap_metrics <- function(theta, theta_star, tol = 1e-10) {
  supp_hat <- which(abs(theta) > tol)
  supp_star <- which(abs(theta_star) > tol)
  tp <- length(intersect(supp_hat, supp_star))
  fp <- length(setdiff(supp_hat, supp_star))
  fn <- length(setdiff(supp_star, supp_hat))
  precision <- if ((tp + fp) == 0) NA_real_ else tp / (tp + fp)
  recall <- if ((tp + fn) == 0) NA_real_ else tp / (tp + fn)
  f1 <- if (is.na(precision) || is.na(recall) || precision + recall == 0) {
    NA_real_
  } else {
    2 * precision * recall / (precision + recall)
  }
  c(tp = tp, fp = fp, fn = fn, precision = precision, recall = recall, f1 = f1)
}

jaccard_distance_bool <- function(a, b) {
  u <- sum(a | b)
  if (u == 0) return(0)
  1 - sum(a & b) / u
}

effective_support_width <- function(freq_row) {
  total <- sum(freq_row)
  if (!is.finite(total) || total <= 0) return(NA_real_)
  q <- freq_row / total
  q <- q[q > 0]
  exp(-sum(q * log(q)))
}

fd_right_trajectory <- function(X, y, theta_star, s, eta, T_max, K,
                                model = c("linear", "logistic"),
                                theta_init = NULL, record_iters = NULL,
                                record_every = 5L, record_dense_until = 50L,
                                radius = Inf, tol = 1e-10) {
  model <- match.arg(model)
  p <- ncol(X)
  grad_func <- get_samplewise_gradient(model)
  theta <- if (is.null(theta_init)) rep(0, p) else as.numeric(theta_init)

  if (is.null(record_iters)) {
    record_iters <- make_record_iters(T_max, record_every, record_dense_until)
  } else {
    record_iters <- sort(unique(as.integer(record_iters)))
    record_iters <- record_iters[record_iters >= 0L & record_iters <= T_max]
  }
  n_record <- length(record_iters)

  theta_mat <- matrix(NA_real_, nrow = n_record, ncol = p)
  support_mat <- matrix(FALSE, nrow = n_record, ncol = p)
  metrics <- data.frame(
    iter = record_iters,
    l2_error = NA_real_,
    tp = NA_real_, fp = NA_real_, fn = NA_real_,
    precision = NA_real_, recall = NA_real_, f1 = NA_real_
  )

  store_state <- function(row_id, theta_val) {
    theta_mat[row_id, ] <<- theta_val
    support_mat[row_id, ] <<- abs(theta_val) > tol
    sm <- support_overlap_metrics(theta_val, theta_star, tol = tol)
    metrics$l2_error[row_id] <<- l2_norm(theta_val - theta_star)
    metrics$tp[row_id] <<- sm["tp"]
    metrics$fp[row_id] <<- sm["fp"]
    metrics$fn[row_id] <<- sm["fn"]
    metrics$precision[row_id] <<- sm["precision"]
    metrics$recall[row_id] <<- sm["recall"]
    metrics$f1[row_id] <<- sm["f1"]
  }

  row_id <- 1L
  store_state(row_id, theta)
  record_set <- rep(FALSE, T_max + 1L)
  record_set[record_iters + 1L] <- TRUE

  for (tt in seq_len(T_max)) {
    theta <- solver_right(
      X = X,
      y = y,
      s = s,
      eta = eta,
      T_max = 1L,
      K = K,
      theta_init = theta,
      grad_func_samplewise = grad_func
    )
    if (is.finite(radius)) {
      nr <- l2_norm(theta)
      if (nr > radius) theta <- theta * (radius / nr)
    }
    if (record_set[tt + 1L]) {
      row_id <- row_id + 1L
      store_state(row_id, theta)
    }
  }

  list(
    record_iters = record_iters,
    theta_mat = theta_mat,
    support_mat = support_mat,
    metrics = metrics,
    theta_final = theta
  )
}

add_support_path_metrics <- function(metrics, support_mat, theta_star, window_records = 5L) {
  n_record <- nrow(support_mat)
  true_support <- which(abs(theta_star) > 0)
  false_support <- setdiff(seq_len(ncol(support_mat)), true_support)
  window_records <- max(1L, as.integer(window_records))

  metrics$jaccard_prev <- NA_real_
  metrics$hamming_prev <- NA_real_
  metrics$window_union_size <- NA_real_
  metrics$window_union_false_size <- NA_real_
  metrics$false_effective_width_current <- NA_real_

  for (i in seq_len(n_record)) {
    if (i > 1L) {
      metrics$jaccard_prev[i] <- jaccard_distance_bool(support_mat[i, ], support_mat[i - 1L, ])
      metrics$hamming_prev[i] <- sum(xor(support_mat[i, ], support_mat[i - 1L, ]))
    }
    lo <- max(1L, i - window_records + 1L)
    win_union <- apply(support_mat[lo:i, , drop = FALSE], 2, any)
    metrics$window_union_size[i] <- sum(win_union)
    metrics$window_union_false_size[i] <- sum(win_union[false_support])
    freq_current <- as.numeric(support_mat[i, false_support])
    metrics$false_effective_width_current[i] <- effective_support_width(freq_current)
  }
  metrics
}

summarize_support_frequency <- function(support_array, record_iters, theta_star) {
  # support_array: reps x n_record x p
  freq <- apply(support_array, c(2, 3), mean)
  p <- ncol(freq)
  true_support <- which(abs(theta_star) > 0)
  noise_support <- setdiff(seq_len(p), true_support)

  max_noise_freq <- if (length(noise_support) > 0) {
    apply(freq[, noise_support, drop = FALSE], 2, max)
  } else {
    numeric(0)
  }
  ordered_noise <- if (length(noise_support) > 0) {
    noise_support[order(max_noise_freq, decreasing = TRUE)]
  } else {
    integer(0)
  }
  coord_order <- c(true_support, ordered_noise)
  freq_ordered <- freq[, coord_order, drop = FALSE]

  heatmap_df <- as.data.frame(freq_ordered)
  colnames(heatmap_df) <- paste0("coord_", seq_along(coord_order))
  heatmap_df$iter <- record_iters
  heatmap_long <- tidyr::pivot_longer(
    heatmap_df,
    cols = starts_with("coord_"),
    names_to = "coord_rank",
    values_to = "selection_frequency"
  )
  heatmap_long$coord_rank <- as.integer(sub("coord_", "", heatmap_long$coord_rank))
  heatmap_long$original_coord <- coord_order[heatmap_long$coord_rank]
  heatmap_long$is_true_support <- heatmap_long$original_coord %in% true_support

  width_all <- apply(freq, 1, effective_support_width)
  width_noise <- if (length(noise_support) > 0) {
    apply(freq[, noise_support, drop = FALSE], 1, effective_support_width)
  } else {
    rep(NA_real_, nrow(freq))
  }
  list(
    freq = freq,
    heatmap_long = heatmap_long,
    coord_order = coord_order,
    width_all = width_all,
    width_noise = width_noise
  )
}

run_support_localization <- function(cfg) {
  model <- match.arg(cfg$model, c("linear", "logistic"))
  out_dir <- file.path(RIGHT_ROOT, "results", "07_schedule_ablation",
                       paste0("E_support_localization_FD_", model))
  dir.create(file.path(out_dir, "raw"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(out_dir, "summary"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(out_dir, "figures"), recursive = TRUE, showWarnings = FALSE)

  theta_star <- make_sparse_theta(
    p = cfg$p,
    s_star = cfg$s_star,
    seed = cfg$theta_seed,
    lower = cfg$theta_lower,
    upper = cfg$theta_upper,
    magnitude = cfg$theta_magnitude,
    random_signs = cfg$random_signs
  )

  K_info <- compute_diag_K(cfg$n, cfg$p, cfg$s_star, cfg$m_min,
                           K1_multiplier = cfg$K1_multiplier,
                           K1_formula = cfg$K1_formula)
  K <- K_info$K1
  if (K < 1L) stop("K is smaller than 1. Increase n or decrease m_min.")
  record_iters <- make_record_iters(cfg$T_max, cfg$record_every, cfg$record_dense_until)
  n_record <- length(record_iters)

  worker <- function(rep_id) {
    dat <- generate_diag_data(
      seed = cfg$seed_base + rep_id,
      model = model,
      n = cfg$n,
      p = cfg$p,
      theta_star = theta_star,
      df_X = cfg$df_X,
      scale_X = cfg$scale_X,
      Sigma_X = diag(cfg$p),
      df_eps = cfg$df_eps,
      scale_eps = cfg$scale_eps
    )
    traj <- fd_right_trajectory(
      X = dat$X,
      y = dat$y,
      theta_star = theta_star,
      s = cfg$s_alg,
      eta = cfg$eta,
      T_max = cfg$T_max,
      K = K,
      model = model,
      record_iters = record_iters,
      radius = cfg$radius
    )
    metrics <- add_support_path_metrics(traj$metrics, traj$support_mat, theta_star,
                                        window_records = cfg$window_records)
    list(
      rep_id = rep_id,
      support_mat = traj$support_mat,
      metrics = cbind(data.frame(model = model, rep_id = rep_id), metrics)
    )
  }

  rep_results <- parallel_lapply(seq_len(cfg$reps), worker, cfg$n_cores,
                                 label = paste0(model, " support-localization reps"))
  support_array <- array(FALSE, dim = c(cfg$reps, n_record, cfg$p))
  metrics_all <- vector("list", cfg$reps)
  for (i in seq_along(rep_results)) {
    rep_id <- rep_results[[i]]$rep_id
    support_array[rep_id, , ] <- rep_results[[i]]$support_mat
    metrics_all[[rep_id]] <- rep_results[[i]]$metrics
  }

  metrics_raw <- dplyr::bind_rows(metrics_all)
  support_summary <- summarize_support_frequency(support_array, record_iters, theta_star)

  curve_summary <- metrics_raw %>%
    group_by(model, iter) %>%
    summarise(
      median_l2 = median(l2_error, na.rm = TRUE),
      q25_l2 = quantile(l2_error, 0.25, na.rm = TRUE),
      q75_l2 = quantile(l2_error, 0.75, na.rm = TRUE),
      mean_precision = mean(precision, na.rm = TRUE),
      mean_recall = mean(recall, na.rm = TRUE),
      mean_f1 = mean(f1, na.rm = TRUE),
      mean_jaccard_prev = mean(jaccard_prev, na.rm = TRUE),
      mean_hamming_prev = mean(hamming_prev, na.rm = TRUE),
      median_window_union_size = median(window_union_size, na.rm = TRUE),
      median_window_union_false_size = median(window_union_false_size, na.rm = TRUE),
      .groups = "drop"
    )
  curve_summary$effective_support_width_all <- support_summary$width_all
  curve_summary$effective_support_width_noise <- support_summary$width_noise

  saveRDS(support_array, file.path(out_dir, "raw", "support_array.rds"))
  utils::write.csv(metrics_raw, file.path(out_dir, "raw", "trajectory_metrics_raw.csv"), row.names = FALSE)
  utils::write.csv(support_summary$heatmap_long, file.path(out_dir, "summary", "support_frequency_heatmap_data.csv"), row.names = FALSE)
  utils::write.csv(curve_summary, file.path(out_dir, "summary", "support_localization_curves.csv"), row.names = FALSE)
  utils::write.csv(data.frame(coord_rank = seq_along(support_summary$coord_order),
                              original_coord = support_summary$coord_order,
                              is_true_support = support_summary$coord_order %in% which(theta_star != 0)),
                   file.path(out_dir, "summary", "coordinate_order.csv"), row.names = FALSE)
  utils::write.csv(data.frame(K1_tar = K_info$K1_tar, K1 = K_info$K1,
                              K2_tar = K_info$K2_tar, K2 = K_info$K2,
                              max_K = K_info$max_K),
                   file.path(out_dir, "summary", "K_values.csv"), row.names = FALSE)
  write_config_csv(cfg, file.path(out_dir, "summary", "config.csv"))

  title_prefix <- if (model == "linear") "Linear" else "Logistic"
  p_heat <- ggplot(support_summary$heatmap_long,
                   aes(x = coord_rank, y = iter, fill = selection_frequency)) +
    geom_tile() +
    geom_vline(xintercept = cfg$s_star + 0.5, linetype = "dashed", linewidth = 0.3) +
    scale_fill_gradient(low = "white", high = "red", limits = c(0, 1)) +
    scale_y_reverse() +
    labs(
      x = "Coordinates: true support first; nuisance coordinates sorted by max frequency",
      y = "Iteration",
      fill = "Selection\nfrequency",
      title = paste0(title_prefix, " FD-RIGHT empirical active-set diagnostic")
    ) +
    theme_bw(base_size = 11)

  ggplot2::ggsave(file.path(out_dir, "figures", paste0(model, "_FD_support_frequency_heatmap.pdf")),
                  p_heat, width = 9, height = 6)
  ggplot2::ggsave(file.path(out_dir, "figures", paste0(model, "_FD_support_frequency_heatmap.png")),
                  p_heat, width = 9, height = 6, dpi = 200)

  noise_long <- support_summary$heatmap_long %>% filter(!is_true_support)
  if (nrow(noise_long) > 0) {
    cap <- if (!is.null(cfg$noise_freq_cap) && is.finite(cfg$noise_freq_cap)) {
      cfg$noise_freq_cap
    } else {
      max(stats::quantile(noise_long$selection_frequency, 0.995, na.rm = TRUE), 0.02)
    }
    noise_long$selection_frequency_capped <- pmin(noise_long$selection_frequency, cap)
    p_noise <- ggplot(noise_long,
                      aes(x = coord_rank - cfg$s_star, y = iter, fill = selection_frequency_capped)) +
      geom_tile() +
      scale_fill_gradient(low = "white", high = "red", limits = c(0, cap)) +
      scale_y_reverse() +
      labs(
        x = "Nuisance coordinates sorted by max frequency",
        y = "Iteration",
        fill = paste0("Frequency\n(capped at ", signif(cap, 3), ")"),
        title = paste0(title_prefix, " nuisance-coordinate selection frequency")
      ) +
      theme_bw(base_size = 11)
    ggplot2::ggsave(file.path(out_dir, "figures", paste0(model, "_FD_nuisance_frequency_heatmap.pdf")),
                    p_noise, width = 9, height = 6)
    ggplot2::ggsave(file.path(out_dir, "figures", paste0(model, "_FD_nuisance_frequency_heatmap.png")),
                    p_noise, width = 9, height = 6, dpi = 200)
  }

  curve_long <- curve_summary %>%
    select(iter, mean_recall, mean_precision, mean_f1,
           mean_jaccard_prev, median_window_union_size,
           median_window_union_false_size,
           effective_support_width_all, effective_support_width_noise) %>%
    tidyr::pivot_longer(cols = -iter, names_to = "metric", values_to = "value")

  p_curve <- ggplot(curve_long, aes(x = iter, y = value)) +
    geom_line(linewidth = 0.8) +
    facet_wrap(~ metric, scales = "free_y", ncol = 2) +
    labs(x = "Iteration", y = "Value", title = paste0(title_prefix, " FD-RIGHT active-set summaries")) +
    theme_bw(base_size = 11)

  ggplot2::ggsave(file.path(out_dir, "figures", paste0(model, "_FD_support_summary_curves.pdf")),
                  p_curve, width = 9, height = 8)
  ggplot2::ggsave(file.path(out_dir, "figures", paste0(model, "_FD_support_summary_curves.png")),
                  p_curve, width = 9, height = 8, dpi = 200)

  p_l2 <- ggplot(curve_summary, aes(x = iter, y = median_l2)) +
    geom_ribbon(aes(ymin = q25_l2, ymax = q75_l2), alpha = 0.20) +
    geom_line(linewidth = 0.8) +
    labs(x = "Iteration", y = "L2 error", title = paste0(title_prefix, " FD-RIGHT error trajectory")) +
    theme_bw(base_size = 11)

  ggplot2::ggsave(file.path(out_dir, "figures", paste0(model, "_FD_l2_error_trajectory.pdf")),
                  p_l2, width = 7, height = 5)
  ggplot2::ggsave(file.path(out_dir, "figures", paste0(model, "_FD_l2_error_trajectory.png")),
                  p_l2, width = 7, height = 5, dpi = 200)

  message("Done. Results written to: ", out_dir)
  invisible(list(raw = metrics_raw, summary = curve_summary, out_dir = out_dir))
}

random_sparse_direction <- function(p, q_support, seed, avoid = integer(0), random_signs = TRUE) {
  stopifnot(q_support >= 1L, q_support <= p)
  old_seed <- if (exists(".Random.seed", envir = .GlobalEnv)) {
    get(".Random.seed", envir = .GlobalEnv)
  } else {
    NULL
  }
  set.seed(seed)
  pool <- setdiff(seq_len(p), avoid)
  if (length(pool) < q_support) stop("Not enough available coordinates for sparse direction.")
  supp <- sample(pool, q_support, replace = FALSE)
  vals <- stats::rnorm(q_support)
  if (!random_signs) vals <- abs(vals)
  vals <- vals / sqrt(sum(vals^2))
  v <- rep(0, p)
  v[supp] <- vals
  if (!is.null(old_seed)) assign(".Random.seed", old_seed, envir = .GlobalEnv)
  v
}

topq_l2_norm <- function(v, q) {
  q <- min(length(v), max(1L, as.integer(q)))
  idx <- order(abs(v), decreasing = TRUE)[seq_len(q)]
  sqrt(sum(v[idx]^2))
}

linear_population_gradient_t <- function(theta, theta_star, df_X = 2.5,
                                         scale_X = 1.0, Sigma_X = NULL) {
  p <- length(theta)
  if (is.null(Sigma_X)) Sigma_X <- diag(p)
  if (df_X <= 2) stop("Closed-form covariance of t design requires df_X > 2.")
  cov_factor <- scale_X^2 * df_X / (df_X - 2)
  as.numeric(cov_factor * Sigma_X %*% (theta - theta_star))
}

logistic_population_gradient_mc <- function(theta, theta_star, p,
                                            df_X = 2.5, scale_X = 1.0,
                                            Sigma_X = NULL, M_ref = 50000L,
                                            chunk_size = 2000L, seed = 12345L) {
  if (is.null(Sigma_X)) Sigma_X <- diag(p)
  set.seed(seed)
  out <- rep(0, p)
  done <- 0L
  while (done < M_ref) {
    m <- min(chunk_size, M_ref - done)
    X <- generate_X_t(n = m, p = p, df_X = df_X, Sigma_X = Sigma_X, scale_X = scale_X)
    diff_prob <- safe_plogis(as.numeric(X %*% theta)) -
      safe_plogis(as.numeric(X %*% theta_star))
    out <- out + colSums(X * diff_prob)
    done <- done + m
  }
  out / M_ref
}

population_gradient_diag <- function(theta, theta_star, model,
                                     df_X = 2.5, scale_X = 1.0, Sigma_X = NULL,
                                     M_ref = 50000L, chunk_size = 2000L,
                                     seed = 12345L) {
  model <- match.arg(model, c("linear", "logistic"))
  if (model == "linear") {
    return(linear_population_gradient_t(theta, theta_star, df_X = df_X,
                                        scale_X = scale_X, Sigma_X = Sigma_X))
  }
  logistic_population_gradient_mc(theta, theta_star, p = length(theta),
                                  df_X = df_X, scale_X = scale_X, Sigma_X = Sigma_X,
                                  M_ref = M_ref, chunk_size = chunk_size, seed = seed)
}

make_candidate_table <- function(cfg, theta_star, n_candidates, r_grid) {
  avoid <- if (isTRUE(cfg$avoid_true_support_in_directions)) which(theta_star != 0) else integer(0)
  rows <- list()
  idx <- 0L
  for (radius_id in seq_along(r_grid)) {
    r <- r_grid[radius_id]
    for (cand_id in seq_len(n_candidates)) {
      v <- random_sparse_direction(
        p = cfg$p,
        q_support = cfg$direction_sparsity,
        seed = cfg$seed_base + 100000L * radius_id + cand_id,
        avoid = avoid,
        random_signs = TRUE
      )
      idx <- idx + 1L
      rows[[idx]] <- list(radius_id = radius_id, radius = r, cand_id = cand_id, v = v,
                          theta = theta_star + r * v)
    }
  }
  rows
}

precompute_reference_gradients <- function(candidates, cfg, theta_star, model) {
  out <- vector("list", length(candidates))
  for (i in seq_along(candidates)) {
    if (i %% 20L == 0L) message("Reference gradient ", i, " / ", length(candidates))
    out[[i]] <- population_gradient_diag(
      theta = candidates[[i]]$theta,
      theta_star = theta_star,
      model = model,
      df_X = cfg$df_X,
      scale_X = cfg$scale_X,
      Sigma_X = diag(cfg$p),
      M_ref = cfg$M_ref,
      chunk_size = cfg$ref_chunk_size,
      seed = cfg$seed_base + 200000L + i
    )
  }
  out
}

run_gradient_complexity_pointwise <- function(cfg) {
  model <- match.arg(cfg$model, c("linear", "logistic"))
  out_dir <- file.path(RIGHT_ROOT, "results", "07_schedule_ablation",
                       paste0("E_gradient_complexity_pointwise_", model))
  dir.create(file.path(out_dir, "raw"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(out_dir, "summary"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(out_dir, "figures"), recursive = TRUE, showWarnings = FALSE)

  theta_star <- make_sparse_theta(
    p = cfg$p,
    s_star = cfg$s_star,
    seed = cfg$theta_seed,
    lower = cfg$theta_lower,
    upper = cfg$theta_upper,
    magnitude = cfg$theta_magnitude,
    random_signs = cfg$random_signs
  )

  K_schemes <- make_K_schemes(cfg$n_grad, cfg$p, cfg$s_star, cfg$m_min,
                              K1_multiplier = cfg$K1_multiplier,
                              K2_multipliers = cfg$K2_multipliers,
                              K1_formula = cfg$K1_formula)
  q_top <- min(cfg$p, 2 * cfg$s_alg + cfg$s_star)

  candidates <- make_candidate_table(cfg, theta_star, cfg$n_directions, cfg$r_grid)
  ref_gradients <- precompute_reference_gradients(candidates, cfg, theta_star, model)
  grad_func <- get_samplewise_gradient(model)

  worker <- function(rep_id) {
    dat <- generate_diag_data(
      seed = cfg$seed_base + rep_id,
      model = model,
      n = cfg$n_grad,
      p = cfg$p,
      theta_star = theta_star,
      df_X = cfg$df_X,
      scale_X = cfg$scale_X,
      Sigma_X = diag(cfg$p),
      df_eps = cfg$df_eps,
      scale_eps = cfg$scale_eps
    )
    rows <- vector("list", length(candidates) * nrow(K_schemes))
    ctr <- 0L
    for (ci in seq_along(candidates)) {
      cand <- candidates[[ci]]
      g_ref <- ref_gradients[[ci]]
      for (kk in seq_len(nrow(K_schemes))) {
        K <- K_schemes$K[kk]
        g_hat <- robust_grad_MoM(
          theta = cand$theta,
          X = dat$X,
          y = dat$y,
          K = K,
          grad_func_samplewise = grad_func
        )
        err <- as.numeric(g_hat - g_ref)
        ctr <- ctr + 1L
        rows[[ctr]] <- data.frame(
          model = model,
          rep_id = rep_id,
          radius_id = cand$radius_id,
          radius = cand$radius,
          cand_id = cand$cand_id,
          K_name = K_schemes$K_name[kk],
          regime = K_schemes$regime[kk],
          multiplier = K_schemes$multiplier[kk],
          K = K,
          block_size = K_schemes$block_size[kk],
          q_top = q_top,
          topq_error = topq_l2_norm(err, q_top),
          l2_error = l2_norm(err),
          linf_error = max(abs(err)),
          ref_grad_norm = l2_norm(g_ref)
        )
      }
    }
    dplyr::bind_rows(rows[seq_len(ctr)])
  }

  raw_list <- parallel_lapply(seq_len(cfg$reps), worker, cfg$n_cores,
                              label = paste0(model, " pointwise gradient reps"))
  raw <- dplyr::bind_rows(raw_list)

  summary_tbl <- raw %>%
    group_by(model, radius, K_name, regime, multiplier, K, block_size, q_top) %>%
    summarise(
      median_topq_error = median(topq_error, na.rm = TRUE),
      q25_topq_error = quantile(topq_error, 0.25, na.rm = TRUE),
      q75_topq_error = quantile(topq_error, 0.75, na.rm = TRUE),
      mean_topq_error = mean(topq_error, na.rm = TRUE),
      median_l2_error = median(l2_error, na.rm = TRUE),
      median_linf_error = median(linf_error, na.rm = TRUE),
      median_ref_grad_norm = median(ref_grad_norm, na.rm = TRUE),
      .groups = "drop"
    )

  utils::write.csv(raw, file.path(out_dir, "raw", paste0(model, "_pointwise_gradient_raw.csv")), row.names = FALSE)
  utils::write.csv(summary_tbl, file.path(out_dir, "summary", paste0(model, "_pointwise_gradient_summary.csv")), row.names = FALSE)
  utils::write.csv(K_schemes, file.path(out_dir, "summary", "K_schemes.csv"), row.names = FALSE)
  write_config_csv(cfg, file.path(out_dir, "summary", "config.csv"))

  title_prefix <- if (model == "linear") "Linear" else "Logistic"
  p1 <- ggplot(summary_tbl, aes(x = radius, y = median_topq_error,
                                group = K_name, color = K_name)) +
    geom_ribbon(aes(ymin = q25_topq_error, ymax = q75_topq_error, fill = K_name),
                alpha = 0.12, color = NA) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2) +
    scale_x_log10(breaks = cfg$r_grid) +
    scale_y_log10() +
    labs(
      x = "Localization radius r in theta(r) = theta* + r v",
      y = paste0("Top-", q_top, " pointwise gradient error"),
      color = "K scheme",
      fill = "K scheme",
      title = paste0(title_prefix, " pointwise gradient diagnostic with K2 constants")
    ) +
    theme_bw(base_size = 11)

  ggplot2::ggsave(file.path(out_dir, "figures", paste0(model, "_pointwise_gradient_topq_error.pdf")),
                  p1, width = 8.5, height = 5)
  ggplot2::ggsave(file.path(out_dir, "figures", paste0(model, "_pointwise_gradient_topq_error.png")),
                  p1, width = 8.5, height = 5, dpi = 200)

  message("Done. Results written to: ", out_dir)
  invisible(list(raw = raw, summary = summary_tbl, out_dir = out_dir))
}

run_gradient_complexity_uniform <- function(cfg) {
  model <- match.arg(cfg$model, c("linear", "logistic"))
  out_dir <- file.path(RIGHT_ROOT, "results", "07_schedule_ablation",
                       paste0("E_gradient_complexity_uniform_", model))
  dir.create(file.path(out_dir, "raw"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(out_dir, "summary"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(out_dir, "figures"), recursive = TRUE, showWarnings = FALSE)

  theta_star <- make_sparse_theta(
    p = cfg$p,
    s_star = cfg$s_star,
    seed = cfg$theta_seed,
    lower = cfg$theta_lower,
    upper = cfg$theta_upper,
    magnitude = cfg$theta_magnitude,
    random_signs = cfg$random_signs
  )

  K_schemes <- make_K_schemes(cfg$n_grad, cfg$p, cfg$s_star, cfg$m_min,
                              K1_multiplier = cfg$K1_multiplier,
                              K2_multipliers = cfg$K2_multipliers,
                              K1_formula = cfg$K1_formula)
  q_top <- min(cfg$p, 2 * cfg$s_alg + cfg$s_star)
  max_candidates <- max(cfg$candidate_grid)
  candidates <- make_candidate_table(cfg, theta_star, max_candidates, cfg$r_grid)
  ref_gradients <- precompute_reference_gradients(candidates, cfg, theta_star, model)
  grad_func <- get_samplewise_gradient(model)

  worker <- function(rep_id) {
    dat <- generate_diag_data(
      seed = cfg$seed_base + rep_id,
      model = model,
      n = cfg$n_grad,
      p = cfg$p,
      theta_star = theta_star,
      df_X = cfg$df_X,
      scale_X = cfg$scale_X,
      Sigma_X = diag(cfg$p),
      df_eps = cfg$df_eps,
      scale_eps = cfg$scale_eps
    )

    err_rows <- vector("list", length(candidates) * nrow(K_schemes))
    ctr <- 0L
    for (ci in seq_along(candidates)) {
      cand <- candidates[[ci]]
      g_ref <- ref_gradients[[ci]]
      for (kk in seq_len(nrow(K_schemes))) {
        K <- K_schemes$K[kk]
        g_hat <- robust_grad_MoM(
          theta = cand$theta,
          X = dat$X,
          y = dat$y,
          K = K,
          grad_func_samplewise = grad_func
        )
        err <- as.numeric(g_hat - g_ref)
        ctr <- ctr + 1L
        err_rows[[ctr]] <- data.frame(
          model = model,
          rep_id = rep_id,
          radius_id = cand$radius_id,
          radius = cand$radius,
          cand_id = cand$cand_id,
          K_name = K_schemes$K_name[kk],
          regime = K_schemes$regime[kk],
          multiplier = K_schemes$multiplier[kk],
          K = K,
          block_size = K_schemes$block_size[kk],
          q_top = q_top,
          topq_error = topq_l2_norm(err, q_top),
          l2_error = l2_norm(err),
          linf_error = max(abs(err)),
          ref_grad_norm = l2_norm(g_ref)
        )
      }
    }
    err_raw <- dplyr::bind_rows(err_rows[seq_len(ctr)])

    # Prefix-max diagnostic: M = 1 is fixed-iterate; larger M approximates uniform control over more candidates.
    summ_rows <- list()
    ctr2 <- 0L
    for (rad in cfg$r_grid) {
      for (kname in unique(err_raw$K_name)) {
        sub <- err_raw %>% filter(radius == rad, K_name == kname) %>% arrange(cand_id)
        for (M in cfg$candidate_grid) {
          ss <- sub %>% filter(cand_id <= M)
          ctr2 <- ctr2 + 1L
          summ_rows[[ctr2]] <- data.frame(
            model = model,
            rep_id = rep_id,
            radius = rad,
            candidate_count = M,
            K_name = kname,
            regime = unique(ss$regime),
            multiplier = unique(ss$multiplier),
            K = unique(ss$K),
            block_size = unique(ss$block_size),
            q_top = q_top,
            max_topq_error = max(ss$topq_error, na.rm = TRUE),
            q90_topq_error = as.numeric(stats::quantile(ss$topq_error, 0.90, na.rm = TRUE)),
            median_topq_error = median(ss$topq_error, na.rm = TRUE)
          )
        }
      }
    }
    list(raw = err_raw, summary = dplyr::bind_rows(summ_rows[seq_len(ctr2)]))
  }

  rep_results <- parallel_lapply(seq_len(cfg$reps), worker, cfg$n_cores,
                                 label = paste0(model, " uniform gradient reps"))
  raw <- dplyr::bind_rows(lapply(rep_results, `[[`, "raw"))
  uniform_raw <- dplyr::bind_rows(lapply(rep_results, `[[`, "summary"))

  uniform_summary <- uniform_raw %>%
    group_by(model, radius, candidate_count, K_name, regime, multiplier, K, block_size, q_top) %>%
    summarise(
      median_max_topq_error = median(max_topq_error, na.rm = TRUE),
      q25_max_topq_error = quantile(max_topq_error, 0.25, na.rm = TRUE),
      q75_max_topq_error = quantile(max_topq_error, 0.75, na.rm = TRUE),
      median_q90_topq_error = median(q90_topq_error, na.rm = TRUE),
      median_median_topq_error = median(median_topq_error, na.rm = TRUE),
      .groups = "drop"
    )

  utils::write.csv(raw, file.path(out_dir, "raw", paste0(model, "_uniform_candidate_errors_raw.csv")), row.names = FALSE)
  utils::write.csv(uniform_raw, file.path(out_dir, "raw", paste0(model, "_uniform_prefix_raw.csv")), row.names = FALSE)
  utils::write.csv(uniform_summary, file.path(out_dir, "summary", paste0(model, "_uniform_gradient_summary.csv")), row.names = FALSE)
  utils::write.csv(K_schemes, file.path(out_dir, "summary", "K_schemes.csv"), row.names = FALSE)
  write_config_csv(cfg, file.path(out_dir, "summary", "config.csv"))

  title_prefix <- if (model == "linear") "Linear" else "Logistic"
  p1 <- ggplot(uniform_summary,
               aes(x = candidate_count, y = median_max_topq_error,
                   group = K_name, color = K_name)) +
    geom_ribbon(aes(ymin = q25_max_topq_error, ymax = q75_max_topq_error, fill = K_name),
                alpha = 0.12, color = NA) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2) +
    scale_x_log10(breaks = cfg$candidate_grid) +
    scale_y_log10() +
    facet_wrap(~ radius, scales = "free_y", labeller = label_both) +
    labs(
      x = "Number of candidate iterates M in prefix maximum",
      y = paste0("Median max top-", q_top, " error over M candidates"),
      color = "K scheme",
      fill = "K scheme",
      title = paste0(title_prefix, " fixed-vs-uniform gradient complexity diagnostic")
    ) +
    theme_bw(base_size = 11)

  ggplot2::ggsave(file.path(out_dir, "figures", paste0(model, "_uniform_gradient_max_error_vs_M.pdf")),
                  p1, width = 10, height = 6)
  ggplot2::ggsave(file.path(out_dir, "figures", paste0(model, "_uniform_gradient_max_error_vs_M.png")),
                  p1, width = 10, height = 6, dpi = 200)

  p2 <- ggplot(uniform_summary,
               aes(x = candidate_count, y = median_q90_topq_error,
                   group = K_name, color = K_name)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2) +
    scale_x_log10(breaks = cfg$candidate_grid) +
    scale_y_log10() +
    facet_wrap(~ radius, scales = "free_y", labeller = label_both) +
    labs(
      x = "Number of candidate iterates M",
      y = paste0("Median 90% top-", q_top, " error over M candidates"),
      color = "K scheme",
      title = paste0(title_prefix, " candidate-class size diagnostic")
    ) +
    theme_bw(base_size = 11)

  ggplot2::ggsave(file.path(out_dir, "figures", paste0(model, "_uniform_gradient_q90_error_vs_M.pdf")),
                  p2, width = 10, height = 6)
  ggplot2::ggsave(file.path(out_dir, "figures", paste0(model, "_uniform_gradient_q90_error_vs_M.png")),
                  p2, width = 10, height = 6, dpi = 200)

  message("Done. Results written to: ", out_dir)
  invisible(list(raw = raw, summary = uniform_summary, out_dir = out_dir))
}
