

#' Run a Single Replication for the Rate vs. Moment Experiment
#'
#' This function encapsulates the logic for a single run of the first experiment,
#' which investigates the convergence rate of the RIGHT estimator under
#' t-distributed noise with varying degrees of freedom.
#'
#' @param n integer. Sample size.
#' @param p integer. Dimension.
#' @param theta_star numeric vector. The true p-dimensional parameter vector.
#' @param delta numeric. The parameter controlling the noise's tail behavior,
#'   where the degrees of freedom of the t-distribution is 1 + delta.
#' @param s integer. The sparsity level to be used by the RIGHT solver.
#' @param eta numeric. The step-size for the RIGHT solver.
#' @param T_max integer. The maximum number of iterations for the RIGHT solver.
#' @param c_K numeric. The constant used to calculate the number of blocks K.
#' @param seed integer. A random seed for reproducibility of this single run.
#'
#' @return A single-row data.frame containing the input parameters and the
#'   resulting performance metrics (l2_error, tpr, fpr).
#'
run_rate_vs_moment_single <- function(n, p, theta_star, delta,
                                      s, eta, T_max, c_K, seed) {
  
  # Set seed for this specific replication
  set.seed(seed)
  
  # --- 1. Generate Data ---
  df_eps <- 1 + delta
  
  sim_data <- generate_linear_data(
    n = n, p = p, theta_star = theta_star,
    generator_X = generate_X_gaussian,
    generator_epsilon = generate_epsilon_t,
    df_eps = df_eps
  )
  
  # --- 2. Calculate K for the RIGHT solver ---
  # K = c * log(p * log(n / log(p)))
  # We ensure K is at least 1 and is an integer.
  K <- max(1, round(c_K * log(p * log(n / log(p)))))
  
  # --- 3. Run the RIGHT Solver ---
  theta_hat <- solver_right(
    X = sim_data$X,
    y = sim_data$y,
    s = s,
    eta = eta,
    T_max = T_max,
    K = K,
    grad_func_samplewise = grad_linear_regression_samplewise
  )
  
  # --- 4. Calculate Performance Metrics ---
  
  # L2 Estimation Error
  l2_error <- sqrt(sum((theta_hat - theta_star)^2))
  
  # Support Recovery Metrics (TPR and FPR)
  support_true <- which(theta_star != 0)
  support_hat <- which(theta_hat != 0)
  s_star <- length(support_true)
  
  true_positives <- length(intersect(support_true, support_hat))
  false_positives <- length(setdiff(support_hat, support_true))
  
  tpr <- if (s_star > 0) true_positives / s_star else 0
  fpr <- if (p > s_star) false_positives / (p - s_star) else 0
  
  # --- 5. Return Results in a Tidy Format ---
  # A single-row data.frame is easy to combine later using rbind.
  result_df <- data.frame(
    n = n,
    p = p,
    delta = delta,
    c_K = c_K,
    K = K,
    l2_error = l2_error,
    tpr = tpr,
    fpr = fpr
  )
  
  return(result_df)
}

