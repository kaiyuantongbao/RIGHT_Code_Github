#' Hard Thresholding Operator (for Vectors and Matrices)
#'
#' Projects an object onto the set of s-sparse objects.
#' - If the input is a vector, it keeps the s elements with the largest
#'   absolute values.
#' - If the input is a matrix, it performs row-wise thresholding, keeping the s
#'   rows with the largest L2 norms.
#'
#' @param v numeric vector or matrix. The input to be projected.
#' @param s integer. The desired sparsity level (number of non-zero elements/rows).
#'
#' @return A new s-sparse object with the same dimensions as v.
#'
hard_threshold <- function(v, s) {
  if (is.matrix(v)) {
    # --- Matrix Case: Row-wise Hard Thresholding ---
    p <- nrow(v)
    if (s <= 0) {
      return(matrix(0, nrow = p, ncol = ncol(v)))
    }
    if (s >= p) {
      return(v)
    }

    # Calculate the L2 norm of each row
    row_norms <- sqrt(rowSums(v^2))

    # Find the indices of the s rows with the largest norms
    indices_to_keep <- order(row_norms, decreasing = TRUE)[1:s]

    # Create a new zero matrix
    v_sparse <- matrix(0, nrow = p, ncol = ncol(v))
    # Copy the entire rows from the original matrix
    v_sparse[indices_to_keep, ] <- v[indices_to_keep, ]

    return(v_sparse)

  } else {
    # --- Vector Case: Element-wise Hard Thresholding (Original Logic) ---
    p <- length(v)
    if (s <= 0) {
      return(rep(0, p))
    }
    if (s >= p) {
      return(v)
    }

    # Find the indices of the s largest elements in absolute value
    indices_to_keep <- order(abs(v), decreasing = TRUE)[1:s]

    # Create a new zero vector
    v_sparse <- rep(0, p)
    # Copy the values from the original vector
    v_sparse[indices_to_keep] <- v[indices_to_keep]

    return(v_sparse)
  }
}


grad_linear_regression <- function(theta, X, y) {
  n <- nrow(X)
  # Calculate the residuals: X*theta - y
  residuals <- X %*% theta - y
  # Calculate the gradient using the matrix form
  gradient <- (t(X) %*% residuals) / n
  return(gradient)
}


#' Sample-wise Gradient for the Linear Regression Model
#'
#' Computes the gradient for each individual sample.
#' The gradient for sample i is grad_i = (x_i^T * theta - y_i) * x_i.
#'
#' @param theta numeric vector. The current p-dimensional parameter estimate.
#' @param X matrix. The n x p design matrix.
#' @param y numeric vector. The n-dimensional response vector.
#'
#' @return An n x p matrix where each row is the gradient for one sample.
#' each row i of sample_gradients corresponds to grad_i.
grad_linear_regression_samplewise <- function(theta, X, y) {
  # Calculate residuals for all samples at once: (X*theta - y)
  residuals <- X %*% theta - y # This is an n x 1 vector

  # Scale each row of X (which is x_i^T) by its corresponding residual.
  # R's recycling rule makes this efficient.
  sample_gradients <- X * as.vector(residuals)

  return(sample_gradients)
}




#' Robust Gradient Estimator using Median-of-Means (MoM)
#'
#' Computes a robust estimate of the gradient by partitioning the data into K
#' blocks, computing the mean gradient for each block, and then taking the
#' element-wise median of these block gradients.
#'
#' @param theta numeric vector or matrix. The current parameter estimate.
#' @param X matrix. The n x p design matrix.
#' @param y numeric vector or matrix. The response.
#' @param K integer. The number of blocks to partition the data into.
#' @param grad_func_samplewise function. A function that computes the sample-wise
#'   gradients, returning an n x p matrix (or n x p x m array for multi-response).
#'
#' @return A robust gradient estimate with the same dimensions as theta.
#'
robust_grad_MoM <- function(theta, X, y, K, grad_func_samplewise) {
  n <- nrow(X)
  sample_grads <- grad_func_samplewise(theta, X, y)
  
  is_multi_response <- is.array(sample_grads) && length(dim(sample_grads)) == 3
  
  block_size <- floor(n / K)
  
  if (is_multi_response) {
    p <- dim(sample_grads)[2]
    m <- dim(sample_grads)[3]
    block_grads <- array(NA, dim = c(K, p, m))
  } else {
    p <- dim(sample_grads)[2]
    block_grads <- matrix(NA, nrow = K, ncol = p)
  }
  
  for (k in 1:K) {
    indices <- ((k - 1) * block_size + 1):(k * block_size)
    
    if (is_multi_response) {
      # Correct 3D array slicing
      grads_in_block <- sample_grads[indices, , , drop = FALSE]
      # apply mean over the first dimension (samples in the block)
      block_grads[k, , ] <- apply(grads_in_block, c(2, 3), mean)
    } else {
      # Standard 2D matrix slicing
      grads_in_block <- sample_grads[indices, , drop = FALSE]
      block_grads[k, ] <- colMeans(grads_in_block)
    }
  }
  
  if (is_multi_response) {
    # apply median over the first dimension (blocks)
    robust_gradient <- apply(block_grads, c(2, 3), median)
  } else {
    robust_gradient <- apply(block_grads, 2, median)
  }
  
  return(robust_gradient)
}




# IHT solver
#' Iterative Hard Thresholding (IHT) Solver
#'
#' A general implementation of the IHT algorithm for sparse estimation.
#' This function is generic and can solve different models by accepting a
#' corresponding gradient function as an argument.
#'
#' @param X matrix. The n x p design matrix.
#' @param y numeric vector or matrix. The response. For multi-response, this is
#'   an n x m matrix.
#' @param s integer. The user given sparsity level.
#' @param eta numeric. The step-size (learning rate).
#' @param T_max integer. The total number of iterations.
#' @param theta_init numeric vector or matrix. The initial guess for the parameters.
#'   If NULL, a zero vector/matrix is used.
#' @param grad_func function. The function to compute the gradient. It must
#'   accept `theta`, `X`, and `y` as its first three arguments.
#'
#' @return The estimated sparse parameter vector or matrix after T_max iterations.
#'
solver_iht <- function(X, y, s, eta, T_max, theta_init = NULL, grad_func) {

  # --- Step 1: Initialization ---
  p <- ncol(X)
  is_multi_response <- is.matrix(y)
  m <- if (is_multi_response) ncol(y) else 1

  if (is.null(theta_init)) {
    # Default initialization is a zero vector or matrix
    theta_current <- if (is_multi_response) matrix(0, p, m) else rep(0, p)
  } else {
    theta_current <- theta_init
  }

  # --- Step 2: Main Iteration Loop ---
  # The loop runs exactly T_max times.
  for (t in 1:T_max) {
    # Step 2a: Calculate the gradient at the current estimate
    gradient <- grad_func(theta_current, X, y)

    # Step 2b: Perform the gradient descent step
    theta_updated <- theta_current - eta * gradient

    # Step 2c: Apply the hard thresholding operator
    theta_current <- hard_threshold(theta_updated, s)
  }

  # --- Step 3: Return the final estimate ---
  return(theta_current)
}



################################################################################################################
#' Robust Iterative Gradient Descent with Hard Thresholding (RIGHT)
#'
#' An implementation of the RIGHT algorithm which uses a Median-of-Means
#' estimator for the gradient to achieve robustness.
#'
#' @param X matrix. The n x p design matrix.
#' @param y numeric vector or matrix. The response.
#' @param s integer. The target sparsity level.
#' @param eta numeric. The step-size (learning rate).
#' @param T_max integer. The total number of iterations.
#' @param K integer. The number of blocks for the MoM estimator.
#' @param theta_init numeric vector or matrix. The initial parameter guess.
#' @param grad_func_samplewise function. The function to compute sample-wise gradients.
#'
#' @return The estimated sparse parameter vector or matrix.
#'
solver_right <- function(X, y, s, eta, T_max, K,
                         theta_init = NULL, grad_func_samplewise) {

  # --- Step 1: Initialization ---
  p <- ncol(X)
  is_multi_response <- is.matrix(y)
  m <- if (is_multi_response) ncol(y) else 1

  if (is.null(theta_init)) {
    theta_current <- if (is_multi_response) matrix(0, p, m) else rep(0, p)
  } else {
    theta_current <- theta_init
  }

  # --- Step 2: Main Iteration Loop ---
  for (t in 1:T_max) {
    # Step 2a: Calculate the ROBUST gradient at the current estimate
    robust_gradient <- robust_grad_MoM(theta_current, X, y, K, grad_func_samplewise)

    # Step 2b: Perform the gradient descent step
    theta_updated <- theta_current - eta * robust_gradient

    # Step 2c: Apply the hard thresholding operator
    theta_current <- hard_threshold(theta_updated, s)
  }

  # --- Step 3: Return the final estimate ---
  return(theta_current)
}



############################################################################################
# Huber loss method
#' Compute the Mean Huber Loss (Corrected Vector/Matrix Version)
#'
#' Calculates the average Huber loss. Handles both single-response (vector y)
#' and multi-response (matrix Y) cases robustly.
#'
#' @param y numeric vector or matrix. The response.
#' @param X matrix. The n x p design matrix.
#' @param theta numeric vector or matrix. The parameter.
#' @param tau numeric. The threshold parameter for the Huber loss.
#'
#' @return A single numeric value representing the mean Huber loss.
#'
compute_huber_loss <- function(y, X, theta, tau) {
  # Calculate residuals
  residuals <- y - X %*% theta
  abs_residuals <- abs(residuals)
  
  # Apply the Huber loss formula element-wise
  loss_quadratic <- 0.5 * residuals^2
  loss_linear <- tau * abs_residuals - 0.5 * tau^2
  huber_loss_per_element <- ifelse(abs_residuals <= tau, loss_quadratic, loss_linear)
  
  # Correctly normalize by the number of SAMPLES (n), not total elements.
  # The loss is (1/n) * sum_{i,k} l_tau(r_ik)
  n_samples <- nrow(X)
  total_loss <- sum(huber_loss_per_element)
  
  return(total_loss / n_samples)
}


# --- File: R/solvers.R (Final Corrected grad_huber_loss) ---

#' Gradient of the Mean Huber Loss (Corrected Vector/Matrix Version)
#'
#' Computes the gradient of the mean Huber loss. This version correctly
#' handles both single-response (vector) and multi-response (matrix) cases
#' by using index assignment to preserve matrix dimensions.
#'
#' @param y numeric vector or matrix. The response.
#' @param X matrix. The n x p design matrix.
#' @param theta numeric vector or matrix. The parameter.
#' @param tau numeric. The threshold parameter for the Huber loss.
#'
#' @return A vector or matrix representing the gradient, with same dims as theta.
#'
grad_huber_loss <- function(y, X, theta, tau) {
  n_samples <- nrow(X)
  
  # Calculate residuals
  residuals <- y - X %*% theta
  
  # --- Correctly apply the psi function while preserving dimensions ---
  # Start with a copy of the residuals
  psi_residuals <- residuals
  
  # Clip the values greater than tau
  psi_residuals[residuals > tau] <- tau
  
  # Clip the values less than -tau
  psi_residuals[residuals < -tau] <- -tau
  # --- End of correction ---
  
  # Calculate the gradient using the matrix form.
  gradient <- -(t(X) %*% psi_residuals) / n_samples
  
  return(gradient)
}

#' Soft-Thresholding Operator
#'
#' Applies the soft-thresholding operator element-wise to a vector.
#' This is the proximal operator for the L1 norm.
#'
#' @param v numeric vector. The input vector.
#' @param lambda numeric. The threshold value (must be non-negative).
#'
#' @return A new vector after applying soft-thresholding.
#'
soft_threshold <- function(v, lambda) {
  if (lambda < 0) {
    stop("'lambda' must be non-negative.")
  }
  
  # sign(v) handles the direction
  # pmax(abs(v) - lambda, 0) handles the shrinkage and setting to zero
  v_shrunk <- sign(v) * pmax(abs(v) - lambda, 0)
  
  return(v_shrunk)
}



#' Block Soft-Thresholding Operator for Group Lasso
#'
#' Applies row-wise soft-thresholding to a matrix. This is the proximal
#' operator for the Group Lasso penalty sum(lambda * ||row||_2).
#'
#' @param V matrix. The input matrix (e.g., p x m).
#' @param lambda numeric. The threshold value.
#'
#' @return A new matrix of the same size as V after row-wise shrinkage.
#'
block_soft_threshold <- function(V, lambda) {
  if (lambda < 0) {
    stop("'lambda' must be non-negative.")
  }
  
  # Calculate the L2 norm for each row
  row_norms <- sqrt(rowSums(V^2))
  
  # Calculate the shrinkage factor for each row.
  # The factor is max(0, 1 - lambda / ||row_norm||).
  # We add a small epsilon to avoid division by zero if a row is all zeros.
  shrinkage_factor <- pmax(0, 1 - lambda / (row_norms + 1e-12))
  
  # Apply the shrinkage factor to each row.
  # R's recycling rule multiplies each row of V by the corresponding
  # element in shrinkage_factor.
  V_shrunk <- V * shrinkage_factor
  
  return(V_shrunk)
}


# --- File: R/solvers.R (continued) ---

#' Solver for Huber-Lasso Regression using LAMM
#'
#' Implements the Local Adaptive Majorize-Minimization (LAMM) algorithm to solve
#' the L1-penalized Huber regression problem. The step-size is determined
#' adaptively using a backtracking line search.
#'
#' @param X matrix. The n x p design matrix.
#' @param y numeric vector. The n-dimensional response vector.
#' @param lambda numeric. The L1 regularization parameter.
#' @param tau numeric. The threshold parameter for the Huber loss.
#' @param T_max integer. The total number of outer iterations.
#' @param beta_init numeric vector. The initial guess for the parameters. If NULL,
#'   a zero vector is used.
#' @param phi_init numeric. The initial guess for the curvature parameter.
#' @param gamma_u numeric. The factor by which phi is increased during line search.
#'
#' @return The estimated sparse parameter vector after T_max iterations.
#'
solver_huber_lasso <- function(X, y, lambda, tau, T_max,
                               beta_init = NULL,
                               phi_init = 1e-4,
                               gamma_u = 2) {

  # --- Step 1: Initialization ---
  p <- ncol(X)
  if (is.null(beta_init)) {
    beta_current <- rep(0, p)
  } else {
    beta_current <- beta_init
  }

  phi_k <- phi_init

  # --- Step 2: Main Iteration Loop ---
  # The loop runs exactly T_max times.
  for (t in 1:T_max) {
    
    # --- Step 2a: Pre-calculate values for the current iteration ---
    loss_current <- compute_huber_loss(y, X, beta_current, tau)
    grad_current <- grad_huber_loss(y, X, beta_current, tau)

    # --- Step 2b: Backtracking Line Search to find a valid phi_k ---
    while (TRUE) {
      # Calculate the candidate solution with the current phi_k
      step_size <- 1 / phi_k
      beta_candidate <- soft_threshold(beta_current - step_size * grad_current,
                                       lambda * step_size)

      # Calculate the loss at the candidate solution
      loss_candidate <- compute_huber_loss(y, X, beta_candidate, tau)

      # Calculate the value of the surrogate function g_k at the candidate
      diff_beta <- beta_candidate - beta_current
      surrogate_val <- loss_current + sum(grad_current * diff_beta) +
                       (phi_k / 2) * sum(diff_beta^2)

      # Check the majorization condition
      if (surrogate_val >= loss_candidate) {
        # Condition is met, exit the line search loop
        break
      } else {
        # Condition not met, increase phi_k and try again
        phi_k <- phi_k * gamma_u
      }
    } # End of while loop (line search)

    # --- Step 2c: Update the parameter estimate ---
    beta_current <- beta_candidate
    
    # --- Step 2d: Prepare phi for the next iteration (optional but good practice) ---
    # This corresponds to the phi.new = max(phi.min, phi.tmp/2) logic
    # It prevents phi from growing indefinitely over iterations.
    phi_k <- max(phi_init, phi_k / gamma_u)

  } # End of for loop (main iteration)

  # --- Step 3: Return the final estimate ---
  return(beta_current)
}


# --- File: R/solvers.R (Final Huber Solver) ---

#' Solver for Huber-GroupLasso Regression using LAMM
#'
#' Implements the LAMM algorithm to solve the Group-Lasso penalized Huber
#' regression problem for multi-response models.
#'
#' @param X matrix. The n x p design matrix.
#' @param Y matrix. The n x m response matrix.
#' @param lambda numeric. The Group Lasso regularization parameter.
#' @param tau numeric. The threshold parameter for the Huber loss.
#' @param T_max integer. The total number of outer iterations.
#' @param Theta_init matrix. The initial p x m guess for the parameters. If NULL,
#'   a zero matrix is used.
#' @param phi_init numeric. The initial guess for the curvature parameter.
#' @param gamma_u numeric. The factor by which phi is increased during line search.
#'
#' @return The estimated sparse p x m parameter matrix after T_max iterations.
#'
solver_huber_group_lasso <- function(X, Y, lambda, tau, T_max,
                                     Theta_init = NULL,
                                     phi_init = 1e-4,
                                     gamma_u = 2) {

  # --- Step 1: Initialization ---
  p <- ncol(X)
  m <- ncol(Y)
  if (is.null(Theta_init)) {
    Theta_current <- matrix(0, nrow = p, ncol = m)
  } else {
    Theta_current <- Theta_init
  }

  phi_k <- phi_init

  # --- Step 2: Main Iteration Loop ---
  for (t in 1:T_max) {
    
    # --- Step 2a: Pre-calculate values for the current iteration ---
    loss_current <- compute_huber_loss(Y, X, Theta_current, tau)
    grad_current <- grad_huber_loss(Y, X, Theta_current, tau)

    # --- Step 2b: Backtracking Line Search to find a valid phi_k ---
    while (TRUE) {
      # Calculate the candidate solution with the current phi_k
      step_size <- 1 / phi_k
      Theta_tilde <- Theta_current - step_size * grad_current
      Theta_candidate <- block_soft_threshold(Theta_tilde, lambda * step_size)

      # Calculate the loss at the candidate solution
      loss_candidate <- compute_huber_loss(Y, X, Theta_candidate, tau)

      # Calculate the value of the surrogate function g_k at the candidate
      diff_Theta <- Theta_candidate - Theta_current
      # Note: sum(grad * diff) is the vectorized form of the Frobenius inner product
      surrogate_val <- loss_current + sum(grad_current * diff_Theta) +
                       (phi_k / 2) * sum(diff_Theta^2)

      # Check the majorization condition
      if (surrogate_val >= loss_candidate) {
        break # Condition met, exit the line search
      } else {
        phi_k <- phi_k * gamma_u # Condition not met, increase phi_k
      }
    } # End of while loop (line search)

    # --- Step 2c: Update the parameter estimate ---
    Theta_current <- Theta_candidate
    
    # --- Step 2d: Prepare phi for the next iteration ---
    phi_k <- max(phi_init, phi_k / gamma_u)

  } # End of for loop (main iteration)

  # --- Step 3: Return the final estimate ---
  return(Theta_current)
}



# Danzig Selector
# (Robust Dantzig Selector Components) ---

#' Element-wise Truncation Operator
#'
#' Applies element-wise truncation to a vector or matrix. For each element a,
#' the truncated value is sign(a) * min(|a|, tau).
#'
#' @param A numeric vector or matrix. The input data.
#' @param tau numeric. The non-negative truncation threshold.
#'
#' @return A new vector or matrix with the same dimensions as A after truncation.
#'
truncate_operator <- function(A, tau) {
  if (tau < 0) {
    stop("'tau' must be non-negative.")
  }
  
  # This is equivalent to the psi_tau function in Huber's loss derivative.
  # We use index assignment to preserve matrix dimensions correctly.
  A_trunc <- A
  A_trunc[A > tau] <- tau
  A_trunc[A < -tau] <- -tau
  
  return(A_trunc)
}

#' Compute Robust Covariance Estimators (Final Version)
#'
#' Calculates robust covariances. This version is optimized for clarity and
#' robust performance across both low-dimensional (p < n) and high-dimensional
#' (p >= n) settings.
#'
#' @param X matrix. The n x p design matrix.
#' @param y numeric vector. The n-dimensional response vector.
#' @param tau_x numeric. The truncation threshold for the feature covariance matrix.
#' @param tau_yx numeric. The truncation threshold for the response-feature
#'   covariance vector.
#'
#' @return A list containing Sigma_x and sigma_yx.
#'
compute_robust_covariances <- function(X, y, tau_x, tau_yx) {
  n <- nrow(X)
  p <- ncol(X)
  
  # --- Step 1: Vectorized computation for sigma_yx (retained) ---
  # This part is already optimal.
  yx_matrix <- as.vector(y)  * X
  trunc_yx_matrix <- truncate_operator(yx_matrix, tau_yx)
  sigma_yx <- colMeans(trunc_yx_matrix)
  
  # --- Step 2: Optimized loop-based computation for Sigma_x ---
  # This approach is memory-efficient and clear, suitable for p >= n.
  sum_trunc_xxT <- matrix(0, nrow = p, ncol = p)
  
  for (i in 1:n) {
    # tcrossprod is faster for outer products.
    # The result is truncated and added in a single, clear step.
    sum_trunc_xxT <- sum_trunc_xxT + 
      truncate_operator(tcrossprod(X[i, ]), tau_x)
  }
  
  Sigma_x <- sum_trunc_xxT / n
  
  return(
    list(
      Sigma_x = Sigma_x,
      sigma_yx = sigma_yx
    )
  )
}


# Load the lpSolve library for solving the linear programming problem.
library(lpSolve)

#' Solver for the Robust Dantzig Selector
#'
#' Solves the robust Dantzig selector problem by first computing robust
#' covariance estimates and then reformulating the problem as a linear program (LP).
#'
#' @param X matrix. The n x p design matrix.
#' @param y numeric vector. The n-dimensional response vector.
#' @param R numeric. The radius parameter for the infinity-norm constraint.
#' @param tau_x numeric. The truncation threshold for the feature covariance.
#' @param tau_yx numeric. The truncation threshold for the response-feature covariance.
#'
#' @return The estimated sparse parameter vector theta.
#'
solver_dantzig_robust <- function(X, y, R, tau_x, tau_yx) {
  p <- ncol(X)
  
  # --- Step 1: Compute robust covariance estimates ---
  robust_covs <- compute_robust_covariances(X, y, tau_x, tau_yx)
  Sigma_x <- robust_covs$Sigma_x
  sigma_yx <- robust_covs$sigma_yx
  
  # --- Step 2: Formulate the Linear Programming problem ---
  # The decision variable is z = (theta_pos, theta_neg), a 2p vector.
  
  # 2a: Objective function: min sum(z_i), so c is a vector of 1s.
  lp_objective <- rep(1, 2 * p)
  
  # 2b: Constraint matrix A for A*z <= b
  # First p constraints:  Sigma_x * theta - sigma_yx <= R
  # Second p constraints: -Sigma_x * theta + sigma_yx <= R
  A_constraint <- rbind(
    cbind(Sigma_x, -Sigma_x),
    cbind(-Sigma_x, Sigma_x)
  )
  
  # 2c: Constraint vector b
  b_constraint <- c(
    sigma_yx + R,
    -sigma_yx + R
  )
  
  # 2d: Constraint directions (all are "<=")
  lp_direction <- rep("<=", 2 * p)
  
  # --- Step 3: Solve the LP using lpSolve ---
  # The lp() function assumes non-negative variables by default, which is
  # exactly what we need for z = (theta_pos, theta_neg).
  lp_result <- lp(
    direction = "min",
    objective.in = lp_objective,
    const.mat = A_constraint,
    const.dir = lp_direction,
    const.rhs = b_constraint
  )
  
  # --- Step 4: Extract the solution and reconstruct theta ---
  if (lp_result$status != 0) {
    warning("LP solver did not converge. Status code: ", lp_result$status)
    return(rep(NA, p))
  }
  
  solution_z <- lp_result$solution
  theta_pos <- solution_z[1:p]
  theta_neg <- solution_z[(p + 1):(2 * p)]
  
  theta_hat <- theta_pos - theta_neg
  
  return(theta_hat)
}



# muti-response Dantzig Selector

#' Compute Robust Covariances for Multi-Response Models
#'
#' Calculates robust covariances for the multi-response setting.
#' Sigma_x is the feature covariance matrix, and Sigma_yx is the p x m
#' feature-response covariance matrix.
#'
#' @param X matrix. The n x p design matrix.
#' @param Y matrix. The n x m response matrix.
#' @param tau_x numeric. Truncation threshold for Sigma_x.
#' @param tau_yx numeric. Truncation threshold for Sigma_yx.
#'
#' @return A list containing Sigma_x and Sigma_yx.
#'
compute_robust_covariances_multi <- function(X, Y, tau_x, tau_yx) {
  n <- nrow(X)
  p <- ncol(X)
  m <- ncol(Y)
  
  # --- Step 1: Computation for Sigma_x (same as before) ---
  sum_trunc_xxT <- matrix(0, nrow = p, ncol = p)
  for (i in 1:n) {
    sum_trunc_xxT <- sum_trunc_xxT + 
      truncate_operator(tcrossprod(X[i, ]), tau_x)
  }
  Sigma_x <- sum_trunc_xxT / n
  
  # --- Step 2: Loop-based computation for Sigma_yx ---
  sum_trunc_xyT <- matrix(0, nrow = p, ncol = m)
  for (i in 1:n) {
    # Outer product of x_i (p x 1) and y_i^T (1 x m) results in a p x m matrix
    xyT_i <- tcrossprod(X[i, ], Y[i, ])
    
    # Truncate and add to the running sum
    sum_trunc_xyT <- sum_trunc_xyT + truncate_operator(xyT_i, tau_yx)
  }
  Sigma_yx <- sum_trunc_xyT / n
  
  return(
    list(
      Sigma_x = Sigma_x,
      Sigma_yx = Sigma_yx
    )
  )
}



# Load the CVXR library for solving the SOCP problem.
#library(CVXR)

#' Solver for the Robust Multi-Response Dantzig Selector
#'
#' Solves the robust multi-response Dantzig selector problem using CVXR.
#' This problem is a Second-Order Cone Program (SOCP).
#'
#' @param X matrix. The n x p design matrix.
#' @param Y matrix. The n x m response matrix.
#' @param R numeric. The radius parameter for the element-wise infinity-norm constraint.
#' @param tau_x numeric. The truncation threshold for the feature covariance.
#' @param tau_yx numeric. The truncation threshold for the feature-response covariance.
#'
#' @return The estimated sparse p x m parameter matrix Theta.
#'
solver_dantzig_robust_multi <- function(X, Y, R, tau_x, tau_yx) {
  p <- ncol(X)
  m <- ncol(Y)
  
  # --- Step 1: Compute robust covariance estimates ---
  robust_covs <- compute_robust_covariances_multi(X, Y, tau_x, tau_yx)
  Sigma_x <- robust_covs$Sigma_x
  Sigma_yx <- robust_covs$Sigma_yx
  
  # --- Step 2: Define the optimization problem using CVXR ---
  
  # 2a: Define the optimization variable: a p x m matrix Theta
  Theta_var <- Variable(p, m)
  
  # 2b: Define the objective function: min sum of L2 norms of rows
  # p_norm(Theta_var, 2, axis = 1) calculates the L2 norm of each row.
  objective <- Minimize(sum(p_norm(Theta_var, 2, axis = 1)))
  
  # 2c: Define the constraints
  # The infinity norm is the maximum of the absolute values of all elements.
  residual_matrix <- Sigma_x %*% Theta_var - Sigma_yx
  constraint <- list(max_entries(abs(residual_matrix)) <= R)
  
  # 2d: Assemble the problem
  problem <- Problem(objective, constraint)
  
  # --- Step 3: Solve the SOCP problem ---
  # CVXR will automatically convert this to an SOCP and call a suitable solver (e.g., ECOS).
  result <- solve(problem)
  
  # --- Step 4: Extract and return the solution ---
  if (result$status != "optimal") {
    warning("CVXR solver did not find an optimal solution. Status: ", result$status)
    return(matrix(NA, nrow = p, ncol = m))
  }
  
  Theta_hat <- result$getValue(Theta_var)
  
  return(Theta_hat)
}


# shrinkage method
# Load the glmnet library for solving the Lasso problem.
library(glmnet)

#' Solver using the Shrinkage Method
#'
#' Implements the robust estimation method described by Zhu and Zhou (2021),
#' which first truncates/shrinks the data (X and y) and then solves a
#' standard L1-penalized GLM (Lasso) on the shrunk data.
#'
#' @param X matrix. The n x p design matrix.
#' @param y numeric vector. The n-dimensional response vector.
#' @param tau_x numeric. The truncation threshold for the features.
#' @param tau_y numeric. The truncation threshold for the response.
#' @param lambda numeric. The L1 regularization parameter for glmnet.
#' @param family character. The model family, either "gaussian" for linear
#'   regression or "binomial" for logistic regression.
#' @param ... Additional arguments to be passed to the `glmnet` function.
#'
#' @return The estimated sparse parameter vector theta.
#'
solver_shrinkage_method <- function(X, y, tau_x, tau_y, lambda,
                                    family = c("gaussian", "binomial"), ...) {
  
  # Match the family argument to ensure it's one of the allowed options
  family <- match.arg(family)
  
  # --- Step 1: Shrink/Truncate the data ---
  # We use the previously defined truncate_operator.
  X_shrunk <- truncate_operator(X, tau_x)
  y_shrunk <- truncate_operator(y, tau_y)
  
  # --- Step 2: Solve the standard Lasso problem on shrunk data ---
  # glmnet is highly optimized for this task.
  # `standardize = FALSE` is important because our data is already processed.
  # `intercept = FALSE` matches the setup of our other solvers.
  # `thresh` is the convergence threshold for glmnet's coordinate descent.
  fit <- glmnet(
    x = X_shrunk,
    y = y_shrunk,
    family = family,
    lambda = lambda,
    standardize = FALSE,
    intercept = FALSE,
    thresh = 1e-7, # A standard convergence threshold
    ...
  )
  
  # --- Step 3: Extract and return the coefficient vector ---
  # The result `fit$beta` is a sparse matrix object from the Matrix package.
  # We convert it to a standard dense vector for consistency.
  theta_hat <- as.vector(fit$beta)
  
  return(theta_hat)
}

# --- File: R/solvers.R (New Multi-Response Shrinkage Solver) ---

#' Solver for Multi-Response Shrinkage Method
#'
#' Adapts the shrinkage method for multi-response linear regression. It truncates
#' the data matrices (X and Y) and then solves a Group Lasso problem using glmnet.
#'
#' @param X matrix. The n x p design matrix.
#' @param Y matrix. The n x m response matrix.
#' @param tau_x numeric. The truncation threshold for the features.
#' @param tau_y numeric. The truncation threshold for the response matrix.
#' @param lambda numeric. The L1/Group Lasso regularization parameter for glmnet.
#' @param ... Additional arguments to be passed to the `glmnet` function.
#'
#' @return The estimated p x m sparse parameter matrix Theta.
#'
solver_shrinkage_method_multi <- function(X, Y, tau_x, tau_y, lambda, ...) {
  
  # --- Step 1: Shrink/Truncate the data matrices ---
  X_shrunk <- truncate_operator(X, tau_x)
  Y_shrunk <- truncate_operator(Y, tau_y)
  
  # --- Step 2: Solve the Group Lasso problem on shrunk data ---
  # Use family = "mgaussian" for multi-response linear regression.
  fit <- glmnet(
    x = X_shrunk,
    y = Y_shrunk,
    family = "mgaussian",
    lambda = lambda,
    standardize = FALSE,
    intercept = FALSE,
    thresh = 1e-7,
    ...
  )
  
  # --- Step 3: Correctly extract and reconstruct the coefficient matrix ---
  # For "mgaussian", fit$beta is a list of sparse matrices (one for each response).
  # We need to extract each one and combine them as columns.
  
  # Get the list of sparse coefficient vectors for the specified lambda
  # Note: glmnet might return coefficients for a sequence of lambdas, so we
  # need to ensure we select the one we passed.
  # The `coef` function is a safer way to extract coefficients for a specific lambda.
  coef_list <- coef(fit, s = lambda)
  
  # The intercept is the first element of the list, beta coefficients are the second.
  # For multi-response, coef_list is a list where each element is a sparse matrix for a response.
  # We need to combine them.
  Theta_hat <- do.call(cbind, lapply(coef_list, function(sparse_matrix) {
    # Each element is a sparse matrix, first row is intercept
    as.vector(sparse_matrix)[-1] 
  }))
  
  return(Theta_hat)
}







# --- File: R/solvers.R (continued) ---

#' Gradient for the Logistic Regression Model
#'
#' Computes the gradient of the empirical negative log-likelihood (cross-entropy)
#' loss for logistic regression.
#' The gradient is grad = (1/n) * X^T * (p(X*theta) - y), where p is the
#' sigmoid function.
#'
#' @param theta numeric vector. The current p-dimensional parameter estimate.
#' @param X matrix. The n x p design matrix.
#' @param y numeric vector. The n-dimensional binary response vector (0s and 1s).
#'
#' @return A p-dimensional vector representing the gradient.
#'
grad_logistic_regression <- function(theta, X, y) {
  n <- nrow(X)
  
  # Calculate the linear predictor: X*theta
  linear_predictor <- X %*% theta
  
  # Calculate the probabilities using the logistic (sigmoid) function.
  # plogis() is R's built-in, numerically stable version of the sigmoid.
  probabilities <- plogis(linear_predictor)
  
  # Calculate the difference between predicted probabilities and actual outcomes
  diff <- probabilities - y
  
  # Calculate the gradient using the matrix form
  gradient <- (t(X) %*% diff) / n
  
  return(gradient)
}




#' Sample-wise Gradient for the Logistic Regression Model
#'
#' Computes the gradient for each individual sample in a logistic regression model.
#' The gradient for sample i is grad_i = (p_i(theta) - y_i) * x_i.
#'
#' @param theta numeric vector. The current p-dimensional parameter estimate.
#' @param X matrix. The n x p design matrix.
#' @param y numeric vector. The n-dimensional binary response vector (0s and 1s).
#'
#' @return An n x p matrix where each row is the gradient for one sample.
#'
grad_logistic_regression_samplewise <- function(theta, X, y) {
  
  # Calculate the linear predictor: X*theta
  linear_predictor <- X %*% theta
  
  # Calculate the probabilities using the logistic (sigmoid) function
  probabilities <- plogis(linear_predictor)
  
  # Calculate the difference for each sample
  diff <- probabilities - y # This is an n x 1 vector
  
  # Scale each row of X (which is x_i^T) by its corresponding difference.
  # R's recycling rule makes this efficient.
  sample_gradients <- X * as.vector(diff)
  
  return(sample_gradients)
}




#' Gradient for the Multi-Response Linear Regression Model
#'
#' Computes the gradient of the empirical least squares loss for a multi-response
#' linear model.
#' The loss is L_n(Theta) = (1 / 2n) * ||Y - X*Theta||_F^2.
#' The gradient is grad = (1 / n) * X^T * (X*Theta - Y).
#'
#' @param Theta matrix. The current p x m parameter estimate.
#' @param X matrix. The n x p design matrix.
#' @param Y matrix. The n x m response matrix.
#'
#' @return A p x m matrix representing the gradient.
#'
grad_multi_response <- function(Theta, X, Y) {
  n <- nrow(X)
  
  # Calculate the residual matrix: X*Theta - Y
  Residuals <- X %*% Theta - Y
  
  # Calculate the gradient using the matrix form
  # (p x n) %*% (n x m) -> p x m
  Gradient <- (t(X) %*% Residuals) / n
  
  return(Gradient)
}



#' Sample-wise Gradient for the Multi-Response Linear Regression Model
#'
#' Computes the gradient for each individual sample in a multi-response model.
#' The gradient for sample i is a p x m matrix: grad_i = x_i * (x_i^T*Theta - y_i^T).
#'
#' @param Theta matrix. The current p x m parameter estimate.
#' @param X matrix. The n x p design matrix.
#' @param Y matrix. The n x m response matrix.
#'
#' @return An n x p x m array where each slice [i, , ] is the gradient matrix
#'   for the i-th sample.
#'
grad_multi_response_samplewise <- function(Theta, X, Y) {
  n <- nrow(X)
  p <- ncol(X)
  m <- ncol(Y)
  
  # Calculate the n x m residual matrix
  Residuals <- X %*% Theta - Y
  
  # Initialize the 3D array to store the results
  sample_gradients <- array(NA, dim = c(n, p, m))
  
  # Loop through each sample to compute the outer product
  for (i in 1:n) {
    # The gradient for sample i is the outer product of
    # x_i (p-dim vector) and the i-th row of the residual matrix (m-dim vector).
    sample_gradients[i, , ] <- tcrossprod(X[i, ], Residuals[i, ])
  }
  
  return(sample_gradients)
}

