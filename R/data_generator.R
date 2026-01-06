#' --- File: R/data_generator.R ---

#' Generate Data for Sparse Linear Regression
#'
#' This is the main function to generate a dataset (y, X) for the sparse
#' linear model y = X * theta_star + epsilon. It uses a flexible design where
#' the specific random generators for features (X) and noise (epsilon) are
#' passed as function arguments.
#'
#' @param n integer. The number of samples.
#' @param p integer. The number of features (dimension).
#' @param theta_star numeric vector. The true p-dimensional coefficient vector.
#'   The function infers the sparsity level s_star from this vector.
#' @param generator_X function. A function to generate the n x p feature matrix X.
#'   It must accept `n` and `p` as its first two arguments.
#' @param generator_epsilon function. A function to generate the n-dimensional
#'   noise vector. It must accept `n` as its first argument.
#' @param ... Additional arguments to be passed down to `generator_X` and
#'   `generator_epsilon`. For example, `df` for a t-distribution.
#'
#' @return A list containing the following elements:
#'   \item{y}{The n-dimensional response vector.}
#'   \item{X}{The n x p feature matrix.}
#'   \item{theta_star}{The true p-dimensional coefficient vector.}
#'   \item{epsilon}{The n-dimensional noise vector.}
#'
generate_linear_data <- function(n, p, theta_star,
                                 generator_X, generator_epsilon, ...) {

  # Generate features and noise using the provided generators ---
  # The '...' argument captures any additional named arguments and passes them
  # to the generator functions.
  args_list <- list(...)

  X <- do.call(generator_X, c(list(n = n, p = p, theta_star = theta_star), args_list))
  epsilon <- do.call(generator_epsilon, c(list(n = n), args_list))

  # --- Step 3: Generate the response variable ---
  # Using matrix multiplication for efficiency.
  y <- X %*% theta_star + epsilon

  # --- Step 4: Return the complete dataset as a list ---
  return(
    list(
      y = y,
      X = X,
      theta_star = theta_star,
      epsilon = epsilon
    )
  )
}




#' Generate Data for Sparse Logistic Regression
#'
#' This function generates a dataset (y, X) for the sparse logistic regression
#' model. The binary response y is generated from a Bernoulli distribution with
#' probability p = 1 / (1 + exp(-X * theta_star)).
#'
#' @param n integer. The number of samples.
#' @param p integer. The number of features (dimension).
#' @param theta_star numeric vector. The true p-dimensional coefficient vector.
#' @param generator_X function. A function to generate the n x p feature matrix X.
#'   It must accept `n`, `p`, and `theta_star` as arguments.
#' @param ... Additional arguments to be passed down to `generator_X`.
#'
#' @return A list containing the following elements:
#'   \item{y}{The n-dimensional binary response vector (0s and 1s).}
#'   \item{X}{The n x p feature matrix.}
#'   \item{theta_star}{The true p-dimensional coefficient vector.}
#'
generate_logistic_data <- function(n, p, theta_star,
                                   generator_X, ...) {

  # --- Step 1: Generate features using the provided generator ---
  args_list <- list(...)
  X <- do.call(generator_X, c(list(n = n, p = p, theta_star = theta_star), args_list))

  # --- Step 2: Generate the response variable using the logistic model ---
  # Calculate the linear predictor
  linear_predictor <- X %*% theta_star
  # Transform to probabilities using the logistic (sigmoid) function
  probabilities <- plogis(linear_predictor)
  # Generate binary outcomes from a Bernoulli distribution
  y <- rbinom(n, size = 1, prob = probabilities)

  # --- Step 3: Return the complete dataset as a list ---
  return(
    list(
      y = y,
      X = X,
      theta_star = theta_star
    )
  )
}






#' Generate Data for Sparse Multi-Response Regression
#'
#' Generates a dataset (Y, X) for the sparse multi-response linear model
#' Y = X * Theta_star + E.
#'
#' @param n integer. The number of samples.
#' @param p integer. The number of features.
#' @param m integer. The number of responses.
#' @param Theta_star matrix. The true p x m coefficient matrix.
#' @param generator_X function. A function to generate the n x p feature matrix X.
#' @param generator_E function. A function to generate the n x m noise matrix E.
#' @param ... Additional arguments for the generator functions.
#'
#' @return A list containing Y, X, Theta_star, and E.
#'
generate_multi_response_data <- function(n, p, m, Theta_star,
                                         generator_X, generator_E, ...) {

  # --- Step 1: Validate input parameters ---
  if (!is.matrix(Theta_star) || nrow(Theta_star) != p || ncol(Theta_star) != m) {
    stop("'Theta_star' must be a p x m matrix.")
  }

  # --- Step 2: Generate features and noise ---
  args_list <- list(...)
  # Note: we pass Theta_star to generator_X for the mixed case
  X <- do.call(generator_X, c(list(n = n, p = p, theta_star = Theta_star), args_list))
  E <- do.call(generator_E, c(list(n = n, m = m), args_list))

  # --- Step 3: Generate the response matrix ---
  Y <- X %*% Theta_star + E

  # --- Step 4: Return the complete dataset ---
  return(
    list(
      Y = Y,
      X = X,
      Theta_star = Theta_star,
      E = E
    )
  )
}



library(MASS)
generate_X_gaussian <- function(n, p, Sigma_X = NULL, ...) {
  # If no covariance matrix is provided, default to the identity matrix I_p.
  if (is.null(Sigma_X)) {
    Sigma_X <- diag(p)
  }
  # Validate dimensions of Sigma_X
  if (!is.matrix(Sigma_X) || nrow(Sigma_X) != p || ncol(Sigma_X) != p) {
    stop("'Sigma_X' must be a p x p matrix.")
  }
  X <- mvrnorm(n = n, mu = rep(0, p), Sigma = Sigma_X)

  return(X)
}

library(mvtnorm)
generate_X_t <- function(n, p, df_X, Sigma_X = NULL, scale_X = 1, ...) {

  # Validate degrees of freedom
  if (df_X <= 0) {
    stop("'df_X' must be positive.")
  }

  # If no scale matrix is provided, default to the identity matrix I_p.
  if (is.null(Sigma_X)) {
    Sigma_X <- diag(p)
  }

  # Validate dimensions of Sigma_X
  if (!is.matrix(Sigma_X) || nrow(Sigma_X) != p || ncol(Sigma_X) != p) {
    stop("'Sigma_X' must be a p x p matrix.")
  }

  # Generate data using the rmvt function from the mvtnorm package.
  # The mean vector is implicitly zero.
  X <- rmvt(n = n, sigma = Sigma_X, df = df_X)

  # Apply the final scaling factor
  X <- X * scale_X

  return(X)
}






#' Generator for mixed-distribution features (X)
#' @param Sigma_S_X matrix. The covariance matrix for the Gaussian part on the
#'   support set. Must be s x s, where s is the sparsity.
#' @param Sigma_Sc_X matrix. The scale matrix for the t-distributed part outside
#'   the support. Must be (p-s) x (p-s).
#' @param df_X numeric. The degrees of freedom for the t-part.
#' @param scale_X numeric. A scaling factor for the t-part. Defaults to 1.
#' @param ... Catches and ignores any other arguments.
#'
#' @return An n x p matrix of features.
#'
generate_X_mixed <- function(n, p, theta_star,
                             Sigma_S_X, Sigma_Sc_X, df_X,
                             scale_X = 1, ...) {
  # --- Identify support sets (works for both vector and matrix theta_star) ---
  if (is.matrix(theta_star)) {
    # For matrix, support is defined by non-zero rows
    support_indices <- which(rowSums(abs(theta_star)) != 0)
  } else {
    # For vector, support is defined by non-zero elements
    support_indices <- which(theta_star != 0)
  }
  if (is.matrix(theta_star)) {
    # For matrix, non-support is defined by zero rows
    non_support_indices <- which(rowSums(abs(theta_star)) == 0)
  } else {
    # For vector, non-support is defined by zero elements
    non_support_indices <- which(theta_star == 0)
  }
  s_star <- length(support_indices)
  p_minus_s <- p - s_star
  # --- Initialize the final X matrix ---
  X <- matrix(0.0, nrow = n, ncol = p)

  # --- Generate and fill in the parts ---
  # Generate Gaussian part for the support columns
  if (s_star > 0) {
    X[, support_indices] <- mvrnorm(n, mu = rep(0, s_star), Sigma = Sigma_S_X)
  }

  # Generate t-distributed part for the non-support columns
  if (p_minus_s > 0) {
    X_Sc <- rmvt(n, sigma = Sigma_Sc_X, df = df_X)
    X[, non_support_indices] <- X_Sc * scale_X
  }

  return(X)
}








#' Generator for Gaussian noise (epsilon)
generate_epsilon_gaussian <- function(n, sigma2_eps = 1, ...) {
  epsilon <- rnorm(n = n, mean = 0, sd = sqrt(sigma2_eps))
  return(epsilon)
}




#' Generator for t-distributed noise (epsilon)
generate_epsilon_t <- function(n, df_eps, scale_eps = 1, ...) {
  # Generate noise using the standard rt function.
  epsilon <- rt(n = n, df = df_eps)
  # Apply the final scaling factor
  epsilon <- epsilon * scale_eps
  return(epsilon)
}


#' Generator for multivariate Gaussian noise matrix (E)
#'
#' Generates an n x m matrix where each row is a sample from N_m(0, Sigma_E).
#'
#' @param n integer. The number of samples.
#' @param m integer. The number of responses.
#' @param Sigma_E matrix. The m x m covariance matrix for the noise.
#'   Defaults to the identity matrix.
#' @param ... Catches and ignores any other arguments.
#'
#' @return An n x m matrix of noise.
#'
generate_E_gaussian <- function(n, m, Sigma_E = NULL, ...) {
  if (is.null(Sigma_E)) {
    Sigma_E <- diag(m)
  }
  if (!is.matrix(Sigma_E) || nrow(Sigma_E) != m || ncol(Sigma_E) != m) {
    stop("'Sigma_E' must be an m x m matrix.")
  }
  E <- mvrnorm(n = n, mu = rep(0, m), Sigma = Sigma_E)
  return(E)
}


#' Generator for multivariate t-distributed noise matrix (E)
#' @param Sigma_E matrix. The m x m scale matrix. Defaults to identity.
#' @param scale_E numeric. A final scaling factor. Defaults to 1.
#' @param ... Catches and ignores any other arguments.
#'
#' @return An n x m matrix of noise.
#'
generate_E_t <- function(n, m, df_E, Sigma_E = NULL, scale_E = 1, ...) {
  if (df_E <= 0) {
    stop("'df_E' must be positive.")
  }
  if (is.null(Sigma_E)) {
    Sigma_E <- diag(m)
  }
  if (!is.matrix(Sigma_E) || nrow(Sigma_E) != m || ncol(Sigma_E) != m) {
    stop("'Sigma_E' must be an m x m matrix.")
  }
  E <- rmvt(n = n, sigma = Sigma_E, df = df_E)
  E <- E * scale_E
  return(E)
}