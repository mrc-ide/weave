log_likelihood <- function(z, mu, space_inv, time_inv, y, n_sites, n_times, dist = "poisson", size = NULL) {
  # Reshape (z - mu) into matrix Z and compute Z %*% time_inv %*% t(Z)
  Z <- matrix(z - mu, nrow = n_sites, ncol = n_times, byrow = TRUE)
  quad_term <- -0.5 * sum((space_inv %*% Z) * (Z %*% time_inv))  # Optimized using element-wise multiplication

  # Distribution term
  w <- !is.na(y)

  if (dist == "poisson") {
    # Poisson term using dpois
    lambda <- exp(z[w])  # Mean parameter for Poisson
    poisson_term <- sum(dpois(y[w], lambda = lambda, log = TRUE))
    likelihood_term <- poisson_term
  } else if (dist == "nbinom") {
    if (is.null(size)) {
      stop("Size parameter must be provided for negative binomial distribution.")
    }
    # Negative Binomial term using dnbinom
    mu_i <- exp(z[w])  # Mean parameter for Negative Binomial
    nb_term <- sum(dnbinom(y[w], size = size, mu = mu_i, log = TRUE))
    likelihood_term <- nb_term
  } else {
    stop("Unsupported distribution")
  }

  as.numeric(quad_term + likelihood_term)
}

log_likelihood_gradient <- function(z, mu, space_inv, time_inv, y, n_sites, n_times, dist = "poisson", size = NULL) {
  # Reshape (z - mu) into matrix Z
  Z <- matrix(z - mu, nrow = n_sites, ncol = n_times, byrow = TRUE)

  # Compute the gradient of the quadratic term using matrix operations
  grad_quad_term_matrix <- -2 * space_inv %*% Z %*% time_inv
  grad_quad_term <- as.vector(t(grad_quad_term_matrix))

  # Initialize the gradient of the likelihood term
  grad_likelihood_term <- numeric(length(y))
  w <- !is.na(y)

  if (dist == "poisson") {
    # Poisson term (only for non-missing indices)
    grad_likelihood_term[w] <- y[w] - exp(z[w])
  } else if (dist == "nbinom") {
    if (is.null(size)) {
      stop("Size parameter must be provided for negative binomial distribution.")
    }
    # Negative Binomial term
    mu <- exp(z[w])  # Mean parameter
    grad_likelihood_term[w] <- (y[w] - mu) * size / (size + mu)
  } else {
    stop("Unsupported distribution")
  }

  as.numeric(grad_quad_term + grad_likelihood_term)
}

log_likelihood_hessian <- function(z, space_inv, time_inv, y, n_sites, n_times, dist = "poisson", size = NULL) {
  # Hessian of the quadratic term
  # The Hessian is -Sigma_inv, represented using the Kronecker product
  Hess_quad <- -kronecker(space_inv, time_inv)

  # Initialize the Hessian of the likelihood term
  Hess_likelihood <- matrix(0, nrow = length(z), ncol = length(z))
  w <- !is.na(y)

  if (dist == "poisson") {
    # Poisson distribution
    mu <- exp(z[w])  # Mean parameter
    tmp <- -mu  # Second derivative of Poisson log-likelihood
    Hess_likelihood[w, w] <- diag(tmp)
  } else if (dist == "nbinom") {
    if (is.null(size)) {
      stop("Size parameter must be provided for the negative binomial distribution.")
    }
    # Negative Binomial distribution
    mu <- exp(z[w])    # Mean parameter
    r <- size          # Dispersion parameter
    numerator <- r * (r + y[w])
    denominator <- (r + mu)^2
    tmp <- -mu * numerator / denominator  # Second derivative
    Hess_likelihood[w, w] <- diag(tmp)
  } else {
    stop("Unsupported distribution")
  }

  # Total Hessian
  Hess <- Hess_quad + Hess_likelihood

  # Return the Hessian matrix
  return(Hess)
}



