log_likelihood <- function(z, dist_k_inv, time_k_inv, y, mu, n, nt, w) {
  Z <- matrix(z, nrow = n, ncol = nt, byrow = TRUE)
  quad_term <- -0.5 * sum((dist_k_inv %*% Z) * (Z %*% time_k_inv))  # Optimized using element-wise multiplication

  # Distribution term
  w <- !is.na(y)

  # Poisson term using dpois
  lambda <- exp(z[w] + mu[w])  # Mean parameter for Poisson
  likelihood_term <- sum(dpois(y[w], lambda = lambda, log = TRUE))

  as.numeric(quad_term + likelihood_term)
}

log_likelihood_gradient <- function(z, dist_k_inv, time_k_inv, y, mu, n, nt, w) {
  Z <- matrix(z, nrow = n, ncol = nt, byrow = TRUE)

  # Compute the gradient of the quadratic term using matrix operations
  grad_quad_term_matrix <- -2 * dist_k_inv %*% Z %*% time_k_inv
  grad_quad_term <- as.vector(t(grad_quad_term_matrix))

  # Initialize the gradient of the likelihood term
  grad_likelihood_term <- numeric(length(y))
  w <- !is.na(y)

  grad_likelihood_term[w] <- y[w] - exp(z[w] + mu[w])

  as.numeric(grad_quad_term + grad_likelihood_term)
}

log_likelihood_hessian <- function(z, dist_k_inv, time_k_inv, y, n, nt) {
  # Hessian of the quadratic term
  # The Hessian is -Sigma_inv, represented using the Kronecker product
  Hess_quad <- -kronecker(dist_k_inv, time_k_inv)

  # Initialize the Hessian of the likelihood term
  Hess_likelihood <- matrix(0, nrow = length(z), ncol = length(z))
  w <- !is.na(y)

  # Poisson distribution
  mu <- exp(z[w])  # Mean parameter
  tmp <- -mu  # Second derivative of Poisson log-likelihood
  Hess_likelihood[w, w] <- diag(tmp)

  # Total Hessian
  Hess <- Hess_quad + Hess_likelihood

  # Return the Hessian matrix
  return(Hess)
}
