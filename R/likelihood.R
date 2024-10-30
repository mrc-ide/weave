log_likelihood <- function(z, mu, space_inv, time_inv, y, n_sites, n_times) {
  # Reshape (z - mu) into matrix Z and compute Z %*% time_inv %*% t(Z)
  Z <- matrix(z - mu, nrow = n_sites, ncol = n_times, byrow = TRUE)
  quad_term <- -0.5 * sum((space_inv %*% Z) * (Z %*% time_inv))  # Optimised using element-wise multiplication

  # Poisson term
  w <- !is.na(y)
  poisson_term <- sum(y[w] * z[w] - exp(z[w]))

  as.numeric(quad_term + poisson_term)
}

log_likelihood_gradient <- function(z, mu, space_inv, time_inv, y, n_sites, n_times) {
  # Reshape (z - mu) into matrix Z
  Z <- matrix(z - mu, nrow = n_sites, ncol = n_times, byrow = TRUE)

  # Compute the gradient of the quadratic term using matrix operations
  grad_quad_term_matrix <- -2 * space_inv %*% Z %*% time_inv
  grad_quad_term <- as.vector(t(grad_quad_term_matrix))

  # Poisson term (only for non-missing indices)
  grad_poisson_term <- numeric(length(y))
  w <- !is.na(y)
  grad_poisson_term[w] <- y[w] - exp(z[w])

  as.numeric(grad_quad_term + grad_poisson_term)
}

# Modified Hessian function
log_likelihood_hessian <- function(z, space_inv, time_inv, y, n_sites, n_times) {
  # Hessian of the quadratic term
  # Since the Hessian is -Sigma_inv, we can represent it using Kronecker product
  Hess_quad <- -kronecker(space_inv, time_inv)

  # Poisson term
  w <- !is.na(y)
  tmp <- rep(0, length(z))
  tmp[w] <- exp(z[w])
  Hess_poisson <- -diag(tmp)

  # Total Hessian
  Hess <- Hess_quad + Hess_poisson

  # Return the Hessian matrix
  as.matrix(Hess)
}
