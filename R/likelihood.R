# Modified log-likelihood function
log_likelihood <- function(z, mu, space_inv, time_inv, y, n_sites, n_times) {
  # Reshape (z - mu) into a matrix Z
  Z <- matrix(z - mu, nrow = n_sites, ncol = n_times, byrow = TRUE)

  # Compute the quadratic term using the correct order
  quad_term <- -0.5 * sum(diag(space_inv %*% Z %*% time_inv %*% t(Z)))

  # Poisson term
  w <- !is.na(y)
  poisson_term <- sum(y[w] * z[w] - exp(z[w]))

  as.numeric(quad_term + poisson_term)
}

log_likelihood_gradient <- function(z, mu, space_inv, time_inv, y, n_sites, n_times) {
  # Reshape (z - mu) into a matrix Z
  Z <- matrix(z - mu, nrow = n_sites, ncol = n_times, byrow = TRUE)

  # Compute the gradient of the quadratic term
  grad_quad_term_matrix <- -2 * space_inv %*% Z %*% time_inv

  # Flatten the gradient matrix
  grad_quad_term <- as.vector(t(grad_quad_term_matrix))

  # Poisson term
  w <- !is.na(y)
  grad_poisson_term <- rep(0, length(y))
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
