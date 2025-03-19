log_likelihood_hyperparameters <- function(par, f, dist_matrix, time_matrix, n, nt, dist_nugget, time_nugget) {

  dist_theta <- par[1]
  periodic_scale <- par[2]
  time_theta <- par[3]

  dist_k <- rbf_kernel(dist_matrix, theta = dist_theta) + diag(x = dist_nugget, nrow = n)
  period_k <- periodic_kernel(time_matrix, periodic_scale = periodic_scale)
  long_term_k <- rbf_kernel(time_matrix, theta = time_theta)
  time_k <- period_k + long_term_k + diag(x = t_nugget, nrow = nt)
  full_k <- kronecker(dist_k, time_k)

  mvtnorm::dmvnorm(f, rep(0, n * nt), full_k, log = TRUE)
}

log_likelihood_hyperparameters_fast <- function(par, f, dist_matrix, time_matrix, n, nt, dist_nugget, time_nugget) {

  dist_theta    <- par[1]
  periodic_scale <- par[2]
  time_theta    <- par[3]

  dist_k <- rbf_kernel(dist_matrix, theta = dist_theta) + diag(x = dist_nugget, nrow = n)
  period_k <- periodic_kernel(time_matrix, periodic_scale = periodic_scale)
  long_term_k <- rbf_kernel(time_matrix, theta = time_theta)
  time_k <- period_k + long_term_k + diag(x = t_nugget, nrow = nt)

  # Cholesky + log det
  dist_chol   <- chol(dist_k)
  time_chol   <- chol(time_k)
  logdet_dist <- 2 * sum(log(diag(dist_chol)))
  logdet_time <- 2 * sum(log(diag(time_chol)))
  logdet_full <- nt * logdet_dist + n * logdet_time

  # Inverse pieces
  dist_inv <- chol2inv(dist_chol)  # n x n
  time_k_inv <- chol2inv(time_chol)  # nt x nt

  # Quad form: reshape f into X (nt x n)
  X <- matrix(f, nrow = nt, ncol = n)
  # f^T Sigma^{-1} f = trace( dist_inv %*% (X^T time_k_inv X) )
  quad_form <- sum( dist_inv * (t(X) %*% time_k_inv %*% X) )

  # final log-likelihood
  loglik <- -0.5 * (quad_form + logdet_full + n * nt * log(2 * pi))
  return(loglik)
}
