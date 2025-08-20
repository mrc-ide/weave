#' Quick multivariate normal samples over two dimensions
#'
#' This is equivalent to estimating the full spatio-temporal covariance matrix
#' and sampling from the multivariate normal distribution:
#' full_k <- kronecker(dist_k, time_k)
#' f  <- mvrnorm(1, rep(0, n * nt), full_k)
#'
#' @param space Space kernel matrix
#' @param time  Time kernel matrix
quick_mvnorm <- function(space, time) {
  n_sites <- nrow(space)
  n_times <- nrow(time)

  # Cholesky factors: L_s lower (so L_s %*% t(L_s) = space), L_t upper
  L_s <- t(chol(space))  # lower-tri
  L_t <- chol(time)      # upper-tri

  # i.i.d. standard normals arranged as [sites x times]
  W <- matrix(rnorm(n_sites * n_times), nrow = n_sites, ncol = n_times)

  # Apply separable transforms: Z has cov(time ⊗ space)
  Z <- L_s %*% W %*% L_t

  # Flatten with times varying fastest
  as.vector(t(Z))
}
