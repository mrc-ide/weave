#' Quick multivariate normal draw over two dimensions
#'
#' Draws one sample from a zero-mean Gaussian with separable space-time
#' covariance, without ever forming the full matrix. This is equivalent to
#' forming the full spatio-temporal covariance and drawing from the
#' multivariate normal distribution:
#'
#' ```r
#' full_k <- kronecker(space, time)
#' f <- MASS::mvrnorm(1, rep(0, nrow(space) * nrow(time)), full_k)
#' ```
#'
#' @param space Space kernel matrix.
#' @param time Time kernel matrix.
#'
#' @return A numeric vector of length `nrow(space) * nrow(time)`, ordered
#'   site by site with time varying fastest (matching
#'   `kronecker(space, time)`).
#' @export
quick_mvnorm <- function(space, time) {
  n_sites <- nrow(space)
  n_times <- nrow(time)

  # Cholesky factors: L_s lower (so L_s %*% t(L_s) = space), L_t upper
  L_s <- t(chol(space))  # lower-tri
  L_t <- chol(time)      # upper-tri

  # i.i.d. standard normals arranged as [sites x times]
  W <- matrix(stats::rnorm(n_sites * n_times), nrow = n_sites, ncol = n_times)

  # Apply separable transforms. vec(Z) has covariance (time ⊗ space); the
  # transpose-flatten below reorders to times-fastest, giving covariance
  # (space ⊗ time) to match kronecker(space, time).
  Z <- L_s %*% W %*% L_t

  # Flatten with times varying fastest
  as.vector(t(Z))
}


#' Quick multivariate normal draw over two dimensions (Cholesky precomputed)
#'
#' As [quick_mvnorm()], but taking precomputed Cholesky factors so repeated
#' draws (e.g. the perturbation draws in [gp_predict()]) skip the
#' factorisation cost.
#'
#' @param space_chol Upper-triangular Cholesky factor of the space kernel
#'   matrix, as returned by [chol()]. Passing the lower-triangular factor
#'   gives silently wrong draws.
#' @param time_chol Upper-triangular Cholesky factor of the time kernel
#'   matrix, as returned by [chol()].
#'
#' @return A numeric vector of length `nrow(space_chol) * nrow(time_chol)`,
#'   ordered site by site with time varying fastest (matching
#'   `kronecker(space, time)`).
#' @export
quick_mvnorm_chol <- function(space_chol, time_chol) {
  n_sites <- nrow(space_chol)
  n_times <- nrow(time_chol)

  # i.i.d. standard normals arranged as [sites x times]
  W <- matrix(stats::rnorm(n_sites * n_times), nrow = n_sites, ncol = n_times)

  # Apply separable transforms. vec(Z) has covariance (time ⊗ space); the
  # transpose-flatten below reorders to times-fastest, giving covariance
  # (space ⊗ time) to match kronecker(space, time).
  Z <- crossprod(space_chol, W) %*% time_chol

  # Flatten with times varying fastest
  as.vector(t(Z))
}
