#' Add a small ridge to a square matrix
#'
#' In plain terms: this adds a tiny value to the diagonal so the matrix is
#' better-behaved numerically (e.g., invertible and Cholesky-able).
#'
#' Technically: returns \eqn{X + \lambda I}, which improves condition number and
#' ensures positive definiteness when \eqn{\lambda > 0}.
#'
#' @param x A square numeric matrix.
#' @param lambda Non-negative ridge (diagonal) value to add. Default `1e-5`.
#'
#' @return A matrix the same size as `x` with `lambda` added to the diagonal.
regularise <- function(x, lambda = 1e-5) {
  x + lambda * diag(nrow(x))
}

#' Fast Kronecker–product matrix–vector multiply (times vary fastest)
#'
#' In plain terms: multiplies a big covariance `K = space ⊗ time` by a vector
#' without ever forming `K`, using a reshape–multiply–reshape trick.
#'
#' Technically: for \eqn{v = \mathrm{vec}(X^\top)} with `times` varying fastest,
#' computes \eqn{(space \otimes time)\,v = \mathrm{vec}\!\big((space\,X\,time^\top)^\top\big)}.
#'
#' @param v Numeric vector of length `nrow(space) * nrow(time)`, ordered with
#'   times varying fastest within site.
#' @param space Spatial kernel matrix (size \eqn{n \times n}).
#' @param time Temporal kernel matrix (size \eqn{nt \times nt}).
#'
#' @return A numeric vector the same length as `v`.
kron_mv <- function(v, space, time) {
  n_sites <- nrow(space); n_times <- nrow(time)
  # Reconstruct X (sites × times) from vec(t(X)) = v
  X <- t(matrix(v, nrow = n_times, ncol = n_sites))
  # Apply K to X: vec(t(space %*% X %*% t(time))) = (space ⊗ time) vec(t(X))
  Y <- space %*% X %*% t(time)
  as.vector(t(Y))
}

#' Fill observed values into a full vector
#'
#' Put the observed entries back into their full-length vector,
#' filling missing positions with zeros.
#'
#' @param x_obs Numeric vector of observed values (length \eqn{m}).
#' @param obs_idx Integer indices (length \eqn{m}) of observed positions in the
#'   full vector.
#' @param N Total length of the full vector.
#'
#' @return A numeric vector of length `N` with `x_obs` scattered at `obs_idx`.
fill_vector <- function(x_obs, obs_idx, N){
  v <- numeric(N)
  v[obs_idx] <- x_obs
  v
}

#' Observed-system matvec: (S K S^T + diag(noise)) v
#'
#' Takes an observed-length vector and applies the GP covariance
#' plus a per-observation noise (nugget), all without building any big matrices.
#'
#' Technically: for \eqn{K = \mathrm{space}\,\otimes\,\mathrm{time}}, returns
#' \deqn{S\,K\,S^{\mathsf T}\,v \;+\; \operatorname{diag}(\sigma^2)\,v,}
#' i.e., the observed block of the GP plus a diagonal nugget. Implemented
#' matrix-free as \code{kron_mv(with_nas(v, obs_idx, N), space_mat, time_mat)[obs_idx] + noise_var * v},
#' where \eqn{S^{\mathsf T}} “scatters’’ into the full vector and \eqn{\sigma^2}
#' denotes the per-observation noise.
#'
#' @param v Numeric vector of length \eqn{m} (observed entries).
#' @param obs_idx Integer indices of observed entries in the full vector.
#' @param N Total length of the full vector.
#' @param space_mat Spatial kernel matrix.
#' @param time_mat Temporal kernel matrix.
#' @param noise_var Scalar or length-\eqn{m} numeric nugget on the observed scale.
#'
#' @return A numeric vector of length \eqn{m}, equal to \eqn{(S K S^\top + D)v}.
Amv <- function(v, obs_idx, N, space_mat, time_mat, noise_var) {
  #browser()
  kron_mv(fill_vector(v, obs_idx, N), space_mat, time_mat)[obs_idx] + noise_var * v
}

#' Diagonal (Jacobi) preconditioner application
#'
#' Divides by an approximation to the diagonal of the system,
#' which makes the iterative solver converge faster.
#'
#' Technically: applies \eqn{M^{-1} v \approx v / \mathrm{diag}(A)}, where
#' \eqn{A = S K S^\top + \mathrm{diag}(\text{noise})} and
#' \eqn{\mathrm{diag}(K) = \mathrm{diag}(space) \otimes \mathrm{diag}(time)}.
#'
#' @param v Numeric vector to precondition (length \eqn{m}).
#' @param kdiag_full Vector \eqn{\mathrm{diag}(K)} of length \eqn{N}
#'   (typically from `as.vector(kronecker(diag(space), diag(time)))`).
#' @param obs_idx Integer indices of observed entries in the full vector.
#' @param noise_var Scalar or length-\eqn{m} numeric nugget to add to the diagonal.
#'
#' @return A numeric vector of length \eqn{m}: elementwise `v / (diagA + 1e-12)`.
M_inv <- function(v, kdiag_full, obs_idx, noise_var) {
  v / (kdiag_full[obs_idx] + noise_var + 1e-12)
}

#' Preconditioned Conjugate Gradient (PCG) solver for the observed system
#'
#' In plain terms: solves the big linear system that gives the GP weights using
#' only matrix–vector products—no huge matrices, no explicit inverse.
#'
#' Technically: solves \eqn{(S K S^\top + \mathrm{diag}(\text{noise}))\,x = b}
#' by PCG, using `Amv` for matrix–vector products and `M_inv` as a Jacobi
#' preconditioner. Stops when the relative residual falls below `tol` or after
#' `maxit` iterations (issues a warning on `maxit`).
#'
#' @param b Right-hand side vector (observed length \eqn{m}).
#' @param obs_idx Integer indices of observed entries in the full vector.
#' @param N Total length of the full vector.
#' @param space_mat Spatial kernel matrix.
#' @param time_mat Temporal kernel matrix.
#' @param noise_var Scalar or length-\eqn{m} nugget on the observed scale.
#' @param kdiag_full Vector \eqn{\mathrm{diag}(K)} of length \eqn{N}.
#' @param tol Relative residual tolerance for convergence (default `1e-8`).
#' @param maxit Maximum number of iterations (default `10000`).
#'
#' @return Numeric solution vector `x` of length \eqn{m}.
pcg <- function(b, obs_idx, N, space_mat, time_mat, noise_var, kdiag_full, tol = 1e-8, maxit = 10000) {
  x <- numeric(length(b))
  r <- b - Amv(x, obs_idx, N, space_mat, time_mat, noise_var)
  z <- M_inv(r, kdiag_full, obs_idx, noise_var)
  p <- z
  rz_old <- sum(r * z)
  for (it in seq_len(maxit)) {
    if(it == maxit){
      warning("maxit reached")
    }
    Ap <- Amv(p, obs_idx, N, space_mat, time_mat, noise_var)
    alpha <- rz_old / sum(p * Ap)
    x <- x + alpha * p
    r <- r - alpha * Ap
    if (sqrt(sum(r * r)) <= tol * sqrt(sum(b * b))) break
    z <- M_inv(r, kdiag_full, obs_idx, noise_var)
    rz_new <- sum(r * z)
    beta <- rz_new / rz_old
    p <- z + beta * p
    rz_old <- rz_new
  }
  x
}
