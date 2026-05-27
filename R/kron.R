# =============================================================================
# Kronecker linear algebra for separable spatio-temporal kernels
#
# We model the GP covariance as a Kronecker product K = K_space (x) K_time of
# an n x n spatial kernel and an nt x nt temporal kernel. The full covariance
# never has to be formed: every operation we need can be implemented as a
# small matrix multiply on a reshaped vector.
#
# This file collects all the Kronecker primitives used by the deterministic
# fit() and the PG-Gibbs sampler:
#
#   kron_mv()         -- matrix-vector multiply K v
#   kron_eigen()      -- cached eigendecompositions of K_space and K_time
#   kron_solve_eigen()-- (K + s2 I)^{-1} v using the eigen cache
#   kron_quad()       -- v' K^{-1} v
#
# Convention: latent fields are stored as a length-(n*nt) vector with the
# TIME dimension varying fastest within site. That is, if X is the n x nt
# matrix of latent values then the stored vector is as.vector(t(X)).
# =============================================================================


#' Cholesky with automatic relative jitter and retry on numerical failure
#'
#' `chol()` errors when a matrix is numerically indefinite. Inside the
#' MCMC sampler we cannot tolerate a crash whenever theta wanders into a
#' regime where K is near-singular -- it kills the entire chain. This
#' wrapper adds a tiny relative jitter (`jitter * mean(diag(A))` per
#' diagonal entry) and retries with the jitter doubled if the
#' factorisation still fails, up to `max_jitter`. Returns the lower-
#' triangular factor L such that `L L^T = A + j * mean(diag(A)) * I` for
#' the smallest successful j.
#'
#' For well-conditioned matrices the perturbation is invisible (jitter of
#' 1e-6 on diag entries of order 1 changes downstream quantities by less
#' than machine precision in any practical metric).
#'
#' @param A A square symmetric numeric matrix.
#' @param jitter Initial relative jitter (default 1e-6).
#' @param max_jitter Maximum jitter to try before giving up (default 1e-2).
#' @return Lower-triangular numeric matrix L with `L %*% t(L) = A + j*I_scaled`.
safe_chol <- function(A, jitter = 1e-6, max_jitter = 1e-2) {
  diag_scale <- mean(diag(A))
  if (!is.finite(diag_scale) || diag_scale <= 0) diag_scale <- 1
  n <- nrow(A)
  j <- jitter
  repeat {
    Ar <- A + (j * diag_scale) * diag(n)
    out <- tryCatch(t(chol(Ar)), error = function(e) NULL)
    if (!is.null(out)) return(out)
    if (j > max_jitter) {
      stop("safe_chol: Cholesky failed even with jitter = ",
           sprintf("%.1e", j),
           ". The matrix is severely ill-conditioned.")
    }
    j <- j * 10
  }
}


#' Symmetric eigendecomposition with automatic relative jitter
#'
#' Mirror of [safe_chol()] for symmetric eigendecomposition. Adds a tiny
#' relative jitter before calling `eigen(_, symmetric = TRUE)`, retrying
#' with bigger jitter if the call errors (rare but possible on extreme
#' inputs). Also clamps the returned eigenvalues at `eig_floor` so that
#' downstream code dividing by them stays finite.
#'
#' @param A A square symmetric numeric matrix.
#' @param jitter Initial relative jitter (default 1e-6).
#' @param max_jitter Maximum jitter to try (default 1e-2).
#' @param eig_floor Lower bound applied to eigenvalues (default 1e-12).
#' @return A list with `$values` (clamped) and `$vectors`, same shape
#'   as [base::eigen()].
safe_eigen <- function(A, jitter = 1e-6, max_jitter = 1e-2,
                       eig_floor = 1e-12) {
  diag_scale <- mean(diag(A))
  if (!is.finite(diag_scale) || diag_scale <= 0) diag_scale <- 1
  n <- nrow(A)
  j <- jitter
  repeat {
    Ar <- A + (j * diag_scale) * diag(n)
    out <- tryCatch(eigen(Ar, symmetric = TRUE), error = function(e) NULL)
    if (!is.null(out)) {
      out$values <- pmax(out$values, eig_floor)
      return(out)
    }
    if (j > max_jitter) {
      stop("safe_eigen: eigen() failed even with jitter = ",
           sprintf("%.1e", j),
           ". The matrix is severely ill-conditioned.")
    }
    j <- j * 10
  }
}


#' Fast Kronecker-product matrix-vector multiply (times vary fastest)
#'
#' In plain terms: multiplies a big covariance K = K_space (x) K_time by a
#' vector without ever forming K, using a reshape-multiply-reshape trick.
#'
#' Technically: for v = vec(t(X)) with times varying fastest,
#'   (K_space (x) K_time) v = vec(t(K_space %*% X %*% t(K_time)))
#'
#' @param v Numeric vector of length nrow(space) * nrow(time).
#' @param space Spatial kernel matrix (n x n).
#' @param time Temporal kernel matrix (nt x nt).
#' @return Numeric vector of the same length as v.
kron_mv <- function(v, space, time) {
  n_sites <- nrow(space)
  n_times <- nrow(time)
  X <- t(matrix(v, nrow = n_times, ncol = n_sites))
  Y <- space %*% X %*% t(time)
  as.vector(t(Y))
}


#' Cache eigendecompositions of the spatial and temporal kernels
#'
#' The Kronecker product factors through the eigendecomposition:
#'   K_space = U_s diag(L_s) U_s'
#'   K_time  = U_t diag(L_t) U_t'
#'   K       = (U_s (x) U_t) diag(L_s (x) L_t) (U_s (x) U_t)'
#'
#' Once we have these we can compute `(K + s2 I)^-1 v`, `K^-1 v`, and
#' log|K| at O(n * nt) cost per application -- which is the only reason
#' we can afford to slice on theta inside the Gibbs sampler.
#'
#' Returned object has fields:
#'   U_s, L_s, U_t, L_t -- eigenvectors / eigenvalues
#'   lam_full           -- outer product L_s x L_t, flattened with times
#'                         varying fastest, matching kron_mv's convention
#'   log_det            -- log|K| = nt * sum(log L_s) + n * sum(log L_t)
#'
#' Eigenvalues are clamped at a tiny positive number to keep the inverse
#' well-defined when the kernels are numerically singular.
#'
#' @param space Spatial kernel matrix (n x n, symmetric PD).
#' @param time Temporal kernel matrix (nt x nt, symmetric PD).
#' @param eig_floor Lower bound for eigenvalues; defaults to 1e-12.
#' @return A list with the cached components above.
kron_eigen <- function(space, time, eig_floor = 1e-12) {
  es <- safe_eigen(space, eig_floor = eig_floor)
  et <- safe_eigen(time,  eig_floor = eig_floor)

  # safe_eigen already clamped the eigenvalues; pmax here is a no-op
  # belt-and-braces.
  L_s <- pmax(es$values, eig_floor)
  L_t <- pmax(et$values, eig_floor)

  # lam_full[(i-1)*nt + j] = L_s[i] * L_t[j]; this is the vec(t(.)) order
  # used everywhere else in the package. outer() is column-major, so we
  # transpose before flattening.
  lam_full <- as.vector(t(outer(L_s, L_t)))

  list(
    U_s     = es$vectors,
    L_s     = L_s,
    U_t     = et$vectors,
    L_t     = L_t,
    lam_full = lam_full,
    log_det = nrow(time) * sum(log(L_s)) + nrow(space) * sum(log(L_t))
  )
}


#' Solve (K + s2 I) x = v using the eigen cache
#'
#' In the (U_s (x) U_t) basis the system is diagonal with entries
#' (L_s_i * L_t_j + s2), so the whole solve is a reshape, a divide, and
#' another reshape. With s2 = 0 this gives `K^-1 v` directly.
#'
#' Used both as the preconditioner inside PCG and to evaluate the GP
#' log-density in the slice sampler.
#'
#' @param v Numeric vector of length n * nt.
#' @param ke Eigen cache from [kron_eigen()].
#' @param sigma2 Non-negative scalar nugget added to every eigenvalue.
#' @return Numeric vector `x = (K + sigma2 I)^-1 v`.
kron_solve_eigen <- function(v, ke, sigma2 = 0) {
  n_sites <- nrow(ke$U_s)
  n_times <- nrow(ke$U_t)

  # Rotate v -> w in the eigenbasis.
  X <- t(matrix(v, nrow = n_times, ncol = n_sites))   # X is n x nt
  W <- t(ke$U_s) %*% X %*% ke$U_t                     # W is n x nt
  w <- as.vector(t(W))                                # back to vec(t(.))

  # Diagonal solve in eigen-space.
  z <- w / (ke$lam_full + sigma2)

  # Rotate back.
  Z <- t(matrix(z, nrow = n_times, ncol = n_sites))
  Y <- ke$U_s %*% Z %*% t(ke$U_t)
  as.vector(t(Y))
}


#' Quadratic form `v' K^-1 v` via the eigen cache
#'
#' Cheap by-product of kron_solve_eigen: we already have w = (U_s' (x) U_t') v
#' in the rotated basis, and `v' K^-1 v = sum(w^2 / lam_full)`.
#'
#' @param v Numeric vector of length n * nt.
#' @param ke Eigen cache from [kron_eigen()].
#' @return Scalar `v' K^-1 v`.
kron_quad <- function(v, ke) {
  n_sites <- nrow(ke$U_s)
  n_times <- nrow(ke$U_t)

  X <- t(matrix(v, nrow = n_times, ncol = n_sites))
  W <- t(ke$U_s) %*% X %*% ke$U_t
  w <- as.vector(t(W))

  sum(w * w / ke$lam_full)
}


#' Cache Cholesky factors of the spatial and temporal kernels
#'
#' Returns lower-triangular `L_s`, `L_t` with `L_s L_s' = K_space`,
#' `L_t L_t' = K_time`, plus the Kronecker log-determinant
#'    log|K| = nt * 2*sum(log(diag(L_s))) + n * 2*sum(log(diag(L_t))).
#'
#' Used inside the theta-slice sampler in place of the full eigen cache.
#' For an n x n matrix Cholesky is ~3-7x cheaper than eigen at the same n
#' (with a typical optimised BLAS), and the slice only needs log|K| and
#' v' K^-1 v -- both available cheaply from the Cholesky factor via
#' [kron_quad_chol()]. The eigen cache is still needed once per sweep
#' for the f-block preconditioner, so this is a *complement* to
#' [kron_eigen()] rather than a replacement.
#'
#' @param space Spatial kernel matrix (n x n, symmetric PD).
#' @param time Temporal kernel matrix (nt x nt, symmetric PD).
#' @return A list with `L_s`, `L_t` (lower-triangular factors) and
#'   `log_det` (scalar log-determinant of the Kronecker product).
kron_chol <- function(space, time) {
  # safe_chol returns lower-triangular L with L L' = A (after a tiny
  # relative jitter on the diagonal if needed for numerical stability).
  L_s <- safe_chol(space)
  L_t <- safe_chol(time)
  list(
    L_s     = L_s,
    L_t     = L_t,
    log_det = nrow(time)  * 2 * sum(log(diag(L_s))) +
              nrow(space) * 2 * sum(log(diag(L_t)))
  )
}


#' Quadratic form `v' K^-1 v` via Kronecker Cholesky factors
#'
#' Given a Cholesky cache from [kron_chol()], computes `v' K^-1 v` via two
#' triangular solves on the reshaped n x nt matrix:
#'
#'   `g = (L_s ⊗ L_t)^-1 v`,    `v' K^-1 v = ||g||^2`.
#'
#' Using the standard reshape identity
#'   `(A ⊗ B) vec(t(X)) = vec(t(A X B^T))`
#' with `A = L_s^-1`, `B = L_t^-1`:
#'   `G1 = L_s^-1 F`    (apply L_s^-1 on the left via forwardsolve)
#'   `G2 = G1 L_t^-T`   (apply L_t^-T on the right, equivalent to
#'                       forwardsolve(L_t, t(G1)) transposed)
#'   `||g||^2 = sum(G2^2)`
#' (we never need to flatten G2; the row-major vec doesn't affect the norm.)
#'
#' Same O(n^2 nt + n nt^2) cost as [kron_quad()] via eigen, but a much
#' smaller leading constant.
#'
#' @param v Numeric vector of length n * nt.
#' @param kc Cholesky cache from [kron_chol()].
#' @return Scalar `v' K^-1 v`.
kron_quad_chol <- function(v, kc) {
  n_sites <- nrow(kc$L_s)
  n_times <- nrow(kc$L_t)

  F_mat <- t(matrix(v, nrow = n_times, ncol = n_sites))   # n x nt
  G1    <- forwardsolve(kc$L_s, F_mat)                    # n x nt
  G2    <- t(forwardsolve(kc$L_t, t(G1)))                 # n x nt
  sum(G2 * G2)
}


#' Add a small ridge to a square matrix for numerical conditioning
#'
#' Adds `lambda * mean(diag(x))` to every diagonal entry. The default 1e-6
#' is a tiny relative jitter that does not perturb the structure of well-
#' conditioned matrices but rescues Cholesky on nearly-singular ones. The
#' multiplier is **relative**: changing the scale of the kernel (e.g. multiplying
#' all entries by 100) does not change how much jitter is added in relative
#' terms.
#'
#' @param x A square numeric matrix.
#' @param lambda Non-negative relative ridge value. Default 1e-6.
#' @return A matrix the same size as `x` with the ridge added to its diagonal.
#' @export
regularise <- function(x, lambda = 1e-6) {
  x + lambda * mean(diag(x)) * diag(nrow(x))
}


#' Fill observed values into a full-length vector
#'
#' Scatter the m-vector x_obs into a length-N vector at positions obs_idx,
#' filling unobserved positions with zero. The transpose of this operation
#' (selection) is just `x_full[obs_idx]`.
#'
#' @param x_obs Numeric vector of length m.
#' @param obs_idx Integer indices of observed positions in the full vector.
#' @param N Total length of the full vector.
#' @return Numeric vector of length N.
fill_vector <- function(x_obs, obs_idx, N) {
  v <- numeric(N)
  v[obs_idx] <- x_obs
  v
}
