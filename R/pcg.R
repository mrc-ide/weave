# =============================================================================
# Preconditioned conjugate gradient (PCG) solver and the "observed system"
# matrix-vector multiply that both fit() and fit_bayes() use.
#
# This is the only linear-system solver in the package. The observed-data
# covariance system is
#
#     A x = b      with      A = S K S' + D
#
# where K = K_space (x) K_time is the GP covariance, S selects the observed
# cells out of the n*nt full grid, and D is a diagonal noise / nugget term
# (per-observation; in the PG sampler D = diag(1/omega)).
#
# We never form A. Each PCG iteration only needs a routine that computes A v
# and a preconditioner that approximates A^{-1}. Both are passed in as
# closures so the same solver serves the deterministic fit() (with the
# Kron-eigen preconditioner) and the PG sampler's f-block.
# =============================================================================


#' Observed-system matrix-vector multiply: (S K S' + diag(noise_var)) v
#'
#' Scatter v from observed length m into the full length-N grid, apply K
#' via the reshape trick, select observed entries back out, then add the
#' diagonal nugget.
#'
#' @param v Numeric vector of length m (observed entries).
#' @param obs_idx Integer indices of observed entries in the full grid.
#' @param N Total length of the full grid (n * nt).
#' @param space_mat Spatial kernel matrix.
#' @param time_mat Temporal kernel matrix.
#' @param noise_var Scalar or length-m numeric diagonal nugget added to S K S' v.
#' @return Numeric vector of length m equal to (S K S' + diag(noise_var)) v.
Amv <- function(v, obs_idx, N, space_mat, time_mat, noise_var) {
  kron_mv(fill_vector(v, obs_idx, N), space_mat, time_mat)[obs_idx] +
    noise_var * v
}


#' Generic preconditioned conjugate gradient solver
#'
#' Solves A x = b given closures for A v and for the preconditioner `M^-1 v`.
#' Stops on relative residual or maxit (and only warns when convergence was
#' not actually achieved at exit).
#'
#' @param b Right-hand side (numeric vector).
#' @param Amv_fun Closure: function(v) returning A v.
#' @param Minv_fun Closure: function(v) returning `M^-1 v` (the preconditioner
#'   applied to v). Pass `identity` to fall back to plain CG.
#' @param tol Relative residual tolerance, default 1e-8.
#' @param maxit Maximum number of iterations, default 500.
#' @param verbose If TRUE, prints the residual every 25 iterations.
#' @return A list with `x` (solution), `iters` (iterations used), `converged`
#'   (logical), `rel_resid` (final relative residual).
pcg <- function(b, Amv_fun, Minv_fun = identity,
                tol = 1e-8, maxit = 500, verbose = FALSE) {
  bnorm <- sqrt(sum(b * b))
  if (bnorm == 0) {
    return(list(x = numeric(length(b)), iters = 0L,
                converged = TRUE, rel_resid = 0))
  }

  x <- numeric(length(b))
  r <- b - Amv_fun(x)
  z <- Minv_fun(r)
  p <- z
  rz_old <- sum(r * z)

  converged <- FALSE
  rel_resid <- sqrt(sum(r * r)) / bnorm

  for (it in seq_len(maxit)) {
    Ap     <- Amv_fun(p)
    pAp    <- sum(p * Ap)
    if (pAp <= 0) {
      # Loss of positive-definiteness -- A or M is misbehaving. Bail with
      # the current iterate rather than dividing by zero.
      warning("pcg: non-positive curvature encountered at iter ", it,
              "; returning current iterate.")
      break
    }
    alpha  <- rz_old / pAp
    x      <- x + alpha * p
    r      <- r - alpha * Ap
    rel_resid <- sqrt(sum(r * r)) / bnorm

    if (verbose && it %% 25 == 0) {
      message(sprintf("pcg: iter %d rel_resid %.3e", it, rel_resid))
    }

    if (rel_resid <= tol) {
      converged <- TRUE
      break
    }

    z      <- Minv_fun(r)
    rz_new <- sum(r * z)
    beta   <- rz_new / rz_old
    p      <- z + beta * p
    rz_old <- rz_new
  }

  if (!converged) {
    warning(sprintf(
      "pcg: failed to converge in %d iters (rel_resid %.3e, tol %.1e)",
      maxit, rel_resid, tol
    ))
  }

  list(x = x, iters = it, converged = converged, rel_resid = rel_resid)
}


#' Build the Kron-eigen preconditioner closure for the observed system
#'
#' The preconditioner approximates A = S K S' + D by replacing the
#' heteroscedastic D with a scalar nugget s2 (typically mean(diag(D))) and
#' inverting (K + s2 I) exactly via the eigen cache, restricted to the
#' observed positions. This is the right thing to use whenever K dominates
#' D; when D dominates K, falls back gracefully toward Jacobi.
#'
#' @param ke Eigen cache from [kron_eigen()].
#' @param sigma2 Scalar nugget proxy for the diagonal of D.
#' @param obs_idx Integer indices of observed cells in the full grid.
#' @param N Total length of the full grid.
#' @return A closure of one argument suitable as Minv_fun in [pcg()].
kron_eigen_preconditioner <- function(ke, sigma2, obs_idx, N) {
  function(v) {
    kron_solve_eigen(fill_vector(v, obs_idx, N), ke, sigma2 = sigma2)[obs_idx]
  }
}


#' Build a Jacobi (diagonal) preconditioner closure
#'
#' Used as a fallback when the Kron-eigen preconditioner underperforms (for
#' instance when omega is highly heteroscedastic and the scalar s2 proxy is
#' a poor approximation). Cheap, robust, but slow to converge.
#'
#' @param kdiag_full Vector diag(K) of length N.
#' @param obs_idx Integer indices of observed entries in the full vector.
#' @param noise_var Scalar or length-m numeric nugget on the observed scale.
#' @return A closure of one argument suitable as Minv_fun in [pcg()].
jacobi_preconditioner <- function(kdiag_full, obs_idx, noise_var) {
  diagA <- kdiag_full[obs_idx] + noise_var
  function(v) v / diagA
}
