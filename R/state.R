#' Kronecker diadiagonalg
#'
#' A helper: Kronecker diagonal without allocating a big dense kronecker()
#'
#' @param space_diag Space matrix diagonal
#' @param time_diag Time matrix diagonal
#' @param n Number of sites
#' @param nt Number of times
#'
#' @returns Kronecker diagonal
kdiag_from_factors <- function(space_diag, time_diag, n, nt) {
  # time varies fastest (sites × times)
  rep(space_diag, each = nt) * rep(time_diag, times = n)
}

#' Build state object
#'
#' This convenience wrapper assembles the ingredients the Gaussian process (GP)
#' model needs so downstream functions can work with a single, tidy bundle.
#'
#' It builds the spatial and temporal kernels, indexes the observed counts,
#' prepares working responses and variances for the Poisson log-scale
#' approximation, and returns the matrices and closures required for GP
#' inference.
#'
#' @param obs_data Observation data
#' @param coordinates Site coordinates
#' @param hyperparameters Vector of hyperparameters
#' @param n N sites
#' @param nt N times
#' @param period Period
#'
#' @returns List of pre-computaed GP elements
#' @export
gp_build_state <- function(obs_data, coordinates, hyperparameters, n, nt, period = 52) {

  time_mat  <- time_kernel(
    times = 1:nt,
    periodic_scale = hyperparameters[2],
    long_term_scale = hyperparameters[3],
    period = period
  )

  space_mat <- space_kernel(
    coordinates = coordinates,
    length_scale = hyperparameters[1]
  )

  obs_idx   <- which(!is.na(obs_data$y_obs))
  N         <- n * nt


  y_work    <- obs_data$f_infer[obs_idx]
  lam_hat   <- exp(obs_data$mu_infer[obs_idx])
  if(length(hyperparameters) == 3){
    # Adjustment for fitting Poisson on the log scale:
    ## approximate count noise so the GP learns the underlying signal correctly.
    noise_var <- lam_hat / (lam_hat + 1)^2
  }
  if(length(hyperparameters) == 4){
    # For NB assumption (if implemented)
    # NB working likelihood on the log-scale (size = k):
    noise_var <- (lam_hat + lam_hat^2 / hyperparameters[4]) / (lam_hat + 1)^2
  }
  # (optional: stabilise)
  noise_var <- pmax(noise_var, 1e-8)

  # Preconditioner diag(K) = diag(space) ⊗ diag(time)
  kdiag_full <- kdiag_from_factors(diag(space_mat), diag(time_mat), n, nt)

  # Optional: precompute Choleskys for fast prior draws
  Rt <- chol(time_mat)   # nt × nt
  Rs <- chol(space_mat)  #  n ×  n

  # A-solve closure for (S K S^T + D) y = rhs
  A_solve <- function(rhs_obs, tol = 1e-6) {
    pcg(rhs_obs, obs_idx, N, space_mat, time_mat, noise_var, kdiag_full, tol = tol)
  }

  list(
    time_mat = time_mat,
    space_mat = space_mat,
    Rt = Rt,
    Rs = Rs,
    obs_idx = obs_idx,
    y_work = y_work,
    noise_var = noise_var,
    kdiag_full = kdiag_full,
    A_solve = A_solve,
    mu_infer = obs_data$mu_infer,
    id = obs_data$id,
    t = obs_data$t,
    N = N,
    n = n,
    nt = nt,
    hyperparameters = hyperparameters
  )
}
