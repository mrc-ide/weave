# =============================================================================
# Fast point estimator: linearised-Gaussian surrogate for the NB-GP model.
#
# This file used to contain the linear-algebra primitives too; those now live
# in R/kron.R and R/pcg.R. What is left here is the public-facing fit()
# function: a fast deterministic estimator that maps observed counts onto
# posterior mean / approximate variance under a *Gaussianised* likelihood,
# good for quick exploration and as a warm-start for fit_bayes().
#
# Model used here:
#
#   y_obs ~ NB / Poisson           (the truth)
#   z = log(1 + y_obs)              (working response, observed cells only)
#   z = mu + f + epsilon            (linearised observation model)
#   f ~ N(0, K_space (x) K_time)
#   epsilon_i ~ N(0, lambda_i / (lambda_i + 1)^2)
#
# The variance of epsilon is the delta-method approximation to Var[log(1+Y)]
# under Poisson(lambda) with lambda hat = exp(mu_i). It is heteroscedastic by
# construction and gives the solver a sensible per-cell weight.
#
# This is the *surrogate* model; the true Poisson/NB likelihood is used by
# fit_bayes(). The surrogate is biased for low counts -- prefer fit_bayes()
# for any uncertainty quantification you actually want to trust.
# =============================================================================


#' Fast point estimator for the spatio-temporal latent field
#'
#' Solves the linearised-Gaussian system via matrix-free PCG with the
#' Kronecker-eigen preconditioner. Computes a stochastic (Hutchinson)
#' estimate of the posterior variance for credible ribbons.
#'
#' @param obs_data A data frame with columns `id`, `t`, `lat`, `lon`, and a
#'   count column (`y_obs` if present, else `n`). NAs in the count column
#'   denote missing observations.
#' @param coordinates **Ignored** -- kept for backward compatibility.
#'   Coordinates are now read from `obs_data`. Will be removed in a future
#'   version.
#' @param hyperparameters Either a named list
#'   `list(length_scale, periodic_scale, long_term_scale)`
#'   or (legacy) a length-3 numeric vector in that order.
#' @param n,nt **Ignored** -- inferred from `obs_data`. Kept for backward
#'   compatibility.
#' @param period Periodic kernel period; default 52 (weekly seasonality).
#' @param n_var_probes Number of Rademacher probes used in the Hutchinson
#'   variance estimator. Default 30.
#' @param pcg_tol Relative residual tolerance for PCG. Default 1e-6.
#' @param pcg_maxit Maximum PCG iterations. Default 500.
#'
#' @return The input `obs_data` augmented with the standard fit columns:
#'   z_est, tausq, z_min, z_max, lambda_est, lambda_min, lambda_max,
#'   pred_Q2.5, pred_Q25, data_Q50, pred_Q75, pred_Q97.5, surprisal, etc.
#'
#' @export
fit <- function(obs_data, coordinates = NULL, hyperparameters,
                n = NULL, nt = NULL, period = 52,
                n_var_probes = 30, pcg_tol = 1e-6, pcg_maxit = 500) {

  # ---- Hyperparameter unpacking (named list preferred; legacy 3-vector OK) --
  if (is.list(hyperparameters)) {
    needed <- c("length_scale", "periodic_scale", "long_term_scale")
    miss   <- setdiff(needed, names(hyperparameters))
    if (length(miss) > 0) stop("hyperparameters missing: ", paste(miss, collapse = ", "))
    length_scale    <- hyperparameters$length_scale
    periodic_scale  <- hyperparameters$periodic_scale
    long_term_scale <- hyperparameters$long_term_scale
  } else if (is.numeric(hyperparameters) && length(hyperparameters) == 3) {
    warning("Passing `hyperparameters` as a numeric vector is deprecated; ",
            "use list(length_scale, periodic_scale, long_term_scale).")
    length_scale    <- hyperparameters[1]
    periodic_scale  <- hyperparameters[2]
    long_term_scale <- hyperparameters[3]
  } else {
    stop("hyperparameters must be a named list or a length-3 numeric vector.")
  }

  # ---- Design ---------------------------------------------------------------
  des <- build_design(obs_data)
  N   <- des$N
  n   <- des$n
  nt  <- des$nt

  # ---- Kernels --------------------------------------------------------------
  time_mat <- time_kernel(
    times          = seq_len(nt),
    periodic_scale = periodic_scale,
    long_term_scale = long_term_scale,
    period         = period
  )
  space_mat <- space_kernel(
    coordinates  = des$coords,
    length_scale = length_scale
  )

  ke <- kron_eigen(space_mat, time_mat)

  # ---- Working response + heteroscedastic nugget ----------------------------
  # Per-site offset from the design (already log-mean of observed counts).
  mu_per_obs <- des$mu_init[des$site_idx_obs]

  # Working response on the log scale.
  y_log <- log1p(des$y_obs)
  y_centred <- y_log - mu_per_obs

  # Delta-method variance of log(1+Y) at lambda_hat = exp(mu_init_per_obs).
  lam_hat   <- exp(mu_per_obs)
  noise_var <- lam_hat / (lam_hat + 1)^2

  # ---- Solve for alpha via PCG ---------------------------------------------
  obs_idx <- des$obs_idx
  Amv_fun <- function(v) Amv(v, obs_idx, N, space_mat, time_mat, noise_var)
  Minv    <- kron_eigen_preconditioner(ke, sigma2 = mean(noise_var),
                                       obs_idx, N)
  alpha_res <- pcg(y_centred, Amv_fun, Minv,
                   tol = pcg_tol, maxit = pcg_maxit)
  alpha <- alpha_res$x

  # ---- Posterior mean of f at every cell -----------------------------------
  f_hat <- kron_mv(fill_vector(alpha, obs_idx, N), space_mat, time_mat)

  # ---- Posterior variance via Hutchinson ------------------------------------
  # We want diag(Sigma_post) where
  #     Sigma_post = K - K S' (S K S' + D)^-1 S K
  # For a Rademacher probe z (in R^N), the unbiased estimator is
  #     diag(Sigma_post) ~= E[ z * (Sigma_post z) ]
  # Each (Sigma_post z) costs one PCG solve.
  tausq <- hutchinson_diag(
    obs_idx     = obs_idx,
    N           = N,
    space_mat   = space_mat,
    time_mat    = time_mat,
    noise_var   = noise_var,
    ke          = ke,
    n_probes    = n_var_probes,
    pcg_tol     = pcg_tol,
    pcg_maxit   = pcg_maxit
  )
  # Tiny numerical floor: the linearised surrogate variance shouldn't be negative
  # but the stochastic estimator can drift slightly below zero on a few cells.
  tausq <- pmax(tausq, 1e-8)

  # ---- Assemble the result --------------------------------------------------
  # Restore the full-length z (NA → estimated, observed → estimated too).
  obs_data <- obs_data |>
    dplyr::arrange(.data$id, .data$t)
  obs_data$z_est <- f_hat + des$mu_init[des$site_idx_full]
  obs_data$tausq <- tausq

  obs_data <- obs_data |>
    dplyr::mutate(
      # Convert log-scale credible intervals to count-scale via lognormal.
      z_min       = .data$z_est - 1.96 * sqrt(.data$tausq),
      z_max       = .data$z_est + 1.96 * sqrt(.data$tausq),
      lambda_est  = stats::qlnorm(0.5,  meanlog = .data$z_est, sdlog = sqrt(.data$tausq)),
      lambda_min  = stats::qlnorm(0.025, meanlog = .data$z_est, sdlog = sqrt(.data$tausq)),
      lambda_max  = stats::qlnorm(0.975, meanlog = .data$z_est, sdlog = sqrt(.data$tausq))
    )

  count_col <- if ("y_obs" %in% names(obs_data)) "y_obs" else "n"

  # Per-site empirical overdispersion (used only to widen ribbons; defensive
  # clamp keeps it positive when var(y) <= mean(y)).
  obs_data <- obs_data |>
    dplyr::mutate(
      mean_lambda    = mean(exp(.data$z_est)),
      var_y          = stats::var(.data[[count_col]], na.rm = TRUE),
      overdispersion = pmax(
        1e-6,
        (.data$mean_lambda^2) / pmax(.data$var_y - .data$mean_lambda, 1e-6)
      ),
      .by = "id"
    ) |>
    dplyr::mutate(
      negbin_size = .data$lambda_est / .data$overdispersion,
      negbin_prob = .data$negbin_size / (.data$negbin_size + .data$lambda_est),
      pred_Q2.5  = stats::qnbinom(0.025, size = .data$negbin_size, prob = .data$negbin_prob),
      pred_Q25   = stats::qnbinom(0.25,  size = .data$negbin_size, prob = .data$negbin_prob),
      data_Q50   = stats::qnbinom(0.5,   size = .data$negbin_size, prob = .data$negbin_prob),
      pred_Q75   = stats::qnbinom(0.75,  size = .data$negbin_size, prob = .data$negbin_prob),
      pred_Q97.5 = stats::qnbinom(0.975, size = .data$negbin_size, prob = .data$negbin_prob)
    ) |>
    dplyr::mutate(
      log_p      = stats::dnbinom(.data[[count_col]], size = .data$negbin_size,
                                  prob = .data$negbin_prob, log = TRUE),
      y_mode     = floor((.data$negbin_size - 1) * (1 - .data$negbin_prob) / .data$negbin_prob),
      log_p_mode = stats::dnbinom(.data$y_mode, size = .data$negbin_size,
                                  prob = .data$negbin_prob, log = TRUE),
      surprisal  = pmax(0, pmin(1, 1 - exp(.data$log_p - .data$log_p_mode)))
    )

  return(obs_data)
}


#' Hutchinson stochastic diagonal of the linearised posterior covariance
#'
#' Estimates diag(Sigma_post) where
#'    Sigma_post = K - K S' (S K S' + D)^-1 S K
#' using `n_probes` Rademacher probes z and the identity
#'    diag(Sigma_post) ~= E[ z .* (Sigma_post z) ]
#' Each probe costs one matrix-free PCG solve.
#'
#' Internal helper to fit(); not exported.
#'
#' @noRd
hutchinson_diag <- function(obs_idx, N, space_mat, time_mat, noise_var,
                            ke, n_probes = 30, pcg_tol = 1e-6,
                            pcg_maxit = 500) {

  Minv    <- kron_eigen_preconditioner(ke, sigma2 = mean(noise_var),
                                       obs_idx, N)
  Amv_fun <- function(v) Amv(v, obs_idx, N, space_mat, time_mat, noise_var)

  accum <- numeric(N)
  for (m in seq_len(n_probes)) {
    z <- sample(c(-1, 1), N, replace = TRUE)

    # K z (full-length)
    Kz <- kron_mv(z, space_mat, time_mat)

    # alpha_z = (S K S' + D)^-1 S K z
    pcg_res <- pcg(Kz[obs_idx], Amv_fun, Minv,
                   tol = pcg_tol, maxit = pcg_maxit)
    # K S' alpha_z
    KSta <- kron_mv(fill_vector(pcg_res$x, obs_idx, N), space_mat, time_mat)

    # Sigma_post z = K z - K S' (...)^-1 S K z
    Spz <- Kz - KSta

    accum <- accum + z * Spz
  }
  accum / n_probes
}
