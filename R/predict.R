# =============================================================================
# GP prediction by preconditioned conjugate gradient (PCG).
#
# Companion to infer_kernel_params(): given estimated kernel hyperparameters,
# predict the latent rate (and a count prediction interval) by conditioning on
# the OBSERVED cells only -- so missing cells are filled by genuine GP
# interpolation and their predictive uncertainty widens, rather than being
# mean-imputed as in a completed-grid smoother.
#
# Two pieces, both matrix-free (reusing pcg()/kron_mv()):
#   * posterior MEAN of the field   -- a single PCG solve,
#       f_hat = K S^T (S K S^T + nu I)^{-1} g_obs,
#     so it is exactly smooth and independent of the number of draws.
#   * posterior VARIANCE of the field -- estimated from `n_draws` perturbation
#     draws (one PCG solve each; this is the expensive part).
#
# The latent log-rate posterior (mean Z, variance Vz) is then turned into a
# Negative-Binomial count prediction interval via the law of total variance and
# a lognormal moment-match.
# =============================================================================

# -----------------------------------------------------------------------------
# One exact posterior draw of the (standardised) latent field, conditioning on
# the observed cells. Perturbation sampler (Papandreou & Yuille 2010):
#   u ~ N(0, K);  e ~ N(0, nu);  solve (S K S^T + nu I) a = (S u + e) - g_obs;
#   f = u - K S^T a.
# -----------------------------------------------------------------------------
gp_field_draw <- function(
  g_obs,
  obs_idx,
  N,
  space_mat,
  time_mat,
  noise_var,
  kdiag,
  Rs_chol,
  Rt_chol,
  tol = 1e-6
) {
  u <- quick_mvnorm_chol(Rs_chol, Rt_chol)
  eps <- stats::rnorm(length(obs_idx), sd = sqrt(noise_var))
  alpha <- pcg(
    (u[obs_idx] + eps) - g_obs,
    obs_idx,
    N,
    space_mat,
    time_mat,
    noise_var,
    kdiag,
    tol = tol
  )
  u - kron_mv(fill_vector(alpha, obs_idx, N), space_mat, time_mat)
}


#' Predict the latent rate and a count prediction interval (PCG)
#'
#' Given kernel hyperparameters (e.g. from [infer_kernel_params()]), predicts the
#' latent rate \eqn{\lambda = e^{\mu_s + f_{st}}} at every site-by-time cell by
#' conditioning a separable Gaussian process on the observed counts only. The
#' posterior mean is obtained from a single matrix-free PCG solve (so it is
#' smooth and deterministic); the posterior variance is estimated from
#' `n_draws` perturbation draws. The latent-rate posterior is then combined with
#' Negative-Binomial observation noise (law of total variance, lognormal
#' moment-match) to give a 95% count prediction interval.
#'
#' Because the fit conditions on the observed set, missing cells are filled by
#' genuine GP interpolation and their prediction interval can widen over gaps --
#' unlike a completed-grid smoother that mean-imputes the gaps.
#'
#' @param obs_data Data frame with `id` (site), `t` (time) and a count column
#'   named by `value` (`NA` where missing).
#' @param coordinates Site coordinates (data frame with `id`, `lon`, `lat`).
#' @param hyperparameters A list with elements `length_scale`, `periodic_scale`,
#'   `long_term_scale`, `nugget_ratio` and `sigma2` -- the value returned by
#'   [infer_kernel_params()].
#' @param nt Number of time points.
#' @param period Period of the seasonal cycle.
#' @param n_draws Number of posterior draws used to estimate the variance
#'   (the prediction interval). Controls only the interval, not the mean. Use
#'   `0` to return the smooth posterior-mean rate only (one solve, no interval).
#' @param r Negative-Binomial dispersion for the count interval. If `NULL`
#'   (default) it is estimated by method of moments from the observed counts.
#' @param value Name of the count column (default `"y_obs"`).
#' @param standardise Logical; standardise the plug-in field per site (default
#'   `TRUE`), matching [infer_kernel_params()].
#' @param pcg_tol Convergence tolerance for the PCG solves.
#' @param progress Logical; show a progress bar over the posterior-draw loop
#'   (the expensive part). Defaults to `TRUE`, but the bar is drawn only in an
#'   interactive UTF-8 / truecolor terminal -- it stays silent in scripts,
#'   knitr, logs and CI. Set `FALSE` to disable it entirely.
#'
#' @return A data frame with one row per cell and columns `id`, `t`, `rate`
#'   (posterior point estimate of \eqn{\lambda}), and -- when `n_draws >= 1` --
#'   `lower` and `upper` (the 95% count prediction interval). The dispersion `r`
#'   used and `n_draws` are attached as attributes.
#'
#' @seealso [infer_kernel_params()]
#' @export
gp_predict <- function(
  obs_data,
  coordinates,
  hyperparameters,
  nt,
  period,
  n_draws = 100,
  r = NULL,
  value = "y_obs",
  standardise = TRUE,
  pcg_tol = 1e-6,
  progress = TRUE
) {
  hp <- hyperparameters
  need <- c(
    "length_scale",
    "periodic_scale",
    "long_term_scale",
    "nugget_ratio",
    "sigma2"
  )
  if (!all(need %in% names(hp))) {
    stop(
      "`hyperparameters` must contain: ",
      paste(need, collapse = ", "),
      call. = FALSE
    )
  }

  ids <- sort(unique(obs_data$id))
  times <- sort(unique(obs_data$t))
  n <- length(ids)
  if (length(times) != nt) {
    stop(
      "`nt` does not match the number of unique times in `obs_data`.",
      call. = FALSE
    )
  }
  coordinates <- coordinates[match(ids, coordinates$id), , drop = FALSE]
  N <- n * nt

  # --- plug-in field with per-site centring (and optional scaling) -----------
  M <- matrix(NA_real_, n, nt)
  M[cbind(
    match(obs_data$id, ids),
    match(obs_data$t, times)
  )] <- log1p(obs_data[[value]])
  row_mean <- rowMeans(M, na.rm = TRUE)
  row_mean[!is.finite(row_mean)] <- 0
  Mc <- M - row_mean
  if (standardise) {
    row_sd <- apply(Mc, 1, stats::sd, na.rm = TRUE)
    row_sd[!is.finite(row_sd) | row_sd == 0] <- 1
  } else {
    row_sd <- rep(1, n)
  }
  G <- Mc / row_sd
  G[is.na(G)] <- 0

  # --- separable kernel (sigma^2 folded into space) + scalar nugget ----------
  space_mat <- hp$sigma2 *
    space_kernel(coordinates, length_scale = hp$length_scale)
  time_mat <- time_kernel(
    times,
    periodic_scale = hp$periodic_scale,
    long_term_scale = hp$long_term_scale,
    period = period
  )
  noise_var <- hp$sigma2 * hp$nugget_ratio
  Rs_chol <- chol(space_mat)
  Rt_chol <- chol(time_mat)
  kdiag <- kdiag_from_factors(diag(space_mat), diag(time_mat), n, nt)

  # --- observed cells in the full (site-major, time-fastest) vector layout ---
  ok <- !is.na(obs_data[[value]])
  obs_grid <- matrix(FALSE, n, nt)
  obs_grid[cbind(
    match(obs_data$id[ok], ids),
    match(obs_data$t[ok], times)
  )] <- TRUE
  obs_idx <- which(as.vector(t(obs_grid)))
  g_obs <- as.vector(t(G))[obs_idx]

  # --- posterior MEAN of the log-rate (single PCG solve, deterministic) ------
  alpha0 <- pcg(
    g_obs,
    obs_idx,
    N,
    space_mat,
    time_mat,
    noise_var,
    kdiag,
    tol = pcg_tol
  )
  f_mean <- kron_mv(fill_vector(alpha0, obs_idx, N), space_mat, time_mat)
  Zmat <- row_mean + row_sd * t(matrix(f_mean, nrow = nt, ncol = n)) # n x nt

  out <- data.frame(
    id = rep(ids, each = nt),
    t = rep(times, times = n),
    rate = as.vector(t(exp(Zmat)))
  )

  if (n_draws >= 1) {
    # --- posterior VARIANCE of the log-rate from draws -----------------------
    Fd <- matrix(NA_real_, nrow = N, ncol = n_draws)
    show_bar <- isTRUE(progress) && interactive()
    if (show_bar) pb <- make_curve_bar(total = n_draws)
    for (i in seq_len(n_draws)) {
      Fd[, i] <- gp_field_draw(
        g_obs,
        obs_idx,
        N,
        space_mat,
        time_mat,
        noise_var,
        kdiag,
        Rs_chol,
        Rt_chol,
        tol = pcg_tol
      )
      if (show_bar) pb$tick()
    }
    if (show_bar) pb$done()
    vstd <- apply(Fd, 1, stats::var)
    Vzmat <- (row_sd^2) * t(matrix(vstd, nrow = nt, ncol = n))

    # --- dispersion r (method of moments) if not supplied --------------------
    if (is.null(r)) {
      lam_obs <- exp(Zmat)[cbind(
        match(obs_data$id[ok], ids),
        match(obs_data$t[ok], times)
      )]
      yk <- obs_data[[value]][ok]
      den <- sum((yk - lam_obs)^2 - lam_obs)
      r <- if (den > 0) max(sum(lam_obs^2) / den, 0.1) else Inf
    }

    # --- count prediction interval (law of total variance + lognormal) -------
    Elam <- exp(Zmat + Vzmat / 2)
    Elam2 <- exp(2 * Zmat + 2 * Vzmat)
    Vlam <- exp(2 * Zmat + Vzmat) * (exp(Vzmat) - 1)
    nb_var <- if (is.finite(r)) Elam2 / r else 0
    pred_var <- Elam + nb_var + Vlam
    mean_safe <- pmax(Elam, 1e-8)
    ss <- log(1 + pred_var / mean_safe^2)
    mln <- log(mean_safe) - ss / 2

    out$lower <- as.vector(t(stats::qlnorm(0.025, mln, sqrt(ss))))
    out$upper <- as.vector(t(stats::qlnorm(0.975, mln, sqrt(ss))))
    attr(out, "r") <- r
  }
  attr(out, "n_draws") <- n_draws
  out
}
