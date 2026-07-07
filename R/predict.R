# =============================================================================
# GP prediction by conjugate gradient (CG).
#
# Companion to infer_kernel_params(): given estimated kernel hyperparameters,
# predict the latent rate (and a count prediction interval) by conditioning on
# the OBSERVED cells only -- so missing cells are filled by genuine GP
# interpolation and their predictive uncertainty widens, rather than being
# mean-imputed as in a completed-grid smoother.
#
# Two pieces, both matrix-free (reusing cg()/kron_mv()):
#   * posterior MEAN of the field   -- a single CG solve,
#       f_hat = K S^T (S K S^T + nu I)^{-1} g_obs,
#     so it is exactly smooth and independent of the number of draws.
#   * posterior VARIANCE of the field -- an exact closed-form "no gaps" part
#     plus a missing-data correction estimated from `n_draws` paired
#     perturbation draws (one CG solve each; this is the expensive part).
#
# The latent log-rate posterior (mean Z, variance Vz) is then turned into a
# Negative-Binomial count prediction interval via the law of total variance and
# a lognormal moment-match.
# =============================================================================

#' Posterior variance of the latent field (exact part + control-variate draws)
#'
#' Estimates the per-cell posterior variance of the GP field conditioned on the
#' observed cells, splitting it into an exact term and a small Monte-Carlo
#' correction:
#' \deqn{\operatorname{Var}(f) = V_{\mathrm{complete}} +
#'   \mathbb{E}\!\left[d_{\mathrm{obs}}^2 - d_{\mathrm{complete}}^2\right].}
#'
#' \eqn{V_{\mathrm{complete}} = \operatorname{diag}(K - K(K+\nu I)^{-1}K)} is the
#' posterior variance had *every* cell been observed -- available in closed form
#' through the Kronecker eigendecomposition, no draws needed. Missing cells
#' change the variance only locally, so the draws are spent purely on that
#' correction: each zero-mean perturbation draw \eqn{d_{\mathrm{obs}}}
#' (Papandreou & Yuille 2010, one CG solve) is paired with an exact
#' complete-grid twin \eqn{d_{\mathrm{complete}}} built from the *same* random
#' numbers \eqn{u, e}, so their difference is nearly noise-free away from gaps
#' (a control variate). Modest draw counts therefore give variances that plain
#' Monte-Carlo would need hundreds of draws to match, and cells far from any
#' gap are essentially exact.
#'
#' The draws are independent, so they run through
#' [future.apply::future_lapply()]: serial under the default
#' [future::plan()], parallel when the caller selects a multi-worker plan.
#' `future.seed = TRUE` gives every draw its own pre-generated
#' L'Ecuyer-CMRG stream, so results are reproducible under [set.seed()] and
#' identical for every backend and worker count.
#'
#' @param obs_idx Integer indices of observed cells in the full vector.
#' @param N Total number of cells (`n * nt`).
#' @param space_mat Spatial kernel matrix (variance-scaled).
#' @param time_mat Temporal kernel matrix.
#' @param noise_var Scalar observation-noise variance \eqn{\nu}.
#' @param Rs_chol,Rt_chol Upper Cholesky factors of `space_mat` / `time_mat`.
#' @param n_draws Number of paired draws for the missing-data correction.
#' @param tol CG tolerance for the draw solves.
#' @param progress_bar Optional progress bar (from `make_curve_bar()`); ticked
#'   once per draw. Only effective under a single-worker plan (parallel
#'   workers cannot tick a bar in the calling session).
#'
#' @return Numeric vector of length `N`: the posterior variance of the field.
#' @keywords internal
gp_posterior_var <- function(obs_idx, N, space_mat, time_mat, noise_var,
                             Rs_chol, Rt_chol, n_draws, tol = 1e-3,
                             progress_bar = NULL) {
  n <- nrow(space_mat)
  nt <- nrow(time_mat)
  eig_s <- eig_sym(space_mat)
  eig_t <- eig_sym(time_mat)
  Us <- eig_s$vectors
  Ut <- eig_t$vectors
  UsT <- t(Us)
  sqrt_nv <- sqrt(noise_var)
  d_eig <- outer(eig_t$values, eig_s$values) # nt x n eigenvalues of K
  shrink <- d_eig / (d_eig + noise_var)

  # exact complete-grid posterior variance: diag(K - K (K + nu I)^{-1} K),
  # per cell (time b, site a) = sum_ij Ut[b,j]^2 Us[a,i]^2 * d_ij nu / (d_ij + nu)
  v_complete <- as.vector((Ut^2) %*% (noise_var * shrink) %*% t(Us^2))

  one_draw <- function(k) {
    # zero-mean perturbation draw conditioned on the observed cells:
    #   u ~ N(0, K);  e ~ N(0, nu I);  d = u - K S^T (S K S^T + nu I)^{-1} (Su + Se)
    u <- quick_mvnorm_chol(Rs_chol, Rt_chol)
    e <- stats::rnorm(N, sd = sqrt_nv)
    alpha <- cg(u[obs_idx] + e[obs_idx], obs_idx, N, space_mat, time_mat,
                noise_var, tol = tol)
    d_obs <- u - kron_mv(fill_vector(alpha, obs_idx, N), space_mat, time_mat)

    # its complete-grid twin, exact in the eigenbasis, sharing the same u and e
    w <- matrix(u + e, nrow = nt, ncol = n)
    coef <- crossprod(Ut, w) %*% Us
    d_complete <- u - as.vector(Ut %*% (coef * shrink) %*% UsT)

    if (!is.null(progress_bar)) progress_bar$tick()
    d_obs^2 - d_complete^2
  }

  # future.seed = TRUE: one CMRG stream per draw, so output is reproducible
  # under set.seed() and identical for every plan/worker count.
  # future.stdout = NA: don't sink stdout, so the (sequential-only) progress
  # bar renders live rather than after the last draw.
  contrib <- future.apply::future_lapply(
    seq_len(n_draws), one_draw,
    future.seed = TRUE, future.stdout = NA
  )
  ss_corr <- Reduce(`+`, contrib)
  pmax(v_complete + ss_corr / n_draws, 0)
}


#' Predict the latent rate and a count prediction interval (CG)
#'
#' Given kernel hyperparameters (e.g. from [infer_kernel_params()]), predicts the
#' latent rate \eqn{\lambda = e^{\mu_s + f_{st}}} at every site-by-time cell by
#' conditioning a separable Gaussian process on the observed counts only. The
#' posterior mean is obtained from a single matrix-free CG solve (so it is
#' smooth and deterministic). The posterior variance splits into an exact
#' closed-form "no gaps" part (via the Kronecker eigendecomposition) plus a
#' missing-data correction estimated from `n_draws` paired perturbation draws
#' (a control variate; see [gp_posterior_var()]), so modest draw counts give
#' tight intervals. The latent-rate posterior is then combined with
#' Negative-Binomial observation noise (law of total variance, lognormal
#' moment-match) to give a 95% count prediction interval.
#'
#' Because the fit conditions on the observed set, missing cells are filled by
#' genuine GP interpolation and their prediction interval can widen over gaps --
#' unlike a completed-grid smoother that mean-imputes the gaps.
#'
#' Predictions are made at every `(id, t)` cell present in `obs_data`. To
#' predict at time points with no data anywhere (e.g. a future week), append
#' rows with `NA` counts for those times (every site) and increase `nt`
#' accordingly; the interval widens with distance from the data. Treat
#' extrapolation beyond the observed range with the usual caution.
#'
#' @param obs_data Data frame with `id` (site), `t` (time) and a count column
#'   named by `value` (`NA` where missing). `t` is a numeric time index whose
#'   *differences* encode real elapsed time, so gaps and uneven spacing between
#'   time points are modelled as genuine time distances (use e.g. weeks or days
#'   since a reference). Must use the same `t` encoding as [infer_kernel_params()].
#' @param coordinates Site coordinates (data frame with `id`, `lon`, `lat`).
#' @param hyperparameters A list with elements `length_scale`, `periodic_scale`,
#'   `long_term_scale`, `nugget_ratio` and `sigma2` -- the value returned by
#'   [infer_kernel_params()].
#' @param nt Number of time points.
#' @param period Period of the seasonal cycle, in the same units as `t`.
#' @param n_draws Number of paired posterior draws used to estimate the
#'   missing-data correction to the variance (the prediction interval).
#'   Controls only the interval, not the mean. Use `0` to return the smooth
#'   posterior-mean rate only (one solve, no interval). Because most of the
#'   variance is computed exactly and the draws only estimate the gap
#'   correction, modest values (25--100) already give tight intervals.
#' @param r Negative-Binomial dispersion for the count interval. If `NULL`
#'   (default) it is estimated by method of moments from the observed counts.
#' @param value Name of the count column (default `"y_obs"`).
#' @param standardise Logical; standardise the plug-in field per site (default
#'   `TRUE`), matching [infer_kernel_params()].
#' @param cg_tol Convergence tolerance for the single posterior-mean CG solve
#'   (the deterministic part of the prediction).
#' @param cg_draw_tol Convergence tolerance for the `n_draws` perturbation-draw
#'   CG solves. Deliberately looser than `cg_tol`: the draws only feed a
#'   Monte-Carlo variance whose own relative error is
#'   \eqn{\approx 1/\sqrt{2\,(n_{draws}-1)}} (about 5% at 200 draws), so
#'   solver error below that is wasted work. At the default `1e-3` the
#'   posterior sd typically changes by well under 1% relative to a tight solve,
#'   while the draw loop needs roughly half the CG iterations.
#' @param progress Logical; show a progress bar over the posterior-draw loop
#'   (the expensive part). Defaults to `TRUE`, but the bar is drawn only in an
#'   interactive UTF-8 / truecolor terminal -- it stays silent in scripts,
#'   knitr, logs, and CI. Under a multi-worker [future::plan()] the bar is
#'   suppressed (workers cannot tick it). Set `FALSE` to disable it entirely.
#'
#' @section Parallel execution:
#' The `n_draws` perturbation draws are independent CG solves and run through
#' the future framework. By default (no [future::plan()] set) they run
#' serially, exactly as before. To spread them across CPU cores, set a plan
#' before calling and reset it after:
#'
#' ```r
#' future::plan(future::multisession, workers = 4)
#' pred <- gp_predict(...)
#' future::plan(future::sequential)
#' ```
#'
#' Results are identical for every backend and worker count, and reproducible
#' under [set.seed()] (each draw gets its own pre-generated L'Ecuyer-CMRG
#' stream). Parallelism pays off when `n * n_draws` is large: each worker
#' costs about a second to start and receives the kernel matrices once. For
#' very large site counts you may need to raise
#' `options(future.globals.maxSize = ...)` (the matrices shipped to workers
#' are ~40 MB at 1000 sites x 260 weeks; the default cap is 500 MiB). If R
#' uses a multithreaded BLAS (e.g. OpenBLAS/MKL), cap its threads inside
#' workers to avoid oversubscription; R's shipped BLAS is single-threaded, so
#' by default there is nothing to do. See
#' `vignette("parallel", package = "weave")` for a walkthrough.
#'
#' @return A data frame with one row per site-by-time cell (site-week) and
#'   columns `id`, `t`, `rate`
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
  cg_tol = 1e-6,
  cg_draw_tol = 1e-3,
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
  coord_idx <- match(ids, coordinates$id)
  if (anyNA(coord_idx)) {
    stop(
      "`coordinates` is missing a row for one or more site `id`s in `obs_data`.",
      call. = FALSE
    )
  }
  coordinates <- coordinates[coord_idx, , drop = FALSE]
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
  # Build the temporal kernel on the actual `t` values (sorted), matching
  # infer_kernel_params(): the hyperparameters were estimated on that axis, so
  # the kernel here must use it too. Using the real `t` values means gaps and
  # uneven spacing between time points are modelled as genuine time distances.
  # `t` should therefore be a numeric encoding actual time position (e.g. weeks
  # or days since a reference), and `period` must be in those same units.
  time_mat <- time_kernel(
    times,
    periodic_scale = hp$periodic_scale,
    long_term_scale = hp$long_term_scale,
    period = period
  )
  noise_var <- hp$sigma2 * hp$nugget_ratio
  Rs_chol <- chol(space_mat)
  Rt_chol <- chol(time_mat)

  # --- observed cells in the full (site-major, time-fastest) vector layout ---
  ok <- !is.na(obs_data[[value]])
  obs_grid <- matrix(FALSE, n, nt)
  obs_grid[cbind(
    match(obs_data$id[ok], ids),
    match(obs_data$t[ok], times)
  )] <- TRUE
  obs_idx <- which(as.vector(t(obs_grid)))
  g_obs <- as.vector(t(G))[obs_idx]

  # --- posterior MEAN of the log-rate (single CG solve, deterministic) -------
  alpha0 <- cg(
    g_obs,
    obs_idx,
    N,
    space_mat,
    time_mat,
    noise_var,
    tol = cg_tol
  )
  f_mean <- kron_mv(fill_vector(alpha0, obs_idx, N), space_mat, time_mat)
  Zmat <- row_mean + row_sd * t(matrix(f_mean, nrow = nt, ncol = n)) # n x nt

  out <- data.frame(
    id = rep(ids, each = nt),
    t = rep(times, times = n),
    rate = as.vector(t(exp(Zmat)))
  )

  if (n_draws >= 1) {
    # --- posterior VARIANCE of the log-rate: exact part + draw correction ----
    # The bar only makes sense under a single-worker plan: sequential futures
    # run in this process (the tick closure fires), parallel workers do not.
    show_bar <- isTRUE(progress) && ansi_tty() && future::nbrOfWorkers() == 1L
    pb <- if (show_bar) make_curve_bar(total = n_draws) else NULL
    vstd <- gp_posterior_var(
      obs_idx,
      N,
      space_mat,
      time_mat,
      noise_var,
      Rs_chol,
      Rt_chol,
      n_draws,
      tol = cg_draw_tol,
      progress_bar = pb
    )
    if (show_bar) pb$done()
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
