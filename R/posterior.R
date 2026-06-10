# =============================================================================
# Posterior-sample storage and post-processing utilities for fit_bayes().
#
# The latent field f has length n * nt -- for n=1000, nt=150 that is 150k
# numbers per sweep. Storing every sweep is a memory bomb (5000 sweeps =
# ~6 GB at 8 bytes/double), so by default we track only the Welford running
# mean and variance of f and keep a thinned sample of full f vectors for
# posterior-predictive plotting.
#
# The Welford recursion is numerically stable and O(N) per update:
#     n_k     = n_{k-1} + 1
#     delta   = x - mean_{k-1}
#     mean_k  = mean_{k-1} + delta / n_k
#     M2_k    = M2_{k-1} + delta * (x - mean_k)
#     var_k   = M2_k / (n_k - 1)
# =============================================================================


welford_new <- function(N) {
  list(n = 0L, mean = numeric(N), M2 = numeric(N))
}

welford_update <- function(state, x) {
  state$n    <- state$n + 1L
  delta      <- x - state$mean
  state$mean <- state$mean + delta / state$n
  state$M2   <- state$M2 + delta * (x - state$mean)
  state
}

welford_var <- function(state) {
  if (state$n < 2L) rep(NA_real_, length(state$mean)) else state$M2 / (state$n - 1L)
}


# Pairwise pooling of two Welford summaries:
#   n_AB    = n_A + n_B
#   mean_AB = (n_A mean_A + n_B mean_B) / n_AB
#   M2_AB   = M2_A + M2_B + delta^2 * n_A * n_B / n_AB,   delta = mean_B - mean_A
# (Chan, Golub & LeVeque 1979; Welford parallel-pooling identity.)
welford_pair <- function(A, B) {
  if (A$n == 0L) return(B)
  if (B$n == 0L) return(A)
  n_tot <- A$n + B$n
  delta <- B$mean - A$mean
  mean_tot <- A$mean + delta * B$n / n_tot
  M2_tot   <- A$M2 + B$M2 + delta^2 * A$n * B$n / n_tot
  list(n = n_tot, mean = mean_tot, M2 = M2_tot)
}

# Reduce a list of Welford summaries to a single pooled one.
welford_combine <- function(summaries) {
  Reduce(welford_pair, summaries)
}


# -----------------------------------------------------------------------------
# Posterior predictive draws
#
# Given a weave_bayes fit and a vector of cells, returns posterior draws of
# y_rep at those cells: for each retained MCMC sample s and each cell c,
# sample y_rep ~ NB(size = r^{(s)}, mu = exp(f^{(s)}_c + mu^{(s)}_{site(c)})).
# -----------------------------------------------------------------------------
#' Posterior predictive draws
#'
#' S3 generic. Dispatches on the fit class -- `weave_bayes` for `fit_bayes()`
#' output, `weave_pfn` for `fit_pfn()` output. Both methods return the same
#' shape so callers don't have to branch on the sampler.
#'
#' @param object A fit object (`weave_bayes` or `weave_pfn`).
#' @param ... Method-specific arguments (e.g. `cells` for `weave_bayes`).
#' @return A matrix with one row per retained posterior sample and one column
#'   per cell.
#' @export
posterior_predict <- function(object, ...) UseMethod("posterior_predict")


#' @rdname posterior_predict
#' @param cells Integer vector of cell indices (positions in the full
#'   times-vary-fastest grid). Default is all cells.
#' @export
posterior_predict.weave_bayes <- function(object, cells = NULL, ...) {
  fit <- object
  if (is.null(fit$f_samples)) {
    stop("This fit was run with store_f = 'summary'; posterior predict needs",
         " store_f = 'thin' or 'all'. Re-fit with a different store_f setting.")
  }
  S <- nrow(fit$f_samples)
  if (is.null(cells)) cells <- seq_len(ncol(fit$f_samples))
  out <- matrix(NA_real_, nrow = S, ncol = length(cells))
  for (s in seq_len(S)) {
    psi <- fit$f_samples[s, cells] +
      fit$mu_samples[s, fit$design$site_idx_full[cells]]
    out[s, ] <- stats::rnbinom(length(cells),
                               size = fit$r_samples[s],
                               mu   = exp(psi))
  }
  out
}
