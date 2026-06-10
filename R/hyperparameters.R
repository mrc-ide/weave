# =============================================================================
# Quick, exact kernel-hyperparameter estimation.
#
# Goal: a fast, deterministic estimate of the separable-GP kernel
# hyperparameters that is more defensible than the cross-validated /
# working-Gaussian approach in fit_hyperparameters.R, but does NOT require an
# MCMC or the matrix-free PCG sampler. Some pragmatic approximation is fine.
#
# Idea: form a cheap *plug-in* latent field g directly from the counts
# (per-site-centred log1p(y), see build_plugin_field), then fit the GP by
# maximising the EXACT marginal likelihood of g under
#
#     g  ~  N(0,  sigma^2 * ( R_space(theta) (x) R_time(theta) + eta * I ) ),
#
# where R_space and R_time are *correlation* kernels (unit diagonal), (x) is the
# Kronecker product, sigma^2 is a global variance and eta is a noise-to-signal
# ratio (a nugget). The nugget is the key defensible improvement over a naive
# correlogram / plug-in fit: it lets the model attribute per-cell observation
# noise (and the log1p attenuation) to eta instead of forcing the length scales
# to shrink. sigma^2 is profiled out analytically.
#
# Everything is O(n^3 + nt^3) via the Kronecker eigendecomposition -- we never
# form, invert, or factorise the (n*nt)-square covariance. Because the noise is
# a scalar multiple of the identity, K = sigma^2 ((R_s (x) R_t) + eta I) shares
# the eigenvectors U_s (x) U_t of the kernel, so the nugget only shifts the
# eigenvalues:  eigenvalue_ij = sigma^2 * (a_i * b_j + eta).
#
# Caveat (honest): this conditions on a noisy plug-in field rather than
# integrating over the latent field, so it is an approximation to full Bayesian
# hyperparameter inference. The nugget mitigates the resulting attenuation but
# does not remove it entirely. It is intended as a fast initialiser / standalone
# estimate, not a posterior.
#
# Vector layout throughout: sites x times with time varying fastest, reshaped to
# an n_sites x n_times matrix with  F <- t(matrix(g, nrow = nt, ncol = n)).
# =============================================================================


#' Default priors for the kernel hyperparameters
#'
#' Weakly-informative log-normal priors (i.e. Normal priors on the log scale) on
#' `length_scale`, `periodic_scale`, `long_term_scale` and the noise-to-signal
#' ratio `nugget_ratio`. They act as mild regularisation on an otherwise
#' maximum-likelihood fit, keeping weakly-identified parameters (notably
#' `long_term_scale`) away from the boundary.
#'
#' @return A named list of `list(meanlog, sdlog)` priors.
#' @export
default_kernel_priors <- function() {
  list(
    length_scale    = list(meanlog = log(1),    sdlog = 2),
    periodic_scale  = list(meanlog = log(1),    sdlog = 2),
    long_term_scale = list(meanlog = log(100),  sdlog = 2),
    nugget_ratio    = list(meanlog = log(0.1),  sdlog = 2)
  )
}

log_prior_kernel <- function(log_pars, priors) {
  stats::dnorm(log_pars[1], priors$length_scale$meanlog,    priors$length_scale$sdlog,    log = TRUE) +
    stats::dnorm(log_pars[2], priors$periodic_scale$meanlog,  priors$periodic_scale$sdlog,  log = TRUE) +
    stats::dnorm(log_pars[3], priors$long_term_scale$meanlog, priors$long_term_scale$sdlog, log = TRUE) +
    stats::dnorm(log_pars[4], priors$nugget_ratio$meanlog,    priors$nugget_ratio$sdlog,    log = TRUE)
}


# -----------------------------------------------------------------------------
# Symmetric eigendecomposition with a floor on the eigenvalues (the kernels from
# space_kernel()/time_kernel() already carry a tiny diagonal nugget, so this
# only guards against negative eigenvalues from numerical noise).
# -----------------------------------------------------------------------------
eig_sym <- function(K, floor = 1e-12) {
  e <- eigen(K, symmetric = TRUE)
  e$values <- pmax(e$values, floor)
  e
}


#' Build a plug-in latent field from observed counts
#'
#' Forms a cheap estimate of the latent log-intensity field as the per-site
#' centred (and optionally scaled) `log1p` of the observed counts. Missing cells
#' are mean-imputed (0 after centring) and therefore contribute nothing to the
#' marginal likelihood.
#'
#' Per-site centring removes the site intercept `mu_s`; per-site scaling
#' homogenises per-site variances so a single global `sigma^2` and the
#' *correlation* kernels apply.
#'
#' @param obs_data Data frame with `id` (site), `t` (time) and the count column
#'   named by `value`.
#' @param n Number of sites.
#' @param nt Number of time points.
#' @param value Name of the count column (default `"y_obs"`).
#' @param standardise Logical; scale each site to unit variance after centring
#'   (default `TRUE`).
#'
#' @return A numeric vector of length `n * nt`, ordered sites x times (time
#'   fastest).
#' @export
build_plugin_field <- function(obs_data, n, nt, value = "y_obs", standardise = TRUE) {
  ids   <- sort(unique(obs_data$id))
  times <- sort(unique(obs_data$t))
  if (length(ids) != n || length(times) != nt) {
    stop("`n`/`nt` do not match the unique ids/times in `obs_data`.", call. = FALSE)
  }

  M <- matrix(NA_real_, nrow = n, ncol = nt)
  M[cbind(match(obs_data$id, ids), match(obs_data$t, times))] <- log1p(obs_data[[value]])

  row_mean <- rowMeans(M, na.rm = TRUE)
  row_mean[!is.finite(row_mean)] <- 0
  M <- M - row_mean
  if (standardise) {
    row_sd <- apply(M, 1, stats::sd, na.rm = TRUE)
    row_sd[!is.finite(row_sd) | row_sd == 0] <- 1
    M <- M / row_sd
  }
  M[is.na(M)] <- 0

  as.vector(t(M))   # sites x times, time fastest
}


#' Exact separable-GP marginal log-likelihood of a field, with a nugget
#'
#' Evaluates the exact log-density of `g ~ N(0, sigma^2 ((R_s (x) R_t) + eta I))`
#' given the eigendecompositions of the spatial and temporal correlation
#' kernels, via the Kronecker log-determinant and quadratic-form identities. The
#' global variance `sigma^2` is profiled out (concentrated log-likelihood) and
#' the profiled value is attached as `attr(., "sigma2")`.
#'
#' @param g Plug-in field, length `n * nt`, ordered sites x times (time fastest).
#' @param n Number of sites.
#' @param nt Number of time points.
#' @param eig_s Eigendecomposition (`eigen` object) of the spatial correlation
#'   kernel.
#' @param eig_t Eigendecomposition (`eigen` object) of the temporal correlation
#'   kernel.
#' @param eta Noise-to-signal ratio (nugget), a non-negative scalar.
#'
#' @return The concentrated log-likelihood (numeric scalar), with the profiled
#'   `sigma2` attached as an attribute.
#' @export
gp_marginal_loglik <- function(g, n, nt, eig_s, eig_t, eta) {
  F_mat <- t(matrix(g, nrow = nt, ncol = n))               # n x nt, time fastest
  G <- crossprod(eig_s$vectors, F_mat) %*% eig_t$vectors   # U_s^T F U_t
  d <- outer(eig_s$values, eig_t$values) + eta             # a_i * b_j + eta

  quad_unit <- sum(G^2 / d)                                # g^T (R + eta I)^{-1} g
  log_det   <- sum(log(d))                                 # log|R + eta I|

  N <- n * nt
  sigma2 <- quad_unit / N
  ll <- -0.5 * (N * log(2 * pi) + log_det + N * log(sigma2) + N)
  attr(ll, "sigma2") <- sigma2
  ll
}


# -----------------------------------------------------------------------------
# Exact log-posterior of (theta, eta) given the plug-in field g, parameterised
# on the log scale as (log length_scale, log periodic_scale, log long_term_scale,
# log eta). Builds the spatial/temporal correlation kernels with the existing
# space_kernel()/time_kernel() builders and eigendecomposes them.
# -----------------------------------------------------------------------------
kernel_log_posterior <- function(log_pars, g, n, nt, coordinates, times, period,
                                 priors) {
  K_space <- space_kernel(coordinates, length_scale = exp(log_pars[1]))
  K_time  <- time_kernel(times, periodic_scale = exp(log_pars[2]),
                         long_term_scale = exp(log_pars[3]), period = period)

  eig_s <- eig_sym(K_space)
  eig_t <- eig_sym(K_time)

  ll <- gp_marginal_loglik(g, n, nt, eig_s, eig_t, eta = exp(log_pars[4]))
  lp <- as.numeric(ll) + log_prior_kernel(log_pars, priors)
  attr(lp, "sigma2") <- attr(ll, "sigma2")
  lp
}


#' Quick exact-marginal-likelihood estimate of the kernel hyperparameters
#'
#' Estimates `(length_scale, periodic_scale, long_term_scale)` plus a
#' noise-to-signal ratio by maximising the exact GP marginal likelihood of a
#' plug-in latent field, using the Kronecker eigendecomposition so the full
#' `(n * nt)`-square covariance is never formed. Fast and deterministic -- no
#' MCMC, no iterative solver.
#'
#' This is the recommended replacement for the cross-validated [fit()] estimator
#' when a quick hyperparameter estimate is wanted (e.g. as a starting point for
#' a downstream sampler, or as a standalone summary).
#'
#' @param obs_data Data frame with `id` (site), `t` (time) and the count column
#'   named by `value`.
#' @param coordinates Site coordinates (data frame with `lon`, `lat`), ordered
#'   to match `sort(unique(obs_data$id))`.
#' @param nt Number of time points.
#' @param period Period of the seasonal cycle.
#' @param value Name of the count column (default `"y_obs"`).
#' @param standardise Logical; standardise the plug-in field per site (default
#'   `TRUE`).
#' @param priors Log-normal priors, see [default_kernel_priors()].
#' @param start Named/length-4 starting values on the natural scale
#'   (`length_scale`, `periodic_scale`, `long_term_scale`, `nugget_ratio`).
#'
#' @return A list with `length_scale`, `periodic_scale`, `long_term_scale`,
#'   `nugget_ratio`, the profiled `sigma2`, the maximised `log_posterior`, and
#'   `convergence` (the `optim` code; 0 = success).
#' @export
infer_kernel_params <- function(obs_data, coordinates, nt, period,
                                value = "y_obs", standardise = TRUE,
                                priors = default_kernel_priors(),
                                start = c(length_scale = 1, periodic_scale = 1,
                                          long_term_scale = 100, nugget_ratio = 0.1)) {
  n <- length(unique(obs_data$id))
  coordinates <- coordinates[match(sort(unique(obs_data$id)), coordinates$id), , drop = FALSE]

  g <- build_plugin_field(obs_data, n, nt, value = value, standardise = standardise)
  times <- seq_len(nt)
  log_start <- log(as.numeric(start))

  neglp <- function(log_pars) {
    -as.numeric(kernel_log_posterior(log_pars, g, n, nt, coordinates, times,
                                     period, priors))
  }
  opt <- stats::optim(log_start, neglp, method = "Nelder-Mead",
                      control = list(reltol = 1e-9, maxit = 2000))
  lp <- kernel_log_posterior(opt$par, g, n, nt, coordinates, times, period, priors)

  list(
    length_scale    = exp(opt$par[1]),
    periodic_scale  = exp(opt$par[2]),
    long_term_scale = exp(opt$par[3]),
    nugget_ratio    = exp(opt$par[4]),
    sigma2          = attr(lp, "sigma2"),
    log_posterior   = as.numeric(lp),
    convergence     = opt$convergence
  )
}
