# =============================================================================
# Quick, exact kernel-hyperparameter estimation.
#
# Goal: a fast, deterministic estimate of the separable-GP kernel
# hyperparameters that is more defensible than the cross-validated /
# working-Gaussian approach in fit_hyperparameters.R, but does NOT require an
# MCMC or the matrix-free CG sampler. Some pragmatic approximation is fine.
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


# -----------------------------------------------------------------------------
# One M-step: maximise the exact marginal posterior of a COMPLETE plug-in field
# `g` over (theta, eta) by Nelder-Mead, returning the estimate list. Factored out
# so the initial fit and each refinement re-fit share one code path.
# -----------------------------------------------------------------------------
fit_kernel_field <- function(g, n, nt, coordinates, times, period, priors, log_start) {
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


# -----------------------------------------------------------------------------
# Refinement E-step: complete the plug-in field by replacing the MISSING cells
# with the GP posterior (conditional) mean given the observed cells under the
# current hyperparameters `hp`. Observed cells keep their plug-in values. This is
# exactly the posterior-mean solve gp_predict() does, reusing the same matrix-free
# CG machinery (cg()/kron_mv()), so it costs one CG solve -- not a dense
# (n*nt)-square factorisation. sigma2 cancels in the conditional mean, so the
# fill is robust to the profiled-variance estimate.
# -----------------------------------------------------------------------------
complete_field_cond_mean <- function(obs_data, coordinates, n, nt, period,
                                     value, standardise, hp) {
  N     <- n * nt
  ids   <- sort(unique(obs_data$id))
  times <- sort(unique(obs_data$t))

  # plug-in field: observed cells carry their standardised values, gaps are 0
  g_vec <- build_plugin_field(obs_data, n, nt, value = value,
                              standardise = standardise)

  ok <- !is.na(obs_data[[value]])
  obs_grid <- matrix(FALSE, n, nt)
  obs_grid[cbind(match(obs_data$id[ok], ids), match(obs_data$t[ok], times))] <- TRUE
  obs_idx  <- which(as.vector(t(obs_grid)))
  miss_idx <- setdiff(seq_len(N), obs_idx)
  if (length(miss_idx) == 0L) return(g_vec)

  space_mat <- hp$sigma2 * space_kernel(coordinates, length_scale = hp$length_scale)
  time_mat  <- time_kernel(times, periodic_scale = hp$periodic_scale,
                           long_term_scale = hp$long_term_scale, period = period)
  noise_var <- hp$sigma2 * hp$nugget_ratio

  alpha  <- cg(g_vec[obs_idx], obs_idx, N, space_mat, time_mat, noise_var)
  f_mean <- kron_mv(fill_vector(alpha, obs_idx, N), space_mat, time_mat)

  g_vec[miss_idx] <- f_mean[miss_idx]
  g_vec
}


#' Quick exact-marginal-likelihood estimate of the kernel hyperparameters
#'
#' Estimates `(length_scale, periodic_scale, long_term_scale)` plus a
#' noise-to-signal ratio by maximising the exact GP marginal likelihood of a
#' plug-in latent field, using the Kronecker eigendecomposition so the full
#' `(n * nt)`-square covariance is never formed. Fast and deterministic -- no
#' MCMC, no iterative solver.
#'
#' This is the recommended quick estimator when a fast hyperparameter estimate is
#' wanted (e.g. as a starting point for a downstream sampler, or as a standalone
#' summary).
#'
#' Set `refine` to enable an EM-style refinement that removes the bias missing
#' cells introduce. Each pass refits after replacing the gaps with the GP
#' posterior (conditional) mean under the current estimate -- a correlation-aware
#' fill, not the flat mean-imputation -- using the same matrix-free CG solve as
#' [gp_predict()]. The expensive observed-cell solve runs only once per pass (not
#' inside the optimiser), so it stays cheap, and it typically converges in 2-3
#' passes to the estimate you would get with no missing data at all. It does not
#' remove the intrinsic plug-in attenuation (conditioning on a noisy field rather
#' than integrating the latent field out), only the part caused by the gaps.
#'
#' @param obs_data Data frame with `id` (site), `t` (time) and the count column
#'   named by `value`. `t` is a numeric time index whose *differences* encode
#'   real elapsed time, so gaps and uneven spacing between time points are
#'   modelled as genuine time distances (use e.g. weeks or days since a
#'   reference). [gp_predict()] must be given the same `t` encoding.
#' @param coordinates Site coordinates (data frame with `lon`, `lat`), ordered
#'   to match `sort(unique(obs_data$id))`.
#' @param nt Number of time points.
#' @param period Period of the seasonal cycle, in the same units as `t`.
#' @param value Name of the count column (default `"y_obs"`).
#' @param standardise Logical; standardise the plug-in field per site (default
#'   `TRUE`).
#' @param priors Log-normal priors, see [default_kernel_priors()].
#' @param start Named/length-4 starting values on the natural scale
#'   (`length_scale`, `periodic_scale`, `long_term_scale`, `nugget_ratio`).
#' @param n_sites Optional integer. If supplied and smaller than the number of
#'   sites, the hyperparameters are estimated from a random subsample of this
#'   many sites. The kernel hyperparameters are shared, population-level
#'   quantities, so a representative site subsample estimates the same length
#'   scales at a fraction of the \eqn{O(n^3)} cost -- useful for very large site
#'   counts. Default `NULL` uses all sites. The subsample is drawn from the
#'   current RNG state, so set a seed beforehand (e.g. [set.seed()]) for a
#'   reproducible estimate. Note: this subsamples *sites* only, not time points
#'   (the temporal kernel needs the full series to resolve the periodic and
#'   long-term scales).
#' @param refine Logical; if `TRUE`, run `refine_iter` EM-style refinement passes
#'   that re-fit after filling the gaps with the GP conditional mean (see
#'   Details). Default `FALSE` (the fast single-pass estimate). Recommended when
#'   missingness is non-trivial.
#' @param refine_iter Number of refinement passes when `refine = TRUE` (default
#'   `3`). Ignored when `refine = FALSE`.
#'
#' @return A list with `length_scale`, `periodic_scale`, `long_term_scale`,
#'   `nugget_ratio`, the profiled `sigma2`, the maximised `log_posterior`, and
#'   `convergence` (the `optim` code; 0 = success).
#' @export
infer_kernel_params <- function(obs_data, coordinates, nt, period,
                                value = "y_obs", standardise = TRUE,
                                priors = default_kernel_priors(),
                                start = c(length_scale = 1, periodic_scale = 1,
                                          long_term_scale = 100, nugget_ratio = 0.1),
                                n_sites = NULL,
                                refine = FALSE,
                                refine_iter = 3L) {
  if (!is.numeric(refine_iter) || length(refine_iter) != 1 || refine_iter < 0) {
    stop("`refine_iter` must be a single non-negative integer.", call. = FALSE)
  }
  refine_iter <- as.integer(refine_iter)
  if (!is.null(n_sites)) {
    site_ids <- sort(unique(obs_data$id))
    if (n_sites < length(site_ids)) {
      keep <- sample(site_ids, n_sites)
      obs_data <- obs_data[obs_data$id %in% keep, , drop = FALSE]
    }
  }

  n <- length(unique(obs_data$id))
  coord_idx <- match(sort(unique(obs_data$id)), coordinates$id)
  if (anyNA(coord_idx)) {
    stop(
      "`coordinates` has no row for every site `id` in `obs_data`.",
      call. = FALSE
    )
  }
  coordinates <- coordinates[coord_idx, , drop = FALSE]

  g <- build_plugin_field(obs_data, n, nt, value = value, standardise = standardise)
  # Temporal kernel axis: the actual `t` values (sorted), so gaps and uneven
  # spacing between time points become genuine time distances rather than being
  # collapsed to a unit-spaced index. This matches the column order of the
  # plug-in field `g` (build_plugin_field() also sorts on unique `t`) and the
  # axis gp_predict() uses, so the estimated hyperparameters transfer correctly.
  times <- sort(unique(obs_data$t))

  est <- fit_kernel_field(g, n, nt, coordinates, times, period, priors,
                          log(as.numeric(start)))

  # EM-style refinement: refit on a grid whose gaps are filled with the GP
  # conditional mean under the current estimate, warm-starting each pass.
  if (isTRUE(refine)) {
    for (i in seq_len(refine_iter)) {
      g <- complete_field_cond_mean(obs_data, coordinates, n, nt, period,
                                    value, standardise, hp = est)
      warm <- log(c(est$length_scale, est$periodic_scale,
                    est$long_term_scale, est$nugget_ratio))
      est <- fit_kernel_field(g, n, nt, coordinates, times, period, priors, warm)
    }
  }

  est
}
