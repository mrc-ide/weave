# =============================================================================
# Post-fit diagnostics helper.
#
# `diagnose_bayes()` produces a single per-parameter table summarising
# whether each of the four kernel/dispersion parameters has been well
# inferred. Two automatic flags:
#
#   - is_poorly_mixed: Rhat > threshold (default 1.10). Multiple chains
#     are not visiting the same posterior region.
#
#   - is_prior_driven: a Kolmogorov-Smirnov test cannot distinguish the
#     posterior from prior samples. The data has not informed this
#     parameter beyond the prior -- e.g. periodic_scale on data with no
#     real seasonality. This is *not* the same as a sampler bug; it just
#     tells you the posterior interval is just the prior interval.
#
# Both checks rely only on the trace data already in the fit object, the
# prior `sample()` closures (already built by bayes_priors()), and base
# R / coda. No new heavy dependencies.
# =============================================================================


#' Per-parameter diagnostic summary of a `weave_bayes` fit
#'
#' Returns a data frame with one row per kernel hyperparameter plus `r`,
#' summarising posterior quantiles, Gelman-Rubin Rhat, effective sample
#' size, and two automatic flags: `is_poorly_mixed` (chains not merged
#' on this parameter) and `is_prior_driven` (data did not move the
#' posterior away from the prior).
#'
#' @param fit A `weave_bayes` object returned by [fit_bayes()].
#' @param rhat_threshold Rhat above which `is_poorly_mixed = TRUE`.
#'   Default 1.10.
#' @param ks_pvalue_threshold KS test p-value above which we *cannot*
#'   reject "posterior == prior" and `is_prior_driven = TRUE`. Default
#'   0.05.
#' @param prior_samples Number of samples to draw from each prior for the
#'   KS test. Default 5000.
#' @return A data frame with columns `param`, `median`, `q025`, `q975`,
#'   `rhat`, `ess`, `prior_overlap_p`, `is_poorly_mixed`,
#'   `is_prior_driven`, `status`.
#' @export
diagnose_bayes <- function(fit,
                           rhat_threshold      = 1.10,
                           ks_pvalue_threshold = 0.05,
                           prior_samples       = 5000) {
  if (!inherits(fit, "weave_bayes")) {
    stop("`fit` must be a weave_bayes object from fit_bayes().")
  }
  post <- (fit$burnin + 1L):fit$n_sweeps

  # Posterior samples pooled across chains, per parameter.
  posts <- list(
    length_scale    = as.vector(fit$theta_trace[post, "length_scale",    , drop = TRUE]),
    periodic_scale  = as.vector(fit$theta_trace[post, "periodic_scale",  , drop = TRUE]),
    long_term_scale = as.vector(fit$theta_trace[post, "long_term_scale", , drop = TRUE]),
    r               = as.vector(fit$r_trace[post, , drop = TRUE])
  )

  # Prior samples via the existing closures. r's prior is in `fit$priors$r`,
  # theta priors live in `fit$priors$theta[[name]]`.
  priors_smp <- list(
    length_scale    = replicate(prior_samples,
                                fit$priors$theta$length_scale$sample()),
    periodic_scale  = replicate(prior_samples,
                                fit$priors$theta$periodic_scale$sample()),
    long_term_scale = replicate(prior_samples,
                                fit$priors$theta$long_term_scale$sample()),
    r               = replicate(prior_samples, fit$priors$r$sample())
  )

  # Per-chain mcmc lists for Rhat / ESS via coda.
  rhat <- ess <- stats::setNames(rep(NA_real_, 4), names(posts))
  if (fit$n_chains > 1 && requireNamespace("coda", quietly = TRUE)) {
    make_mcmc <- function(name) {
      coda::as.mcmc.list(lapply(seq_len(fit$n_chains), function(ch) {
        values <- if (name == "r") fit$r_trace[post, ch] else
          fit$theta_trace[post, name, ch]
        coda::mcmc(values, start = fit$burnin + 1L)
      }))
    }
    for (nm in names(posts)) {
      mc <- make_mcmc(nm)
      rhat[nm] <- coda::gelman.diag(mc, autoburnin = FALSE)$psrf[1]
      ess[nm]  <- sum(coda::effectiveSize(mc))
    }
  } else if (requireNamespace("coda", quietly = TRUE)) {
    for (nm in names(posts)) {
      mc <- coda::mcmc(posts[[nm]])
      ess[nm] <- as.numeric(coda::effectiveSize(mc))
    }
  }

  # KS test of posterior vs prior. Warning suppression: ks.test() warns on
  # ties, which are common with discrete posteriors; we don't care.
  ks_p <- vapply(names(posts), function(nm) {
    suppressWarnings(
      stats::ks.test(posts[[nm]], priors_smp[[nm]])$p.value
    )
  }, numeric(1))

  qfun <- function(x) stats::quantile(x, c(0.025, 0.5, 0.975))
  rows <- lapply(names(posts), function(nm) {
    q <- qfun(posts[[nm]])
    poorly_mixed <- isTRUE(rhat[nm] > rhat_threshold)
    prior_driven <- isTRUE(ks_p[nm] > ks_pvalue_threshold)
    status <- if (poorly_mixed && prior_driven) "poorly-mixed; prior-driven"
              else if (poorly_mixed) "poorly-mixed"
              else if (prior_driven) "prior-driven"
              else "ok"
    data.frame(
      param            = nm,
      median           = unname(q[2]),
      q025             = unname(q[1]),
      q975             = unname(q[3]),
      rhat             = unname(rhat[nm]),
      ess              = unname(ess[nm]),
      prior_overlap_p  = unname(ks_p[nm]),
      is_poorly_mixed  = poorly_mixed,
      is_prior_driven  = prior_driven,
      status           = status,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  class(out) <- c("weave_diagnosis", "data.frame")
  attr(out, "rhat_threshold")      <- rhat_threshold
  attr(out, "ks_pvalue_threshold") <- ks_pvalue_threshold
  out
}


#' @export
print.weave_diagnosis <- function(x, ...) {
  cat("<weave diagnosis>\n")
  # Single-row formatting for quick scan.
  flag <- ifelse(x$status == "ok", "  ok",
          ifelse(x$status == "poorly-mixed",  "  POORLY MIXED",
          ifelse(x$status == "prior-driven",  "  PRIOR-DRIVEN",
                                              "  POORLY MIXED + PRIOR-DRIVEN")))
  fmt <- function(v, w, n = 2) formatC(v, width = w, digits = n, format = "g")
  for (i in seq_len(nrow(x))) {
    cat(sprintf("  %-16s med=%s [%s, %s]  Rhat=%s  ESS=%s%s\n",
                x$param[i],
                fmt(x$median[i],  8),
                fmt(x$q025[i],    8),
                fmt(x$q975[i],    8),
                fmt(x$rhat[i],    5),
                fmt(x$ess[i],     6, 4),
                flag[i]))
  }
  rt <- attr(x, "rhat_threshold")
  pt <- attr(x, "ks_pvalue_threshold")
  cat(sprintf("\n  thresholds: Rhat > %.2f -> poorly mixed; ", rt))
  cat(sprintf("KS p > %.2f -> prior-driven\n", pt))
  invisible(x)
}
