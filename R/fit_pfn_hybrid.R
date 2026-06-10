# =============================================================================
# fit_pfn_hybrid() -- amortised (theta, mu_s, r) + exact f
#
# The pure PFN samples every component from the network's diagonal-Gaussian
# heads. That works for the scalar / per-site outputs (theta, log_r, mu_s)
# but breaks badly on f, because the per-cell diagonal Gaussian cannot
# represent a strongly correlated 200-dim GP-shaped posterior.
#
# Hybrid: keep the NN's draws for (theta, r, mu_s), but draw f *exactly*
# conditional on those, using the existing Polya-Gamma machinery from
# R/pg_sampler.R. The conditional p(f | theta, mu_s, r, y) is Gaussian after
# PG augmentation, so two-block (omega, f) Gibbs at fixed theta converges in
# tens of iterations -- no slow theta-mixing, no MH, no slice.
#
# Per posterior draw:
#
#   1. sample (theta, mu_s, log_r)   from the PFN's predicted Gaussian heads
#   2. build space / time kernels at this theta; eigen-cache once
#   3. inner Gibbs:
#         for k = 1 .. n_inner:
#           omega | (f, theta, mu_s, r, y)        -- PG draws
#           f     | (omega, theta, mu_s, r, y)    -- exact Gaussian via pg_draw_f
#         return the last f
#
# Output layout mirrors fit_pfn() / fit_bayes() so downstream code
# (posterior_predict, summaries, plots) is class-agnostic.
# =============================================================================


# Internal: one inner (omega, f) Gibbs trajectory at fixed (theta, mu_s, r).
# `f_warm` is the starting f (typically the PFN's f-mean for the dataset, but
# can be zero). Returns the final f after `n_inner` sweeps.
.inner_omega_f <- function(design, theta, r, mu_s, f_warm,
                            n_inner, period, control) {
  space_mat <- space_kernel(coordinates  = design$coords,
                            length_scale = theta$length_scale)
  time_mat  <- time_kernel(times           = seq_len(design$nt),
                           periodic_scale  = theta$periodic_scale,
                           long_term_scale = theta$long_term_scale,
                           period          = period)
  ke <- kron_eigen(space_mat, time_mat)

  state <- list(
    f         = f_warm,
    mu        = mu_s,
    r         = r,
    space_mat = space_mat,
    time_mat  = time_mat,
    ke        = ke
  )
  pcg_iters_trace <- integer(n_inner)
  for (k in seq_len(n_inner)) {
    psi <- state$f[design$obs_idx] + mu_s[design$site_idx_obs]
    state$omega <- pg_draw_omega(psi, design$y_obs, r)
    f_res <- pg_draw_f(state, design, control)
    state$f <- f_res$f
    pcg_iters_trace[k] <- f_res$pcg_iters
  }
  list(f = state$f, pcg_iters = pcg_iters_trace)
}


#' Hybrid PFN inference: amortised (theta, r, mu_s) + exact f
#'
#' For each of `n_post` posterior draws, sample `(theta, r, mu_s)` from the
#' PFN's predicted diagonal-Gaussian heads (which SBC has shown are well-
#' calibrated), then draw `f` from its exact conditional via a short
#' Polya-Gamma Gibbs trajectory at the sampled hyperparameters. The
#' output is the joint posterior `(theta, r, mu_s, f)`, with each component
#' drawn from a correct conditional.
#'
#' Why this works: the PFN's per-component heads are trained against the
#' true simulated parameters, so they learn the marginals `p(theta | y)`,
#' `p(r | y)`, `p(mu_s | y)` regardless of how well the f-head does its
#' (impossible-for-diagonal-Gaussian) job. The PG augmentation gives an
#' exact Gaussian conditional `p(f | theta, mu_s, r, y)`, and the existing
#' [pg_draw_f()] block samples from it. Composing the two yields an
#' amortised posterior with no diagonal-Gaussian-on-f approximation.
#'
#' @param obs_data Tidy data frame; same shape `fit_bayes()` expects.
#' @param weights Path to a PFN checkpoint produced by [pfn_train()].
#' @param n_post Number of posterior draws (default 500).
#' @param n_inner Number of inner Gibbs sweeps per draw (default 20). With
#'   theta held fixed, the (omega, f) chain mixes essentially instantly --
#'   20 is conservative; 10 is usually enough.
#' @param device Torch device.
#' @param seed Optional RNG seed.
#' @param period Periodic kernel period. Defaults to the model spec.
#' @param pcg_tol,pcg_maxit PCG control for the inner f draws.
#' @param verbose Show a progress bar.
#'
#' @return A `weave_pfn_hybrid` S3 object (which also inherits from
#'   `weave_pfn`, so [posterior_predict()] works unchanged) with elements:
#'   - `theta_samples`: matrix `(n_post, 3)` native-scale draws
#'   - `r_samples`:     numeric length `n_post`
#'   - `mu_samples`:    matrix `(n_post, n)`
#'   - `f_samples`:     matrix `(n_post, N)`, time-fastest
#'   - `f_mean`, `f_var`: posterior summaries
#'   - `pcg_iters_mean`: mean inner-PCG iterations (per draw)
#'   - `design`, `n_post`, `n_inner`, `weights_path`
#' @export
fit_pfn_hybrid <- function(obs_data, weights,
                           n_post   = 500L,
                           n_inner  = 20L,
                           device   = "cpu",
                           seed     = NULL,
                           period   = NULL,
                           pcg_tol  = 1e-4,
                           pcg_maxit = 500L,
                           verbose  = TRUE) {
  if (!requireNamespace("torch", quietly = TRUE)) {
    stop("`torch` is required for fit_pfn_hybrid().")
  }
  n_post  <- as.integer(n_post)
  n_inner <- as.integer(n_inner)
  if (!is.finite(n_post)  || n_post  < 1L) stop("`n_post` must be >= 1.")
  if (!is.finite(n_inner) || n_inner < 1L) stop("`n_inner` must be >= 1.")

  design <- build_design(obs_data)
  model  <- pfn_load_model(weights, device = device)
  .check_geometry(design, model)
  if (is.null(period)) period <- model$spec$period

  tt        <- .design_to_tensors(design, device = device)
  posterior <- torch::with_no_grad({ model(tt$y, tt$mask) })

  if (!is.null(seed)) {
    set.seed(seed)
    torch::torch_manual_seed(seed)
  }

  # ---- Pull Gaussian-head parameters into plain R arrays --------------------
  theta_mean_log <- as.numeric(posterior$theta_mean[1, ])
  theta_sd_log   <- as.numeric(posterior$theta_logvar[1, ]$exp()$sqrt())
  mu_mean_v      <- as.numeric(posterior$mu_mean[1, ])
  mu_sd_v        <- as.numeric(posterior$mu_logvar[1, ]$exp()$sqrt())

  # Warm-start f at the PFN's predicted f-mean. This is biased (it's the
  # broken head's mean) but pulls each inner chain into the basin faster
  # than starting at zero would. 20 inner sweeps still washes out the
  # bias; the warm start is just a convergence aid.
  f_mean_mat <- as.matrix(posterior$f_mean[1, , ])      # n x nt
  f_warm     <- as.vector(t(f_mean_mat))                # length N, time-fast

  # ---- Storage -------------------------------------------------------------
  theta_out <- matrix(NA_real_, n_post, 3L,
                      dimnames = list(NULL, c("length_scale",
                                              "periodic_scale",
                                              "long_term_scale")))
  r_out     <- numeric(n_post)
  mu_out    <- matrix(NA_real_, n_post, design$n)
  f_out     <- matrix(NA_real_, n_post, design$N)
  pcg_mean  <- numeric(n_post)

  control <- list(pcg_tol = pcg_tol, pcg_maxit = pcg_maxit)

  pb <- if (verbose) progress::progress_bar$new(
    format = sprintf("  hybrid PFN [:bar] :percent  draw :current/:total"),
    total = n_post, clear = FALSE, width = 70) else NULL

  for (s in seq_len(n_post)) {
    log4 <- stats::rnorm(4L, mean = theta_mean_log, sd = theta_sd_log)
    theta_s <- list(
      length_scale    = exp(log4[1]),
      periodic_scale  = exp(log4[2]),
      long_term_scale = exp(log4[3])
    )
    r_s    <- exp(log4[4])
    mu_s_s <- stats::rnorm(design$n, mean = mu_mean_v, sd = mu_sd_v)

    inner <- .inner_omega_f(design, theta_s, r_s, mu_s_s, f_warm,
                            n_inner = n_inner, period = period,
                            control = control)

    theta_out[s, ] <- c(theta_s$length_scale, theta_s$periodic_scale,
                        theta_s$long_term_scale)
    r_out[s]       <- r_s
    mu_out[s, ]    <- mu_s_s
    f_out[s, ]     <- inner$f
    pcg_mean[s]    <- mean(inner$pcg_iters)

    if (!is.null(pb)) pb$tick()
  }

  # ---- Summaries ----------------------------------------------------------
  f_mean <- colMeans(f_out)
  f_var  <- apply(f_out, 2, stats::var)

  result <- list(
    design          = design,
    n_post          = n_post,
    n_inner         = n_inner,
    weights_path    = normalizePath(weights, mustWork = TRUE),
    theta_samples   = theta_out,
    r_samples       = r_out,
    mu_samples      = mu_out,
    f_samples       = f_out,
    f_mean          = f_mean,
    f_var           = f_var,
    pcg_iters_mean  = pcg_mean,
    posterior       = lapply(posterior, function(t) as.array(t$cpu()))
  )
  # Inherit from weave_pfn so posterior_predict.weave_pfn works unchanged.
  class(result) <- c("weave_pfn_hybrid", "weave_pfn", "list")
  result
}


#' @export
print.weave_pfn_hybrid <- function(x, ...) {
  cat("<weave_pfn_hybrid fit>\n")
  cat(sprintf("  weights     : %s\n", x$weights_path))
  cat(sprintf("  n_post      : %d (with %d inner Gibbs sweeps each)\n",
              x$n_post, x$n_inner))
  cat(sprintf("  sites x time: %d x %d  (N = %d, observed = %d)\n",
              x$design$n, x$design$nt, x$design$N,
              length(x$design$obs_idx)))
  cat(sprintf("  inner PCG   : %.1f iters / sweep (mean)\n",
              mean(x$pcg_iters_mean)))
  invisible(x)
}
