# =============================================================================
# fit_pfn() -- PFN inference, the drop-in alternative to fit_bayes()
#
# Trade-off vs fit_bayes():
#   * 1 forward pass + n_post Gaussian draws, no MCMC.
#   * No chains, no Rhat -- the network is deterministic given the data.
#   * Calibration is empirical (SBC), not asymptotic. Validate against
#     fit_bayes on small geometries before trusting at scale.
#
# Output shape mirrors weave_bayes objects so downstream code (plots,
# posterior_predict, summaries) can treat both interchangeably.
# =============================================================================


#' Load a trained PFN model from a checkpoint
#'
#' Rebuilds the architecture from the spec saved by [pfn_train()] and
#' restores trained weights. The model is moved to `device` and put in
#' `eval()` mode (so any future dropout/batchnorm would behave correctly --
#' the current architecture has none, but this is the right default).
#'
#' @param weights Path to a checkpoint `.pt` produced by [pfn_train()].
#' @param device Torch device string.
#' @return An `nn_module` ready for inference.
#' @export
pfn_load_model <- function(weights, device = "cpu") {
  if (!requireNamespace("torch", quietly = TRUE)) {
    stop("`torch` is required for pfn_load_model().")
  }
  if (!file.exists(weights)) {
    stop("Checkpoint not found: ", weights)
  }
  ckpt <- torch::torch_load(weights)
  if (is.null(ckpt$spec)) {
    stop("Checkpoint has no `spec` field. It was produced by an older ",
         "version of pfn_train() that didn't persist the model ",
         "construction args. Re-train with the current version.")
  }
  # torch_save() strips the `data.frame` class on round-trip; restore it
  # before passing back to nn_pfn(), which expects a data frame for coords.
  if (!is.data.frame(ckpt$spec$coords)) {
    ckpt$spec$coords <- as.data.frame(ckpt$spec$coords)
  }
  model <- do.call(nn_pfn, ckpt$spec)
  model$load_state_dict(ckpt$model)
  model$to(device = device)
  model$eval()
  model
}


# Internal: pack a build_design() output into (1, n, nt) tensors matching
# what the model was trained on. y_obs lives at obs_idx in a length-N
# (time-fastest) vector; zero-fill missing cells and build a corresponding
# mask.
.design_to_tensors <- function(design, device = "cpu") {
  n  <- design$n
  nt <- design$nt
  N  <- design$N

  y_flat <- numeric(N)               # 0 at missing
  y_flat[design$obs_idx] <- design$y_obs
  mask_flat <- numeric(N)
  mask_flat[design$obs_idx] <- 1

  # Time-fastest flat -> (n, nt) matrix with mat[i, j] = cell (site i, time j).
  y_mat    <- t(matrix(y_flat,    nrow = nt, ncol = n))
  mask_mat <- t(matrix(mask_flat, nrow = nt, ncol = n))

  list(
    y    = torch::torch_tensor(array(y_mat,    dim = c(1L, n, nt)),
                               dtype = torch::torch_float(), device = device),
    mask = torch::torch_tensor(array(mask_mat, dim = c(1L, n, nt)),
                               dtype = torch::torch_float(), device = device)
  )
}


# Internal: assert that the design's coords match the model's geometry.
# The model has baked in the spatial distance matrix as a buffer; running
# inference on a different geometry would silently return nonsense, so we
# refuse outright.
.check_geometry <- function(design, model, tol = 1e-6) {
  spec_coords <- model$spec$coords
  if (nrow(spec_coords) != nrow(design$coords)) {
    stop("Geometry mismatch: model trained on ", nrow(spec_coords),
         " sites, obs_data has ", nrow(design$coords),
         ". A PFN is tied to a fixed geometry; retrain or use a different ",
         "checkpoint.")
  }
  d_spec   <- as.matrix(dist(spec_coords[, c("lon", "lat")]))
  d_design <- as.matrix(dist(design$coords[, c("lon", "lat")]))
  if (max(abs(d_spec - d_design)) > tol) {
    stop("Geometry mismatch: pairwise distances between the model's ",
         "trained coords and obs_data's coords differ by up to ",
         signif(max(abs(d_spec - d_design)), 3), ". ",
         "A PFN is tied to a fixed geometry; retrain or supply matching ",
         "coords.")
  }
  invisible(TRUE)
}


#' PFN inference: full Bayesian posterior in one forward pass
#'
#' Drop-in alternative to [fit_bayes()]: same `obs_data` input, an
#' analogous `weave_pfn` object out, the same `posterior_predict()`
#' method, and -- if the checkpoint is well-trained and calibrated --
#' a posterior that matches what `fit_bayes()` would have produced after
#' thousands of sweeps.
#'
#' @param obs_data Tidy data frame with columns `id`, `t`, `lat`, `lon`
#'   and a count column (`y_obs` or `n`). Same shape `fit_bayes()` expects.
#' @param weights Path to a PFN checkpoint produced by [pfn_train()].
#' @param n_post Number of posterior draws to return (default 1000). Cheap
#'   -- each draw is a Gaussian sample from the predicted posterior.
#' @param device Torch device (default `"cpu"`).
#' @param seed Optional RNG seed for the posterior draws (set both R and
#'   torch seeds for reproducibility).
#'
#' @return A `weave_pfn` S3 object with elements mirroring the relevant
#'   fields of `weave_bayes`:
#'   - `design`, `n_post`
#'   - `theta_samples`: matrix `(n_post, 3)` natural-scale draws
#'   - `r_samples`:     numeric `n_post`
#'   - `mu_samples`:    matrix `(n_post, n)`
#'   - `f_samples`:     matrix `(n_post, N)`, time-fastest order
#'   - `f_mean`, `f_var`: per-cell posterior mean / variance (length N)
#'   - `posterior`: list with the raw model output (for diagnostics)
#' @export
fit_pfn <- function(obs_data, weights,
                    n_post = 1000L,
                    device = "cpu",
                    seed   = NULL) {
  if (!requireNamespace("torch", quietly = TRUE)) {
    stop("`torch` is required for fit_pfn().")
  }
  n_post <- as.integer(n_post)
  if (!is.finite(n_post) || n_post < 1L) {
    stop("`n_post` must be a positive integer; got ", n_post, ".")
  }

  design <- build_design(obs_data)

  model <- pfn_load_model(weights, device = device)
  .check_geometry(design, model)

  tt <- .design_to_tensors(design, device = device)

  # Single forward pass; no grad needed.
  posterior <- torch::with_no_grad({ model(tt$y, tt$mask) })

  if (!is.null(seed)) {
    set.seed(seed)
    torch::torch_manual_seed(seed)
  }

  # ---- Sample from each Gaussian head ---------------------------------------
  # theta head outputs log(length, periodic, long_term, r); the first three
  # are theta, the fourth is log r.
  theta_mean_log <- as.numeric(posterior$theta_mean[1, ])
  theta_sd_log   <- as.numeric(posterior$theta_logvar[1, ]$exp()$sqrt())
  log4 <- matrix(stats::rnorm(n_post * 4L,
                              mean = rep(theta_mean_log, each = n_post),
                              sd   = rep(theta_sd_log,   each = n_post)),
                 nrow = n_post, ncol = 4L)
  colnames(log4) <- c("length_scale", "periodic_scale",
                      "long_term_scale", "log_r")
  theta_samples <- exp(log4[, 1:3, drop = FALSE])
  r_samples     <- exp(log4[, 4])

  # mu head: (n,) means + log-vars.
  mu_mean_v <- as.numeric(posterior$mu_mean[1, ])
  mu_sd_v   <- as.numeric(posterior$mu_logvar[1, ]$exp()$sqrt())
  mu_samples <- matrix(stats::rnorm(n_post * design$n,
                                    mean = rep(mu_mean_v, each = n_post),
                                    sd   = rep(mu_sd_v,   each = n_post)),
                       nrow = n_post, ncol = design$n)

  # f head: (n, nt) means + log-vars per dataset. Flatten with time fastest
  # to match fit_bayes()'s `f_samples` convention.
  f_mean_mat <- as.matrix(posterior$f_mean[1, , ])    # n x nt
  f_sd_mat   <- as.matrix(posterior$f_logvar[1, , ]$exp()$sqrt())
  f_mean_flat <- as.vector(t(f_mean_mat))             # length N, time fastest
  f_sd_flat   <- as.vector(t(f_sd_mat))
  N <- design$N
  f_samples <- matrix(stats::rnorm(n_post * N,
                                   mean = rep(f_mean_flat, each = n_post),
                                   sd   = rep(f_sd_flat,   each = n_post)),
                      nrow = n_post, ncol = N)

  result <- list(
    design        = design,
    n_post        = n_post,
    weights_path  = normalizePath(weights, mustWork = TRUE),
    theta_samples = theta_samples,
    r_samples     = r_samples,
    mu_samples    = mu_samples,
    f_samples     = f_samples,
    f_mean        = f_mean_flat,
    f_var         = f_sd_flat^2,
    posterior     = lapply(posterior, function(t) as.array(t$cpu()))
  )
  class(result) <- c("weave_pfn", "list")
  result
}


#' @export
print.weave_pfn <- function(x, ...) {
  cat("<weave_pfn fit>\n")
  cat(sprintf("  weights     : %s\n", x$weights_path))
  cat(sprintf("  n_post      : %d\n", x$n_post))
  cat(sprintf("  sites x time: %d x %d  (N = %d, observed = %d)\n",
              x$design$n, x$design$nt, x$design$N,
              length(x$design$obs_idx)))
  invisible(x)
}


#' @export
summary.weave_pfn <- function(object, ...) {
  quant <- function(v) stats::quantile(v, c(0.025, 0.5, 0.975))
  out <- rbind(
    length_scale    = quant(object$theta_samples[, "length_scale"]),
    periodic_scale  = quant(object$theta_samples[, "periodic_scale"]),
    long_term_scale = quant(object$theta_samples[, "long_term_scale"]),
    r               = quant(object$r_samples)
  )
  colnames(out) <- c("2.5%", "50%", "97.5%")
  print(out)
  invisible(out)
}


#' Posterior predictive draws from a weave_pfn fit
#'
#' Returns an `(n_post, N)` matrix of `y_rep` draws, one per cell in the
#' full site x time grid. Same column layout as `posterior_predict` on a
#' `weave_bayes` object (R/posterior.R), so downstream plotting and
#' coverage checks work without branching on the fit class.
#'
#' @param object A `weave_pfn` fit.
#' @param ... Ignored.
#' @return Numeric matrix `(n_post, N)` of NB(r, mu) draws with time
#'   varying fastest within site.
#' @export
posterior_predict.weave_pfn <- function(object, ...) {
  # lambda^(s)_{i,t} = exp(f^(s)_{i,t} + mu^(s)_i)
  # mu_samples is (n_post, n); broadcast across the time axis via
  # site_idx_full which maps each cell to its site.
  site_idx <- object$design$site_idx_full
  mu_full  <- object$mu_samples[, site_idx, drop = FALSE]     # (n_post, N)
  lam      <- exp(object$f_samples + mu_full)
  # NB sampling per cell with per-sample r.
  matrix(stats::rnbinom(length(lam),
                        size = rep(object$r_samples, ncol(lam)),
                        mu   = as.vector(lam)),
         nrow = nrow(lam), ncol = ncol(lam))
}
