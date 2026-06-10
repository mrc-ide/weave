# =============================================================================
# Simulation-Based Calibration (SBC) for the PFN
#
# SBC is the canonical sanity check for amortised Bayesian inference: if the
# network's predicted posterior is well-calibrated, then for each parameter
# the QUANTILE of the truth under the predicted posterior CDF should be
# Uniform(0, 1) across many fresh synthetic datasets.
#
# Because the current PFN heads are diagonal Gaussians, we can compute these
# quantiles analytically -- no posterior sampling needed:
#
#     q_i = Phi( (target_i - mean_i) / sd_i )
#
# Under perfect calibration:
#   - Histogram of q_i is flat
#   - Kolmogorov-Smirnov test against Uniform(0, 1) does not reject
#
# Common failure modes and what they look like:
#   - U-shape    : posterior too narrow (overconfident)
#   - Inverted-U : posterior too wide (underconfident)
#   - Skew       : posterior mean is biased
# =============================================================================


#' Run simulation-based calibration on a trained PFN
#'
#' Simulates `n_datasets` fresh synthetic datasets using the same prior as
#' training, runs the model forward on each, and reports the rank quantile
#' of the truth under the predicted posterior for every scalar component
#' the model outputs.
#'
#' For Gaussian heads (the current architecture) the rank quantile is the
#' standard-normal CDF of the residual / sd; we compute it analytically.
#' KS p-values < 0.01 are flagged as miscalibrated.
#'
#' This is fast: simulating 200 datasets at the walkthrough scale takes a
#' few seconds, and the model forward over the whole batch is cheap. Run it
#' after every training run.
#'
#' @param model A trained `nn_module` from [nn_pfn()] (or
#'   [pfn_load_model()]).
#' @param n_datasets Number of fresh datasets to simulate (default 200).
#' @param coords,nt,period Geometry / time grid; defaults pulled from
#'   `model$spec` so the network is tested on the same grid it was trained
#'   on (which is what SBC requires).
#' @param priors Optional [pfn_priors()] (defaults to fresh priors
#'   constructed from `coords` / `nt`).
#' @param device Torch device.
#' @param seed Optional RNG seed.
#'
#' @return A `pfn_sbc` S3 object with components:
#'   - `rank_theta`: matrix `(n_datasets, 4)` of rank quantiles for
#'     `(length_scale, periodic_scale, long_term_scale, log_r)`
#'   - `rank_mu`:    matrix `(n_datasets, n)` for site offsets
#'   - `rank_f`:     array  `(n_datasets, n, nt)` for the latent field
#'   - `ks`: named list of KS p-values against Uniform(0, 1)
#'   - `batch`: the simulated batch (for downstream inspection)
#' @export
pfn_sbc <- function(model,
                    n_datasets = 200L,
                    coords     = NULL,
                    nt         = NULL,
                    period     = NULL,
                    priors     = NULL,
                    device     = "cpu",
                    seed       = NULL) {
  if (!requireNamespace("torch", quietly = TRUE)) {
    stop("`torch` is required for pfn_sbc().")
  }

  if (is.null(coords)) coords <- model$spec$coords
  if (is.null(nt))     nt     <- model$spec$nt
  if (is.null(period)) period <- model$spec$period
  if (!is.data.frame(coords)) coords <- as.data.frame(coords)

  if (!is.null(seed)) {
    set.seed(seed)
    torch::torch_manual_seed(seed)
  }

  batch <- pfn_simulate_batch(coords, nt = nt, n_datasets = n_datasets,
                              priors = priors, period = period,
                              verbose = FALSE)
  tt <- pfn_to_tensors(batch, device = device)

  model$eval()
  on.exit(model$train(), add = TRUE)
  out <- torch::with_no_grad({ model(tt$y, tt$mask) })

  # Pull predicted Gaussian parameters into R arrays on CPU.
  theta_mean <- as.matrix(out$theta_mean$cpu())          # (B, 4)
  theta_sd   <- as.matrix(out$theta_logvar$cpu()$exp()$sqrt())
  mu_mean    <- as.matrix(out$mu_mean$cpu())              # (B, n)
  mu_sd      <- as.matrix(out$mu_logvar$cpu()$exp()$sqrt())
  f_mean     <- as.array (out$f_mean$cpu())               # (B, n, nt)
  f_sd       <- as.array (out$f_logvar$cpu()$exp()$sqrt())

  target_log4 <- cbind(batch$log_theta, log_r = batch$log_r)

  rank_theta <- stats::pnorm((target_log4 - theta_mean) / theta_sd)
  colnames(rank_theta) <- c("length_scale", "periodic_scale",
                            "long_term_scale", "log_r")
  rank_mu    <- stats::pnorm((batch$mu_s   - mu_mean)    / mu_sd)
  rank_f     <- stats::pnorm((batch$f      - f_mean)     / f_sd)

  # KS tests. Suppress the "ties" warning -- with float ranks there are
  # essentially no ties and the warning is noise.
  ks_p <- function(v) suppressWarnings(stats::ks.test(v, "punif")$p.value)
  ks <- list(
    length_scale    = ks_p(rank_theta[, "length_scale"]),
    periodic_scale  = ks_p(rank_theta[, "periodic_scale"]),
    long_term_scale = ks_p(rank_theta[, "long_term_scale"]),
    log_r           = ks_p(rank_theta[, "log_r"]),
    mu_s            = ks_p(as.vector(rank_mu)),
    f               = ks_p(as.vector(rank_f))
  )

  result <- list(
    rank_theta = rank_theta,
    rank_mu    = rank_mu,
    rank_f     = rank_f,
    ks         = ks,
    n_datasets = n_datasets,
    batch      = batch
  )
  class(result) <- c("pfn_sbc", "list")
  result
}


#' @export
print.pfn_sbc <- function(x, ...) {
  cat(sprintf("<pfn_sbc: %d synthetic datasets>\n", x$n_datasets))
  cat("KS p-value vs Uniform(0, 1) (lower => miscalibrated):\n")
  for (nm in names(x$ks)) {
    p   <- x$ks[[nm]]
    tag <- if (is.na(p)) "" else if (p < 0.01) "  *** miscalibrated"
           else if (p < 0.05) "  *  borderline" else ""
    cat(sprintf("  %-18s  %.3f%s\n", nm, p, tag))
  }
  invisible(x)
}


#' @export
plot.pfn_sbc <- function(x, ...) {
  oldpar <- graphics::par(mfrow = c(2, 3), mar = c(4, 4, 3, 1))
  on.exit(graphics::par(oldpar))

  hist_one <- function(v, main) {
    graphics::hist(v, breaks = 20, freq = FALSE, col = "grey80",
                   border = "white", xlim = c(0, 1),
                   main = main, xlab = "rank quantile")
    graphics::abline(h = 1, col = "red", lwd = 2, lty = 2)
  }
  hist_one(x$rank_theta[, "length_scale"],
           sprintf("length_scale (p=%.3f)",    x$ks$length_scale))
  hist_one(x$rank_theta[, "periodic_scale"],
           sprintf("periodic_scale (p=%.3f)",  x$ks$periodic_scale))
  hist_one(x$rank_theta[, "long_term_scale"],
           sprintf("long_term_scale (p=%.3f)", x$ks$long_term_scale))
  hist_one(x$rank_theta[, "log_r"],
           sprintf("log_r (p=%.3f)",           x$ks$log_r))
  hist_one(as.vector(x$rank_mu),
           sprintf("mu_s pooled (p=%.3f)",     x$ks$mu_s))
  hist_one(as.vector(x$rank_f),
           sprintf("f pooled (p=%.3f)",        x$ks$f))
  invisible(x)
}
