# =============================================================================
# PFN training-data simulator
#
# Draws synthetic datasets from the NB-GP generative model that fit_bayes()
# targets, for training a Prior-Fitted Network (see implementation/
# pfn_walkthrough.R for the end-to-end pipeline).
#
# Why this lives in the package (whereas implementation/simulation.R does
# not): the PFN cannot train without this simulator -- it's part of the
# package's inference toolchain, not an example. The implementation/
# simulator is kept separate because it ships demonstration scaffolding
# (plots, infer-columns) that the package itself doesn't need.
#
# Convention: f is stored with TIME varying fastest within site, matching
# kron_mv() and quick_mvnorm() (R/kron.R, R/sample.R). For an n x nt matrix
# Z with Z[i, j] = GP value at site i, time j:
#   flat   = as.vector(t(Z))            # length n*nt
#   Z_back = t(matrix(flat, nrow=nt))   # round-trip
# =============================================================================


#' Hyperprior over (theta, r, mu_s) for PFN training-data simulation
#'
#' Mirrors [bayes_priors()] for `theta` and `r` so the PFN trains under the
#' same prior shapes that `fit_bayes()` uses at inference time. The `mu_s`
#' prior has an extra hierarchy: each synthetic dataset draws a per-dataset
#' offset `m0 ~ Uniform(m0_range)` (since at training time we have no data
#' to anchor `m0` against, unlike `bayes_priors()`), then `mu_s ~ N(m0, mu_sd)`.
#'
#' `m0_range` defaults cover ~5 to ~500 counts/cell, which spans the routine
#' HF metrics this package is built for. Widen it if your metric range is
#' larger.
#'
#' @param coords Data frame with columns `lat`, `lon` (the fixed geometry).
#' @param nt Number of timepoints in the (fixed) grid.
#' @param m0_range Length-2 numeric, range for the per-dataset log-mean
#'   offset `m0` (default `c(log(5), log(500))`).
#' @param mu_sd Standard deviation of `mu_s | m0` (default 1.5).
#'
#' @return A list with element `$sample` containing zero-argument samplers
#'   for `length_scale`, `periodic_scale`, `long_term_scale`, `r`, `m0`,
#'   plus the scalar `mu_sd` for drawing `mu_s | m0` and metadata
#'   (`median_dist`, `nt`, `m0_range`, `mu_sd`).
#' @export
pfn_priors <- function(coords, nt,
                       m0_range = c(log(5), log(500)),
                       mu_sd    = 1.5) {
  d <- get_spatial_distance(coords[, c("lon", "lat")])
  diag(d) <- NA
  median_dist <- stats::median(d, na.rm = TRUE)
  if (!is.finite(median_dist) || median_dist <= 0) {
    stop("pfn_priors(): median pairwise distance is non-positive. Check ",
         "that `coords` has at least 2 distinct (lat, lon) rows.")
  }

  list(
    median_dist = median_dist,
    nt          = nt,
    m0_range    = m0_range,
    mu_sd       = mu_sd,
    sample = list(
      length_scale    = function() stats::rlnorm(1, log(median_dist), 1),
      periodic_scale  = function() stats::rlnorm(1, 0, 1),
      long_term_scale = function() stats::rlnorm(1, log(nt), 1),
      r               = function() stats::rgamma(1, shape = 2, rate = 0.1),
      m0              = function() stats::runif(1, m0_range[1], m0_range[2])
    )
  )
}


# Internal: clustered 0/1 missingness mirror of generate_clustered_binary()
# in implementation/simulation.R. Kept in-package so the PFN simulator is
# self-contained. Returns a length-n integer vector where 1 = missing.
clustered_missing <- function(n, p_one, p_switch) {
  out <- integer(n)
  out[1L] <- stats::rbinom(1, 1, p_one)
  if (n >= 2L) {
    for (i in seq.int(2L, n)) {
      if (stats::runif(1) < p_switch) {
        out[i] <- stats::rbinom(1, 1, p_one)
      } else {
        out[i] <- out[i - 1L]
      }
    }
  }
  out
}


#' Draw one synthetic dataset under the NB-GP generative model
#'
#' Returns an `n x nt` count matrix `y` with a clustered missingness mask,
#' plus the full ground truth (`f`, `mu_s`, `theta`, `r`) needed for the PFN
#' supervised loss. Internal helper for [pfn_simulate_batch()].
#'
#' @param coords Fixed-geometry coords data frame (`lat`, `lon`).
#' @param nt Number of timepoints.
#' @param priors A [pfn_priors()] object.
#' @param period Periodic kernel period (default 52).
#' @param missingness Named list with `p_one`, `p_switch` for the clustered
#'   missingness pattern.
#' @param psi_clamp Symmetric clamp on the log-rate `psi = mu_s + f` to keep
#'   extreme prior draws from producing un-representable counts. Default 20
#'   (i.e. `exp(psi) <= 4.85e8`, well within `.Machine$integer.max`).
#'
#' @return A list with `y`, `mask` (both `n x nt`, integer), `f` (`n x nt`,
#'   numeric), `mu_s` (length `n`), `log_theta` (length 3), `log_r` (scalar).
#' @keywords internal
pfn_simulate_one <- function(coords, nt, priors,
                             period      = 52,
                             missingness = list(p_one = 0.2, p_switch = 0.3),
                             psi_clamp   = 20) {
  n <- nrow(coords)

  theta <- list(
    length_scale    = priors$sample$length_scale(),
    periodic_scale  = priors$sample$periodic_scale(),
    long_term_scale = priors$sample$long_term_scale()
  )
  r    <- priors$sample$r()
  m0   <- priors$sample$m0()
  mu_s <- stats::rnorm(n, mean = m0, sd = priors$mu_sd)

  space_k <- space_kernel(coords, length_scale = theta$length_scale)
  time_k  <- time_kernel(seq_len(nt),
                         periodic_scale  = theta$periodic_scale,
                         long_term_scale = theta$long_term_scale,
                         period          = period)

  # quick_mvnorm returns the flat (time-fastest) vector; recover the n x nt
  # matrix with Z[i, j] = GP value at site i, time j.
  f_flat <- quick_mvnorm(space_k, time_k)
  f_mat  <- t(matrix(f_flat, nrow = nt, ncol = n))    # n x nt

  # psi_mat[i, j] = mu_s[i] + f_mat[i, j]; R recycles mu_s down columns
  # because matrices are column-major. Clamp to avoid pathological draws.
  psi <- f_mat + mu_s
  psi <- pmin(pmax(psi, -psi_clamp), psi_clamp)
  lam <- exp(psi)

  y_mat <- matrix(stats::rnbinom(n * nt, size = r, mu = as.vector(lam)),
                  nrow = n, ncol = nt)

  # Clustered missingness, independent per site (matches observed_data() in
  # implementation/simulation.R).
  mask_mat <- matrix(1L, nrow = n, ncol = nt)
  for (i in seq_len(n)) {
    miss <- clustered_missing(nt, p_one = missingness$p_one,
                              p_switch = missingness$p_switch)
    mask_mat[i, ] <- 1L - miss
  }

  list(
    y         = y_mat,                 # n x nt, integer
    mask      = mask_mat,              # n x nt, integer (1 = observed)
    f         = f_mat,                 # n x nt, numeric
    mu_s      = mu_s,                  # length n
    log_theta = log(c(theta$length_scale,
                      theta$periodic_scale,
                      theta$long_term_scale)),
    log_r     = log(r)
  )
}


#' Draw a batch of synthetic datasets for PFN training
#'
#' Wraps [pfn_simulate_one()] to produce `n_datasets` synthetic datasets and
#' stacks them into batched arrays of shape `(B, n, nt)` (or `(B, n)`,
#' `(B, 3)`, `B` for the lower-dimensional outputs). Results are plain R
#' arrays; convert to torch tensors with [pfn_to_tensors()] only when
#' actually training.
#'
#' Optionally caches the entire batch to disk via `qs::qsave()`. Subsequent
#' calls with the same `cache_path` short-circuit the simulation -- useful
#' during model iteration when you don't want to pay simulator cost on
#' every restart.
#'
#' Parallelism: if `parallel = TRUE` and `future.apply` is installed, draws
#' are dispatched via `future_lapply()`. The caller controls the plan
#' (`future::plan(future::multisession, workers = 8)`); with no plan set
#' you get the sequential fallback.
#'
#' @param coords Fixed-geometry coords data frame; must have `lat`, `lon`
#'   (and any other columns are ignored).
#' @param nt Number of timepoints.
#' @param n_datasets Batch size B.
#' @param priors A [pfn_priors()] object; if NULL, one is constructed from
#'   `coords` and `nt` with default ranges.
#' @param period Periodic kernel period (default 52).
#' @param missingness Named list with `p_one`, `p_switch` for clustered
#'   missingness (default p_one = 0.2, p_switch = 0.3).
#' @param cache_path If non-NULL, a path; results are saved/loaded via
#'   `saveRDS()` / `readRDS()`. The parent directory is created if needed.
#' @param parallel If TRUE and `future.apply` is installed, run draws in
#'   parallel; caller sets the future plan.
#' @param verbose Show progress bar (serial path only).
#'
#' @return A list with elements:
#'   - `y`         : integer array `(B, n, nt)` (zeros at missing cells)
#'   - `mask`      : integer array `(B, n, nt)`, 1 = observed
#'   - `f`         : numeric array `(B, n, nt)`, ground-truth latent field
#'   - `mu_s`      : numeric matrix `(B, n)`, ground-truth site offsets
#'   - `log_theta` : numeric matrix `(B, 3)`, log(length, periodic, long_term)
#'   - `log_r`     : numeric vector length `B`, log(r)
#'   - `coords`    : the input coords
#'   - `n`, `nt`, `period`, `B` : scalars
#' @export
pfn_simulate_batch <- function(coords, nt, n_datasets,
                               priors      = NULL,
                               period      = 52,
                               missingness = list(p_one = 0.2,
                                                  p_switch = 0.3),
                               cache_path  = NULL,
                               parallel    = FALSE,
                               verbose     = TRUE) {

  # Cache hit -- shortcut.
  if (!is.null(cache_path) && file.exists(cache_path)) {
    if (verbose) message("[pfn_simulate_batch] loading cached batch from ",
                         cache_path)
    cached <- readRDS(cache_path)
    # Sanity-check that the cached batch matches the requested shape;
    # otherwise the user has changed coords / nt / n and the cache is stale.
    if (cached$n != nrow(coords) || cached$nt != nt ||
        cached$B != n_datasets) {
      warning("[pfn_simulate_batch] cached batch shape (n=", cached$n,
              ", nt=", cached$nt, ", B=", cached$B, ") does not match ",
              "requested (n=", nrow(coords), ", nt=", nt, ", B=",
              n_datasets, "); re-simulating.")
    } else {
      return(cached)
    }
  }

  if (is.null(priors)) priors <- pfn_priors(coords, nt)
  n <- nrow(coords)
  B <- as.integer(n_datasets)

  one <- function(.) {
    pfn_simulate_one(coords, nt, priors, period = period,
                     missingness = missingness)
  }

  use_parallel <- isTRUE(parallel) &&
                  requireNamespace("future.apply", quietly = TRUE)
  if (isTRUE(parallel) && !use_parallel) {
    warning("`parallel = TRUE` requested but future.apply is not installed; ",
            "falling back to serial.")
  }

  draws <- if (use_parallel) {
    future.apply::future_lapply(
      seq_len(B), one,
      future.seed     = TRUE,
      future.packages = "weave"   # workers need internal pfn_simulate_one
    )
  } else {
    pb <- if (verbose) {
      progress::progress_bar$new(
        format = "  simulating [:bar] :percent eta::eta",
        total = B, clear = FALSE, width = 60)
    } else NULL
    out <- vector("list", B)
    for (i in seq_len(B)) {
      out[[i]] <- one(i)
      if (!is.null(pb)) pb$tick()
    }
    out
  }

  # Stack -- pre-allocate and fill, much cheaper than abind on large B.
  batch <- list(
    y         = array(0L,        dim = c(B, n, nt)),
    mask      = array(0L,        dim = c(B, n, nt)),
    f         = array(0,         dim = c(B, n, nt)),
    mu_s      = matrix(0,        nrow = B, ncol = n),
    log_theta = matrix(0,        nrow = B, ncol = 3,
                       dimnames = list(NULL, c("length_scale",
                                               "periodic_scale",
                                               "long_term_scale"))),
    log_r     = numeric(B),
    coords    = coords,
    n         = n,
    nt        = nt,
    period    = period,
    B         = B,
    priors    = priors
  )

  for (i in seq_len(B)) {
    d <- draws[[i]]
    batch$y[i, , ]    <- d$y
    batch$mask[i, , ] <- d$mask
    batch$f[i, , ]    <- d$f
    batch$mu_s[i, ]   <- d$mu_s
    batch$log_theta[i, ] <- d$log_theta
    batch$log_r[i]    <- d$log_r
  }

  if (!is.null(cache_path)) {
    dir.create(dirname(cache_path), recursive = TRUE, showWarnings = FALSE)
    if (verbose) message("[pfn_simulate_batch] caching to ", cache_path)
    saveRDS(batch, cache_path)
  }

  batch
}


#' Convert a PFN batch (R arrays) to torch tensors on the requested device
#'
#' Lazy conversion -- keeps the simulator itself torch-free so it can run
#' in environments without libtorch (CI, plain R sessions), and lets us
#' move whole batches to GPU with one call when that becomes relevant.
#'
#' @param batch A list returned by [pfn_simulate_batch()].
#' @param device Torch device string (default `"cpu"`). Pass `"cuda"` once
#'   you're on a GPU box.
#'
#' @return A list of torch tensors mirroring `batch`'s numeric components.
#'   Non-tensor metadata (`coords`, `n`, `nt`, `period`, `B`, `priors`) is
#'   copied through unchanged.
#' @export
pfn_to_tensors <- function(batch, device = "cpu") {
  if (!requireNamespace("torch", quietly = TRUE)) {
    stop("`torch` is required for pfn_to_tensors(). ",
         "Install with install.packages('torch') and torch::install_torch().")
  }
  tt <- function(x, dtype) {
    torch::torch_tensor(x, dtype = dtype, device = device)
  }
  list(
    y         = tt(batch$y,         torch::torch_float()),
    mask      = tt(batch$mask,      torch::torch_float()),
    f         = tt(batch$f,         torch::torch_float()),
    mu_s      = tt(batch$mu_s,      torch::torch_float()),
    log_theta = tt(batch$log_theta, torch::torch_float()),
    log_r     = tt(batch$log_r,     torch::torch_float()),
    coords    = batch$coords,
    n         = batch$n,
    nt        = batch$nt,
    period    = batch$period,
    B         = batch$B,
    priors    = batch$priors
  )
}
