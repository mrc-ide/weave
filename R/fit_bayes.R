# =============================================================================
# Public Bayesian fit -- the orchestration around the PG-Gibbs sampler.
#
# Workflow:
#   1. build_design()                  -- pack the data into compact form.
#   2. bayes_priors()                  -- assemble prior densities and samplers.
#   3. for each chain: initialise state (overdispersed from priors), run
#      burn-in + sampling, store per-chain traces.
#   4. pool chain outputs into a single weave_bayes object.
#
# Multi-chain support: `n_chains` > 1 runs independent chains serially with
# inits drawn from the priors so chains start at different points. Per-chain
# traces are stored with the chain as the last array dimension, which
# matches the layout coda::as.mcmc.list expects (and which the summary
# method uses to report Gelman-Rubin Rhat when n_chains > 1).
# =============================================================================


#' Build a default set of priors for fit_bayes()
#'
#' Each prior entry exposes two closures:
#'   * `logp(x)` -- log prior density at `x`
#'   * `sample()` -- one draw from the prior (used to overdisperse chain inits)
#'
#' Replace any single entry to swap in a custom prior; the only contract is
#' that `logp` and `sample` are consistent.
#'
#' Defaults are weakly informative:
#'   * length_scale     ~ LogNormal(meanlog = log(median_dist), sdlog = 1)
#'   * periodic_scale   ~ LogNormal(meanlog = 0,                sdlog = 1)
#'   * long_term_scale  ~ LogNormal(meanlog = log(n_times),     sdlog = 1)
#'   * mu_s             ~ Normal(m0 = log(median count + 1), v0 = 4)
#'   * r                ~ Gamma(shape = 2, rate = 0.1)            (mean 20)
#'
#' All three length scales are positive; their LogNormal priors keep the
#' slice samplers from wandering into negative territory.
#'
#' @param design A design list from [build_design()] (used to anchor priors
#'   to the data scale).
#' @export
bayes_priors <- function(design) {
  d <- get_spatial_distance(design$coords)
  diag(d) <- NA
  median_dist <- stats::median(d, na.rm = TRUE)

  median_y <- stats::median(design$y_obs, na.rm = TRUE)
  m0_mu    <- log(max(median_y + 1, 1))

  list(
    theta = list(
      length_scale    = list(
        logp   = function(x) stats::dlnorm(x, log(median_dist), 1, log = TRUE),
        sample = function() stats::rlnorm(1, log(median_dist), 1)
      ),
      periodic_scale  = list(
        logp   = function(x) stats::dlnorm(x, 0, 1, log = TRUE),
        sample = function() stats::rlnorm(1, 0, 1)
      ),
      long_term_scale = list(
        logp   = function(x) stats::dlnorm(x, log(design$nt), 1, log = TRUE),
        sample = function() stats::rlnorm(1, log(design$nt), 1)
      )
    ),
    mu = list(
      m0 = m0_mu,
      v0 = 4
    ),
    r = list(
      logp   = function(x) stats::dgamma(x, shape = 2, rate = 0.1, log = TRUE),
      sample = function() stats::rgamma(1, shape = 2, rate = 0.1)
    )
  )
}


# Internal: validate the user-facing fit_bayes() inputs. Called at the very
# top of fit_bayes() so problems surface with a clear, single message
# instead of an obscure crash inside the inner loop.
validate_fit_bayes_inputs <- function(obs_data, n_sweeps, burnin, n_chains,
                                      fix, period, slice_widths, n_thin_f) {

  # obs_data must be a data frame / tibble with the required structure.
  if (!is.data.frame(obs_data)) {
    stop("`obs_data` must be a data frame.")
  }
  required <- c("id", "t", "lat", "lon")
  miss <- setdiff(required, names(obs_data))
  if (length(miss) > 0) {
    stop("`obs_data` is missing required column(s): ",
         paste(miss, collapse = ", "), ".")
  }
  count_col <- if ("y_obs" %in% names(obs_data)) "y_obs" else "n"
  if (!count_col %in% names(obs_data)) {
    stop("`obs_data` must contain a count column named `y_obs` or `n`.")
  }
  y <- obs_data[[count_col]]
  # Counts must be non-negative integers (NA allowed for missingness).
  y_present <- y[!is.na(y)]
  if (length(y_present) > 0 &&
      (any(y_present < 0) || any(y_present != floor(y_present)))) {
    stop("`", count_col, "` must contain non-negative integers ",
         "(NA for missing).")
  }

  # Scalar integer checks for sampling-loop dimensions.
  for (nm in c("n_sweeps", "burnin", "n_chains")) {
    v <- get(nm)
    if (!is.numeric(v) || length(v) != 1L || !is.finite(v) || v < 0 ||
        v != floor(v)) {
      stop("`", nm, "` must be a single non-negative integer.")
    }
  }
  if (n_chains < 1L) stop("`n_chains` must be >= 1.")
  if (burnin >= n_sweeps) {
    stop("`burnin` (", burnin, ") must be < `n_sweeps` (", n_sweeps, ").")
  }
  if (!is.numeric(n_thin_f) || length(n_thin_f) != 1L || n_thin_f < 1) {
    stop("`n_thin_f` must be a positive integer.")
  }

  # period must be a positive integer; we'll separately check period <= nt
  # once design has been built.
  if (!is.numeric(period) || length(period) != 1L || !is.finite(period) ||
      period <= 0 || period != floor(period)) {
    stop("`period` must be a positive integer.")
  }

  # fix$r, if supplied, must be a finite positive scalar.
  if (!is.null(fix$r)) {
    if (!is.numeric(fix$r) || length(fix$r) != 1L || !is.finite(fix$r) ||
        fix$r <= 0) {
      stop("`fix$r` must be a single positive finite number; got ",
           format(fix$r), ".")
    }
  }

  # slice_widths must have exactly the three expected components and all
  # be positive scalars.
  needed_sw <- c("length_scale", "periodic_scale", "long_term_scale")
  if (!is.list(slice_widths) || !all(needed_sw %in% names(slice_widths))) {
    stop("`slice_widths` must be a named list with entries ",
         paste(needed_sw, collapse = ", "), ".")
  }
  for (nm in needed_sw) {
    v <- slice_widths[[nm]]
    if (!is.numeric(v) || length(v) != 1L || !is.finite(v) || v <= 0) {
      stop("`slice_widths$", nm, "` must be a single positive number.")
    }
  }

  invisible(TRUE)
}


# Internal: run a single PG-Gibbs chain. Returns the per-chain storage that
# fit_bayes() pools across chains.
run_one_chain <- function(design, priors, init_state, n_sweeps, burnin,
                          period, pcg_tol, pcg_maxit, slice_widths,
                          store_f, n_thin_f, chain_id, n_chains, verbose) {

  state <- init_state
  state$space_mat <- space_kernel(
    coordinates  = design$coords,
    length_scale = state$theta$length_scale
  )
  state$time_mat <- time_kernel(
    times           = seq_len(design$nt),
    periodic_scale  = state$theta$periodic_scale,
    long_term_scale = state$theta$long_term_scale,
    period          = period
  )
  state$ke <- kron_eigen(state$space_mat, state$time_mat)

  control <- list(
    pcg_tol      = pcg_tol,
    pcg_maxit    = pcg_maxit,
    slice_widths = slice_widths,
    mh_state     = list(log_step = log(0.3),
                        gamma    = 1.0,
                        iter     = 0L,
                        accepts  = 0L,
                        attempts = 0L)
  )

  theta_trace    <- matrix(NA_real_, nrow = n_sweeps, ncol = 3,
                           dimnames = list(NULL, c("length_scale",
                                                   "periodic_scale",
                                                   "long_term_scale")))
  mu_trace       <- matrix(NA_real_, nrow = n_sweeps, ncol = design$n)
  r_trace        <- numeric(n_sweeps)
  logpost_trace  <- numeric(n_sweeps)
  pcg_iter_trace <- integer(n_sweeps)
  # Health counters for the f-block PCG. We track these per chain so
  # summary.weave_bayes() can flag a sampler that has been quietly
  # struggling rather than burying it in 5000 warnings.
  pcg_fails      <- 0L
  jacobi_uses    <- 0L

  f_summary <- welford_new(design$N)

  n_post <- n_sweeps - burnin
  thin_idx <- integer(0)
  if (store_f != "summary") {
    thin_idx <- if (store_f == "all") {
      seq(burnin + 1L, n_sweeps)
    } else {
      round(seq(burnin + 1L, n_sweeps, length.out = min(n_thin_f, n_post)))
    }
  }
  f_samples <- if (length(thin_idx) > 0)
    matrix(NA_real_, nrow = length(thin_idx), ncol = design$N) else NULL
  thin_ptr  <- 1L

  pb <- if (verbose) progress::progress_bar$new(
    format = sprintf("  chain %d/%d [:bar] :percent eta::eta",
                     chain_id, n_chains),
    total = n_sweeps, clear = FALSE, width = 60) else NULL

  for (sw in seq_len(n_sweeps)) {
    adapt_r <- (sw <= burnin) && !state$r_fixed
    step <- pg_sweep(state, design, priors, control, adapt_r = adapt_r)
    state   <- step$state
    control <- step$control

    theta_trace[sw, ] <- c(state$theta$length_scale,
                           state$theta$periodic_scale,
                           state$theta$long_term_scale)
    mu_trace[sw, ]     <- state$mu
    r_trace[sw]        <- state$r
    logpost_trace[sw]  <- state$log_post
    pcg_iter_trace[sw] <- step$diag$pcg_iters
    if (isFALSE(step$diag$pcg_converged)) pcg_fails   <- pcg_fails   + 1L
    if (isTRUE(step$diag$used_jacobi))    jacobi_uses <- jacobi_uses + 1L

    if (sw > burnin) {
      f_summary <- welford_update(f_summary, state$f)
      if (length(thin_idx) > 0 && thin_ptr <= length(thin_idx) &&
          sw == thin_idx[thin_ptr]) {
        f_samples[thin_ptr, ] <- state$f
        thin_ptr <- thin_ptr + 1L
      }
    }

    if (!is.null(pb)) pb$tick()
  }

  list(
    theta_trace    = theta_trace,
    mu_trace       = mu_trace,
    r_trace        = r_trace,
    logpost_trace  = logpost_trace,
    pcg_iter_trace = pcg_iter_trace,
    pcg_fails      = pcg_fails,
    jacobi_uses    = jacobi_uses,
    f_summary      = f_summary,
    f_samples      = f_samples,
    thin_idx       = thin_idx,
    mh_acceptance  = if (control$mh_state$attempts > 0)
      control$mh_state$accepts / control$mh_state$attempts else NA_real_
  )
}


# Build a per-chain initial state. If n_chains > 1 and no explicit init was
# provided, draw theta and r from the priors so chains start at different
# points (overdispersion is what makes the Gelman-Rubin diagnostic useful).
init_state_for_chain <- function(chain_id, n_chains, design, priors,
                                 init, fix) {
  use_overdispersed <- (n_chains > 1) && is.null(init$theta) && chain_id > 1

  theta_init <- if (!is.null(init$theta)) {
    init$theta
  } else if (use_overdispersed) {
    list(
      length_scale    = priors$theta$length_scale$sample(),
      periodic_scale  = priors$theta$periodic_scale$sample(),
      long_term_scale = priors$theta$long_term_scale$sample()
    )
  } else {
    list(
      length_scale    = stats::median(stats::dist(design$coords[, c("lat","lon")])),
      periodic_scale  = 1,
      long_term_scale = design$nt
    )
  }

  r_init <- if (!is.null(init$r)) init$r else if (use_overdispersed) {
    priors$r$sample()
  } else 10

  state <- list(
    f       = numeric(design$N),
    mu      = if (!is.null(init$mu)) init$mu else design$mu_init,
    theta   = theta_init,
    r       = r_init,
    r_fixed = FALSE
  )
  if (!is.null(fix$r)) {
    state$r       <- fix$r
    state$r_fixed <- TRUE
  }
  state
}


#' Full Bayesian fit via PG-Gibbs sampling
#'
#' Samples from the posterior over (f, mu, theta, r) under the NB-GP model
#' using Pólya-Gamma data augmentation. The latent-field block is an exact
#' Gaussian draw via PCG; theta is updated by univariate slice on the log
#' scale; r by adaptive random-walk MH; no Laplace approximation anywhere.
#'
#' For a 1000-site, 150-week dataset on a laptop, expect ~2-5 seconds per
#' sweep per chain with the default settings. Multiple chains run serially.
#'
#' @param obs_data Data frame with columns `id`, `t`, `lat`, `lon`, and a
#'   count column (`y_obs` if present, else `n`). NA in the count column
#'   denotes missing.
#' @param n_sweeps Total number of sweeps per chain.
#' @param burnin Number of sweeps to discard per chain. r-MH adapts during
#'   burn-in.
#' @param n_chains Number of independent chains (default 1). Chains use
#'   overdispersed initial values drawn from the priors so Gelman-Rubin
#'   Rhat is informative; the first chain still starts at deterministic
#'   "reasonable" inits to give a fast warm-up.
#' @param thin Keep every `thin`-th post-burnin sweep in stored f samples.
#' @param priors Prior list from [bayes_priors()]; pass a custom list to
#'   override individual entries.
#' @param init Optional init list (`f`, `mu`, `theta`, `r`); any missing
#'   entry is filled with sensible defaults. Applied identically to every
#'   chain (disables overdispersion if `theta` is set).
#' @param fix Optional `list(r = <value>)` to fix r (e.g. a large finite
#'   number for the Poisson limit).
#' @param period Periodic kernel period (default 52).
#' @param pcg_tol PCG relative-residual tolerance in the f-block. The
#'   default 1e-4 is intentionally loose: in a Gibbs sampler the noise from
#'   the perturbation draws (u, e) dominates anything the solver could
#'   tighten below this. Use 1e-6 if you want tighter f-block convergence
#'   at the cost of more PCG iterations per sweep.
#' @param pcg_maxit PCG iteration cap in the f-block.
#' @param slice_widths Initial slice widths for the three theta components.
#' @param store_f One of `"summary"`, `"thin"`, `"all"`.
#' @param n_thin_f Number of full f samples to retain per chain when
#'   `store_f = "thin"`. Default 100.
#' @param verbose Show a progress bar per chain.
#' @param parallel If TRUE and `n_chains > 1`, run chains in parallel via
#'   `future.apply::future_lapply()`. The user is responsible for the
#'   `future::plan()` (eg `future::plan("multisession", workers = 4)`);
#'   without one the default sequential plan is used and you'll see no
#'   speedup. Requires the `future.apply` package; falls back to serial
#'   with a warning if it's not installed. Default FALSE.
#' @return A `weave_bayes` S3 object. Per-chain traces are stored with the
#'   chain as the last array dimension, eg `theta_trace[sweep, param, chain]`.
#' @export
fit_bayes <- function(obs_data,
                      n_sweeps     = 1000,
                      burnin       = 500,
                      n_chains     = 1L,
                      thin         = 1,
                      priors       = NULL,
                      init         = NULL,
                      fix          = list(),
                      period       = 52,
                      pcg_tol      = 1e-4,
                      pcg_maxit    = 500,
                      slice_widths = list(length_scale = 0.5,
                                          periodic_scale = 0.5,
                                          long_term_scale = 0.5),
                      store_f      = c("summary", "thin", "all"),
                      n_thin_f     = 100L,
                      verbose      = TRUE,
                      parallel     = FALSE) {

  store_f <- match.arg(store_f)
  n_chains <- as.integer(n_chains)

  # Validate user inputs up front and produce a clear error per category,
  # rather than letting a malformed input crash a hot inner loop with a
  # cryptic "non-conformable arguments" message twenty minutes in.
  validate_fit_bayes_inputs(
    obs_data     = obs_data,
    n_sweeps     = n_sweeps,
    burnin       = burnin,
    n_chains     = n_chains,
    fix          = fix,
    period       = period,
    slice_widths = slice_widths,
    n_thin_f     = n_thin_f
  )

  design  <- build_design(obs_data)
  design$period <- period
  if (is.null(priors)) priors <- bayes_priors(design)
  # period vs nt is a check that needs the design to be built first.
  if (period > design$nt) {
    stop("period (", period, ") must be <= nt (", design$nt, ").")
  }

  # Decide on parallel vs serial. We need future.apply available and we
  # need more than one chain for parallel to be meaningful.
  use_parallel <- isTRUE(parallel) && n_chains > 1L
  if (use_parallel && !requireNamespace("future.apply", quietly = TRUE)) {
    warning("`parallel = TRUE` requested but future.apply is not installed; ",
            "falling back to serial execution.")
    use_parallel <- FALSE
  }

  # Per-chain integer seeds derived from the user's RNG state.
  # `sample.int()` consumes one position from the user's RNG (so they retain
  # full control via set.seed before calling fit_bayes), then each chain
  # set.seed()s with its assigned integer before any sampling happens.
  # This gives bit-identical results across serial and parallel execution
  # and across any order the parallel scheduler dispatches chains.
  chain_seeds <- sample.int(.Machine$integer.max, n_chains)

  # Run a single chain. Used both in the serial loop and the future.apply
  # call. Encapsulates "set seed, build init state, call run_one_chain".
  run_chain <- function(ch) {
    set.seed(chain_seeds[ch])
    init_st <- init_state_for_chain(ch, n_chains, design, priors, init, fix)
    run_one_chain(
      design = design, priors = priors, init_state = init_st,
      n_sweeps = n_sweeps, burnin = burnin,
      period = period, pcg_tol = pcg_tol, pcg_maxit = pcg_maxit,
      slice_widths = slice_widths,
      store_f = store_f, n_thin_f = n_thin_f,
      chain_id = ch, n_chains = n_chains,
      # Per-chain progress bars interleave noisily across workers, so
      # silence them when running in parallel.
      verbose = verbose && !use_parallel
    )
  }

  # ---- Run chains -----------------------------------------------------------
  if (use_parallel) {
    if (verbose) {
      message(sprintf("Running %d chains in parallel via future.apply (%s).",
                      n_chains, class(future::plan())[1]))
    }
    chain_results <- future.apply::future_lapply(
      seq_len(n_chains), run_chain,
      future.seed     = NULL,           # we've handled RNG ourselves above
      future.packages = "weave"          # workers need package functions
    )
  } else {
    chain_results <- lapply(seq_len(n_chains), run_chain)
  }

  # ---- Pool per-chain traces into arrays with chain as last dim -------------
  theta_trace <- array(
    NA_real_, dim = c(n_sweeps, 3, n_chains),
    dimnames = list(NULL,
                    c("length_scale", "periodic_scale", "long_term_scale"),
                    paste0("chain", seq_len(n_chains)))
  )
  mu_trace       <- array(NA_real_, dim = c(n_sweeps, design$n, n_chains))
  r_trace        <- matrix(NA_real_, nrow = n_sweeps, ncol = n_chains)
  logpost_trace  <- matrix(NA_real_, nrow = n_sweeps, ncol = n_chains)
  pcg_iter_trace <- matrix(NA_integer_, nrow = n_sweeps, ncol = n_chains)
  mh_acceptance  <- numeric(n_chains)
  pcg_fails      <- integer(n_chains)
  jacobi_uses    <- integer(n_chains)

  for (ch in seq_len(n_chains)) {
    cr <- chain_results[[ch]]
    theta_trace[, , ch]   <- cr$theta_trace
    mu_trace[, , ch]      <- cr$mu_trace
    r_trace[, ch]         <- cr$r_trace
    logpost_trace[, ch]   <- cr$logpost_trace
    pcg_iter_trace[, ch]  <- cr$pcg_iter_trace
    mh_acceptance[ch]     <- cr$mh_acceptance
    # Defensive: an older worker (loaded from an installed package) may
    # not yet include these fields. Default to 0 in that case.
    pcg_fails[ch]         <- if (is.null(cr$pcg_fails))   0L else cr$pcg_fails
    jacobi_uses[ch]       <- if (is.null(cr$jacobi_uses)) 0L else cr$jacobi_uses
  }

  # ---- Pool Welford summaries across chains (post-burnin samples only) -----
  pooled <- welford_combine(lapply(chain_results, function(cr) cr$f_summary))

  # ---- Concatenate thinned f / mu / r samples across chains -----------------
  f_samples <- if (!is.null(chain_results[[1]]$f_samples)) {
    do.call(rbind, lapply(chain_results, function(cr) cr$f_samples))
  } else NULL

  mu_samples <- if (!is.null(f_samples)) {
    do.call(rbind, lapply(seq_len(n_chains), function(ch) {
      # mu_trace[, , ch] is n_sweeps x n_sites; index rows by thin_idx and
      # drop the chain dimension. Without drop = TRUE the slice would stay
      # 3-D and rbind would mis-stack the chains.
      mu_trace[chain_results[[ch]]$thin_idx, , ch]
    }))
  } else NULL

  r_samples <- if (!is.null(f_samples)) {
    unlist(lapply(seq_len(n_chains), function(ch) {
      r_trace[chain_results[[ch]]$thin_idx, ch]
    }))
  } else NULL

  result <- list(
    design          = design,
    priors          = priors,
    n_sweeps        = n_sweeps,
    burnin          = burnin,
    n_chains        = n_chains,
    theta_trace     = theta_trace,
    mu_trace        = mu_trace,
    r_trace         = r_trace,
    logpost_trace   = logpost_trace,
    pcg_iter_trace  = pcg_iter_trace,
    pcg_fails       = pcg_fails,
    jacobi_uses     = jacobi_uses,
    f_mean          = pooled$mean,
    f_var           = welford_var(pooled),
    f_samples       = f_samples,
    mu_samples      = mu_samples,
    r_samples       = r_samples,
    mh_acceptance   = mh_acceptance
  )
  class(result) <- c("weave_bayes", "list")
  result
}


#' @export
print.weave_bayes <- function(x, ...) {
  cat("<weave_bayes fit>\n")
  cat(sprintf("  chains      : %d\n", x$n_chains))
  cat(sprintf("  sweeps      : %d (burnin %d) per chain\n",
              x$n_sweeps, x$burnin))
  cat(sprintf("  sites x time: %d x %d  (N = %d, observed = %d)\n",
              x$design$n, x$design$nt, x$design$N,
              length(x$design$obs_idx)))
  if (length(x$mh_acceptance) == 1L) {
    cat(sprintf("  r-MH accept : %.3f\n", x$mh_acceptance))
  } else {
    cat(sprintf("  r-MH accept : %s\n",
                paste(sprintf("%.2f", x$mh_acceptance), collapse = " / ")))
  }
  cat(sprintf("  median PCG  : %d iters\n",
              as.integer(stats::median(x$pcg_iter_trace))))
  invisible(x)
}


#' @export
summary.weave_bayes <- function(object, ...) {
  post <- (object$burnin + 1L):object$n_sweeps

  # Pool post-burnin draws across chains for quantile summary.
  ls_all <- as.vector(object$theta_trace[post, "length_scale",    , drop = TRUE])
  ps_all <- as.vector(object$theta_trace[post, "periodic_scale",  , drop = TRUE])
  lt_all <- as.vector(object$theta_trace[post, "long_term_scale", , drop = TRUE])
  r_all  <- as.vector(object$r_trace[post, , drop = TRUE])

  quant <- function(x) stats::quantile(x, c(0.025, 0.5, 0.975))
  out <- rbind(
    length_scale    = quant(ls_all),
    periodic_scale  = quant(ps_all),
    long_term_scale = quant(lt_all),
    r               = quant(r_all)
  )
  colnames(out) <- c("2.5%", "50%", "97.5%")

  # Gelman-Rubin Rhat if we have multiple chains and coda is available.
  if (object$n_chains > 1 && requireNamespace("coda", quietly = TRUE)) {
    make_mcmc <- function(arr) {
      coda::as.mcmc.list(lapply(seq_len(object$n_chains), function(ch) {
        coda::mcmc(arr[post, ch], start = object$burnin + 1L)
      }))
    }
    make_theta_mcmc <- function(name) {
      coda::as.mcmc.list(lapply(seq_len(object$n_chains), function(ch) {
        coda::mcmc(object$theta_trace[post, name, ch],
                   start = object$burnin + 1L)
      }))
    }
    rhats <- c(
      length_scale    = coda::gelman.diag(make_theta_mcmc("length_scale"),
                                          autoburnin = FALSE)$psrf[1],
      periodic_scale  = coda::gelman.diag(make_theta_mcmc("periodic_scale"),
                                          autoburnin = FALSE)$psrf[1],
      long_term_scale = coda::gelman.diag(make_theta_mcmc("long_term_scale"),
                                          autoburnin = FALSE)$psrf[1],
      r               = coda::gelman.diag(make_mcmc(object$r_trace),
                                          autoburnin = FALSE)$psrf[1]
    )
    out <- cbind(out, Rhat = rhats)
  }

  print(out)
  if (length(object$mh_acceptance) == 1L) {
    cat(sprintf("\nMH acceptance for r : %.3f\n", object$mh_acceptance))
  } else {
    cat(sprintf("\nMH acceptance for r (per chain): %s\n",
                paste(sprintf("%.2f", object$mh_acceptance), collapse = ", ")))
  }

  # PCG-block health summary. Show a one-line summary always; warn if
  # more than 10% of sweeps were unhealthy (non-converged PCG or fallback
  # to Jacobi). This makes a quietly-struggling sampler hard to miss.
  if (!is.null(object$pcg_fails)) {
    total_sweeps <- object$n_sweeps * object$n_chains
    n_fail   <- sum(object$pcg_fails)
    n_jacobi <- sum(object$jacobi_uses)
    cat(sprintf("PCG: %d non-converged, %d Jacobi fallback (of %d sweeps).\n",
                n_fail, n_jacobi, total_sweeps))
    if (n_fail + n_jacobi > 0.10 * total_sweeps) {
      warning("More than 10% of sweeps had a PCG failure or Jacobi ",
              "fallback. Consider increasing pcg_maxit or relaxing pcg_tol.")
    }
  }
  invisible(out)
}
