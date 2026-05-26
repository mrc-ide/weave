# =============================================================================
# Pólya-Gamma augmented Gibbs sampler -- block updates
#
# This file implements the five conditional draws that make up one sweep of
# the PG-Gibbs sampler for the NB-GP model:
#
#   y_i  ~ NB(r, mu_i),    log(mu_i) = psi_i = f_{cell(i)} + mu_{site(i)}
#   f    ~ N(0, K_space(theta_s) (x) K_time(theta_t))
#
# Augmentation. The PG identity (Polson, Scott & Windle 2013) writes
#
#   P(y | r, psi) ∝ ∫ exp(kappa * eta - (omega/2) eta^2) p_PG(omega; y+r, 0) d omega
#
# with eta = psi - log(r) and kappa = (y - r)/2. Conditional on omega the
# NB likelihood is GAUSSIAN in eta -- and therefore Gaussian in psi, which
# is what the GP latent field needs.
#
# After augmentation one sweep looks like:
#
#   1. omega | f, mu, r, y       -- independent PG draws (pgdraw, vectorised)
#   2. f     | omega, mu, theta  -- Gaussian; perturbation sampler via PCG
#   3. mu    | omega, f, y, r    -- closed-form Normal per site
#   4. theta | f                 -- univariate slice on log scales
#   5. r     | y, psi            -- random-walk MH on log r
#
# Conventions:
#   - omega lives at observed cells only (length m).
#   - f lives on the full grid (length N = n*nt), times varying fastest.
#   - psi = f[obs_idx] + mu[site_idx_obs] cached as design$N-sized? No: length m.
#   - y_star_i = (y_i - r) / (2 omega_i) + log(r)   -- pseudo-obs of psi_i.
#   - All slice / MH steps work on the LOG scale of positive parameters
#     (length scales, r) for numerical stability and better mixing.
# =============================================================================


# -----------------------------------------------------------------------------
# Block 1: omega | f, mu, r, y
# -----------------------------------------------------------------------------
# For each observed cell:
#     eta_i  = psi_i - log(r)        (logit of NB success probability)
#     omega_i ~ PG(y_i + r, eta_i)
# BayesLogit::rpg accepts real-valued shape (h) -- needed because r is
# real-valued in the NB likelihood.
# -----------------------------------------------------------------------------
pg_draw_omega <- function(psi, y, r) {
  BayesLogit::rpg(num = length(y), h = y + r, z = psi - log(r))
}


# -----------------------------------------------------------------------------
# Block 2: f | omega, mu, theta -- the heavy block
#
# Conditional on (omega, mu, theta, y, r), the latent field is Gaussian:
#     Sigma_post^-1 = K^-1 + S' Omega S
#     mean_post     = Sigma_post * S' Omega (y_star - mu_per_obs)
# We never form Sigma_post. We use the Bhattacharya-Papandreou perturbation
# sampler, which produces an exact draw from N(mean_post, Sigma_post) via:
#
#     u  ~ N(0, K)                                   -- prior draw
#     e  ~ N(0, Omega^-1) (length m, diagonal)
#     y_p = (y_star - mu_per_obs) - S u - e
#     gamma = (S K S' + Omega^-1)^-1 y_p             -- one PCG solve
#     f  = u + K S' gamma
#
# Correctness: see eg Bhattacharya, Chakraborty & Mallick 2016. Cost: one
# kron_mv-based PCG solve, dominated by the matrix-free matvec.
# -----------------------------------------------------------------------------
pg_draw_f <- function(state, design, control) {
  obs_idx <- design$obs_idx
  N       <- design$N
  m       <- length(obs_idx)
  omega   <- state$omega

  # Pseudo-observations of psi (length m).
  y_star <- (design$y_obs - state$r) / (2 * omega) + log(state$r)
  mu_per_obs <- state$mu[design$site_idx_obs]

  # 1. Prior draw u ~ N(0, K). quick_mvnorm uses Cholesky factors of the
  #    small kernels; for n ~ 1000 this is the most expensive single line
  #    after the slice block.
  u <- quick_mvnorm(state$space_mat, state$time_mat)

  # 2. Observation noise e ~ N(0, Omega^-1).
  e <- stats::rnorm(m) / sqrt(omega)

  # 3. Right-hand side.
  noise_var <- 1 / omega
  b <- (y_star - mu_per_obs) - u[obs_idx] - e

  # 4. PCG solve with Kron-eigen preconditioner.
  Amv_fun <- function(v) Amv(v, obs_idx, N, state$space_mat, state$time_mat,
                             noise_var)
  Minv    <- kron_eigen_preconditioner(state$ke,
                                       sigma2 = mean(noise_var),
                                       obs_idx, N)
  pcg_res <- pcg(b, Amv_fun, Minv,
                 tol = control$pcg_tol, maxit = control$pcg_maxit)

  # 5. Sample.
  f <- u + kron_mv(fill_vector(pcg_res$x, obs_idx, N),
                   state$space_mat, state$time_mat)

  list(f = f, pcg_iters = pcg_res$iters, pcg_converged = pcg_res$converged)
}


# -----------------------------------------------------------------------------
# Block 3: mu_s | omega, f, y, r -- closed-form Normal per site
#
# The pseudo-observation y_star_i is N(psi_i, 1/omega_i) ; subtracting f gives
# y_star_i - f_{cell(i)} ~ N(mu_{site(i)}, 1/omega_i). Combined with a
# Normal(m_0, v_0) prior, each site has independent Normal posterior:
#     prec_post = 1/v_0 + sum_{i in site s} omega_i
#     mean_post = (m_0 / v_0 + sum_{i in site s} omega_i * (y_star_i - f_i)) / prec_post
# Implemented with two tapply() / accumarray-style sums.
# -----------------------------------------------------------------------------
pg_draw_mu <- function(state, design, mu_prior) {
  obs_idx     <- design$obs_idx
  n           <- design$n
  omega       <- state$omega
  y_star      <- (design$y_obs - state$r) / (2 * omega) + log(state$r)
  resid       <- y_star - state$f[obs_idx]
  site_idx_obs <- design$site_idx_obs

  # Per-site accumulators.
  sum_omega <- tapply(omega,            site_idx_obs, sum, default = 0)
  sum_oy    <- tapply(omega * resid,    site_idx_obs, sum, default = 0)

  # tapply returns a named vector indexed by the unique site ids that appear
  # in site_idx_obs; pad to length n.
  sum_omega_full <- numeric(n)
  sum_oy_full    <- numeric(n)
  present        <- as.integer(names(sum_omega))
  sum_omega_full[present] <- sum_omega
  sum_oy_full[present]    <- sum_oy

  prec_post <- 1 / mu_prior$v0 + sum_omega_full
  mean_post <- (mu_prior$m0 / mu_prior$v0 + sum_oy_full) / prec_post

  stats::rnorm(n, mean = mean_post, sd = 1 / sqrt(prec_post))
}


# -----------------------------------------------------------------------------
# Block 4: theta | f -- univariate slice on each of log(length_scale),
#                       log(periodic_scale), log(long_term_scale).
#
# Target log-density (up to constants in theta):
#     log p(theta | f) = log prior(theta) - 0.5 log|K(theta)| - 0.5 f' K(theta)^-1 f
#
# Every slice evaluation rebuilds the relevant kernel and its eigendecomp.
# We accept the cost (n x n eigen ~ 1s for n=1000) once per *accepted* slice
# proposal; rejected proposals re-use the cached log-density at the original
# theta. Sample on the log scale: the support is (0, infty) and the
# posterior tends to be heavy-tailed.
# -----------------------------------------------------------------------------
pg_draw_theta <- function(state, design, theta_priors, slice_widths) {

  # Build a closure that, given a candidate theta value for one component,
  # rebuilds the kernels and evaluates the log target.
  build_kernels <- function(theta) {
    space_mat <- space_kernel(
      coordinates  = design$coords,
      length_scale = theta$length_scale
    )
    time_mat <- time_kernel(
      times           = seq_len(design$nt),
      periodic_scale  = theta$periodic_scale,
      long_term_scale = theta$long_term_scale,
      period          = design$period
    )
    ke <- kron_eigen(space_mat, time_mat)
    list(space_mat = space_mat, time_mat = time_mat, ke = ke)
  }
  log_target <- function(theta, mats) {
    -0.5 * mats$ke$log_det -
      0.5 * kron_quad(state$f, mats$ke) +
      theta_priors$length_scale$logp(theta$length_scale) +
      theta_priors$periodic_scale$logp(theta$periodic_scale) +
      theta_priors$long_term_scale$logp(theta$long_term_scale)
  }

  # Current cached values.
  theta <- state$theta
  mats  <- list(space_mat = state$space_mat, time_mat = state$time_mat,
                ke = state$ke)
  current_lp <- log_target(theta, mats)

  # Component-wise slice on the log scale.
  for (component in c("length_scale", "periodic_scale", "long_term_scale")) {
    w <- slice_widths[[component]]

    # Stepping-out + shrinkage slice on log(theta[[component]]).
    x_cur  <- log(theta[[component]])
    log_y  <- current_lp + log(stats::runif(1))

    # Stepping-out.
    u <- stats::runif(1)
    L <- x_cur - w * u
    R <- L + w
    # Cap the number of stepping-out steps to keep cost bounded.
    max_steps <- 25L
    steps <- 0L
    repeat {
      th_try <- theta
      th_try[[component]] <- exp(L)
      mats_try <- build_kernels(th_try)
      if (log_target(th_try, mats_try) <= log_y || steps >= max_steps) break
      L <- L - w
      steps <- steps + 1L
    }
    steps <- 0L
    repeat {
      th_try <- theta
      th_try[[component]] <- exp(R)
      mats_try <- build_kernels(th_try)
      if (log_target(th_try, mats_try) <= log_y || steps >= max_steps) break
      R <- R + w
      steps <- steps + 1L
    }

    # Shrinkage.
    repeat {
      x_new <- stats::runif(1, L, R)
      th_try <- theta
      th_try[[component]] <- exp(x_new)
      mats_try <- build_kernels(th_try)
      lp_new  <- log_target(th_try, mats_try)
      if (lp_new > log_y) {
        theta      <- th_try
        mats       <- mats_try
        current_lp <- lp_new
        break
      }
      if (x_new < x_cur) L <- x_new else R <- x_new
    }
  }

  list(
    theta     = theta,
    space_mat = mats$space_mat,
    time_mat  = mats$time_mat,
    ke        = mats$ke,
    log_post  = current_lp
  )
}


# -----------------------------------------------------------------------------
# Block 5: r | y, psi -- random-walk MH on log r.
#
# Likelihood is the full NB log-likelihood at the observed cells:
#     log p(y | psi, r) = sum_i dnbinom(y_i, size = r, mu = exp(psi_i), log = TRUE)
#
# Proposal:  log r' = log r + N(0, s).  Adapted during burn-in via
# Robbins-Monro on log s targeting 0.44 acceptance (standard 1-D RW-MH).
#
# Note the jacobian: working in log r requires adding (log r' - log r) to the
# log acceptance ratio.
# -----------------------------------------------------------------------------
pg_draw_r <- function(state, design, r_prior, mh_state, adapt = TRUE) {
  if (!is.null(state$r_fixed) && isTRUE(state$r_fixed)) {
    return(list(r = state$r, mh_state = mh_state))
  }

  obs_idx <- design$obs_idx
  psi_obs <- state$f[obs_idx] + state$mu[design$site_idx_obs]

  log_r <- log(state$r)
  s     <- exp(mh_state$log_step)
  log_r_prop <- log_r + stats::rnorm(1, 0, s)
  r_prop     <- exp(log_r_prop)

  log_lik_cur  <- sum(stats::dnbinom(design$y_obs, size = state$r,
                                     mu = exp(psi_obs), log = TRUE))
  log_lik_prop <- sum(stats::dnbinom(design$y_obs, size = r_prop,
                                     mu = exp(psi_obs), log = TRUE))

  log_alpha <-
    log_lik_prop - log_lik_cur +
    r_prior$logp(r_prop) - r_prior$logp(state$r) +
    (log_r_prop - log_r)                                # jacobian for log r

  accept <- log(stats::runif(1)) < log_alpha
  new_r  <- if (accept) r_prop else state$r

  # Robbins-Monro adaptation during burn-in: nudge log_step toward target.
  if (adapt) {
    accept_prob <- min(1, exp(log_alpha))
    mh_state$log_step <- mh_state$log_step +
      mh_state$gamma * (accept_prob - 0.44)
    mh_state$gamma    <- 1 / (mh_state$iter + 1)
    mh_state$iter     <- mh_state$iter + 1L
  }

  mh_state$attempts <- mh_state$attempts + 1L
  if (accept) mh_state$accepts <- mh_state$accepts + 1L

  list(r = new_r, mh_state = mh_state, accepted = accept)
}


# -----------------------------------------------------------------------------
# One full PG-Gibbs sweep.
#
# Updates state in canonical order: omega, f, mu, theta, r. Returns the
# updated state together with diagnostics (PCG iters used, r-MH acceptance).
# -----------------------------------------------------------------------------
pg_sweep <- function(state, design, priors, control, adapt_r = TRUE) {

  # 1. omega | rest.
  psi_obs <- state$f[design$obs_idx] + state$mu[design$site_idx_obs]
  state$omega <- pg_draw_omega(psi_obs, design$y_obs, state$r)

  # 2. f | omega, ...
  f_out <- pg_draw_f(state, design, control)
  state$f <- f_out$f

  # 3. mu | omega, f, ...
  state$mu <- pg_draw_mu(state, design, priors$mu)

  # 4. theta | f.
  th_out <- pg_draw_theta(state, design, priors$theta, control$slice_widths)
  state$theta     <- th_out$theta
  state$space_mat <- th_out$space_mat
  state$time_mat  <- th_out$time_mat
  state$ke        <- th_out$ke

  # 5. r | y, psi.
  r_out <- pg_draw_r(state, design, priors$r, control$mh_state, adapt = adapt_r)
  state$r        <- r_out$r
  control$mh_state <- r_out$mh_state

  state$log_post <- th_out$log_post

  list(
    state   = state,
    control = control,
    diag    = list(
      pcg_iters   = f_out$pcg_iters,
      pcg_converged = f_out$pcg_converged,
      r_accepted  = r_out$accepted
    )
  )
}
