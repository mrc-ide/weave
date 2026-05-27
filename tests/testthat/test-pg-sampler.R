test_that("Polya-Gamma draws have the correct first moment", {
  # Polson-Scott-Windle: E[PG(b, c)] = b/(2c) * tanh(c/2)
  skip_if_not_installed("BayesLogit")
  set.seed(31)
  b <- 3.5
  c_ <- 0.8
  draws <- BayesLogit::rpg(num = 5e4, h = b, z = c_)
  analytic <- (b / (2 * c_)) * tanh(c_ / 2)
  expect_equal(mean(draws), analytic, tolerance = 0.02)
})

test_that("pg_draw_f returns a draw from the analytic Gaussian conditional", {
  # Tiny full-data case with fixed (omega, mu, theta, r) so the conditional
  # is just a Gaussian. We check that the empirical mean and covariance from
  # many pg_draw_f calls match the closed-form conditional.
  skip_if_not_installed("BayesLogit")
  set.seed(32)
  n  <- 4; nt <- 3; N <- n * nt
  S  <- crossprod(matrix(rnorm(n * n),  n))  + diag(n)
  T_ <- crossprod(matrix(rnorm(nt * nt), nt)) + diag(nt)
  ke <- kron_eigen(S, T_)

  obs_idx <- seq_len(N)                            # full data
  omega   <- runif(N, 0.5, 1.5)
  y_obs   <- rep(2, N)
  r       <- 5
  mu      <- runif(n, 0.1, 0.3)
  state <- list(
    f = numeric(N), mu = mu, r = r, omega = omega,
    space_mat = S, time_mat = T_, ke = ke
  )
  design <- list(
    obs_idx = obs_idx, N = N, n = n, nt = nt,
    y_obs = y_obs, site_idx_obs = rep(seq_len(n), each = nt)
  )
  control <- list(pcg_tol = 1e-12, pcg_maxit = 500)

  # Analytic conditional.
  y_star  <- (y_obs - r) / (2 * omega) + log(r)
  mu_full <- mu[design$site_idx_obs]
  K       <- kronecker(S, T_)
  Sigma_post <- solve(solve(K) + diag(omega))
  mu_post    <- Sigma_post %*% (omega * (y_star - mu_full))

  M <- 2000
  draws <- matrix(NA_real_, nrow = M, ncol = N)
  for (k in seq_len(M)) {
    draws[k, ] <- pg_draw_f(state, design, control)$f
  }
  emp_mean <- colMeans(draws)
  emp_cov  <- stats::cov(draws)
  scale_mu <- max(abs(mu_post))
  expect_lt(max(abs(emp_mean - as.vector(mu_post))) / scale_mu, 0.05)
  expect_lt(max(abs(emp_cov - Sigma_post)), 0.1)
})

test_that("pg_draw_mu has the correct closed-form conditional", {
  # One site, fixed f and omega -- the conditional on mu_s is closed-form
  # Normal. We sample many times and check the moments.
  skip_if_not_installed("BayesLogit")
  set.seed(33)
  n  <- 1; nt <- 6; N <- n * nt
  omega  <- runif(N, 1, 2)
  y_obs  <- rep(3, N)
  r      <- 4
  f_full <- runif(N, -0.2, 0.2)
  mu_prior <- list(m0 = 0, v0 = 4)

  state <- list(f = f_full, mu = 0, r = r, omega = omega)
  design <- list(obs_idx = seq_len(N), n = n,
                 y_obs = y_obs, site_idx_obs = rep(1L, N))

  y_star <- (y_obs - r) / (2 * omega) + log(r)
  prec_post <- 1 / mu_prior$v0 + sum(omega)
  mean_post <- (mu_prior$m0 / mu_prior$v0 +
                  sum(omega * (y_star - f_full))) / prec_post

  M <- 5000
  draws <- replicate(M, pg_draw_mu(state, design, mu_prior))
  expect_lt(abs(mean(draws) - mean_post) * sqrt(prec_post * M), 4)
  expect_lt(abs(stats::var(draws) - 1 / prec_post) * prec_post, 0.1)
})


test_that("pg_draw_theta is deterministic given identical state and seed", {
  # The Cholesky-in-slice rewrite must preserve determinism: with the same
  # state and the same RNG seed it should return the same theta, same
  # log_post, and a refreshed eigen cache consistent with the new kernels.
  skip_if_not_installed("BayesLogit")
  set.seed(50)
  n  <- 5; nt <- 6; N <- n * nt
  coords <- data.frame(
    id  = factor(seq_len(n)),
    lat = runif(n), lon = runif(n)
  )
  state <- list(
    f         = rnorm(N) * 0.5,
    mu        = rep(0, n),
    theta     = list(length_scale = 1, periodic_scale = 1, long_term_scale = 20),
    r         = 10,
    omega     = rep(1, N),
    r_fixed   = FALSE,
    space_mat = space_kernel(coords, length_scale = 1),
    time_mat  = time_kernel(seq_len(nt), periodic_scale = 1,
                            long_term_scale = 20, period = 52)
  )
  state$ke <- kron_eigen(state$space_mat, state$time_mat)
  design <- list(
    coords = coords, n = n, nt = nt, period = 52, N = N,
    obs_idx = seq_len(N), y_obs = rep(0, N), site_idx_obs = rep(seq_len(n), each = nt)
  )
  priors <- bayes_priors(c(design, list(y_obs = rep(1, N))))
  slice_widths <- list(length_scale = 0.5, periodic_scale = 0.5,
                       long_term_scale = 0.5)

  set.seed(99)
  out1 <- pg_draw_theta(state, design, priors$theta, slice_widths)
  set.seed(99)
  out2 <- pg_draw_theta(state, design, priors$theta, slice_widths)

  expect_equal(out1$theta,    out2$theta)
  expect_equal(out1$log_post, out2$log_post)

  # Eigen cache should reflect the returned kernels: rebuilding from
  # out1$space_mat / out1$time_mat must give the same eigenvalues.
  ke_recheck <- kron_eigen(out1$space_mat, out1$time_mat)
  expect_equal(sort(out1$ke$L_s), sort(ke_recheck$L_s), tolerance = 1e-8)
  expect_equal(sort(out1$ke$L_t), sort(ke_recheck$L_t), tolerance = 1e-8)
})
