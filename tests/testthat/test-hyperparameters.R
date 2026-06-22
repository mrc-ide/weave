test_that("gp_marginal_loglik matches a dense brute-force computation", {
  set.seed(1)
  n <- 4; nt <- 6; N <- n * nt
  coords <- data.frame(id = 1:n, lon = runif(n), lat = runif(n))
  times  <- seq_len(nt)

  Ks <- space_kernel(coords, length_scale = 1.5)
  Kt <- time_kernel(times, periodic_scale = 1, long_term_scale = 80, period = 52)
  eig_s <- eig_sym(Ks)
  eig_t <- eig_sym(Kt)

  g <- rnorm(N)

  for (eta in c(0.05, 0.5, 2)) {
    # Kronecker (sites slow, times fast) + scalar nugget.
    Bracket <- kronecker(Ks, Kt) + eta * diag(N)
    quad    <- as.numeric(t(g) %*% solve(Bracket, g))
    logdet  <- as.numeric(determinant(Bracket, logarithm = TRUE)$modulus)
    sigma2  <- quad / N
    ll_dense <- -0.5 * (N * log(2 * pi) + logdet + N * log(sigma2) + N)

    ll_fast <- gp_marginal_loglik(g, n, nt, eig_s, eig_t, eta = eta)

    expect_equal(as.numeric(ll_fast), ll_dense, tolerance = 1e-8)
    expect_equal(attr(ll_fast, "sigma2"), sigma2, tolerance = 1e-8)
  }
})


test_that("build_plugin_field has the right shape, ordering and NA handling", {
  n <- 3; nt <- 4
  obs <- expand.grid(t = 1:nt, id = 1:n)        # time fastest within site
  obs <- obs[order(obs$id, obs$t), ]
  obs$y_obs <- c(
    1, 2, 3, 4,        # site 1
    10, 20, NA, 40,    # site 2 (one missing)
    5, 5, 5, 5         # site 3 (constant -> sd 0)
  )

  g <- build_plugin_field(obs, n, nt)
  expect_length(g, n * nt)

  M <- t(matrix(g, nrow = nt, ncol = n))         # back to n x nt
  # Missing cell imputed to 0 (the per-site mean after centring).
  expect_equal(M[2, 3], 0)
  # Constant site -> centred to all zeros (sd guard prevents divide-by-zero).
  expect_true(all(M[3, ] == 0))
  # Observed cells of a varying site are finite and centred (mean ~ 0).
  expect_true(all(is.finite(M[1, ])))
  expect_equal(mean(M[1, ]), 0, tolerance = 1e-8)
})


test_that("refine validation and no-op behaviour", {
  set.seed(1)
  n <- 4; nt <- 10; period <- 52
  coords <- data.frame(id = 1:n, lon = runif(n), lat = runif(n))
  obs <- data.frame(id = rep(1:n, each = nt), t = rep(1:nt, n),
                    y_obs = stats::rpois(n * nt, 20))

  expect_error(
    infer_kernel_params(obs, coords, nt = nt, period = period, refine_iter = -1),
    "non-negative"
  )
  # refine = FALSE and refine_iter = 0 (refine = TRUE) both skip the loop.
  a <- infer_kernel_params(obs, coords, nt = nt, period = period, refine = FALSE)
  b <- infer_kernel_params(obs, coords, nt = nt, period = period,
                           refine = TRUE, refine_iter = 0)
  expect_equal(a$long_term_scale, b$long_term_scale)
})


test_that("refine reduces the gap-induced temporal attenuation", {
  # GP truth with a long temporal scale; clustered missingness attenuates it,
  # and refinement should pull the estimate back toward the no-gap fit.
  skip_on_cran()
  set.seed(42)
  n <- 10; nt <- 104; period <- 52
  coords <- data.frame(id = factor(1:n), lon = runif(n, 0, 5), lat = runif(n, 0, 5))
  sk <- space_kernel(coords, length_scale = 0.5)
  tk <- time_kernel(1:nt, periodic_scale = 2, long_term_scale = 80, period = period)
  f  <- quick_mvnorm(sk, tk)
  y  <- stats::rpois(n * nt, exp(3 + f))
  full <- data.frame(id = factor(rep(1:n, each = nt)), t = rep(1:nt, n), y_obs = y)

  # clustered missingness
  miss <- numeric(n * nt); miss[1] <- 0
  for (i in 2:(n * nt)) miss[i] <- if (runif(1) < 0.06) rbinom(1, 1, 0.15) else miss[i - 1]
  gap <- full; gap$y_obs[miss == 1] <- NA

  oracle <- infer_kernel_params(full, coords, nt = nt, period = period)
  naive  <- infer_kernel_params(gap,  coords, nt = nt, period = period)
  refined <- infer_kernel_params(gap, coords, nt = nt, period = period,
                                 refine = TRUE, refine_iter = 3)

  expect_equal(refined$convergence, 0)
  # gaps attenuate the long-term scale downward; refinement corrects it upward
  expect_gt(refined$long_term_scale, naive$long_term_scale)
  # and lands closer to the no-gap oracle than the naive (mean-imputed) fit
  expect_lt(abs(refined$long_term_scale - oracle$long_term_scale),
            abs(naive$long_term_scale - oracle$long_term_scale))
})


test_that("build_plugin_field errors when n/nt disagree with the data", {
  obs <- expand.grid(t = 1:4, id = 1:3)
  obs$y_obs <- 1
  expect_error(build_plugin_field(obs, n = 2, nt = 4), "do not match")
  expect_error(build_plugin_field(obs, n = 3, nt = 5), "do not match")
})


test_that("infer_kernel_params errors when coordinates miss a site", {
  set.seed(1)
  n <- 4; nt <- 6
  coords <- data.frame(id = 1:n, lon = runif(n), lat = runif(n))
  obs <- data.frame(id = rep(1:n, each = nt), t = rep(1:nt, n),
                    y_obs = stats::rpois(n * nt, 5))
  expect_error(
    infer_kernel_params(obs, coords[1:3, ], nt = nt, period = 52),
    "coordinates"
  )
})


test_that("infer_kernel_params uses real t spacing (gaps change the fit)", {
  set.seed(1)
  n <- 4; nt <- 8; period <- 52
  coords <- data.frame(id = 1:n, lon = runif(n), lat = runif(n))
  y <- stats::rpois(n * nt, 20)

  # Same counts, two time encodings. Under the old seq_len(nt) axis these were
  # identical fits; with the real-t axis the gap must shift the likelihood.
  even <- data.frame(id = rep(1:n, each = nt), t = rep(1:nt, n), y_obs = y)
  gap  <- data.frame(id = rep(1:n, each = nt), t = rep(c(1:7, 30), n), y_obs = y)

  e_even <- infer_kernel_params(even, coords, nt = nt, period = period)
  e_gap  <- infer_kernel_params(gap,  coords, nt = nt, period = period)

  expect_false(isTRUE(all.equal(e_even$log_posterior, e_gap$log_posterior)))
})


test_that("infer_kernel_params n_sites subsamples sites and is seed-reproducible", {
  n <- 8; nt <- 12; period <- 6
  coords <- data.frame(id = factor(1:n), lon = runif(n, 0, 5), lat = runif(n, 0, 5))
  obs <- expand.grid(t = seq_len(nt), id = factor(1:n))
  obs$y_obs <- stats::rpois(nrow(obs), lambda = 5)

  # Same seed -> same subsample -> identical estimate.
  set.seed(42); e1 <- infer_kernel_params(obs, coords, nt = nt, period = period, n_sites = 4)
  set.seed(42); e2 <- infer_kernel_params(obs, coords, nt = nt, period = period, n_sites = 4)
  expect_equal(e1$length_scale, e2$length_scale)
  expect_equal(e1$log_posterior, e2$log_posterior)

  # n_sites >= number of sites is a no-op (matches using all sites).
  e_all  <- infer_kernel_params(obs, coords, nt = nt, period = period)
  e_full <- infer_kernel_params(obs, coords, nt = nt, period = period, n_sites = n)
  expect_equal(e_all$log_posterior, e_full$log_posterior)
})


test_that("infer_kernel_params recovers known kernel params and the nugget helps", {
  skip_on_cran()
  set.seed(123)
  n <- 25; nt <- 104; period <- 52
  coords <- data.frame(id = 1:n, lon = runif(n, 0, 5), lat = runif(n, 0, 5))
  true_ls <- 2; true_ps <- 1.2; true_lts <- 150; true_r <- 15

  Ks <- space_kernel(coords, length_scale = true_ls)
  Kt <- time_kernel(1:nt, periodic_scale = true_ps, long_term_scale = true_lts,
                    period = period)
  f  <- quick_mvnorm(Ks, Kt)
  mu <- log(runif(n, 15, 70))
  psi <- f + rep(mu, each = nt)
  y   <- rnbinom(n * nt, size = true_r, mu = exp(psi))
  obs <- data.frame(id = rep(1:n, each = nt), t = rep(1:nt, n), y_obs = y)

  est <- infer_kernel_params(obs, coords, nt = nt, period = period)

  # With the nugget, the spatial length scale should land in a sensible band
  # around the truth (this is a plug-in estimate, so allow generous tolerance).
  expect_gt(est$length_scale, 1)
  expect_lt(est$length_scale, 3.5)
  expect_gt(est$nugget_ratio, 0)          # a real noise component is detected
  expect_equal(est$convergence, 0)

  # Pinning the nugget ~0 distorts the length scale badly -- the fit with the
  # nugget must be closer to the truth than the fit without it.
  priors0 <- default_kernel_priors()
  priors0$nugget_ratio <- list(meanlog = log(1e-8), sdlog = 1e-4)
  est0 <- infer_kernel_params(obs, coords, nt = nt, period = period,
                              priors = priors0,
                              start = c(1, 1, 100, 1e-8))

  expect_lt(abs(est$length_scale  - true_ls),
            abs(est0$length_scale - true_ls))
})
