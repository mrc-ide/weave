library(testthat)
library(weave)

test_that("make_time_folds creates contiguous blocks", {
  df <- data.frame(t = rep(1:6, each = 2))
  folds <- make_time_folds(df, K = 3)
  expect_length(folds, 3)
  expect_true(all(vapply(folds, length, integer(1)) == 4))
  expect_equal(folds[[1]], 1:4)
})

test_that("make_time_folds_interleave alternates times", {
  df <- data.frame(t = rep(1:6, each = 2))
  folds <- make_time_folds_interleave(df, K = 3)
  times_in_folds <- lapply(folds, function(idx) unique(df$t[idx]))
  expect_equal(times_in_folds[[1]], c(1, 4))
  expect_equal(times_in_folds[[2]], c(2, 5))
  expect_equal(times_in_folds[[3]], c(3, 6))
})

test_that("tune_hyperparameters_optim returns list of results", {
  set.seed(1)
  n_sites <- 3
  nt <- 4
  ids <- rep(1:n_sites, each = nt)
  t <- rep(1:nt, times = n_sites)
  y <- rpois(n_sites * nt, lambda = 1)
  mu <- rep(0, n_sites * nt)
  f <- log1p(y) - mu
  obs_data <- data.frame(id = ids, t = t, y_obs = y, mu_infer = mu, f_infer = f)
  coordinates <- data.frame(id = 1:n_sites, lon = rnorm(n_sites), lat = rnorm(n_sites))

  res <- tune_hyperparameters_optim(
    obs_data = obs_data,
    coordinates = coordinates,
    n_sites_sample = 2,
    K_folds = 2,
    init = c(space = 1, t_per = 1, t_long = 1),
    lower = c(space = 0.1, t_per = 0.1, t_long = 0.1),
    upper = c(space = 5, t_per = 5, t_long = 5),
    seed = 1
  )

  expect_named(res, c("best_theta", "best_cv_score", "convergence", "message", "folds_info"))
  expect_length(res$best_theta, 3)
  expect_true(is.finite(res$best_cv_score))
})

test_that("fit requires coordinates for sampled sites", {
  set.seed(2)
  n_sites <- 3
  nt <- 4
  ids <- rep(seq_len(n_sites), each = nt)
  t <- rep(seq_len(nt), times = n_sites)
  y <- rpois(n_sites * nt, lambda = 1)
  mu <- rep(0, n_sites * nt)
  f <- log1p(y) - mu
  obs_data <- data.frame(id = ids, t = t, y_obs = y, mu_infer = mu, f_infer = f)
  coordinates <- data.frame(id = seq_len(n_sites), lon = rnorm(n_sites), lat = rnorm(n_sites))

  res <- fit(
    obs_data = obs_data,
    coordinates = coordinates,
    nt = nt,
    period = 52,
    n_sites = 2
  )
  expect_true(all(c("par", "value", "counts", "convergence", "message") %in% names(res)))

  missing_coords <- coordinates[-1, , drop = FALSE]
  expect_error(
    fit(
      obs_data = obs_data,
      coordinates = missing_coords,
      nt = nt,
      period = 52,
      n_sites = n_sites
    ),
    "Coordinates missing"
  )
})

