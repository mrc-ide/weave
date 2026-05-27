# Each invalid input to fit_bayes() must error with a specific,
# informative message rather than crashing mid-loop.
tiny <- data.frame(
  id = factor(c(1, 1, 2, 2)),
  t  = c(1, 2, 1, 2),
  lat = c(0, 0, 1, 1),
  lon = c(0, 0, 1, 1),
  y_obs = c(3L, 5L, 2L, 1L)
)

test_that("obs_data must be a data frame", {
  expect_error(fit_bayes(list()), "must be a data frame")
})

test_that("obs_data must have id/t/lat/lon columns", {
  bad <- tiny; bad$lat <- NULL
  expect_error(fit_bayes(bad), "missing required column")
})

test_that("obs_data must have a y_obs or n column", {
  bad <- tiny; bad$y_obs <- NULL
  expect_error(fit_bayes(bad), "y_obs.* or .*n")
})

test_that("counts must be non-negative integers", {
  bad <- tiny; bad$y_obs <- c(1.5, 2, 3, 4)
  expect_error(fit_bayes(bad), "non-negative integers")
  bad2 <- tiny; bad2$y_obs <- c(-1L, 0L, 1L, 2L)
  expect_error(fit_bayes(bad2), "non-negative integers")
})

test_that("burnin must be less than n_sweeps", {
  expect_error(
    fit_bayes(tiny, n_sweeps = 10, burnin = 10, period = 2),
    "must be < `n_sweeps`",
    fixed = TRUE
  )
})

test_that("n_chains must be a positive integer", {
  expect_error(fit_bayes(tiny, n_sweeps = 10, burnin = 2, n_chains = 0,
                         period = 2), "n_chains")
  expect_error(fit_bayes(tiny, n_sweeps = 10, burnin = 2, n_chains = -1,
                         period = 2), "non-negative integer")
})

test_that("fix$r must be a positive finite number", {
  expect_error(
    fit_bayes(tiny, n_sweeps = 10, burnin = 2, fix = list(r = 0),
              period = 2),
    "fix$r", fixed = TRUE
  )
  expect_error(
    fit_bayes(tiny, n_sweeps = 10, burnin = 2, fix = list(r = -3),
              period = 2),
    "fix$r", fixed = TRUE
  )
  expect_error(
    fit_bayes(tiny, n_sweeps = 10, burnin = 2, fix = list(r = Inf),
              period = 2),
    "fix$r", fixed = TRUE
  )
})

test_that("period must be <= nt", {
  expect_error(
    fit_bayes(tiny, n_sweeps = 10, burnin = 2, period = 99),
    "must be <= nt", fixed = TRUE
  )
})

test_that("slice_widths must have all three components positive", {
  bad <- list(length_scale = 0.5, periodic_scale = 0.5)
  expect_error(
    fit_bayes(tiny, n_sweeps = 10, burnin = 2, period = 2,
              slice_widths = bad),
    "slice_widths"
  )
  bad2 <- list(length_scale = -1, periodic_scale = 0.5, long_term_scale = 0.5)
  expect_error(
    fit_bayes(tiny, n_sweeps = 10, burnin = 2, period = 2,
              slice_widths = bad2),
    "positive number"
  )
})
