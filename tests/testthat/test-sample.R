test_that("quick_mvnorm and quick_mvnorm_chol agree and have the right length", {
  set.seed(1)
  n <- 4; nt <- 5
  coords <- data.frame(id = 1:n, lon = runif(n), lat = runif(n))
  space <- space_kernel(coords, length_scale = 1.5)
  time  <- time_kernel(seq_len(nt), periodic_scale = 1, long_term_scale = 80,
                       period = 52)

  # Same RNG state -> the precomputed-Cholesky variant must reproduce the other.
  set.seed(99); a <- quick_mvnorm(space, time)
  set.seed(99); b <- quick_mvnorm_chol(chol(space), chol(time))

  expect_length(a, n * nt)
  expect_equal(a, b)
})


test_that("quick_mvnorm draws have the separable Kronecker covariance", {
  set.seed(2)
  n <- 3; nt <- 4
  coords <- data.frame(id = 1:n, lon = runif(n), lat = runif(n))
  space <- space_kernel(coords, length_scale = 1.2)
  time  <- time_kernel(seq_len(nt), periodic_scale = 1, long_term_scale = 50,
                       period = 12)

  draws <- replicate(20000, quick_mvnorm(space, time))   # (n*nt) x reps
  emp   <- stats::cov(t(draws))
  target <- kronecker(space, time)                        # times fastest

  # Monte-Carlo, so allow a generous tolerance on the empirical covariance.
  expect_equal(emp, target, tolerance = 0.05, ignore_attr = TRUE)
})
