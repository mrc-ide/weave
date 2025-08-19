test_that("mu can be scalar or vector", {
  n <- 3
  nt <- 2
  coordinates <- data.frame(
    id = factor(1:n),
    lat = 1:n,
    lon = 1:n
  )
  space_k <- diag(n)
  time_k <- diag(nt)

  set.seed(123)
  out_scalar <- simulate_data(
    n = n, nt = nt,
    coordinates = coordinates,
    space_k = space_k, time_k = time_k,
    mu = 1
  )
  expect_equal(unique(out_scalar$mu), rep(1, n))

  mu_vec <- seq_len(n)
  set.seed(123)
  out_vec <- simulate_data(
    n = n, nt = nt,
    coordinates = coordinates,
    space_k = space_k, time_k = time_k,
    mu = mu_vec
  )
  expect_equal(unique(out_vec$mu), mu_vec)
  expect_equal(out_vec$z, out_vec$mu + out_vec$f)

  expect_error(
    simulate_data(n, nt, coordinates, space_k, time_k, mu = 1:4),
    "`mu` must have length 1 or `n`"
  )
})
