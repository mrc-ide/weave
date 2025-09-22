test_that("bounds returns quantiles for minimal state", {
  set.seed(123)

  obs_data <- data.frame(
    id = 1L,
    t = 1L,
    y_obs = 0,
    f_infer = 0,
    mu_infer = 0
  )

  coordinates <- data.frame(
    id = 1L,
    lon = 0,
    lat = 0
  )

  state <- gp_build_state(
    obs_data = obs_data,
    coordinates = coordinates,
    hyperparameters = c(1, 1, 1),
    n = 1,
    nt = 1,
    period = 52
  )

  result <- bounds(state, n_lambda = 2, n_draw = 2, quantiles = c(0.5))

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 1)
  expect_named(result, c("id", "t", "q0.5"))
  expect_equal(result$id, obs_data$id)
  expect_equal(result$t, obs_data$t)
})
