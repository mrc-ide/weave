test_that("first row of time_cor corresponds to lag 0", {
  set.seed(123)
  nt <- 5
  n <- 3
  data <- data.frame(z_infer = rnorm(nt * n))
  res <- infer_time_kernel_params(data, period = 52, nt = nt, n = n, plot = FALSE)
  expect_identical(res$time_cor$time_distance[1], 0)
})
