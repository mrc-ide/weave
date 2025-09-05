test_that("observed summary calculated correctly", {
  input_data <-
    data.frame(
      admin1 = rep("A", 4),
      admin2 = rep("A", 4),
      t = 1:4,
      n = c(1, 2, NA, 4),
      lat = 1,
      lon = 1
    )

  result <- data_observed_summary(input_data, admin1, admin2)

  exp_mu <- log(mean(input_data$n, na.rm = TRUE))
  exp_sigma <- log((sd(input_data$n, na.rm = TRUE) / mean(input_data$n, na.rm = TRUE))^2 + 1)

  expect_equal(unique(result$observed_mu), exp_mu)
  expect_equal(unique(result$observed_sigmasq), exp_sigma)
})

test_that("initial parameters derived from counts", {
  input_data <-
    data.frame(
      admin1 = rep("A", 4),
      admin2 = rep("A", 4),
      t = 1:4,
      n = c(1, 2, NA, 4),
      lat = 1,
      lon = 1
    ) |> 
    data_observed_summary(admin1, admin2)

  result <- data_initial_par(input_data)

  exp_start <- ifelse(is.na(input_data$n), input_data$observed_mu, log1p(input_data$n))

  expect_equal(result$start_par, exp_start)
})

test_that("index assigned after ordering", {
  input_data <-
    data.frame(
      admin1 = rep("A", 6),
      admin2 = c(rep("A", 4), rep("B", 2)),
      t = c(1, 3, 4, 2, 2, 1),
      n = 1:6,
      lat = 1,
      lon = 1
    )

  result <- data_order_index(input_data, admin1, admin2)

  expect_s3_class(result$id, "factor")
  expect_identical(result$admin2, c(rep("A", 4), rep("B", 2)))
  expect_identical(result$t, c(1:4, 1:2))
  expect_identical(as.integer(result$id), c(rep(1, 4), rep(2, 2)))
})

test_that("data processing pipeline builds structure", {
  input_data <-
    data.frame(
      admin1 = rep("A", 6),
      admin2 = c(rep("A", 4), rep("B", 2)),
      t = c(1:4, 1:2),
      n = c(1, 2, 3, NA, 4, NA),
      lat = c(1, NA, NA, NA, 2, NA),
      lon = c(1, NA, NA, NA, 2, NA)
    )

  result <- data_process(input_data, admin1, admin2)

  expect_identical(
    colnames(result),
    c(
      "admin1", "admin2", "t", "n", "lat", "lon",
      "observed_mu", "observed_sigmasq", "start_par", "id"
    )
  )
  expect_s3_class(result$id, "factor")
  expect_identical(result$admin2, c(rep("A", 4), rep("B", 4)))
  expect_identical(result$t, rep(1:4, times = 2))
  expect_identical(as.integer(result$id), c(rep(1, 4), rep(2, 4)))

  exp_start <- ifelse(is.na(result$n), result$observed_mu, log1p(result$n))
  expect_equal(result$start_par, exp_start)
})
