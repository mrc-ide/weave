test_that("data completed correctly", {
  input_data <-
    data.frame(
      admin1 = rep("A", 6),
      admin2 = c(rep("A", 4), rep("B", 2)),
      t = c(1:4, 1:2),
      n = c(1, 2, 3, NA, 4, NA),
      lat = c(1, NA, NA, NA, 2, NA),
      lon = c(1, NA, NA, NA, 2, NA)
    )

  output_data <- data_complete(input_data, admin1, admin2)

  expect_identical(
    as.data.frame(output_data),
    data.frame(
      admin1 = rep("A", 8),
      admin2 = c(rep("A", 4), rep("B", 4)),
      t = rep(1:4, 2),
      n = c(1, 2, 3, NA, 4, NA, NA, NA),
      lat = c(rep(1, 4), rep(2, 4)),
      lon = c(rep(1, 4), rep(2, 4))
    )
  )
})

test_that("missing data dropped correctly", {
  input_data <-
    data.frame(
      admin1 = rep("A", 6),
      admin2 = c(rep("A", 4), rep("B", 2)),
      t = c(1:4, 1:2),
      n = c(1, 2, 3, NA, NA, NA),
      lat = c(1, NA, NA, NA, 2, NA),
      lon = c(1, NA, NA, NA, 2, NA)
    )

  output_data <- data_missing(input_data, admin1, admin2)

  expect_identical(
    as.data.frame(output_data),
    data.frame(
      admin1 = rep("A", 4),
      admin2 = rep("A", 4),
      t = 1:4,
      n = c(1, 2, 3, NA),
      lat = c(1, NA, NA, NA),
      lon = c(1, NA, NA, NA)
    )
  )
})

test_that("all-zero sites are retained by default and dropped with drop_zero", {
  input_zero <-
    data.frame(
      admin1 = rep("A", 4),
      admin2 = c(rep("A", 2), rep("B", 2)),
      t = c(1:2, 1:2),
      n = c(5, 7, 0, 0),
      lat = c(1, 1, 2, 2),
      lon = c(1, 1, 2, 2)
    )

  kept <- data_missing(input_zero, admin1, admin2)
  expect_true("B" %in% kept$admin2)

  dropped <- data_missing(input_zero, admin1, admin2, drop_zero = TRUE)
  expect_false("B" %in% dropped$admin2)
})

test_that("data_process validates required and protected columns", {
  ok <- data.frame(admin1 = "A", t = 1, n = 1, lat = 1, lon = 1)

  expect_error(
    data_process(ok[, c("admin1", "t", "n")], admin1),
    "must include"
  )
  expect_error(
    data_process(cbind(ok, id = 1), admin1),
    "protected"
  )
})

test_that("data processing pipeline returns a model-ready bundle", {
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

  expect_named(result, c("obs_data", "coordinates", "nt"))

  # obs_data: site keys + id + t + the renamed count column
  expect_identical(
    colnames(result$obs_data),
    c("admin1", "admin2", "id", "t", "y_obs")
  )
  expect_s3_class(result$obs_data$id, "factor")
  expect_identical(result$obs_data$t, rep(1:4, times = 2))

  # coordinates: one row per site, id/lon/lat only
  expect_identical(colnames(result$coordinates), c("id", "lon", "lat"))
  expect_identical(nrow(result$coordinates), length(unique(result$obs_data$id)))

  expect_identical(result$nt, 4L)
})
