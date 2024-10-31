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
