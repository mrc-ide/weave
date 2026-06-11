test_that("index assigned after ordering", {
  input_data <-
    data.frame(
      admin1 = rep("A", 6),
      admin2 = c(rep("A", 4), rep("B", 2)),
      t = as.integer(c(1, 3, 4, 2, 2, 1)),
      n = 1:6,
      lat = 1,
      lon = 1
    )

  result <- data_order_index(input_data, admin1, admin2)

  expect_s3_class(result$id, "factor")
  expect_identical(result$admin2, c(rep("A", 4), rep("B", 2)))
  expect_identical(result$t, c(1:4, 1:2))
  expect_identical(as.integer(result$id), as.integer(c(rep(1, 4), rep(2, 2))))
})
