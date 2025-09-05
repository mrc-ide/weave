test_that("rbf kernel returns correct values", {
  vec <- c(0, 1)
  res_vec <- rbf_kernel(vec, theta = 1)
  expect_equal(res_vec, exp(-vec^2 / 2))

  mat <- matrix(c(0, 1, 1, 0), nrow = 2)
  res_mat <- rbf_kernel(mat, theta = 1)
  expect_equal(dim(res_mat), c(2, 2))
  expect_equal(res_mat[1, 2], exp(-1 / 2))
})

test_that("periodic kernel returns correct values", {
  vec <- c(0, 1)
  res_vec <- periodic_kernel(vec, alpha = 1, period = 2)
  expect_equal(res_vec, c(1, exp(-2)))

  mat <- matrix(c(0, 1, 1, 0), nrow = 2)
  res_mat <- periodic_kernel(mat, alpha = 1, period = 2)
  expect_equal(dim(res_mat), c(2, 2))
  expect_equal(res_mat[1, 2], exp(-2))
})
