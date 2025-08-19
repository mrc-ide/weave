test_that("infer_space_kernel_params requires z_infer", {
  data <- data.frame(id = 1, t = 1, lat = 0, lon = 0)
  expect_error(infer_space_kernel_params(data), "z_infer")
})

test_that("infer_time_kernel_params requires z_infer", {
  data <- data.frame(id = 1, t = 1, lat = 0, lon = 0)
  expect_error(infer_time_kernel_params(data, period = 1), "z_infer")
})

test_that("infer parameter functions derive counts from data", {
  set.seed(1)
  data <- data.frame(
    id = rep(1:2, each = 3),
    t = rep(1:3, times = 2),
    lon = rep(c(0, 1), each = 3),
    lat = rep(c(0, 1), each = 3),
    z_infer = rnorm(6)
  )
  space <- infer_space_kernel_params(data)
  time <- infer_time_kernel_params(data, period = 3)
  expect_named(space, "length_scale")
  expect_named(time, c("periodic_scale", "long_term_scale", "period"))
})
