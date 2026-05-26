test_that("build_design exposes the expected fields", {
  set.seed(21)
  n  <- 4; nt <- 5
  d <- tidyr::expand_grid(
    id = factor(1:n),
    t  = 1:nt
  ) |>
    dplyr::mutate(
      lat   = as.integer(id),
      lon   = -as.integer(id),
      y_obs = ifelse(t == 1, NA_integer_, sample(0:50, dplyr::n(), replace = TRUE))
    )

  des <- build_design(d)

  expect_named(des, c("y_full", "y_obs", "obs_idx",
                      "n", "nt", "N",
                      "site_idx_full", "site_idx_obs",
                      "time_idx_full", "coords", "mu_init"))
  expect_equal(des$n,  n)
  expect_equal(des$nt, nt)
  expect_equal(des$N,  n * nt)
  expect_equal(length(des$obs_idx), sum(!is.na(d$y_obs)))
  expect_equal(des$y_obs, des$y_full[des$obs_idx])
  expect_equal(des$site_idx_obs, des$site_idx_full[des$obs_idx])
  # Times must vary fastest within site.
  expect_true(!is.unsorted(des$site_idx_full))
  for (s in seq_len(des$n)) {
    expect_true(!is.unsorted(des$time_idx_full[des$site_idx_full == s]))
  }
  # One row per site in `coords`.
  expect_equal(nrow(des$coords), des$n)
})

test_that("build_design accepts the `n` (real-data) column convention", {
  d <- tidyr::expand_grid(id = factor(1:3), t = 1:4) |>
    dplyr::mutate(lat = 1, lon = 1, n = c(rep(NA, 4), 1:8))
  des <- build_design(d)
  expect_equal(length(des$obs_idx), 8)
})

test_that("build_design errors loudly on incomplete grids", {
  d <- data.frame(id = factor(c(1, 1, 2)), t = c(1, 2, 1),
                  lat = 0, lon = 0, y_obs = c(1, 2, 3))
  expect_error(build_design(d), "complete site")
})
