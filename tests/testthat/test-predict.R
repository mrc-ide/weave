make_obs <- function(n, nt, missing = integer(0), seed = 2) {
  set.seed(seed)
  y <- stats::rpois(n * nt, 20)
  y[missing] <- NA
  data.frame(
    id = factor(rep(seq_len(n), each = nt), levels = seq_len(n)),
    t  = rep(seq_len(nt), times = n),
    y_obs = y
  )
}
hp_fixed <- list(length_scale = 1.5, periodic_scale = 1, long_term_scale = 80,
                 nugget_ratio = 0.3, sigma2 = 1)


test_that("gp_predict returns the expected shape and respects n_draws", {
  n <- 4; nt <- 6; period <- 52
  coords <- data.frame(id = factor(1:n), lon = runif(n), lat = runif(n))
  obs <- make_obs(n, nt, missing = c(3, 10, 15))

  out0 <- gp_predict(obs, coords, hp_fixed, nt = nt, period = period, n_draws = 0)
  expect_setequal(names(out0), c("id", "t", "rate"))
  expect_equal(nrow(out0), n * nt)
  expect_true(all(out0$rate > 0))
  expect_equal(attr(out0, "n_draws"), 0)

  out1 <- gp_predict(obs, coords, hp_fixed, nt = nt, period = period, n_draws = 20,
                     progress = FALSE)
  expect_true(all(c("lower", "upper") %in% names(out1)))
  expect_true(all(out1$lower <= out1$upper))
  expect_true(attr(out1, "r") > 0) # Poisson data -> r = Inf (valid: no overdispersion)

  # the mean is a deterministic single solve -> unaffected by the draw count
  expect_equal(out0$rate, out1$rate, tolerance = 1e-8)
})


test_that("gp_predict progress bar is cosmetic (same numbers, no error)", {
  n <- 4; nt <- 6; period <- 52
  coords <- data.frame(id = factor(1:n), lon = runif(n), lat = runif(n))
  obs <- make_obs(n, nt, missing = c(3, 10, 15))

  set.seed(7)
  quiet <- gp_predict(obs, coords, hp_fixed, nt = nt, period = period,
                      n_draws = 20, progress = FALSE)

  set.seed(7)
  expect_no_error(
    loud <- gp_predict(obs, coords, hp_fixed, nt = nt, period = period,
                       n_draws = 20, progress = TRUE)
  )

  # The bar is purely cosmetic: identical numeric output under the same seed.
  expect_equal(loud$lower, quiet$lower)
  expect_equal(loud$upper, quiet$upper)
  expect_equal(loud$rate, quiet$rate)
})


test_that("make_curve_bar draws a braille wave and tracks progress", {
  out <- utils::capture.output({
    pb <- make_curve_bar(total = 10, width = 12)
    for (i in 1:10) pb$tick()
    pb$done()
  })
  txt <- paste(out, collapse = "")

  expect_true(grepl("[⠀-⣿]", txt))   # contains braille glyphs
  expect_true(grepl("100%", txt))              # reaches 100%
  expect_true(grepl("10/10", txt))             # final n/n count

  # set() jumps to an arbitrary value without error
  expect_no_error(
    utils::capture.output({
      pb2 <- make_curve_bar(total = 100, width = 8)
      pb2$set(50)
      pb2$done()
    })
  )
})


test_that("gp_predict posterior mean matches a dense GP computation", {
  n <- 4; nt <- 6; period <- 52
  coords <- data.frame(id = factor(1:n), lon = runif(n), lat = runif(n))
  obs <- make_obs(n, nt, missing = c(3, 10, 15))

  out <- gp_predict(obs, coords, hp_fixed, nt = nt, period = period,
                    n_draws = 0, pcg_tol = 1e-11)

  # --- dense replication of the posterior-mean rate ------------------------
  ids <- sort(unique(obs$id)); times <- sort(unique(obs$t)); N <- n * nt
  coordsm <- coords[match(ids, coords$id), , drop = FALSE]
  M <- matrix(NA_real_, n, nt)
  M[cbind(match(obs$id, ids), match(obs$t, times))] <- log1p(obs$y_obs)
  row_mean <- rowMeans(M, na.rm = TRUE); row_mean[!is.finite(row_mean)] <- 0
  Mc <- M - row_mean
  row_sd <- apply(Mc, 1, stats::sd, na.rm = TRUE)
  row_sd[!is.finite(row_sd) | row_sd == 0] <- 1
  G <- Mc / row_sd; G[is.na(G)] <- 0

  space_mat <- hp_fixed$sigma2 * space_kernel(coordsm, length_scale = hp_fixed$length_scale)
  time_mat  <- time_kernel(times, periodic_scale = hp_fixed$periodic_scale,
                           long_term_scale = hp_fixed$long_term_scale, period = period)
  K  <- kronecker(space_mat, time_mat)
  nu <- hp_fixed$sigma2 * hp_fixed$nugget_ratio

  ok <- !is.na(obs$y_obs)
  og <- matrix(FALSE, n, nt)
  og[cbind(match(obs$id[ok], ids), match(obs$t[ok], times))] <- TRUE
  obs_idx <- which(as.vector(t(og)))
  g_obs <- as.vector(t(G))[obs_idx]

  alpha  <- solve(K[obs_idx, obs_idx] + nu * diag(length(obs_idx)), g_obs)
  f_mean <- as.vector(K[, obs_idx] %*% alpha)              # K S^T alpha
  Z <- row_mean + row_sd * t(matrix(f_mean, nrow = nt, ncol = n))
  rate_dense <- as.vector(t(exp(Z)))

  expect_equal(out$rate, rate_dense, tolerance = 1e-5)
})
