test_that("pcg matches solve() with the Kron-eigen preconditioner", {
  set.seed(11)
  n  <- 10
  nt <- 6
  S  <- crossprod(matrix(rnorm(n * n),  n))  + diag(n)
  T_ <- crossprod(matrix(rnorm(nt * nt), nt)) + diag(nt)
  ke <- kron_eigen(S, T_)
  N  <- n * nt
  obs_idx   <- sort(sample.int(N, floor(0.7 * N)))
  m         <- length(obs_idx)
  noise_var <- runif(m, 0.5, 1.5)
  b         <- rnorm(m)

  Amv_fun <- function(v) Amv(v, obs_idx, N, S, T_, noise_var)
  Minv    <- kron_eigen_preconditioner(ke, sigma2 = mean(noise_var),
                                       obs_idx, N)
  res <- pcg(b, Amv_fun, Minv, tol = 1e-10, maxit = 200)
  expect_true(res$converged)
  expect_lt(res$iters, 100)

  K  <- kronecker(S, T_)
  A  <- K[obs_idx, obs_idx] + diag(noise_var)
  expect_equal(res$x, as.vector(solve(A, b)), tolerance = 1e-7)
})

test_that("pcg's maxit warning only fires when convergence was not achieved", {
  set.seed(12)
  n  <- 6
  nt <- 4
  S  <- crossprod(matrix(rnorm(n * n),  n))  + diag(n)
  T_ <- crossprod(matrix(rnorm(nt * nt), nt)) + diag(nt)
  N  <- n * nt
  obs_idx <- seq_len(N)
  noise_var <- rep(1, N)
  b <- rnorm(N)
  Amv_fun <- function(v) Amv(v, obs_idx, N, S, T_, noise_var)

  # Comfortable tol, plenty of iters: no warning.
  expect_no_warning({
    res_ok <- pcg(b, Amv_fun, tol = 1e-6, maxit = 500)
  })
  expect_true(res_ok$converged)

  # Tight tol, tiny maxit: must warn.
  expect_warning(pcg(b, Amv_fun, tol = 1e-30, maxit = 2), "converge")
})
