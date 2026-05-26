test_that("kron_mv matches the explicit Kronecker matvec", {
  set.seed(1)
  n  <- 6
  nt <- 4
  S  <- crossprod(matrix(rnorm(n * n),  n))  + diag(n)
  T_ <- crossprod(matrix(rnorm(nt * nt), nt)) + diag(nt)
  v  <- rnorm(n * nt)

  ref <- as.vector(kronecker(S, T_) %*% v)
  expect_equal(kron_mv(v, S, T_), ref, tolerance = 1e-10)
})

test_that("kron_solve_eigen matches solve(K + sigma^2 I) for several sigma", {
  set.seed(2)
  n  <- 5
  nt <- 4
  S  <- crossprod(matrix(rnorm(n * n),  n))  + diag(n)
  T_ <- crossprod(matrix(rnorm(nt * nt), nt)) + diag(nt)
  v  <- rnorm(n * nt)
  ke <- kron_eigen(S, T_)
  K  <- kronecker(S, T_)

  for (s2 in c(0, 0.1, 1, 10)) {
    ref <- as.vector(solve(K + s2 * diag(n * nt), v))
    got <- kron_solve_eigen(v, ke, sigma2 = s2)
    expect_equal(got, ref, tolerance = 1e-8,
                 info = sprintf("sigma2 = %g", s2))
  }
})

test_that("kron_quad matches t(v) K^-1 v", {
  set.seed(3)
  n  <- 6
  nt <- 5
  S  <- crossprod(matrix(rnorm(n * n),  n))  + diag(n)
  T_ <- crossprod(matrix(rnorm(nt * nt), nt)) + diag(nt)
  v  <- rnorm(n * nt)
  ke <- kron_eigen(S, T_)
  K  <- kronecker(S, T_)

  expect_equal(kron_quad(v, ke),
               as.numeric(t(v) %*% solve(K, v)),
               tolerance = 1e-8)
})

test_that("log_det in the eigen cache matches the explicit determinant", {
  set.seed(4)
  n  <- 5
  nt <- 4
  S  <- crossprod(matrix(rnorm(n * n),  n))  + diag(n)
  T_ <- crossprod(matrix(rnorm(nt * nt), nt)) + diag(nt)
  ke <- kron_eigen(S, T_)
  K  <- kronecker(S, T_)
  expect_equal(ke$log_det,
               as.numeric(determinant(K, logarithm = TRUE)$modulus),
               tolerance = 1e-8)
})

test_that("regularise() uses a relative jitter", {
  set.seed(5)
  X <- crossprod(matrix(rnorm(16), 4))
  Y <- regularise(X, lambda = 1e-4)
  # The added value should equal lambda * mean(diag(X)) on each diagonal,
  # and zero off-diagonal.
  expect_equal(diag(Y) - diag(X),
               rep(1e-4 * mean(diag(X)), 4), tolerance = 1e-12)
  expect_equal(Y[upper.tri(Y)], X[upper.tri(X)], tolerance = 1e-12)
})
