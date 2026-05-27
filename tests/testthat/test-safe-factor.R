test_that("safe_chol succeeds on a well-conditioned matrix", {
  set.seed(1)
  A <- crossprod(matrix(rnorm(36), 6)) + diag(6)
  L <- safe_chol(A)
  # L L^T should reconstruct A up to the tiny jitter we added.
  expect_lt(max(abs(L %*% t(L) - A)), 1e-4)
  # Lower-triangular: above-diagonal entries are zero.
  expect_true(all(L[upper.tri(L)] == 0))
})

test_that("safe_chol succeeds on a near-singular matrix", {
  # An "all-ones plus tiny diagonal" matrix is essentially rank-1 and is
  # the kind of structure that arises when length_scale wanders large
  # enough that K_space saturates at 1. safe_chol must succeed regardless.
  B <- matrix(1, 5, 5)
  diag(B) <- 1 + 1e-14
  L_safe <- safe_chol(B)
  expect_equal(dim(L_safe), c(5, 5))
  # Reconstruction should be close to the jittered matrix.
  recon <- L_safe %*% t(L_safe)
  expect_lt(max(abs(diag(recon) - diag(B))) / mean(diag(B)), 1e-3)
})

test_that("safe_eigen succeeds on a well-conditioned matrix and clamps eigenvalues", {
  set.seed(2)
  A <- crossprod(matrix(rnorm(36), 6)) + diag(6)
  out <- safe_eigen(A)
  expect_true(all(out$values > 0))
  # Reconstruction U diag(L) U^T ≈ A (up to jitter).
  recon <- out$vectors %*% diag(out$values) %*% t(out$vectors)
  expect_lt(max(abs(recon - A)), 1e-3)
})

test_that("safe_eigen clamps tiny eigenvalues at eig_floor", {
  # Rank-deficient: one direction has eigenvalue 0.
  M <- diag(c(1, 1, 1, 0))
  out <- safe_eigen(M, eig_floor = 1e-8)
  expect_true(all(out$values >= 1e-8))
})

test_that("safe_chol errors clearly if max_jitter is exhausted", {
  # Force failure by passing a matrix that no jitter level will fix.
  bad <- matrix(NA_real_, 3, 3)
  expect_error(safe_chol(bad), "safe_chol")
})
