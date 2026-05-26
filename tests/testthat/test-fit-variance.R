test_that("Hutchinson diag matches the exact posterior diagonal of the surrogate", {
  # Tiny full-data linearised problem: we can compute the exact
  # diag(K - K S' (S K S' + D)^-1 S K) for comparison.
  set.seed(41)
  n  <- 4; nt <- 5; N <- n * nt
  coords <- data.frame(
    id  = factor(1:n),
    lat = runif(n), lon = runif(n)
  )
  space_mat <- space_kernel(coords, length_scale = 1, nugget = 1e-4)
  time_mat  <- time_kernel(seq_len(nt), periodic_scale = 1,
                           long_term_scale = 30, period = 52, nugget = 1e-4)
  K  <- kronecker(space_mat, time_mat)

  obs_idx   <- seq_len(N)
  noise_var <- runif(N, 0.3, 0.7)
  Sigma_post <- K - K[, obs_idx] %*%
    solve(K[obs_idx, obs_idx] + diag(noise_var)) %*% K[obs_idx, ]
  diag_exact <- diag(Sigma_post)

  ke <- kron_eigen(space_mat, time_mat)
  diag_hat <- weave:::hutchinson_diag(
    obs_idx   = obs_idx,
    N         = N,
    space_mat = space_mat,
    time_mat  = time_mat,
    noise_var = noise_var,
    ke        = ke,
    n_probes  = 800,
    pcg_tol   = 1e-10,
    pcg_maxit = 500
  )
  # Stochastic estimator: tolerate non-trivial Monte-Carlo error on N=20
  # cells with 800 probes.
  expect_lt(max(abs(diag_hat - diag_exact)) / max(diag_exact), 0.2)
})
