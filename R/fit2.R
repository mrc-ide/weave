# Fast matvec for K = space ⊗ time, with times varying fastest
kron_mv <- function(v, space, time) {
  n_sites <- nrow(space); n_times <- nrow(time)
  # Reconstruct X (sites × times) from vec(t(X)) = v
  X <- t(matrix(v, nrow = n_times, ncol = n_sites))
  # Apply K to X: vec(t(space %*% X %*% t(time))) = (space ⊗ time) vec(t(X))
  Y <- space %*% X %*% t(time)
  as.vector(t(Y))
}

with_nas <- function(x_obs, obs_idx, N){
  v <- numeric(N)
  v[obs_idx] <- x_obs
  v
}

# Define the observed-system matvec: (S K S^T + noise * I) v
Amv <- function(v, obs_idx, N, space_mat, time_mat, noise_var) {
  #browser()
  kron_mv(with_nas(v, obs_idx, N), space_mat, time_mat)[obs_idx] + noise_var * v
}

# Apply M^{-1} v ≈ v / (diag(K_oo) [+ noise]); +1e-12 avoids divide-by-zero.
# (If using per-observation noise_var, prefer: v / (kdiag_full[obs_idx] + noise_var + 1e-12))
M_inv <- function(v, kdiag_full, obs_idx, noise_var) {
  v / (kdiag_full[obs_idx] + noise_var + 1e-12)
}

pcg <- function(b, obs_idx, N, space_mat, time_mat, noise_var, kdiag_full, tol = 1e-8, maxit = 10000) {
  x <- numeric(length(b))
  r <- b - Amv(x, obs_idx, N, space_mat, time_mat, noise_var)
  z <- M_inv(r, kdiag_full, obs_idx, noise_var)
  p <- z
  rz_old <- sum(r * z)
  for (it in seq_len(maxit)) {
    if(it == maxit){
      warning("maxit reached")
    }
    Ap <- Amv(p, obs_idx, N, space_mat, time_mat, noise_var)
    alpha <- rz_old / sum(p * Ap)
    x <- x + alpha * p
    r <- r - alpha * Ap
    if (sqrt(sum(r * r)) <= tol * sqrt(sum(b * b))) break
    z <- M_inv(r, kdiag_full, obs_idx, noise_var)
    rz_new <- sum(r * z)
    beta <- rz_new / rz_old
    p <- z + beta * p
    rz_old <- rz_new
  }
  x
}

fit2 <- function(obs_data, coordinates, hyperparameters, n, nt){
  # Build kernels
  time_mat  <- time_kernel(
    times = 1:nt,
    periodic_scale = hyperparameters[2],
    long_term_scale = hyperparameters[3],
    period = 52
  )
  space_mat <- space_kernel(
    coordinates = coordinates, length_scale = hyperparameters[1]
  )

  obs_idx <- which(!is.na(obs_data$y_obs))  # ensure ordering matches vec stacking
  N <- n * nt

  # Centered observed response (log with eps, minus mu_infer)
  y_obs <- obs_data$f_infer[obs_idx]

  # Heteroscedastic nugget on the log scale:
  ## y_obs = log(Y+1) − mu_infer is noisy (Poisson + log). By the delta method,
  ## Var[log(Y+1) | λ] ≈ λ/(λ+1)^2. Using λ̂_i = exp(mu_infer_i) gives
  ## noise_var_i = λ̂_i/(λ̂_i+1)^2, which stabilises (S K S^T + diag(noise_var))
  ## and avoids overfitting log-count noise.
  lam_hat <- exp(obs_data$mu_infer[obs_idx])
  noise_var <- lam_hat / (lam_hat + 1)^2

  # Jacobi (diagonal) preconditioner:
  ## For A = S K S^T [+ diag(noise)], approximate A^{-1} with 1/diag(A).
  ## Since diag(K) = diag(space) ⊗ diag(time), build it once, then select observed entries.
  kdiag_full <- as.vector(kronecker(diag(space_mat), diag(time_mat)))

  # Solve for α using PCG:
  ## α are the GP conditioning weights solving (S K S^T + diag(noise_var)) α = y_obs.
  ## They map observations to the posterior mean via f_hat = K S^T α.
  ## PCG (Preconditioned Conjugate Gradient) is an iterative solver that needs only
  ## matrix–vector products with A (no full matrix), using M_inv as a diagonal
  ## preconditioner to speed convergence.
  alpha <- pcg(y_obs, obs_idx, N, space_mat, time_mat, noise_var, kdiag_full, tol = 1e-6)

  # Posterior mean at all entries: f_hat = K S^T alpha  (then add mu)
  f_hat <- kron_mv(with_nas(alpha, obs_idx, N), space_mat, time_mat)

  # Write back z_est per row
  obs_data$z_est <- f_hat + obs_data$mu_infer

  # Variance (Hessian diag at ML): use your efficient diag trick, global form
  space_inv <- chol2inv(chol(regularise(space_mat)))
  time_inv  <- chol2inv(chol(regularise(time_mat)))
  kron_diag_inv <- as.vector(kronecker(diag(space_inv), diag(time_inv)))
  hess_diag <- -kron_diag_inv - exp(obs_data$z_est)  # evaluated at f_hat
  obs_data$tausq <- pmax(0, -1 / hess_diag)

  fit_data <- obs_data |>
    dplyr::mutate(
      # Overdisperson estimates, by site
      mean_lambda = mean(exp(z_est)),
      var_y = var(y_obs, na.rm = TRUE),
      overdispersion = (mean_lambda^2) / (var_y - mean_lambda),
      .by = id
    ) |>
    dplyr::mutate(
      # Compute the 95% confidence interval for z
      z_min = z_est - 1.96 * sqrt(tausq),
      z_max = z_est + 1.96 * sqrt(tausq),

      # Convert this uncertainty to the count scale. We cannot simply exponentiate
      # z_min and z_max because exponentiation is nonlinear, and normal confidence
      # intervals don’t transform correctly.
      # Instead, we compute quantiles of the corresponding lognormal distribution.
      lambda_est = qlnorm(0.5, meanlog = z_est, sdlog = sqrt(tausq)),
      lambda_min = qlnorm(0.025, meanlog = z_est, sdlog = sqrt(tausq)),
      lambda_max = qlnorm(0.975, meanlog = z_est, sdlog = sqrt(tausq)),

      # Prediction intervals assuming Negative Binomially distributed data.
      # The Negative Binomial accounts for overdispersion beyond Poisson variability.
      negbin_size = lambda_est / overdispersion,  # Size parameter for NB
      negbin_prob = negbin_size / (negbin_size + lambda_est), # NB probability

      # Compute quantiles of the Negative Binomial distribution as prediction intervals.
      pred_Q2.5 = qnbinom(0.025, size = negbin_size, prob = negbin_prob),
      pred_Q25 = qnbinom(0.25, size = negbin_size, prob = negbin_prob),
      data_Q50 = qnbinom(0.5, size = negbin_size, prob = negbin_prob),
      pred_Q75 = qnbinom(0.75, size = negbin_size, prob = negbin_prob),
      pred_Q97.5 = qnbinom(0.975, size = negbin_size, prob = negbin_prob)
    ) |>
    dplyr::mutate(
      # Surprisal (information theory)
      log_p = dnbinom(y_obs, size = negbin_size, prob = negbin_prob, log = TRUE),
      y_mode = floor((negbin_size - 1) * (1 - negbin_prob) / negbin_prob),
      log_p_mode = dnbinom(y_mode,  size = negbin_size, prob = negbin_prob, log = TRUE),
      surprisal = 1 - exp(log_p - log_p_mode) # 0: most expected, 1: highly surprising
    )

  return(fit_data)
}
