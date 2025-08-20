# Fast matvec for K = space ⊗ time, with times varying fastest
kron_mv <- function(v, space, time) {
  n_sites <- nrow(space); n_times <- nrow(time)
  # Reconstruct X (sites × times) from vec(t(X)) = v
  X <- t(matrix(v, nrow = n_times, ncol = n_sites))
  # Apply K to X: vec(t(space %*% X %*% t(time))) = (space ⊗ time) vec(t(X))
  Y <- space %*% X %*% t(time)
  as.vector(t(Y))
}

# Selection helpers for masking observed entries
make_select <- function(obs_idx, N) {
  # obs_idx: integer positions (1..N) of observed entries in vec order (times fastest)
  list(
    S   = function(x_full) x_full[obs_idx],
    ST  = function(x_obs) { v <- numeric(N); v[obs_idx] <- x_obs; v }
  )
}

# Minimal (preconditioned) Conjugate Gradient for A x = b, where A is given as a function
pcg <- function(Amv, b, M_inv = NULL, tol = 1e-8, maxit = 10000) {
  x <- numeric(length(b))
  r <- b - Amv(x)
  z <- if (is.null(M_inv)) r else M_inv(r)
  p <- z
  rz_old <- sum(r * z)
  for (it in seq_len(maxit)) {
    if(it == maxit){
      warning("maxit reached")
    }
    Ap <- Amv(p)
    alpha <- rz_old / sum(p * Ap)
    x <- x + alpha * p
    r <- r - alpha * Ap
    if (sqrt(sum(r * r)) <= tol * sqrt(sum(b * b))) break
    z <- if (is.null(M_inv)) r else M_inv(r)
    rz_new <- sum(r * z)
    beta <- rz_new / rz_old
    p <- z + beta * p
    rz_old <- rz_new
  }
  x
}


fit2 <- function(obs_data, coordinates, hyperparameters){
  # Build kernels once
  time_mat  <- time_kernel(1:nt, periodic_scale = hyperparameters[2],
                           long_term_scale = hyperparameters[3], period = 52)
  space_mat <- space_kernel(coordinates, length_scale = hyperparameters[1])

  # Indices of observed entries in your vec order (times fastest within site)
  # Suppose obs_data has columns id (1..n), t (1..nt), y_obs
  obs_idx <- which(!is.na(obs_data$y_obs))  # ensure ordering matches vec stacking
  N <- nrow(space_mat) * nrow(time_mat)

  sel <- make_select(obs_idx, N)

  # Centered observed response (log with eps, minus mu_infer)
  y_obs <- log(obs_data$y_obs[obs_idx] + 1e-4) - obs_data$mu_infer[obs_idx]

  # Define the observed-system matvec: (S K S^T + noise * I) v
  A_mv <- function(v) {
    sel$S(kron_mv(sel$ST(v), space_mat, time_mat)) + noise_var * v
  }

  noise_var =
    # (Optional but helpful) diagonal preconditioner from diag(K_oo) + noise
    # diag(K) = kron(diag(space), diag(time)); pick observed entries
    kdiag_full <- as.vector(kronecker(diag(space_mat), diag(time_mat)))
  M_inv <- function(v) v / (kdiag_full[obs_idx] + 1e-12)

  # Solve for alpha using PCG
  alpha <- pcg(A_mv, y_obs, M_inv = M_inv, tol = 1e-6)

  # Posterior mean at all entries: f_hat = K S^T alpha  (then add mu)
  f_hat <- kron_mv(sel$ST(alpha), space_mat, time_mat)

  # Write back z_est per row (times-fastest vec matches your (id,t) order)
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
#
# fit_data2 <- fit2(obs_data, coordinates, hyperparameters)

# fit_plot <- sim_data +
#   geom_ribbon(
#     data = fit_data2,
#     aes(x = t, ymin = pred_Q2.5, ymax = pred_Q97.5, fill = id), alpha = 0.25
#   ) +
#   geom_ribbon(
#     data = fit_data2,
#     aes(x = t, ymin = pred_Q25, ymax = pred_Q75, fill = id, alpha = 0.5)
#   ) +
#   geom_line(
#     data = fit_data2,
#     aes(x = t, y = data_Q50, col = id), linewidth = 1
#   ) +
#   ggtitle("Filling missingness")
# fit_plot
# outlier_pd <- fit_data |>
#   mutate(
#     size = ifelse(surprisal < 0.8, 1, 5)
#   )


