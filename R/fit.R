regularise <- function(x, lambda = 1e-5) {
  x + lambda * diag(nrow(x))
}


llh <- function(f, dist_k_inv, time_k_inv) {
  # Hessian of the quadratic term
  # The Hessian is -Sigma_inv, represented using the Kronecker product
  Hess_quad <- -kronecker(dist_k_inv, time_k_inv)

  # Initialize the Hessian of the likelihood term
  Hess_likelihood <- matrix(0, nrow = length(f), ncol = length(f))

  # Poisson distribution
  mu <- exp(f)  # Mean parameter
  tmp <- -mu  # Second derivative of Poisson log-likelihood
  Hess_likelihood <- diag(tmp)

  # Total Hessian
  Hess <- Hess_quad + Hess_likelihood

  # Return the Hessian matrix
  return(Hess)
}

fit <- function(obs_data, coordinates, hyperparameters, sub_n = 5, noise_var = 1e-3){
  time_matrix <- time_kernel(
    1:nt,
    periodic_scale = hyperparmameters[2],
    long_term_scale = hyperparmameters[3],
    period = 52
  )

  time_inv_reg <- time_matrix |>
    regularise() |>
    solve()

  spatial_dist <- get_spatial_distance(coordinates)
  space_matrix <- space_kernel(coordinates, length_scale = hyperparmameters[1])

  obs_data$z_est <- NA
  obs_data$tausq <- NA


  pb <- progress_bar$new(total = n)


  for(i in 1:n){
    pb$tick()
    dists <- spatial_dist[i,]
    indices <- order(dists)[1:(sub_n + 1)]
    sub_d <- dplyr::filter(obs_data, id %in% indices)

    sub_space_matrix <- space_matrix[indices, indices]

    sub_space_inv_reg <- sub_space_matrix |>
      regularise() |>
      solve()

    # Full Spatiotemporal Kernel (Kronecker Product)
    K_spacetime <- kronecker(sub_space_matrix, time_matrix)

    miss_index <- 1:nrow(sub_d)
    obs_index <- which(!is.na(sub_d$y_obs))

    y_obs <- log(sub_d$y_obs + 0.0001)[obs_index]
    y_obs <- y_obs - sub_d$mu_infer[obs_index]

    # 2) Subset into blocks for observed and missing
    K_obs_obs   <- K_spacetime[obs_index, obs_index]
    K_miss_obs  <- K_spacetime[miss_index, obs_index]

    # 3) GP posterior (Gaussian version)
    K_obs_obs_noisy <- K_obs_obs + diag(noise_var, nrow(K_obs_obs))

    # Solve for alpha
    alpha <- solve(K_obs_obs_noisy, y_obs)

    # Posterior mean for missing points
    post_mean <- K_miss_obs %*% alpha
    post_mean <- post_mean + sub_d$mu_infer

    obs_data[obs_data$id == i, "z_est"] <- post_mean[sub_d$id == i]

    # calculate Hessian at the ML value
    hess <- llh(
      f = post_mean,
      dist_k_inv = sub_space_inv_reg,
      time_k_inv = time_inv_reg
    )
    # Approximate the Hessian as a diagonal matrix
    hess_diag <- diag(hess)
    # Approximate the variances
    sub_tausq <- -1 / hess_diag
    sub_tausq[sub_tausq < 0] <- 0  # Ensure non-negative variances
    # Approximate the Hessian as a diagonal matrix
    hess_diag <- diag(hess)
    # Approximate the variances
    sub_tausq <- -1 / hess_diag
    sub_tausq[sub_tausq < 0] <- 0  # Ensure non-negative variances

    obs_data[obs_data$id == i, "tausq"] <- sub_tausq[sub_d$id == i]
  }

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
