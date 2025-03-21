log_likelihood_hyperparameters <- function(pars, f, n, nt, coordinates) {
  space_k <- space_kernel(coordinates, length_scale = pars[1])
  time_k <- time_kernel(times = 1:nt, periodic_scale = pars[2], long_term_scale = pars[3])

  reg_space <- regularise(space_k)
  reg_time <- regularise(time_k)

  space_inv_reg <- solve(reg_space)
  time_inv_reg <- solve(reg_time)

  f_mat <- matrix(f, nrow = n, ncol = nt, byrow = TRUE)
  quad_term <- -0.5 * sum((space_inv_reg %*% f_mat) * (f_mat %*% time_inv_reg))

  # Compute log-determinants of the regularised space and time kernels
  log_det_space <- as.numeric(determinant(reg_space, logarithm = TRUE)$modulus)
  log_det_time  <- as.numeric(determinant(reg_time, logarithm = TRUE)$modulus)

  # Log-determinant term: accounts for the normalization of the Gaussian density
  log_det_term <- -0.5 * (nt * log_det_space + n * log_det_time)

  # Return the sum of the terms
  quad_term + log_det_term
}
