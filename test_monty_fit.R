library(monty)

ll <- function(length_scale, periodic_scale, long_term_scale, f, n, nt, coordinates) {
  space_k <- space_kernel(coordinates, length_scale = length_scale)
  time_k <- time_kernel(times = 1:nt, periodic_scale = periodic_scale, long_term_scale = long_term_scale)

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

packer <- monty::monty_packer(
  scalar = c("length_scale", "periodic_scale", "long_term_scale"),
  fixed = list(
    "f" = obs_data$z_infer,
    "n" = n,
    "nt" = nt,
    "coordinates" = coordinates
  )
)

likelihood <- monty_model_function(density = ll, packer = packer)
likelihood
prior <- monty_dsl({
  length_scale ~ Gamma(shape = 2 , rate = 2)
  periodic_scale ~ Uniform(0.8, 1.2)
  long_term_scale ~ Gamma(shape = 10 , rate = 0.05)
})
posterior <- likelihood + prior
posterior
vcv <- diag(3) * c(0.5, 0.5, 100)
sampler <- monty_sampler_random_walk(vcv = vcv)
samples <- monty_sample(
  model = posterior,
  sampler = sampler,
  n_steps = 1000,
  initial = c(1, 1, 300),
  n_chains = 4
)
draws <- posterior::as_draws_df(samples)
posterior::summarise_draws(draws)
bayesplot::mcmc_trace(draws)
bayesplot::mcmc_dens(draws)
bayesplot::mcmc_scatter(draws, c("periodic_scale", "long_term_scale"))


par(mfrow = c(1, 3))
plot(density(rgamma(10000, shape = 2, rate = 2)), t = "l")
lines(density(draws$length_scale), col = "red")
plot(density(rgamma(10000, shape = 100, rate = 100)), t = "l")
lines(density(draws$periodic_scale), col = "red")
plot(density(rgamma(10000, shape = 20, rate = 0.1)), t = "l")
lines(density(draws$long_term_scale), col = "red")

#   par = c(20, 1, 200),
#   fn = log_likelihood_hyperparameters,
#   f = obs_data$f_infer,
#   n = n,
#   nt = nt,
