devtools::load_all()
library(progress)
library(ggplot2)
library(dplyr)
library(tidyr)

# True parameters --------------------------------------------------------------
set.seed(321)
# Number of sites
n = 16
# Number of timesteps
nt = 52 * 10
# Site mean case count
site_means = round(runif(n, 10, 100))
# length_scale: determines how quickly correlation decays with distance
#   - higher length_scale: correlation persists over longer distances (smoother spatial variation)
#   - lower length_scale: correlation decays rapidly, indicating localised variation
length_scale <- 3
# periodic_scale: strength of repeating (seasonal or cyclical) patterns
#   - higher periodic_scale: stronger seasonal patterns
#   - lower periodic_scale: weaker seasonal patterns
periodic_scale = 1
# long_term_scale: scale controlling decay rate of long-term temporal correlation
#   - higher long_term_scale: smoother long-term trends (correlation persists longer)
#   - lower long_term_scale: rapid loss of correlation, short-term variation dominates
long_term_scale = 200
# period: duration of the repeating cycle (e.g., 52 weeks for annual seasonality)
period = 52
# p_one: probability that a new cluster begins with a missing value.
#   - higher p_one: more missing data overall.
#   - lower p_one: less missing data overall.
p_one = 0.6
# p_switch: probability of switching between missing and observed states.
#   - lower p_switch: longer sequences (clusters) of missingness or non-missingness.
#   - higher p_switch: shorter clusters, more frequent switching between states.
p_switch = 0.3
# ------------------------------------------------------------------------------

# Simulated data ---------------------------------------------------------------
coordinates <- data.frame(
  id = factor(1:n),
  lat = 1:n,
  lon = 1:n,
  mu = log(site_means)
)

space_k <- space_kernel(
  coordinates = coordinates,
  length_scale = length_scale
)

time_k <- time_kernel(
  times = 1:nt,
  periodic_scale = periodic_scale,
  long_term_scale = long_term_scale,
  period = period
)

true_data <- simulate_data(
  n = n,
  nt = nt,
  coordinates = coordinates,
  space_k = space_k,
  time_k = time_k
)

obs_data <- observed_data2(
  data = true_data,
  p_one = p_one,
  p_switch = p_switch
)

space_pd <- data.frame(
  distance = 1:100
) |>
  mutate(
    k = rbf_kernel(distance, theta = length_scale),
    group = "True"
  )

space_plot <- ggplot(space_pd, aes(x = distance, y = k, colour = group)) +
  geom_line() +
  theme_bw()

time_pd <- data.frame(
  week = 1:(nt)
) |>
  mutate(
    k = periodic_kernel(x = week, alpha = periodic_scale, period = period) *
      rbf_kernel(x = week, theta = long_term_scale),
    group = "True"
  )

time_plot <- ggplot(time_pd, aes(x = week, y = k, colour = group)) +
  geom_line() +
  theme_bw()

hf_labeller <- function(value) {
  paste("HF:", value)
}

sim_data <- ggplot() +
  geom_point(data = true_data, aes(x = t, y = y), size = 0.3, colour = "red") +
  geom_point(data = obs_data, aes(x = t, y = y_obs), size = 0.3, colour = "black") +
  geom_line(data = true_data, aes(x = t, y = lambda, colour = id)) +
  facet_wrap(~ id, scales = "free_y", labeller = labeller(id = hf_labeller)) +
  ylab("Cases") +
  xlab("Week") +
  theme_bw() +
  theme(
    legend.position = "none",
    strip.background = element_rect(fill = "white", colour = "grey50"),
    strip.text = element_text(size = 8, face = "bold"),
    panel.spacing = unit(0.5, "lines")
  ) +
  ggtitle("Simulated data")
# ------------------------------------------------------------------------------

# Infer kernel hyper-parameters ------------------------------------------------
infer_space <- infer_space_kernel_params(obs_data, TRUE)
infer_time <- infer_time_kernel_params(obs_data, 52, TRUE)

hyperparmameters <- c(infer_space$length_scale, infer_time$periodic_scale, infer_time$long_term_scale)
#hyperparmameters <- c(3, 1, 200)
## Options
# Fix these conservatively - simulated sensitvity analyses?
# Fit Bayesian with priors
# ------------------------------------------------------------------------------

# Fit --------------------------------------------------------------------------
sub_n <- 4

time_matrix <- time_kernel(
  1:nt,
  periodic_scale = hyperparmameters[2],
  long_term_scale = hyperparmameters[3],
  period = 52
)

time_inv_reg <- time_matrix |>
  regularise() |>
  solve()

space_matrix <- get_spatial_distance(coordinates)

pars <- rep(NA, nrow(obs_data))
tausq <- rep(NA, nrow(obs_data))
pb <- progress_bar$new(total = n)

for(i in 1:n){
  pb$tick()
  dists <- space_matrix[i,]
  indices <- order(dists)[1:(sub_n + 1)]
  sub_d <- dplyr::filter(obs_data, id %in% indices)

  sub_coordinates <- sub_d |>
    dplyr::select(id, lat, lon) |>
    dplyr::distinct() |>
    dplyr::select(lat, lon)

  sub_space_matrix <- space_kernel(
    coordinates = sub_coordinates,
    length_scale = hyperparmameters[1]
  )

  sub_space_inv_reg <- sub_space_matrix |>
    regularise() |>
    solve()

  # Full Spatiotemporal Kernel (Kronecker Product)
  K_spacetime <- kronecker(sub_space_matrix, time_matrix)

  #miss_index <- which(is.na(sub_d$y_obs))
  miss_index <- 1:nrow(sub_d)
  obs_index <- which(!is.na(sub_d$y_obs))

  y_obs <- log(sub_d$y_obs + 0.0001)[obs_index]
  y_obs <- y_obs - sub_d$mu_infer[obs_index]

  # 2) Subset into blocks for observed and missing
  K_obs_obs   <- K_spacetime[obs_index, obs_index]
  K_obs_miss  <- K_spacetime[obs_index, miss_index]
  K_miss_obs  <- K_spacetime[miss_index, obs_index]
  K_miss_miss <- K_spacetime[miss_index, miss_index]

  # 3) GP posterior (Gaussian version):
  #    (Add noise variance on the diagonal of K_obs_obs if needed)
  noise_var <- 1e-3
  K_obs_obs_noisy <- K_obs_obs + diag(noise_var, nrow(K_obs_obs))

  # Solve for alpha
  alpha <- solve(K_obs_obs_noisy, y_obs)

  # Posterior mean for missing points
  post_mean <- K_miss_obs %*% alpha
  post_mean <- post_mean + sub_d$mu_infer
  # opt <- optim(
  #   par = sub_d$f_infer,
  #   fn = log_likelihood,
  #   gr = log_likelihood_gradient,
  #   mu = sub_d$mu_infer ,
  #   dist_k_inv = sub_space_inv_reg,
  #   time_k_inv = time_inv_reg,
  #   n = sub_n + 1,
  #   nt = nt,
  #   y = sub_d$y_obs,
  #   method = "BFGS",  # BFGS is a gradient-based optimization method
  #   control = list(
  #     fnscale = -1,     # To maximize the function
  #     maxit = 100000,     # Increase the maximum number of iterations
  #     reltol = 0.000001
  #   )
  # )

  pars[obs_data$id == i] <- post_mean[sub_d$id == i]

  # calculate Hessian at the ML value
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

  tausq[obs_data$id == i] <- sub_tausq[sub_d$id == i]
}

mean_lambda <- mean(exp(pars))
var_y <- var(obs_data$y_obs, na.rm = TRUE)
overdispersion <- (mean_lambda^2) / (var_y - mean_lambda)

fit_data <- obs_data |>
  dplyr::mutate(

    # Estimate of z, conditioned on observations
    z_est = pars,
    # Posterior variance of z
    tausq = tausq,

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
  )

sim_data +
  geom_ribbon(
    data = fit_data,
    aes(x = t, ymin = pred_Q2.5, ymax = pred_Q97.5, fill = id), alpha = 0.25
  ) +
  geom_ribbon(
    data = fit_data,
    aes(x = t, ymin = pred_Q25, ymax = pred_Q75, fill = id, alpha = 0.5)
  ) +
  geom_line(
    data = fit_data,
    aes(x = t, y = data_Q50, col = id), linewidth = 1
  )

