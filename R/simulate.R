#' Simulate data from the model
#' Case counts:
#'  y_st ~ Pois(λ_st)
#' Cases mean:
#'  λ_st = exp(z_st)
#' Cases mean (log scale):
#'  z_st = µ_s + f_st
#' Where we have a site specific intercept:
#'  µ_s
#' and a latent GP:
#'  f_st ~ MVN(0,Σ)
#' Spatiotemporal covariance:
#'  Σ = dist_k ⊗ time_k
simulate_data <- function(
    n, nt, site_means,
    theta = 0.1, alpha = 0.65,
    beta = 500, period = 52, p_one = 0, p_switch = 0){

  coordinates <- data.frame(
    id = factor(1:n),
    lat = 1:n,
    lon = 1:n,
    mu = log(site_means)
  )

  dist_k <- space_kernel(coordinates, theta = theta)
  time_k <- time_kernel(times = 1:nt, alpha = alpha, beta = beta)

  output_df <- tidyr::expand_grid(
    id = factor(1:n),
    t = 1:nt
  ) |>
    dplyr::left_join(
      coordinates, by = "id"
    ) |>
    dplyr::mutate(
      f = quick_mvnorm(dist_k, time_k),
      z = mu + f,
      lambda = exp(z),
      y = rpois(n * nt, lambda)
    ) |>
    dplyr::mutate(
      true_y = y,
      missing = generate_clustered_binary(dplyr::n(), p_one, p_switch),
      y = ifelse(missing == 1, NA, y)
    ) |>
    dplyr::select(- missing)

  return(output_df)
}
