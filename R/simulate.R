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
#'
#' @param n Number of sites
#' @param nt Number of time points
#' @param coordinates Data frame of site coordinates
#' @param space_k Spatial kernel matrix
#' @param time_k Temporal kernel matrix
#' @param mu Site-specific intercepts; either scalar or vector of length `n`
#'
#' @return A data frame with simulated counts and latent variables
#'
#' @examples
#' n <- 2; nt <- 3
#' coords <- data.frame(id = factor(1:n), lat = 1:n, lon = 1:n)
#' simulate_data(
#'   n, nt, coords,
#'   space_k = diag(n), time_k = diag(nt),
#'   mu = c(0.1, 0.2)
#' )
simulate_data <- function(
    n, nt,
    coordinates,
    space_k, time_k,
    mu = 0){

  if (length(mu) == 1) {
    mu <- rep(mu, n)
  } else if (length(mu) != n) {
    stop("`mu` must have length 1 or `n`")
  }

  mu_vec <- rep(mu, each = nt)

  output_df <- tidyr::expand_grid(
    id = factor(1:n),
    t = 1:nt
  ) |>
    dplyr::left_join(
      coordinates, by = "id"
    ) |>
    dplyr::mutate(
      f = quick_mvnorm(space_k, time_k),
      mu = mu_vec,
      z = mu + f,
      lambda = exp(z),
      y = stats::rpois(n * nt, lambda)
    )

  return(output_df)
}

observed_data <- function(data, p_one, p_switch){
  data |>
    dplyr::mutate(
      y_obs = y,
      missing = generate_clustered_binary(dplyr::n(), p_one, p_switch),
      y_obs = ifelse(missing == 1, NA, y_obs)
    ) |>
    dplyr::mutate(
      mu_infer = log(mean(y_obs, na.rm = TRUE)),
      .by = (id)
    ) |>
    dplyr::mutate(
      z_infer = ifelse(is.na(y_obs), mu_infer, log(y_obs + 1)),
      f_infer = z_infer - mu_infer
    ) |>
    dplyr::select(id, t, lat, lon, y_obs, mu_infer, z_infer, f_infer)
}

observed_data2 <- function(data, p_one, p_switch){
  data |>
    dplyr::mutate(
      y_obs = y,
      missing = generate_clustered_binary(dplyr::n(), p_one, p_switch),
      y_obs = ifelse(missing == 1, NA, y_obs)
    ) |>
    dplyr::mutate(
      z_infer = log(y_obs + 1),
      mu_infer = mean(z_infer, na.rm = TRUE),
      f_infer = z_infer - mu_infer,
      #f_infer = f_infer - mean(f_infer),
      .by = "id"
    ) |>
    dplyr::select(id, t, lat, lon, y_obs, mu_infer, z_infer, f_infer)
}
