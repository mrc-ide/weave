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
    n, nt,
    coordinates,
    space_k, time_k){

  output_df <- tidyr::expand_grid(
    id = factor(1:n),
    t = 1:nt
  ) |>
    dplyr::left_join(
      coordinates, by = "id"
    ) |>
    dplyr::mutate(
      f = quick_mvnorm(space_k, time_k),
      z = mu + f,
      lambda = exp(z),
      y = rpois(n * nt, lambda)
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
