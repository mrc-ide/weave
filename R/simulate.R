#' Simulate spatiotemporal case counts
#'
#' Generates Poisson counts for `n` sites across `nt` time points from a latent
#' Gaussian process with separable spatial and temporal kernels.
#'
#' @param n Number of spatial locations.
#' @param nt Number of time points.
#' @param coordinates Data frame with columns `id`, `lat` and `lon` describing
#'   site locations.
#' @param space_k Spatial covariance matrix.
#' @param time_k Temporal covariance matrix.
#'
#' @return Tibble containing the latent effect `f`, mean `z`, rate `lambda` and
#'   simulated counts `y` for each site-time combination.
#'
#' @details Randomness comes from multivariate normal and Poisson draws. Use
#'   [set.seed()] before calling for reproducible simulations.
#'
#' @keywords internal
simulate_data <- function(
    n, nt,
    coordinates,
    space_k, time_k) {

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

#' Generate observed data with missingness
#'
#' Introduces missing values into simulated counts and computes simple
#' estimates of the latent components for each site.
#'
#' @param data Output of [simulate_data()].
#' @param p_one Probability that a new missingness cluster starts with a missing
#'   value.
#' @param p_switch Probability of switching between missing and observed
#'   clusters.
#'
#' @return Tibble with observed counts `y_obs` and inferred `mu_infer`,
#'   `z_infer` and `f_infer`.
#'
#' @details Missingness is generated stochastically via
#'   [generate_clustered_binary()]. Set a seed with [set.seed()] to reproduce
#'   results.
#'
#' @keywords internal
observed_data <- function(data, p_one, p_switch) {
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

#' Alternative observed data generation
#'
#' Variant of [observed_data()] that derives `mu_infer` from the mean of the
#' log-transformed observations.
#'
#' @param data Output of [simulate_data()].
#' @param p_one Probability that a new missingness cluster starts with a missing
#'   value.
#' @param p_switch Probability of switching between missing and observed
#'   clusters.
#'
#' @return Tibble with observed counts `y_obs` and inferred `mu_infer`,
#'   `z_infer` and `f_infer`.
#'
#' @details Missingness is generated stochastically via
#'   [generate_clustered_binary()]. Set a seed with [set.seed()] to reproduce
#'   results.
#'
#' @keywords internal
observed_data2 <- function(data, p_one, p_switch) {
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
      .by = "id"
    ) |>
    dplyr::select(id, t, lat, lon, y_obs, mu_infer, z_infer, f_infer)
}
