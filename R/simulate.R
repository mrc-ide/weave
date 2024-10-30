simulate_data <- function(n_sites, max_weeks, mu, sd, space_sigma, periodic_scale, long_term_scale, period){

  df_site_pos <- data.frame(
    site_index = 1:n_sites,
    lon = 1:n_sites,
    lat = 1:n_sites
  )

  output_df <- tidyr::expand_grid(
    site_index = 1:n_sites,
    week = 1:max_weeks
  ) |>
    dplyr::left_join(
      df_site_pos, by = "site_index"
    )

  output_df$lambda_mean <- mu[output_df$site_index]
  output_df$lambda_sd <- sd[output_df$site_index]
  sigma_sq <- log((sd / mu)^2 + 1)
  output_df$sigmasq <- sigma_sq[output_df$site_index]
  output_df$mu <- log(output_df$lambda_mean) - output_df$sigmasq / 2

  space_matrix <- create_spatial_matrix(
    coords = df_site_pos[, c("lat", "lon")],
    sigma_sq = sigma_sq,
    space_sigma = space_sigma
    )

  time_matrix <- create_time_matrix(
    times = 1:max_weeks,
    periodic_scale = periodic_scale,
    long_term_scale = long_term_scale,
    period = period
    )

  # Todo: put chol in directly
  output_df$z <- quick_mvnorm(space_matrix$space, time_matrix$time) + output_df$mu

  output_df$lambda <- exp(output_df$z)
  output_df$n <- rpois(nrow(output_df), output_df$lambda)

  return(output_df)
}
