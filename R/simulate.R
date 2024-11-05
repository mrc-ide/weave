simulate_data <- function(n_sites, max_t, mu, sd, space_sigma, periodic_scale, long_term_scale, period, nugget = 0){

  df_site_pos <- data.frame(
    id = 1:n_sites,
    lat = 1:n_sites,
    lon = 1:n_sites
  )

  output_df <- tidyr::expand_grid(
    id = 1:n_sites,
    t = 1:max_t
  ) |>
    dplyr::left_join(
      df_site_pos, by = "id"
    )

  output_df$lambda_mean <- mu[output_df$id]
  output_df$lambda_sd <- sd[output_df$id]
  sigma_sq <- log((sd / mu)^2 + 1)
  output_df$observed_sigmasq <- sigma_sq[output_df$id]
  output_df$mu <- log(output_df$lambda_mean) - output_df$observed_sigmasq / 2

  #browser()
  space_matrix <- create_spatial_matrix(
    output_df,
    space_sigma = space_sigma,
    nugget = nugget
    )

  time_matrix <- create_time_matrix(
    output_df,
    periodic_scale = periodic_scale,
    long_term_scale = long_term_scale,
    period = period
    )

  # Todo: put chol in directly
  output_df$z <- quick_mvnorm(space_matrix$space, time_matrix$time) + output_df$mu
  output_df$lambda <- exp(output_df$z)
  output_df$n <- rpois(nrow(output_df), output_df$lambda)

  output_df <- output_df |>
    dplyr::rename("site_index" = id)

  return(output_df)
}
