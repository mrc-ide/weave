simulate_data <- function(n_sites, max_t, mu, sigma_time, sigma_space, theta, periodic_scale, long_term_scale, period, nugget = 0){

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

  lambda_mean <- mu[output_df$id]
  lambda_sigma_space <- sigma_space[output_df$id]
  lambda_sd_space <- sqrt(lambda_sigma_space^2)
  output_df$sigma_sq_space <- log((lambda_sigma_space / lambda_mean)^2 + 1)
  lambda_sigma_time <- sigma_time[output_df$id]
  lambda_sd_time<- sqrt(lambda_sigma_time^2)
  output_df$sigma_sq_time <- log((lambda_sigma_time / lambda_mean)^2 + 1)
  lamda_sd_combined <- sqrt(lambda_sigma_space^2 + lambda_sigma_time^2)
  sigma_sq_combined <- log((lamda_sd_combined / lambda_mean)^2 + 1)
  output_df$mu <- log(lambda_mean) - sigma_sq_combined / 2


  space_matrix <- create_spatial_matrix(
    output_df,
    theta = theta,
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
