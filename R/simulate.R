simulate_data <- function(n_sites, max_weeks, mu, sd, space_sigma, periodic_scale, long_term_scale, period){

  output_df <- expand.grid(
    site_index = 1:n_sites,
    week = 1:max_weeks
  )

  output_df$lambda_mean <- mu[output_df$site_index]
  output_df$lambda_sd <- sd[output_df$site_index]
  output_df$sigmasq <- log((output_df$lambda_sd / output_df$lambda_mean)^2 + 1)
  output_df$mu <- log(output_df$lambda_mean) - output_df$sigmasq / 2

  df_site_pos <- data.frame(
    lon = 1:n_sites,
    lat = 1:n_sites
  )

  spatial_distance <- df_site_pos |>
    dist() |>
    as.matrix()
  space <- spatial_distance |>
    rbf_kernel(space_sigma)
  space <- Matrix::nearPD(space)
  space <- as.matrix(space$mat)
  space_chol <- chol(space)

  time_distance <- outer(1:max_weeks,  1:max_weeks, "-")
  time <- time_distance |>
    periodic_kernel(periodic_scale, long_term_scale, period)
  time <- Matrix::nearPD(time)
  time <- as.matrix(time$mat)
  time_chol <- chol(time)

  output_df$z <- quick_mvnorm(space, time) + output_df$mu

  output_df$lambda <- exp(output_df$z)
  output_df$n <-rpois(nrow(output_df), output_df$lambda)

  return(output_df)
}
