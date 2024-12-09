simulate_data <- function(
    n_sites, max_t, site_means,
    site_sds, theta = 0.1, periodic_scale = 0.65,
    long_term_scale = 500, period = 52,
    dist = "poisson", size = NULL){

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

  sigmasq <- log((site_sds / site_means)^2 + 1)
  mu <- log(site_means) - sigmasq / 2

  space <- create_spatial_matrix(df_site_pos[,c("lat", "lon")], sigmasq, theta = theta)
  time <- create_time_matrix(1:max_t, periodic_scale = periodic_scale, long_term_scale = long_term_scale, period = period)

  output_df$mu <- mu[output_df$id]
  output_df$sigmasq <- sigmasq[output_df$id]
  output_df$z <- quick_mvnorm(space, time) + output_df$mu
  output_df$lambda <- exp(output_df$z)
  if (dist == "poisson") {
    output_df$n <- rpois(nrow(output_df), output_df$lambda)
  } else if (dist == "nbinom") {
    if (is.null(size)) {
      stop("Size parameter must be provided for negative binomial distribution.")
    }
    output_df$n <- rnbinom(nrow(output_df), mu = output_df$lambda, size = size)
  } else {
    stop("Unsupported distribution")
  }

  output_df <- output_df |>
    dplyr::rename("site_index" = id)

  return(output_df)
}
