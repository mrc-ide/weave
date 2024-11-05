create_time_matrix <- function(data, periodic_scale, long_term_scale, period, epsilon = 0.01, nugget = 0){
  times <- unique(data$t)
  time_distance <- outer(times,  times, "-")
  time <- time_distance |>
    periodic_kernel(periodic_scale, long_term_scale, period)
  time <- time * (1 - nugget)
  time <- Matrix::nearPD(time)
  time <- as.matrix(time$mat)
  time_chol <- chol(time)
  time_reg <- time + epsilon * diag(max(times))
  time_inv_reg <- MASS::ginv(time_reg)

  return(
    list(
      time_distance = time_distance,
      time = time,
      time_reg = time_reg,
      time_inv_reg = time_inv_reg,
      time_chol = time_chol
    )
  )
}

create_spatial_matrix <- function(data, space_sigma, epsilon = 0.01, nugget = 0){
  coordinates <- data |>
    dplyr::select(id, lat, lon) |>
    dplyr::distinct()

  sigmasq <- data |>
    dplyr::select(id, observed_sigmasq) |>
    dplyr::distinct()

  spatial_distance <- coordinates |>
    dist() |>
    as.matrix()

  space <- spatial_distance |>
    rbf_kernel(space_sigma)

  spcae <- space * (1 - nugget)

  # Create D_space
  d_space <- diag(sqrt(sigmasq$observed_sigmasq))
  # Compute Sigma_space
  space <- d_space %*% space %*% d_space
  space <- Matrix::nearPD(space)
  space <- as.matrix(space$mat)
  space_reg <- space + epsilon * diag(nrow(coordinates))
  space_inv_reg <- MASS::ginv(space_reg)
  space_chol <- chol(space)

  return(
    list(
      spatial_distance = spatial_distance,
      space = space,
      space_reg = space_reg,
      space_inv_reg = space_inv_reg,
      space_chol = space_chol
    )
  )
}
