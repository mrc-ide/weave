create_time_matrix <- function(times, periodic_scale, long_term_scale, period, epsilon = 0.01){
  time_distance <- outer(times,  times, "-")
  time <- time_distance |>
    periodic_kernel(periodic_scale, long_term_scale, period)
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

create_spatial_matrix <- function(coords, sigma_sq, space_sigma, epsilon = 0.01){
  spatial_distance <- coords |>
    dist() |>
    as.matrix()
  space <- spatial_distance |>
    rbf_kernel(space_sigma)

  # Create D_space
  d_space <- diag(sqrt(sigma_sq))
  # Compute Sigma_space
  space <- d_space %*% space %*% d_space
  space <- Matrix::nearPD(space)
  space <- as.matrix(space$mat)
  space_reg <- space + epsilon * diag(nrow(coords))
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
