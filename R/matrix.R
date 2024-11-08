create_time_matrix <- function(times, periodic_scale, long_term_scale, period){
  time_distance <- abs(outer(times,  times, "-"))
  time <- time_distance |>
    periodic_kernel(periodic_scale, long_term_scale, period)
  time <- Matrix::nearPD(time)
  time <- as.matrix(time$mat)


  return(time)
}

create_spatial_matrix <- function(coordinates, sigma2, theta){

  spatial_kernel <- coordinates |>
    dist() |>
    as.matrix() |>
    rbf_kernel(
      theta = theta
    )

  sigma <- outer(1:nrow(coordinates), 1:nrow(coordinates), function(x, y) sqrt(sigma2[x] * sigma2[y]))

  space <- spatial_kernel * sigma
  space <- Matrix::nearPD(space)
  space <- as.matrix(space$mat)

  return(space)
}

regularise <- function(x, lambda = 1e-5) {
  x + lambda * diag(nrow(x))
}
