#' RBF kernel function
rbf_kernel <- function(x, theta) {
  exp(-x^2 / (2 * theta^2))
}

#' Periodic kernel function
periodic_kernel <- function(x, alpha, period) {
  exp(-2 * sin(pi * x / period)^2 / alpha^2)
}

get_spatial_distance <- function(coordinates){
  dist(coordinates[,c("lon", "lat")], diag = TRUE, upper = TRUE) |>
    as.matrix()
}

get_temporal_distance <- function(times){
  dist(times, diag = TRUE, upper = TRUE) |>
    as.matrix()
}

#' Estimate the spatial kernel
space_kernel <- function(coordinates, length_scale, nugget = 1e-9){
  space_matrix <- get_spatial_distance(coordinates)
  rbf_kernel(space_matrix, theta = length_scale) + diag(x = nugget, nrow = nrow(space_matrix))
}

#' Extimate the temporal kernel
time_kernel <- function(times, periodic_scale, long_term_scale, nugget = 1e-9, period = 52){
  time_matrix <- get_temporal_distance(times)
  period_k <- periodic_kernel(x = time_matrix, alpha = periodic_scale, period = period)
  long_term_k <- rbf_kernel(x = time_matrix, theta = long_term_scale)
  period_k * long_term_k + diag(x = nugget, nrow = nrow(time_matrix))
}
