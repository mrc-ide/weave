#' Radial basis function kernel
#'
#' Computes the radial basis function (RBF) kernel for a distance vector
#' or matrix.
#'
#' @param x A numeric vector or matrix of distances.
#' @param theta A positive numeric scalar giving the length-scale parameter.
#'
#' @return A numeric vector or matrix with RBF kernel values.
rbf_kernel <- function(x, theta) {
  exp(-x^2 / (2 * theta^2))
}

#' Periodic kernel
#'
#' Computes a periodic kernel for a distance vector or matrix.
#'
#' @param x A numeric vector or matrix of distances.
#' @param alpha A positive numeric scalar controlling the amplitude.
#' @param period A positive numeric scalar giving the period.
#'
#' @return A numeric vector or matrix of periodic kernel values.
periodic_kernel <- function(x, alpha, period) {
  exp(-2 * sin(pi * x / period)^2 / alpha^2)
}

#' Pairwise spatial distances
#'
#' Computes pairwise Euclidean distances between locations.
#'
#' @param coordinates A data frame with columns `lon` and `lat` in degrees.
#'
#' @return A symmetric matrix of pairwise spatial distances.
get_spatial_distance <- function(coordinates) {
  dist(coordinates[, c("lon", "lat")], diag = TRUE, upper = TRUE) |>
    as.matrix()
}

#' Pairwise temporal distances
#'
#' Computes pairwise distances between time points.
#'
#' @param times A numeric vector of time indices.
#'
#' @return A symmetric matrix of pairwise temporal distances.
get_temporal_distance <- function(times) {
  dist(times, diag = TRUE, upper = TRUE) |>
    as.matrix()
}

#' Estimate the spatial kernel
#'
#' Builds a spatial covariance matrix using an RBF kernel with a nugget
#' term for numerical stability.
#'
#' @param coordinates A data frame with columns `lon` and `lat` in degrees.
#' @param length_scale A positive numeric scalar for the spatial length scale.
#' @param nugget A non-negative numeric scalar added to the diagonal for
#'   numerical stability.
#'
#' @return A positive-definite matrix representing spatial covariance.
space_kernel <- function(coordinates, length_scale, nugget = 1e-9) {
  space_matrix <- get_spatial_distance(coordinates)
  rbf_kernel(space_matrix, theta = length_scale) +
    diag(x = nugget, nrow = nrow(space_matrix))
}

#' Estimate the temporal kernel
#'
#' Builds a temporal covariance matrix by combining periodic and
#' long-term RBF components with a nugget term for numerical stability.
#'
#' @param times A numeric vector of time indices.
#' @param periodic_scale A positive numeric scalar controlling the periodic
#'   variation.
#' @param long_term_scale A positive numeric scalar for the long-term length
#'   scale.
#' @param nugget A non-negative numeric scalar added to the diagonal for
#'   numerical stability.
#' @param period A positive numeric scalar giving the period of the seasonal
#'   component.
#'
#' @return A positive-definite matrix representing temporal covariance.
time_kernel <- function(times, periodic_scale, long_term_scale,
                        nugget = 1e-9, period = 52) {
  time_matrix <- get_temporal_distance(times)
  period_k <- periodic_kernel(x = time_matrix, alpha = periodic_scale, period = period)
  long_term_k <- rbf_kernel(x = time_matrix, theta = long_term_scale)
  period_k * long_term_k + diag(x = nugget, nrow = nrow(time_matrix))
}
