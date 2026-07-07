#' Radial basis function kernel
#'
#' Computes the radial basis function (RBF) kernel for a distance vector
#' or matrix:
#' \deqn{k(d) = \exp\left(-\frac{d^2}{2\theta^2}\right).}
#' This is the correlation form of the kernel (\eqn{k(0) = 1}); any global
#' variance is applied separately.
#'
#' @param x A numeric vector or matrix of distances.
#' @param theta A positive numeric scalar giving the length-scale parameter
#'   (the \eqn{\ell} of textbook presentations): correlation decays with
#'   distance over this scale, so larger values give smoother functions.
#'
#' @return A numeric vector or matrix with RBF kernel values.
#' @export
rbf_kernel <- function(x, theta) {
  exp(-x^2 / (2 * theta^2))
}

#' Periodic kernel
#'
#' Computes a periodic kernel for a distance vector or matrix:
#' \deqn{k(d) = \exp\left(-\frac{2\sin^2(\pi d / p)}{\alpha^2}\right).}
#'
#' @param x A numeric vector or matrix of distances.
#' @param alpha A positive numeric scalar controlling how sharply correlation
#'   falls within each cycle: smaller values allow sharp seasonal peaks;
#'   larger values give a gentler, smoother cycle.
#' @param period A positive numeric scalar giving the period \eqn{p}.
#'
#' @return A numeric vector or matrix of periodic kernel values.
#' @export
periodic_kernel <- function(x, alpha, period) {
  exp(-2 * sin(pi * x / period)^2 / alpha^2)
}

#' Pairwise spatial distances
#'
#' Computes pairwise Euclidean distances between locations, in the units of
#' the coordinates. No great-circle correction is applied: with raw
#' longitude/latitude degrees, one degree of longitude shrinks with latitude,
#' so for large or high-latitude extents project the coordinates first (e.g.
#' to km) and interpret `length_scale` in those units.
#'
#' @param coordinates A data frame with columns `lon` and `lat`.
#'
#' @return A symmetric matrix of pairwise spatial distances.
#' @export
get_spatial_distance <- function(coordinates) {
  stats::dist(coordinates[, c("lon", "lat")], diag = TRUE, upper = TRUE) |>
    as.matrix()
}

#' Pairwise temporal distances
#'
#' Computes pairwise distances between time points.
#'
#' @param times A numeric vector of time indices.
#'
#' @return A symmetric matrix of pairwise temporal distances.
#' @export
get_temporal_distance <- function(times) {
  stats::dist(times, diag = TRUE, upper = TRUE) |>
    as.matrix()
}

#' Build the spatial correlation matrix
#'
#' Builds a spatial correlation matrix using an RBF kernel with a nugget
#' term for numerical stability. Distances are Euclidean in the coordinate
#' units (see [get_spatial_distance()]), so `length_scale` is in those same
#' units.
#'
#' @param coordinates A data frame with columns `lon` and `lat`.
#' @param length_scale A positive numeric scalar for the spatial length-scale.
#' @param nugget A non-negative numeric scalar added to the diagonal for
#'   numerical stability.
#'
#' @return A positive-definite matrix representing spatial correlation.
#' @export
space_kernel <- function(coordinates, length_scale, nugget = 1e-9) {
  space_matrix <- get_spatial_distance(coordinates)
  rbf_kernel(space_matrix, theta = length_scale) +
    diag(x = nugget, nrow = nrow(space_matrix))
}

#' Build the temporal correlation matrix
#'
#' Builds a temporal correlation matrix by combining periodic and
#' long-term RBF components with a nugget term for numerical stability.
#'
#' @param times A numeric vector of time indices.
#' @param periodic_scale A positive numeric scalar controlling how sharply
#'   correlation falls within each seasonal cycle: smaller values allow sharp
#'   seasonal peaks; larger values give a gentler, smoother cycle.
#' @param long_term_scale A positive numeric scalar for the long-term
#'   length-scale.
#' @param nugget A non-negative numeric scalar added to the diagonal for
#'   numerical stability.
#' @param period A positive numeric scalar giving the period of the seasonal
#'   component, in the same units as `times`.
#'
#' @return A positive-definite matrix representing temporal correlation.
#' @export
time_kernel <- function(times, periodic_scale, long_term_scale,
                        nugget = 1e-9, period = 52) {
  time_matrix <- get_temporal_distance(times)
  period_k <- periodic_kernel(x = time_matrix, alpha = periodic_scale, period = period)
  long_term_k <- rbf_kernel(x = time_matrix, theta = long_term_scale)
  period_k * long_term_k + diag(x = nugget, nrow = nrow(time_matrix))
}
