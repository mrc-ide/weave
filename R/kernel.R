#' Radial basis function kernel
#'
#' Computes the radial basis function (RBF) kernel for a distance vector
#' or matrix.
#'
#' @param x A numeric vector or matrix of distances.
#' @param theta A positive numeric scalar giving the length-scale parameter.
#'
#' @return A numeric vector or matrix with RBF kernel values.
#' @export
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
#' @export
periodic_kernel <- function(x, alpha, period) {
  exp(-2 * sin(pi * x / period)^2 / alpha^2)
}

#' Pairwise spatial distances (Euclidean)
#'
#' Computes pairwise Euclidean distances between locations. Note that this
#' treats `lon`/`lat` as Cartesian coordinates -- correct for projected
#' coordinate systems and for synthetic data, but **wrong for real
#' geographic coordinates in degrees** (one degree of longitude is not one
#' degree of latitude). For real lat/lon data pass
#' `distance_fn = haversine_distance` to [space_kernel()].
#'
#' @param coordinates A data frame with columns `lon` and `lat`.
#'
#' @return A symmetric matrix of pairwise spatial distances.
#' @export
get_spatial_distance <- function(coordinates) {
  stats::dist(coordinates[, c("lon", "lat")], diag = TRUE, upper = TRUE) |>
    as.matrix()
}

#' Pairwise great-circle (haversine) distances
#'
#' Treats `lon` and `lat` as degrees on the surface of the Earth and returns
#' pairwise distances in kilometres. Use this whenever sites are spread over
#' a region large enough that the Euclidean approximation breaks down (more
#' than a few hundred kilometres) or when sites span a wide range of
#' latitudes.
#'
#' @param coordinates A data frame with columns `lon` and `lat` in degrees.
#' @param earth_radius_km Earth radius to use; defaults to 6371 km (mean).
#'
#' @return A symmetric matrix of pairwise great-circle distances in km.
#' @export
haversine_distance <- function(coordinates, earth_radius_km = 6371) {
  lat <- coordinates$lat * pi / 180
  lon <- coordinates$lon * pi / 180
  n   <- length(lat)

  # Pairwise differences via outer().
  dlat <- outer(lat, lat, `-`)
  dlon <- outer(lon, lon, `-`)
  a    <- sin(dlat / 2)^2 +
          cos(outer(lat, lat, function(a, b) (a + b) / 2 - (a - b) / 2)) *
          cos(outer(lat, lat, function(a, b) (a + b) / 2 + (a - b) / 2)) *
          sin(dlon / 2)^2
  # The above expands to cos(lat_i) cos(lat_j); rewrite explicitly for clarity.
  a    <- sin(dlat / 2)^2 +
          outer(cos(lat), cos(lat)) * sin(dlon / 2)^2
  c_   <- 2 * atan2(sqrt(a), sqrt(pmax(0, 1 - a)))

  earth_radius_km * c_
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

#' Estimate the spatial kernel
#'
#' Builds a spatial covariance matrix using an RBF kernel with a nugget
#' term for numerical stability.
#'
#' By default uses Euclidean distance on `(lon, lat)`. Pass
#' `distance_fn = haversine_distance` for real geographic coordinates in
#' degrees.
#'
#' @param coordinates A data frame with columns `lon` and `lat`.
#' @param length_scale A positive numeric scalar for the spatial length scale.
#'   Interpreted in the units returned by `distance_fn` (degrees by default,
#'   kilometres with haversine).
#' @param nugget A non-negative numeric scalar added to the diagonal for
#'   numerical stability.
#' @param distance_fn Function taking `coordinates` and returning a symmetric
#'   distance matrix. Defaults to [get_spatial_distance()] (Euclidean).
#'
#' @return A positive-definite matrix representing spatial covariance.
#' @export
space_kernel <- function(coordinates, length_scale, nugget = 1e-9,
                         distance_fn = get_spatial_distance) {
  space_matrix <- distance_fn(coordinates)
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
#' @export
time_kernel <- function(times, periodic_scale, long_term_scale,
                        nugget = 1e-9, period = 52) {
  time_matrix <- get_temporal_distance(times)
  period_k <- periodic_kernel(x = time_matrix, alpha = periodic_scale, period = period)
  long_term_k <- rbf_kernel(x = time_matrix, theta = long_term_scale)
  period_k * long_term_k + diag(x = nugget, nrow = nrow(time_matrix))
}
