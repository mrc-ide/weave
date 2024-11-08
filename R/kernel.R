#' Modified periodic kernel with long-term decay
#'
#' Period kernel refs https://tinyurl.com/52mzsrhh, https://tinyurl.com/yc4dtm68
#'
#' @param t A scalar representing the time point(s) for the first input
#' @param t_prime A scalar or vector representing the time point(s) for the second input
#' @param periodic_scale The length scale parameter for the periodic kernel (in time units)
#' @param long_term_scale The length scale parameter for the long term kernel (in time units)
#' @param period The period of the data, representing the time interval over which periodic behavior repeats (in time units)
#'
#' @export
periodic_kernel <- function(time_distance, periodic_scale = 0.65, long_term_scale = 500, period = 52) {
  periodic_component <- exp(-2 * sin(pi * time_distance / period)^2 / periodic_scale^2)
  long_term_component <- exp(-(1 / long_term_scale) * time_distance)
  kernel_value <- periodic_component * long_term_component
  return(kernel_value)
}

#' RBF kernel function
#'
#' @param x A scalar representing the first input point
#' @param x_prime A scalar or vector representing the second input point(s)
#' @param sigma The length scale parameter, controlling the smoothness of the kernel
#'
#' @export
rbf_kernel <- function(distance, theta) {
  kernel_value <- exp(-distance^2 / (2 * theta^2))
  return(kernel_value)
}

get_spatial_distance <- function(data) {
  data |>
    dplyr::select(id, lat, lon) |>
    dplyr::distinct() |>
    dplyr::arrange(id) |>
    dplyr::select(-id) |>
    dist() |>
    as.matrix()
}


