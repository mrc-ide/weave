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
periodic_kernel <- function(distance, periodic_scale, long_term_scale, period) {
  # Periodic component
  sine_term <- sin(pi * (distance) / period) ^ 2
  periodic_component <- exp(-2 * periodic_scale ^ 2 * sine_term)

  # Long-term decay component
  time_distance <- abs(distance)
  long_term_decay <- exp(-(1 / long_term_scale) * time_distance)

  # Combine both components: periodic and long-term decay
  kernel_value <- periodic_component * long_term_decay

  return(kernel_value)
}

#' RBF kernel function
#'
#' @param x A scalar representing the first input point
#' @param x_prime A scalar or vector representing the second input point(s)
#' @param sigma The length scale parameter, controlling the smoothness of the kernel
#'
#' @export
rbf_kernel <- function(distance, sigma) {
  # Squared distance between x and x_prime
  distance_squared <- (distance)^2

  # Apply the RBF kernel formula
  kernel_value <- exp(-distance_squared / (2 * sigma^2))

  return(kernel_value)
}

#' Exponential kernel function
#'
#' @param x A scalar representing the first input point
#' @param x_prime A scalar or vector representing the second input point(s)
#' @param sigma The length scale parameter, controlling the rate of decay
#'
#' @export
exponential_kernel <- function(distance, sigma) {
  # Absolute distance between x and x_prime
  distance <- abs(distance)

  # Apply the Exponential kernel formula
  kernel_value <- exp(-distance / sigma)

  return(kernel_value)
}
