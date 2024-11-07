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
periodic_kernel <- function(distance, periodic_scale, long_term_scale, period, sigma, nugget) {
  # Periodic component
  sine_term <- sin(pi * (distance) / period) ^ 2
  periodic_component <- exp(-2 * periodic_scale ^ 2 * sine_term)

  # Long-term decay component
  time_distance <- abs(distance)
  long_term_decay <- exp(-(1 / long_term_scale) * time_distance)

  # Combine both components: periodic and long-term decay
  kernel_value <- sigma^2 * periodic_component * long_term_decay
  kernel_value[distance == 0] <- kernel_value[distance == 0] + nugget
  kernel_value <- pmin(kernel_value, 1)

  return(kernel_value)
}

#' RBF kernel function
#'
#' @param x A scalar representing the first input point
#' @param x_prime A scalar or vector representing the second input point(s)
#' @param sigma The length scale parameter, controlling the smoothness of the kernel
#'
#' @export
rbf_kernel <- function(distance, sigma, theta, nugget = 0) {
  kernel_value <- ifelse(distance == 0, sigma^2 + nugget, sigma^2 * exp(-distance^2 / (2 * theta^2)))
  kernel_value <- pmin(kernel_value, 1)
  return(kernel_value)
}


