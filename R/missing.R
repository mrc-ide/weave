#' Generate clustered binary sequence
#'
#' Creates a binary vector of length `n` where consecutive values tend to
#' cluster. Values are freshly drawn using `p_one` when a new cluster starts and
#' otherwise repeat the previous value.
#'
#' @param n Number of draws to generate.
#' @param p_one Probability that a new cluster begins with a one. Also used for
#'   the first draw.
#' @param p_switch Probability of switching to a new cluster at each step.
#'
#' @return Numeric vector of 0s and 1s.
#'
#' @details Randomness is generated using [stats::rbinom()] and
#'   [stats::runif()]. Set a seed via [set.seed()] for reproducible results.
#'
#' @examples
#' set.seed(1)
#' generate_clustered_binary(5, 0.2, 0.1)
#'
#' @export
generate_clustered_binary <- function(n, p_one, p_switch) {
  result <- numeric(n)
  result[1] <- rbinom(1, 1, p_one)
  for (i in 2:n) {
    if (runif(1) < p_switch) {
      result[i] <- rbinom(1, 1, p_one)
    } else {
      result[i] <- result[i - 1]
    }
  }
  return(result)
}

#' Add clustered missingness to counts
#'
#' Applies a clustered missingness pattern to a count column, replacing selected
#' values with `NA`. The original counts are preserved in `true_n`.
#'
#' @param data Data frame containing column `n`.
#' @param p_one Probability that a new missingness cluster starts with a missing
#'   value.
#' @param p_switch Probability of switching between missing and observed
#'   clusters.
#'
#' @return Data frame with modified `n` and a `true_n` column.
#'
#' @details Uses [generate_clustered_binary()] for randomness. Results are
#'   stochastic unless a seed is set via [set.seed()] before calling.
#'
#' @export
add_missingness <- function(data, p_one, p_switch) {
  data |>
    dplyr::mutate(
      true_n = .data$n,
      missing = generate_clustered_binary(dplyr::n(), p_one, p_switch),
      n = ifelse(missing == 1, NA, .data$n)
    ) |>
    dplyr::select(-dplyr::all_of("missing"))
}
