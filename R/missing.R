# generate a series of binary draws with probability p_one then tend to cluster
# with switch rate p_switch
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

add_missingness <- function(data, p_one, p_switch){
  data |>
    dplyr::mutate(
      true_n = n,
      missing = generate_clustered_binary(dplyr::n(), p_one, p_switch),
      n = ifelse(missing == 1, NA, n)
    ) |>
  dplyr::select(- missing)
}
