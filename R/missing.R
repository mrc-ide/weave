# generate a series of binary draws with probability p_one then tend to cluster
# with switch rate p_switch
generate_clustered_binary <- function(n, p_one = 0.5, p_switch = 0.1) {
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
