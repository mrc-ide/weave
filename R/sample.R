#' Quick multivariate normal samples over two dimensions
#'
#' @param space Space kernel matrix
#' @param time  Time kernel matrix
quick_mvnorm <- function(space, time){
  n_sites <- nrow(space)
  n_times <- nrow(time)

  L_s <- chol(space)
  L_t <- chol(time)

  # Generate a matrix of i.i.d. standard normal variables
  W <- matrix(rnorm(n_sites * n_times), nrow = n_sites, ncol = n_times)

  browser()
  # Construct the sample matrix
  Z <- L_s %*% W %*% L_t

  # Flatten the matrix
  z <- as.vector(Z)

  # We need to reorder 'z' so that times vary faster than sites
  z_matrix <- matrix(z, nrow = n_sites, ncol = n_times)
  z_reordered <- as.vector(t(z_matrix))

  # Add the mean
  z <- as.vector(z_reordered)

  return(z)
}
