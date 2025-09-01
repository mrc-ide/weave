infer_space_kernel_params <- function(data, nt, n, plot = FALSE, max_pairs = 1000){

  zmat <- matrix(data$z_infer, nrow = nt, ncol = n, byrow = FALSE)

  # Pairwise complete correlation/covariance
  cov_mat <- cov(zmat, use = "pairwise.complete.obs")
  corr_mat <- cov2cor(cov_mat)
  corr_mat[corr_mat < 0] <- 0

  dist_mat <- get_spatial_distance(unique(data[,c("lon", "lat")]))
  space_cor <- data.frame(distance = as.vector(dist_mat), cor = as.vector(corr_mat))

  # Downsample pairs if needed
  if(nrow(space_cor) > max_pairs){
    space_cor <- space_cor[sample.int(nrow(space_cor), max_pairs), ]
  }

  # Fitting theta to empirical correlations
  fit_sigma <- function(theta, space_cor) {
    predicted_correlations <- rbf_kernel(space_cor$distance, theta)
    sum((predicted_correlations - space_cor$cor)^2, na.rm = TRUE)
  }

  optimal_theta <- optimise(f = fit_sigma, interval = c(0, 100), space_cor = space_cor)$minimum

  if(plot){
    pred <- data.frame(
      distance = seq(0, max(space_cor$distance), length.out = 100)
    ) |>
      dplyr::mutate(
        y = rbf_kernel(distance, optimal_theta)
      )

    cor_plot <- ggplot2::ggplot() +
      ggplot2::geom_point(data = space_cor, ggplot2::aes(x = distance, y = cor), alpha = 0.5) +
      ggplot2::geom_line(data = pred, ggplot2::aes(x = distance, y = y), col = "deeppink") +
      ggplot2::xlab("Spatial distance") +
      ggplot2::ylab("") +
      ggplot2::theme_bw()
    print(cor_plot)
  }

  list(length_scale = optimal_theta)
}

infer_time_kernel_params <- function(data, period, nt, n, plot = FALSE, max_pairs = 1000){

  fmat <- t(matrix(data$z_infer, nrow = nt, ncol = n, byrow = FALSE))
  cov_mat <- cov(fmat, use = "pairwise.complete.obs")
  corr_mat <- cov2cor(cov_mat)

  lags <- 0:(nt - 1)
  mean_corr_by_lag <- numeric(length(lags))

  for (h in lags) {
    idx_i <- 1:(nt - h)
    idx_j <- idx_i + h
    correlations <- mapply(function(i, j) corr_mat[i, j], idx_i, idx_j)
    if(length(correlations) > max_pairs){
      correlations <- sample(correlations, max_pairs)
    }
    mean_corr_by_lag[h + 1] <- mean(correlations, na.rm = TRUE)
  }

  time_cor <- data.frame(time_distance = 1:nt, cor = mean_corr_by_lag)

  fit_sigma <- function(params, period, time_cor) {
    predicted_correlations <- periodic_kernel(x = time_cor$time_distance,
                                              alpha = params[1], period = period) *
      rbf_kernel(x = time_cor$time_distance,
                 theta = params[2])
    sum((predicted_correlations - time_cor$cor)^2, na.rm = TRUE)
  }

  optim_result_time <- optim(
    par = c(1, 200),
    fn = fit_sigma,
    method = "L-BFGS-B",
    lower = c(0.1, 52),
    upper = c(10, 52 * 10000),
    period = period,
    time_cor = time_cor
  )

  if(plot){
    pred <- data.frame(
      distance = seq(0, max(time_cor$time_distance), length.out = 100)
    ) |>
      dplyr::mutate(
        y = periodic_kernel(x = distance,
                            alpha = optim_result_time$par[1],
                            period = period) *
          rbf_kernel(x = distance,
                     theta = optim_result_time$par[2])
      )

    cor_plot <- ggplot2::ggplot() +
      ggplot2::geom_point(data = time_cor, ggplot2::aes(x = time_distance, y = cor), alpha = 0.5) +
      ggplot2::geom_line(data = pred, ggplot2::aes(x = distance, y = y), col = "deeppink") +
      ggplot2::xlab("Temporal distance") +
      ggplot2::ylab("") +
      ggplot2::theme_bw()
    print(cor_plot)
  }

  list(
    periodic_scale = optim_result_time$par[1],
    long_term_scale = optim_result_time$par[2],
    period = period
  )
}
