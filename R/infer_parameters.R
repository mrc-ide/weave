#' Estimate spatial RBF length scale from empirical correlations
#'
#' Computes cross-site correlations (pairwise complete) of `z_infer` across
#' time, pairs them with inter-site distances, and fits an RBF (`exp(-d^2/(2θ^2))`)
#' by least squares to obtain the spatial length scale.
#'
#' @param data Data frame with columns `z_infer`, `lon`, `lat`.
#' @param nt Integer, number of time points (used to reshape).
#' @param n Integer, number of sites (used to reshape).
#' @param plot Logical; if `TRUE`, show empirical vs fitted correlation–distance curve.
#' @param max_pairs Integer; optional downsampling of site pairs for speed.
#'
#' @return A list with `length_scale`.
#' @export
infer_space_kernel_params <- function(data, nt, n, plot = FALSE, max_pairs = 1000){

  zmat <- matrix(data$z_infer, nrow = nt, ncol = n, byrow = FALSE)

  # Pairwise complete correlation/covariance
  cov_mat <- stats::cov(zmat, use = "pairwise.complete.obs")
  corr_mat <- stats::cov2cor(cov_mat)
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

  optimal_theta <- stats::optimise(f = fit_sigma, interval = c(0, 100), space_cor = space_cor)$minimum

  if(plot){
    pred <- data.frame(
      distance = seq(0, max(space_cor$distance), length.out = 100)
    ) |>
      dplyr::mutate(
        y = rbf_kernel(.data$distance, optimal_theta)
      )

    cor_plot <- ggplot2::ggplot() +
      ggplot2::geom_point(data = space_cor, ggplot2::aes(x = .data$distance, y = .data$cor), alpha = 0.5) +
      ggplot2::geom_line(data = pred, ggplot2::aes(x = .data$distance, y = .data$y), col = "deeppink") +
      ggplot2::xlab("Spatial distance") +
      ggplot2::ylab("") +
      ggplot2::theme_bw()
    print(cor_plot)
  }

  list(length_scale = optimal_theta)
}

#' Estimate time-kernel (periodic × RBF) parameters from empirical ACF
#'
#' Forms an across-site time–time correlation matrix of `z_infer` (pairwise
#' complete), averages super-diagonals by lag to get an empirical ACF, then fits
#' a product kernel `k(h)=k_per(h; alpha, period) * k_rbf(h; theta)` by
#' least-squares.
#'
#' @param data Data frame with column `z_infer` and facetting columns `id`, `t`
#'   (times fastest within site).
#' @param period Numeric period used in the periodic kernel.
#' @param nt Integer, number of time points.
#' @param n Integer, number of sites.
#' @param plot Logical; if `TRUE`, show empirical vs fitted correlation–lag curve.
#' @param max_pairs Integer; optional downsampling of per-lag pairs for speed.
#'
#' @return A list with `periodic_scale`, `long_term_scale` and `period`.
#' @export
infer_time_kernel_params <- function(data, period, nt, n, plot = FALSE, max_pairs = 1000){

  fmat <- t(matrix(data$z_infer, nrow = nt, ncol = n, byrow = FALSE))
  cov_mat <- stats::cov(fmat, use = "pairwise.complete.obs")
  corr_mat <- stats::cov2cor(cov_mat)

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

  time_cor <- data.frame(time_distance = lags, cor = mean_corr_by_lag)

  fit_sigma <- function(params, period, time_cor) {
    predicted_correlations <- periodic_kernel(x = time_cor$time_distance,
                                              alpha = params[1], period = period) *
      rbf_kernel(x = time_cor$time_distance,
                 theta = params[2])
    sum((predicted_correlations - time_cor$cor)^2, na.rm = TRUE)
  }

  optim_result_time <- stats::optim(
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
        y = periodic_kernel(x = .data$distance,
                            alpha = optim_result_time$par[1],
                            period = period) *
          rbf_kernel(x = .data$distance,
                     theta = optim_result_time$par[2])
      )

    cor_plot <- ggplot2::ggplot() +
      ggplot2::geom_point(data = time_cor, ggplot2::aes(x = .data$time_distance, y = .data$cor), alpha = 0.5) +
      ggplot2::geom_line(data = pred, ggplot2::aes(x = .data$distance, y = .data$y), col = "deeppink") +
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
