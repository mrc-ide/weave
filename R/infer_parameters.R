infer_space_kernel_params <- function(data, plot = FALSE){

  data$z_t_hat = log(data$y + 1)

  spatial_distance <- get_spatial_distance(unique(data[,c("lon", "lat")]))

  ids <- unique(data$id)
  times <- unique(data$t)

  space_cor <- expand.grid(id1 = ids, id2 = ids, t = times) |>
    dplyr::filter(as.numeric(id1) < as.numeric(id2)) |>
    dplyr::left_join(dplyr::select(data, id, t, z_t_hat), by = c("id1" = "id", "t" = "t")) |>
    dplyr::left_join(dplyr::select(data, id, t, z_t_hat), by = c("id2" = "id", "t" = "t")) |>
    dplyr::filter(!is.na(z_t_hat.x), !is.na(z_t_hat.y)) |>
    dplyr::filter(dplyr::n() > 5, .by = c(id1, id2)) |>
    dplyr::summarise(
      cor = cor(z_t_hat.x, z_t_hat.y),
      .by = c("id1", "id2")
    ) |>
    dplyr::mutate(
      cor = ifelse(cor < 0, 0, cor),
      distance = purrr::map2_dbl(id1, id2, ~ spatial_distance[.x, .y])
    ) |>
    dplyr::filter(!is.na(cor))


  # Fitting theta to empirical correlations
  fit_sigma <- function(theta, space_cor) {
    predicted_correlations <- rbf_kernel(space_cor$distance, theta)
    sum((predicted_correlations - space_cor$cor)^2)
  }

  # Optimise theta to minimise the difference between empirical and predicted correlations
  optimal_theta <- optimise(f = fit_sigma, interval = c(0, 100), space_cor = space_cor)$minimum

  if(plot){
    pred <- data.frame(
      distance = seq(0, max(space_cor$distance), length.out = 100)
    ) |>
      dplyr::mutate(
        y = rbf_kernel(distance, optimal_theta)
      )

    cor_plot <- ggplot() +
      geom_point(data = space_cor, aes(x = distance, y = cor), alpha = 0.5) +
      geom_line(data = pred, aes(x = distance, y = y), col = "deeppink") +
      xlab("Spatial distance") +
      ylab("") +
      theme_bw()
    print(cor_plot)
  }

  return(
    list(
      theta = optimal_theta
    )
  )
}

infer_time_kernel_params <- function(data, period, plot = FALSE){

  data$z_t_hat = log(data$y + 1)

  time_cor <- expand.grid(t1 = unique(data$t), t2 = unique(data$t), id = unique(data$id)) |>
    dplyr::filter(t1 < t2) |>
    dplyr::left_join(dplyr::select(data, id, t, z_t_hat), by = c("t1" = "t", "id" = "id")) |>
    dplyr::left_join(dplyr::select(data, id, t, z_t_hat), by = c("t2" = "t", "id" = "id")) |>
    dplyr::filter(!is.na(z_t_hat.x), !is.na(z_t_hat.y)) |>
    dplyr::filter(dplyr::n() > 5, .by = c(t1, t2)) |>
    dplyr::summarise(
      cor = cor(z_t_hat.x, z_t_hat.y),
      .by = c("t1", "t2")
    ) |>
    dplyr::mutate(
      time_distance = abs(t2 - t1)
    )

  # Fitting sigma to empirical correlations,
  fit_sigma <- function(params, period, time_cor) {
    predicted_correlations <- periodic_kernel(x = time_cor$time_distance, alpha = params[1], period = period) *
      rbf_kernel(x = time_cor$time_distance, theta = params[2])
    sum((predicted_correlations - time_cor$cor)^2)
  }

  # Optimise sigma to minimise the difference between empirical and predicted correlations
  optim_result_time <- optim(
    par = c(1, 1),
    fn = fit_sigma,
    method = "L-BFGS-B",
    lower = c(0.1, 10),
    upper = c(10, 52 * 10000),
    period = period,
    time_cor = time_cor
  )
  if(plot){
    pred <- data.frame(
      distance = seq(0, max(time_cor$time_distance), length.out = 100)
    ) |>
      dplyr::mutate(
        y = periodic_kernel(x = distance, alpha = optim_result_time$par[1], period = period) *
          rbf_kernel(x =distance, theta = optim_result_time$par[2])
      )

    cor_plot <- ggplot() +
      geom_point(data = time_cor, aes(x = time_distance, y = cor), alpha = 0.5) +
      geom_line(data = pred, aes(x = distance, y = y), col = "deeppink") +
      xlab("Temporal distance") +
      ylab("") +
      theme_bw()
    print(cor_plot)
  }

  return(
    list(
      alpha = optim_result_time$par[1],
      beta = optim_result_time$par[2],
      period = period
    )
  )
}
