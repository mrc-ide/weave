infer_space_kernel_params <- function(data, spatial_distance){
  space_cor <- expand.grid(id1 = unique(data$id), id2 = unique(data$id), t = unique(data$t)) |>
    dplyr::filter(id1 < id2) |>
    dplyr::left_join(dplyr::select(data, id, t, n), by = c("id1" = "id", "t" = "t")) |>
    dplyr::left_join(dplyr::select(data, id, t, n), by = c("id2" = "id", "t" = "t")) |>
    dplyr::filter(!is.na(n.x), !is.na(n.y)) |>
    dplyr::summarise(
      cor = cor(n.x, n.y),
      .by = c("id1", "id2")
    ) |>
    dplyr::mutate(
      distance = purrr::map2_dbl(id1, id2, ~ spatial_distance[.x, .y])
    )

  # Fitting sigma to empirical correlations
  fit_sigma <- function(sigma) {
    predicted_correlations <- rbf_kernel(space_cor$distance, sigma)
    sum((predicted_correlations - space_cor$cor)^2)
  }

  # Optimise sigma to minimise the difference between empirical and predicted correlations
  optimal_sigma <- optimise(f = fit_sigma, interval = c(0, 100))$minimum
  return(
    list(
      sigma = optimal_sigma
    )
  )
}

infer_time_kernel_params <- function(data, period){
  time_cor <- expand.grid(t1 = unique(data$t), t2 = unique(data$t), id = unique(data$id)) |>
    dplyr::filter(t2 < t1) |>
    dplyr::left_join(dplyr::select(data, id, t, n), by = c("t1" = "t", "id" = "id")) |>
    dplyr::left_join(dplyr::select(data, id, t, n), by = c("t2" = "t", "id" = "id")) |>
    dplyr::filter(!is.na(n.x), !is.na(n.y)) |>
    dplyr::summarise(
      cor = cor(n.x, n.y),
      .by = c("t1", "t2")
    ) |>
    dplyr::mutate(
      time_distance = t2 - t1
    )

  # Fitting sigma to empirical correlations,
  fit_sigma <- function(params, period) {
    predicted_correlations <- periodic_kernel(
      time_cor$time_distance,
      periodic_scale = params[1],
      long_term_scale = params[2],
      period = period)
    sum((predicted_correlations - time_cor$cor)^2)
  }

  # Optimise sigma to minimise the difference between empirical and predicted correlations
  optim_result_time <- optim(
    par = c(1, 1),
    fn = fit_sigma,
    method = "L-BFGS-B",
    lower = c(0.01, 0.01),
    period = 12
  )

  return(
    list(
      periodic_scale = optim_result_time$par[1],
      long_term_scale = optim_result_time$par[2],
      period = period
    )
  )
}
