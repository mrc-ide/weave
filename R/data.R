data_complete <- function(data, ...){
  site_names <- rlang::enquos(...)

  # Complete all site x time combinations
  data <- data |>
    tidyr::complete(!!!site_names, t, fill = list(n = NA)) |>
    dplyr::group_by(!!!site_names)|>
    dplyr::mutate(
      lat = dplyr::first(lat, na_rm = TRUE),
      lon = dplyr::first(lon, na_rm = TRUE)
    ) |>
    dplyr::ungroup()

  return(data)
}

data_missing <- function(data, ...){
  site_names <- rlang::enquos(...)

  sites_to_drop_n <- data |>
    dplyr::group_by(!!!site_names) |>
    dplyr::filter(
      all(is.na(n))
    ) |>
    dplyr::distinct(!!!site_names)

  sites_to_drop_lat_lon <- data |>
    dplyr::group_by(!!!site_names) |>
    dplyr::filter(
      all(is.na(lat)) | all(is.na(lon))
    ) |>
    dplyr::distinct(!!!site_names)

  sites_to_drop_no_cases <- data |>
    dplyr::group_by(!!!site_names) |>
    dplyr::filter(
      sum(n, na.rm = TRUE) == 0
    ) |>
    dplyr::distinct(!!!site_names)

  sites_to_drop <-
    sites_to_drop_n |>
    dplyr::bind_rows(sites_to_drop_lat_lon) |>
    dplyr::bind_rows(sites_to_drop_no_cases) |>
    dplyr::distinct()


  if(nrow(sites_to_drop) > 0){
    cat("Sites dropped as all data missing, or all counts = 0: ")
    knitr::kable(sites_to_drop, format = "pipe", align = "c") |>
      print()

    data <- data |>
      dplyr::anti_join(
        sites_to_drop,
        by = dplyr::join_by(...)
      )
  }

  return(data)
}

data_observed_summmary <- function(data, ...){
  site_names <- rlang::enquos(...)

  data <- data |>
    dplyr::group_by(!!!site_names) |>
    dplyr::mutate(
      observed_mu = log(mean(n, na.rm = T)),
      observed_sigmasq = log((sd(n, na.rm = T) / mean(n, na.rm = T))^2 + 1)
    ) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      observed_mu = ifelse(is.na(observed_mu), mean(observed_mu, na.rm = T), observed_mu),
      observed_sigmasq = ifelse(is.na(observed_sigmasq), mean(observed_sigmasq, na.rm = T), observed_sigmasq)
    )

  return(data)
}

data_initial_par <- function(data){
  data <- data |>
    dplyr::mutate(
      start_par = log(n + 1),
      start_par = ifelse(is.na(n), observed_mu, start_par)
    )

  return(data)
}

data_order_index <- function(data, ...){
  site_names <- rlang::enquos(...)

  data <- data |>
    dplyr::arrange(
      !!!site_names,
      t
    ) |>
    dplyr::group_by(!!!site_names)|>
    dplyr::mutate(id = dplyr::cur_group_id()) |>
    dplyr::ungroup()

  return(data)
}


data_process <- function(data, ...){
  if(!all(c("t", "n", "lat", "lon") %in% colnames(data))){
    stop("Input data must include the following columns: t, n, lat and lon"
    )
  }

  if("id" %in% colnames(data)){
    stop("The column name 'id' is protected and cannot be used in the input data")
  }

  data <- data |>
    data_complete(...) |>
    data_missing(...) |>
    data_observed_summmary(...) |>
    data_initial_par() |>
    data_order_index() |>
    data_order_index(...)

  return(data)
}
