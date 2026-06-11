#' Complete site-time combinations
#'
#' @description
#' Adds rows for all combinations of site identifiers and time `t`.
#'
#' @param data A data frame containing site identifiers, time `t`,
#'   counts `n`, and coordinates `lat` and `lon`.
#' @param ... Columns identifying sites passed to [dplyr::group_by()]
#'   (unquoted).
#'
#' @return A data frame with missing site-time combinations filled in and
#'   `n` set to `NA`.
data_complete <- function(data, ...){
  site_names <- rlang::ensyms(...)

  # Complete all site x time combinations
  data <- data |>
    tidyr::complete(tidyr::nesting(!!!site_names), .data$t, fill = list(n = NA)) |>
    dplyr::group_by(!!!site_names) |>
    dplyr::mutate(
      lat = dplyr::first(.data$lat, na_rm = TRUE),
      lon = dplyr::first(.data$lon, na_rm = TRUE)
    ) |>
    dplyr::ungroup()

  return(data)
}

#' Drop sites with missing data
#'
#' @description
#' Removes sites that cannot be modelled and reports them. Sites with no
#' observed counts (all `n` missing) or without coordinates (all `lat`/`lon`
#' missing) are always dropped. Sites whose observed counts are all zero are
#' dropped only when `drop_zero = TRUE`.
#'
#' Under the `log1p` per-site-centred rate model an all-zero site is valid
#' low-rate data, not a defect, so it is retained by default; set
#' `drop_zero = TRUE` if all-zero facilities should be treated as reporting
#' artefacts and removed.
#'
#' @param data A data frame containing site identifiers, time `t`,
#'   counts `n`, and coordinates `lat` and `lon`.
#' @param ... Columns identifying sites passed to [dplyr::group_by()]
#'   (unquoted).
#' @param drop_zero Logical; also drop sites whose observed counts sum to zero
#'   (default `FALSE`).
#'
#' @return A data frame with problem sites removed.
data_missing <- function(data, ..., drop_zero = FALSE){
  site_names <- rlang::enquos(...)

  sites_to_drop_n <- data |>
    dplyr::group_by(!!!site_names) |>
    dplyr::filter(
      all(is.na(.data$n))
    ) |>
    dplyr::distinct(!!!site_names)

  sites_to_drop_lat_lon <- data |>
    dplyr::group_by(!!!site_names) |>
    dplyr::filter(
      all(is.na(.data$lat)) | all(is.na(.data$lon))
    ) |>
    dplyr::distinct(!!!site_names)

  sites_to_drop <- sites_to_drop_n |>
    dplyr::bind_rows(sites_to_drop_lat_lon)

  if(drop_zero){
    sites_to_drop_no_cases <- data |>
      dplyr::group_by(!!!site_names) |>
      dplyr::filter(
        sum(.data$n, na.rm = TRUE) == 0
      ) |>
      dplyr::distinct(!!!site_names)

    sites_to_drop <- sites_to_drop |>
      dplyr::bind_rows(sites_to_drop_no_cases)
  }

  sites_to_drop <- sites_to_drop |>
    dplyr::distinct()

  if(nrow(sites_to_drop) > 0){
    msg <- if(drop_zero){
      "Sites dropped (all data missing, missing coordinates, or all counts = 0): "
    } else {
      "Sites dropped (all data missing or missing coordinates): "
    }
    cat(msg)
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

#' Order data and assign identifiers
#'
#' @description
#' Arranges data by site and time and creates a factor `id` per site.
#'
#' @param data A data frame containing site identifiers and time `t`.
#' @param ... Columns identifying sites passed to [dplyr::group_by()]
#'   (unquoted).
#'
#' @return A data frame ordered by site and time with an `id` column.
data_order_index <- function(data, ...){
  site_names <- rlang::enquos(...)

  data <- data |>
    dplyr::arrange(
      !!!site_names,
      t
    ) |>
    dplyr::group_by(!!!site_names) |>
    dplyr::mutate(
      id = dplyr::cur_group_id(),
      id = factor(.data$id)
    ) |>
    dplyr::ungroup()

  return(data)
}


#' Process raw epidemiological data for the GP model
#'
#' @description
#' Validates and prepares input data into the shape consumed by
#' [infer_kernel_params()] and [gp_predict()]: it completes the site-by-time
#' grid, drops sites that cannot be modelled, assigns a factor site `id`, and
#' returns the observations and coordinates as separate, ready-to-use frames.
#'
#' @param data A data frame containing site identifiers, time `t`,
#'   counts `n`, and coordinates `lat` and `lon`.
#' @param ... Columns identifying sites passed to [dplyr::group_by()]
#'   (unquoted).
#' @param drop_zero Logical; passed to [data_missing()] -- also drop sites whose
#'   observed counts sum to zero (default `FALSE`).
#'
#' @return A list with three elements ready to pass to the model functions:
#'   \describe{
#'     \item{`obs_data`}{Observations with the site-identifier columns, the
#'       factor `id`, time `t`, and the count column `y_obs` (`NA` where
#'       missing).}
#'     \item{`coordinates`}{One row per site with `id`, `lon` and `lat`.}
#'     \item{`nt`}{The number of time points.}
#'   }
#' @export
data_process <- function(data, ..., drop_zero = FALSE){
  if(!all(c("t", "n", "lat", "lon") %in% colnames(data))){
    stop("Input data must include the following columns: t, n, lat and lon"
    )
  }

  if("id" %in% colnames(data)){
    stop("The column name 'id' is protected and cannot be used in the input data")
  }

  processed <- data |>
    data_complete(...) |>
    data_missing(..., drop_zero = drop_zero) |>
    data_order_index(...)

  site_names <- rlang::enquos(...)

  obs_data <- processed |>
    dplyr::select(!!!site_names, "id", "t", y_obs = "n") |>
    as.data.frame()

  coordinates <- processed |>
    dplyr::distinct(.data$id, .data$lon, .data$lat) |>
    as.data.frame()

  list(
    obs_data    = obs_data,
    coordinates = coordinates,
    nt          = length(unique(processed$t))
  )
}
