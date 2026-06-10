# =============================================================================
# Data utilities
#
# Two pipelines coexist here:
#
#   data_process(data, ...)  -- cleans a tidy user-supplied data frame:
#                               completes the site x time grid, drops empty
#                               sites, assigns integer id codes. Designed for
#                               real-world inputs with column `n` for counts.
#
#   build_design(data)       -- packs an already-tidy data frame into the
#                               compact list that the fitters consume:
#                               flat y vector, observed-index vector,
#                               coordinates table, etc. Works on either
#                               `y_obs` (simulation output) or `n` (user data).
#
# Both produce data ordered so that TIME varies fastest within site -- the
# convention assumed by kron_mv() and the rest of the linear algebra.
# =============================================================================


data_complete <- function(data, ...){
  site_names <- rlang::enquos(...)

  # Complete all site x time combinations
  data <- data |>
    tidyr::complete(!!!site_names, .data$t, fill = list(n = NA)) |>
    dplyr::group_by(!!!site_names)|>
    dplyr::mutate(
      lat = dplyr::first(.data$lat, na_rm = TRUE),
      lon = dplyr::first(.data$lon, na_rm = TRUE)
    ) |>
    dplyr::ungroup()

  return(data)
}

data_missing <- function(data, ...){
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

  sites_to_drop_no_cases <- data |>
    dplyr::group_by(!!!site_names) |>
    dplyr::filter(
      sum(.data$n, na.rm = TRUE) == 0
    ) |>
    dplyr::distinct(!!!site_names)

  sites_to_drop <- sites_to_drop_n |>
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

# Per-site starting offset on the log scale -- used as the initial value of
# mu_s in the Gibbs sampler and as the centring offset in the deterministic
# fit. Kept robust to all-NA sites by filling with the global mean of the
# computable values.
data_mu_init <- function(data, ...) {
  site_names <- rlang::enquos(...)

  data |>
    dplyr::group_by(!!!site_names) |>
    dplyr::mutate(
      mu_init = log(mean(.data$n, na.rm = TRUE))
    ) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      mu_init = ifelse(
        is.na(.data$mu_init) | !is.finite(.data$mu_init),
        mean(.data$mu_init[is.finite(.data$mu_init)], na.rm = TRUE),
        .data$mu_init
      )
    )
}

data_order_index <- function(data, ...){
  site_names <- rlang::enquos(...)

  data <- data |>
    dplyr::arrange(
      !!!site_names,
      t
    ) |>
    dplyr::group_by(!!!site_names)|>
    dplyr::mutate(
      id = dplyr::cur_group_id(),
      id = factor(.data$id)
      ) |>
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
    data_mu_init(...) |>
    data_order_index(...)

  return(data)
}


#' Pack a tidy data frame into the compact design list the fitter consumes
#'
#' Given a long data frame with columns `id` (factor), `t`, `lat`, `lon`, and
#' a count column (`y_obs` if present, else `n`), returns the fields needed
#' by `fit_bayes()`: a flat length-N vector of counts in the times-vary-
#' fastest order, the indices of the observed (non-NA) cells, the per-site
#' coordinates, and the per-site initial intercept `mu_init`.
#'
#' This is the single source of truth for the data shape passed to the
#' sampler; `fit_bayes()` calls it internally if given a raw data frame.
#'
#' @param data Tidy data frame; must include columns id, t, lat, lon and
#'   one of (y_obs, n).
#' @return A list with elements:
#'   - `y_full`: numeric length n*nt, NA at unobserved cells
#'   - `y_obs`:  numeric length m, observed counts only
#'   - `obs_idx`: integer length m, positions of observed cells in the full grid
#'   - `n`, `nt`: scalar dimensions
#'   - `N`: n * nt
#'   - `site_idx_full`: integer length n*nt, site index per cell
#'   - `site_idx_obs`:  integer length m, site index per observed cell
#'   - `time_idx_full`: integer length n*nt, time index per cell
#'   - `coords`: data frame with one row per site (id, lat, lon)
#'   - `mu_init`: numeric length n, starting site intercept on log scale
#' @export
build_design <- function(data) {
  count_col <- if ("y_obs" %in% names(data)) "y_obs" else "n"
  if (!count_col %in% names(data)) {
    stop("`data` must contain a count column named `y_obs` or `n`.")
  }
  required <- c("id", "t", "lat", "lon")
  missing  <- setdiff(required, names(data))
  if (length(missing) > 0) {
    stop("`data` is missing required column(s): ", paste(missing, collapse = ", "))
  }

  d <- data |>
    dplyr::arrange(.data$id, .data$t)

  ids   <- as.integer(d$id)
  times <- as.integer(d$t)
  y_full <- d[[count_col]]

  n  <- length(unique(ids))
  nt <- length(unique(times))
  N  <- n * nt
  if (nrow(d) != N) {
    stop("Data is not a complete site x time grid: expected ", N,
         " rows, got ", nrow(d), ". Run data_process() first.")
  }

  obs_idx <- which(!is.na(y_full))
  y_obs   <- y_full[obs_idx]

  coords <- d |>
    dplyr::distinct(.data$id, .data$lat, .data$lon) |>
    dplyr::arrange(.data$id)
  if (nrow(coords) != n) {
    stop("Inconsistent (lat, lon) per site: ", nrow(coords),
         " unique coords for ", n, " sites.")
  }

  # Initial site intercept: log mean of observed counts in that site, with a
  # global-mean fallback for sites that happen to be all-NA after subsetting.
  mu_init <- vapply(seq_len(n), function(s) {
    yy <- y_obs[ids[obs_idx] == s]
    if (length(yy) == 0 || all(yy == 0)) NA_real_ else log(mean(yy))
  }, numeric(1))
  if (any(is.na(mu_init))) {
    mu_init[is.na(mu_init)] <- mean(mu_init, na.rm = TRUE)
  }

  list(
    y_full        = y_full,
    y_obs         = y_obs,
    obs_idx       = obs_idx,
    n             = n,
    nt            = nt,
    N             = N,
    site_idx_full = ids,
    site_idx_obs  = ids[obs_idx],
    time_idx_full = times,
    coords        = coords,
    mu_init       = mu_init
  )
}
