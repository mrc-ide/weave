#' Create contiguous time folds
#'
#' Split unique time points into contiguous blocks, returning row indices
#' for each fold.
#'
#' @param df Data frame containing a time column.
#' @param time_col Name of the time column. Defaults to `"t"`.
#' @param K Number of folds.
#'
#' @return A list of integer vectors with row indices for each fold.
#' @examples
#' df <- data.frame(t = rep(1:10, each = 2))
#' make_time_folds(df, K = 2)
#' @export
make_time_folds <- function(df, time_col = "t", K = 5) {
  tt <- sort(unique(df[[time_col]]))
  blocks <- split(tt, cut(seq_along(tt), breaks = K, labels = FALSE))
  lapply(blocks, function(ttk) which(df[[time_col]] %in% ttk))
}

#' Create interleaved time folds
#'
#' Alternate time points across folds so each fold spans the full series.
#'
#' @inheritParams make_time_folds
#'
#' @return A list of integer vectors with row indices for each fold.
#' @examples
#' df <- data.frame(t = rep(1:10, each = 2))
#' make_time_folds_interleave(df, K = 3)
#' @export
make_time_folds_interleave <- function(df, time_col = "t", K = 5) {
  tt <- sort(unique(df[[time_col]]))
  folds_idx <- split(tt, (seq_along(tt) - 1) %% K)
  lapply(folds_idx, function(ttk) which(df[[time_col]] %in% ttk))
}

#' Mean predictive log-likelihood across folds
#'
#' Evaluate a hyperparameter vector by its mean predictive
#' log-likelihood over predefined folds under a Negative Binomial model.
#'
#' @param theta Numeric vector of length three: spatial length scale,
#'   periodic time scale, and long-term time scale.
#' @param obs_data Data frame with columns `y_obs` and `f_infer`.
#' @param coordinates Data frame of coordinates used by [fit()].
#' @param n Number of sites.
#' @param nt Number of time points.
#' @param folds List of integer vectors giving indices for each fold.
#'
#' @return Mean predictive log-likelihood across folds.
#' @keywords internal
fold_score_nb <- function(theta, obs_data, coordinates, n, nt, folds) {
  fold_ll <- numeric(length(folds))
  for (k in seq_along(folds)) {
    idx_hold <- folds[[k]]

    dat_mask <- obs_data
    dat_mask$y_obs[idx_hold] <- NA
    dat_mask$f_infer[idx_hold] <- NA

    fit_k <- fit(
      obs_data = dat_mask,
      coordinates = coordinates,
      hyperparameters = theta,
      n = n, nt = nt
    )

    y_true <- obs_data$y_obs[idx_hold]
    size <- fit_k$negbin_size[idx_hold]
    prob <- fit_k$negbin_prob[idx_hold]

    bad <- is.na(y_true) | is.na(size) | is.na(prob) | size <= 0 | prob <= 0 | prob >= 1
    ll <- rep(NA_real_, length(idx_hold))
    ll[!bad] <- stats::dnbinom(y_true[!bad], size = size[!bad], prob = prob[!bad], log = TRUE)
    ll[bad] <- -1e6

    fold_ll[k] <- mean(ll)
  }
  mean(fold_ll)
}

#' Tune hyperparameters by cross-validation
#'
#' Subsample sites, build time-based folds and optimise the mean predictive
#' log-likelihood using `optim` with `L-BFGS-B` bounds.
#'
#' @param obs_data Data frame with one row per site-time. Must contain at
#'   least `id`, `t`, `y_obs`, `mu_infer` and `f_infer`.
#' @param coordinates Data frame with site-level information including `id`.
#' @param id_col Name of the site identifier column.
#' @param time_col Name of the time column.
#' @param n_sites_sample Number of sites to randomly sample for tuning.
#' @param K_folds Number of cross-validation folds.
#' @param init Initial hyperparameter values.
#' @param lower Lower bounds for hyperparameters.
#' @param upper Upper bounds for hyperparameters.
#' @param seed Random seed for reproducibility.
#'
#' @return A list containing the best hyperparameters and optimisation details.
#' @examples
#' \dontrun{
#' init <- c(space = 3, t_per = 4, t_long = 12)
#' res <- tune_hyperparameters_optim(
#'   obs_data = obs_data,
#'   coordinates = coordinates,
#'   n_sites_sample = 20,
#'   K_folds = 10,
#'   init = init,
#'   lower = c(space = 0.01, t_per = 0.1, t_long = 1),
#'   upper = c(space = 100, t_per = 100, t_long = 24)
#' )
#' res$best_theta
#' }
#' @export
tune_hyperparameters_optim <- function(
    obs_data,
    coordinates,
    id_col = "id",
    time_col = "t",
    n_sites_sample = 40,
    K_folds = 5,
    init = c(space = 20, t_per = 4, t_long = 12),
    lower = c(space = 5, t_per = 1, t_long = 3),
    upper = c(space = 200, t_per = 20, t_long = 52),
    seed = 123
) {
  set.seed(seed)

  ids_all <- unique(obs_data[[id_col]])
  chosen_ids <- sample(ids_all, size = min(n_sites_sample, length(ids_all)), replace = FALSE)
  obs_sub <- obs_data[obs_data[[id_col]] %in% chosen_ids, , drop = FALSE]
  coords_sub <- coordinates[coordinates[[id_col]] %in% chosen_ids, , drop = FALSE]

  obs_sub <- obs_sub[order(obs_sub[[id_col]], obs_sub[[time_col]]), ]
  coords_sub <- coords_sub[order(coords_sub[[id_col]]), ]

  n <- length(unique(obs_sub[[id_col]]))
  nt <- length(unique(obs_sub[[time_col]]))

  folds <- make_time_folds_interleave(obs_sub, time_col = time_col, K = K_folds)

  obj <- function(par) {
    theta <- as.numeric(par)
    -fold_score_nb(theta, obs_data = obs_sub, coordinates = coords_sub, n = n, nt = nt, folds = folds)
  }

  fit_opt <- optim(
    par = init,
    fn = obj,
    method = "L-BFGS-B",
    lower = lower,
    upper = upper,
    control = list(trace = 1, maxit = 100)
  )

  best_theta <- as.numeric(fit_opt$par)
  names(best_theta) <- c("space_length_scale", "time_periodic_scale", "time_long_term_scale")

  list(
    best_theta = best_theta,
    best_cv_score = -fit_opt$value,
    convergence = fit_opt$convergence,
    message = fit_opt$message,
    folds_info = list(n_sites = n, n_times = nt, K_folds = K_folds, chosen_ids = chosen_ids)
  )
}

