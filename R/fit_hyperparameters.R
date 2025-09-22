#' Fit hyperparmeters
#'
#' @param obs_data Data frame of observations with at least `id`, `y_obs`,
#'   `mu_infer`, and `f_infer` columns.
#' @param coordinates Data frame of site coordinates containing an `id` column
#'   that matches the site identifiers in `obs_data`.
#' @param nt Number of time points per site.
#' @param period Period used for the temporal kernel.
#' @param n_sites Number of sites to sample for optimisation.
#' @param mask_prop Proportion of observed points per site to hold out.
#' @param verbose Logical; whether to print objective values during
#'   optimisation.
#' @param par0 Initial hyperparameter values passed to `optim()`.
#' @param lower Lower bounds for the hyperparameters passed to `optim()`.
#' @param upper Upper bounds for the hyperparameters passed to `optim()`.
#'
#' @return A list as returned by `stats::optim()`.
#' @export
fit <- function(obs_data, coordinates, nt, period, n_sites, mask_prop = 0.2, verbose = FALSE, par0 = c(1, 5, 100), lower = c(1e-4, 0.8, 52 * 1.5), upper = c(2, 10, 500)) {
  if (!"id" %in% names(coordinates)) {
    stop("`coordinates` must contain an `id` column.", call. = FALSE)
  }

  ids <- sample(unique(obs_data$id), n_sites)
  missing_ids <- setdiff(ids, coordinates$id)
  if (length(missing_ids)) {
    stop(
      "Coordinates missing for sampled site IDs: ",
      paste(missing_ids, collapse = ", "),
      call. = FALSE
    )
  }

  fitting_data <- obs_data[obs_data$id %in% ids, ]
  fitting_coords <- coordinates[match(ids, coordinates$id), , drop = FALSE]

  # Hold-out mask: choose a fraction of NON-NA y_obs per site; keep at least 1 observed point/site
  split_idx <- split(seq_len(nrow(fitting_data)), fitting_data$id)
  held_idx <- unlist(lapply(split_idx, function(ix) {
    cand <- ix[!is.na(fitting_data$y_obs[ix])]
    n_c  <- length(cand)
    if (n_c <= 1L) return(integer(0))                 # skip sites with 0/1 observed points
    k <- max(1L, floor(n_c * mask_prop))              # desired hold-out count
    k <- min(k, n_c - 1L)                             # ensure at least 1 observed remains
    sample(cand, k)
  }), use.names = FALSE)

  fitting_data$y_comp <- fitting_data$y_obs
  if (length(held_idx)) fitting_data$y_obs[held_idx] <- NA

  fit_f <- function(par) {
    state <- gp_build_state(
      obs_data       = fitting_data,
      coordinates    = fitting_coords,
      hyperparameters = par,
      n              = n_sites,
      nt             = nt,
      period         = period
    )
    lam_hat <- gp_posterior_mean(state)
    # Score ONLY the held-out cells
    ll <- sum(dpois(fitting_data$y_comp[held_idx], lam_hat[held_idx], log = TRUE))
    if(verbose){
      print(ll)
    }
    return(ll)
  }

  optim(
    par = par0,
    fn = fit_f,
    method = "L-BFGS-B",
    lower = lower,
    upper = upper,
    control = list(
      maxit = 100,
      fnscale = -1,
      factr = 1e4,
      pgtol = 1e-8
    )
  )
}
