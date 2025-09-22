#' Fit hyperparameters
#'
#' This helper picks a small group of monitoring sites, hides a few of their
#' observed counts, and tweaks the model settings until those withheld values
#' are predicted well.
#'
#' Under the hood it samples `n_sites` IDs, masks a proportion of their observed
#' counts, and uses `optim()` with an L-BFGS-B search to maximise the Poisson
#' log-likelihood of the held-out data, returning the hyperparameters that score
#' best.
#' @export
fit <- function(obs_data, nt, period, n_sites, mask_prop = 0.2, verbose = FALSE, par0 = c(1, 5, 100), lower = c(1e-4, 0.8, 52 * 1.5), upper = c(2, 10, 500)) {
  ids <- sample(unique(obs_data$id), n_sites)
  fitting_data   <- obs_data[obs_data$id %in% ids, ]
  fitting_coords <- coordinates[coordinates$id %in% ids, ]

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
