#' Posterior mean of the intensity surface
#'
#' Computes the posterior mean of \eqn{\lambda = \exp(z)} where
#' \eqn{z = f + \mu_{\mathrm{infer}}}. The latent field \eqn{f} follows a
#' zero–mean GP with separable covariance \eqn{K = K_{\mathrm{space}} \otimes
#' K_{\mathrm{time}}}. A Gaussian working likelihood on the log scale with
#' heteroscedastic diagonal variance \eqn{D = \mathrm{diag}(\text{noise\_var})}
#' is used. This evaluates
#' \deqn{f_{\text{hat}} = K S^\top (S K S^\top + D)^{-1} y_{\text{work}},}
#' then returns \eqn{\exp(f_{\text{hat}} + \mu_{\mathrm{infer}})}.
#'
#' @param state A sampler state created by \code{gp_build_state()}, containing
#'   at least \code{space_mat}, \code{time_mat}, \code{obs_idx}, \code{N},
#'   \code{y_work}, \code{noise_var}, \code{kdiag_full}, \code{A_solve},
#'   and \code{mu_infer}. The vector layout is sites × times with time varying
#'   fastest.
#' @param tol Convergence tolerance passed to the inner PCG solve.
#'
#' @return A numeric vector of length \code{state$N} giving the posterior mean
#'   intensity \eqn{\lambda} in the same ordering as \code{state$mu_infer}.
#'
#' @export
gp_posterior_mean <- function(state, tol = 1e-6) {
  alpha <- state$A_solve(state$y_work, tol = tol)             # A^{-1} y
  f_hat <- kron_mv(
    fill_vector(alpha, state$obs_idx, state$N),
    state$space_mat, state$time_mat
  )
  z_hat <- f_hat + state$mu_infer
  exp(z_hat)  # lambda_hat
}

#' One posterior draw of the intensity surface
#'
#' Generates a single posterior sample of \eqn{\lambda = \exp(z)} under the
#' Gaussian working model on the log scale with heteroscedastic variance
#' \eqn{D}. It draws a GP prior realisation \eqn{\eta \sim \mathcal N(0, K)},
#' adds working-likelihood noise \eqn{\varepsilon \sim \mathcal N(0, D)} on the
#' observed log scale, solves
#' \eqn{(S K S^\top + D)\alpha_{\sim} = (S\eta + \varepsilon) - y_{\text{work}}},
#' and returns \eqn{\exp(\eta - K S^\top \alpha_{\sim} + \mu_{\mathrm{infer}})}.
#'
#' @param state A sampler state created by \code{gp_build_state()}, containing
#'   at least \code{space_mat}, \code{time_mat}, \code{obs_idx}, \code{N},
#'   \code{y_work}, \code{noise_var}, \code{kdiag_full}, \code{A_solve},
#'   and \code{mu_infer}. The vector layout is sites × times with time varying
#'   fastest.
#' @param tol Convergence tolerance passed to the inner PCG solve.
#'
#' @return A numeric vector of length \code{state$N} containing one posterior
#'   draw of \eqn{\lambda} in the same ordering as \code{state$mu_infer}.
#' @export
gp_draw <- function(state, tol = 1e-6) {
  m   <- length(state$obs_idx)

  # Prior draw using precomputed Choleskys
  eta <- quick_mvnorm_chol(state$Rs, state$Rt)

  # Pseudo-noise on observed log scale (heteroscedastic)
  eps <- stats::rnorm(m, sd = sqrt(state$noise_var))

  # Solve A alpha_tilde = (S eta + eps) - y_work
  rhs <- (eta[state$obs_idx] + eps) - state$y_work
  alpha_tilde <- state$A_solve(rhs, tol = tol)

  # Posterior draw for f (centred log scale): eta - K S^T alpha_tilde
  f_draw <- eta - kron_mv(
    fill_vector(alpha_tilde, state$obs_idx, state$N),
    state$space_mat, state$time_mat
  )

  z_draw <- f_draw + state$mu_infer
  exp(z_draw)
}

#' Compute quantile bounds for GP-based count draws
#'
#' This function generates draws from a Gaussian process state, applies
#' Poisson observation noise, and computes quantile bounds of the resulting
#' counts for each observation.
#'
#' @param state A fitted GP state object used as input to \code{gp_draw()}.
#' @param n_lambda Integer. Number of latent Gaussian process draws.
#' @param n_draw Integer. Number of replicate draws with Poisson noise applied.
#' @param quantiles Numeric vector. Quantiles to compute for the count draws
#'   (values between 0 and 1).
#'
#' @return A data frame with one row per observation. Contains the columns:
#'   \itemize{
#'     \item \code{id} Observation identifier.
#'     \item \code{t} Observation time index.
#'     \item Additional columns for each quantile requested, named
#'       \code{q{quantile}} (e.g. \code{q0.025}, \code{q0.5}).
#'   }
#'
#' @seealso \code{\link{gp_draw}}, \code{\link{quantile}}
#'
#' @export
bounds <- function(state, n_lambda = 30, n_draw = 100,
                   quantiles = c(0.025, 0.25, 0.75, 0.975)) {
  # Lambda draws
  #lam_draws <- replicate(n_lambda, gp_draw(state), simplify = "array")

  lam_draws <- pbapply::pbsapply(
    seq_len(n_lambda),
    function(i) gp_draw(state)
  ) |>
    as.matrix()

  # Add observation noise
  count_draws <- lapply(seq_len(n_draw), function(x) {
    tmp <- lam_draws
    tmp[] <- stats::rpois(length(lam_draws), lambda = lam_draws)
    tmp
  })
  count_draws <- do.call("cbind", count_draws)

  # Apply quantiles: result is (n_obs × n_quantiles)
  qs <- t(apply(
    count_draws,
    1,
    function(x) as.vector(stats::quantile(x, probs = quantiles))
  ))

  qs <- matrix(qs, ncol = length(quantiles))

  # Build output
  out <- data.frame(
    id = state$id,
    t = state$t,
    qs
  )

  # Optional: give nice column names
  colnames(out)[-(1:2)] <- paste0("q", quantiles)

  return(out)
}
