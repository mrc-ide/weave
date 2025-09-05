#' Add a small ridge to a square matrix
#'
#' In plain terms: this adds a tiny value to the diagonal so the matrix is
#' better-behaved numerically (e.g., invertible and Cholesky-able).
#'
#' Technically: returns \eqn{X + \lambda I}, which improves condition number and
#' ensures positive definiteness when \eqn{\lambda > 0}.
#'
#' @param x A square numeric matrix.
#' @param lambda Non-negative ridge (diagonal) value to add. Default `1e-5`.
#'
#' @return A matrix the same size as `x` with `lambda` added to the diagonal.
regularise <- function(x, lambda = 1e-5) {
  x + lambda * diag(nrow(x))
}

#' Hessian of the (log posterior) for Poisson log-Gaussian model
#'
#' Builds the curvature matrix used for uncertainty—combining
#' the GP prior and the Poisson likelihood—at a given latent vector `f`.
#'
#' Technically: returns the Hessian
#' \deqn{H = -\Sigma^{-1} - \mathrm{diag}(\exp(f))}
#' where \eqn{\Sigma^{-1} = \text{dist\_k\_inv} \otimes \text{time\_k\_inv}}.
#' This is the Hessian of the log posterior (negative definite).
#'
#' @param f Numeric vector of latent values (length \eqn{n \times nt}).
#' @param dist_k_inv Inverse spatial kernel (precision) matrix.
#' @param time_k_inv Inverse temporal kernel (precision) matrix.
#'
#' @return A dense square matrix \eqn{H} of size `length(f)` × `length(f)`.
llh <- function(f, dist_k_inv, time_k_inv) {
  # Hessian of the quadratic term
  # The Hessian is -Sigma_inv, represented using the Kronecker product
  Hess_quad <- -kronecker(dist_k_inv, time_k_inv)

  # Initialize the Hessian of the likelihood term
  Hess_likelihood <- matrix(0, nrow = length(f), ncol = length(f))

  # Poisson distribution
  mu <- exp(f)  # Mean parameter
  tmp <- -mu  # Second derivative of Poisson log-likelihood
  Hess_likelihood <- diag(tmp)

  # Total Hessian
  Hess <- Hess_quad + Hess_likelihood

  # Return the Hessian matrix
  return(Hess)
}

#' Fast Kronecker–product matrix–vector multiply (times vary fastest)
#'
#' In plain terms: multiplies a big covariance `K = space ⊗ time` by a vector
#' without ever forming `K`, using a reshape–multiply–reshape trick.
#'
#' Technically: for \eqn{v = \mathrm{vec}(X^\top)} with `times` varying fastest,
#' computes \eqn{(space \otimes time)\,v = \mathrm{vec}\!\big((space\,X\,time^\top)^\top\big)}.
#'
#' @param v Numeric vector of length `nrow(space) * nrow(time)`, ordered with
#'   times varying fastest within site.
#' @param space Spatial kernel matrix (size \eqn{n \times n}).
#' @param time Temporal kernel matrix (size \eqn{nt \times nt}).
#'
#' @return A numeric vector the same length as `v`.
kron_mv <- function(v, space, time) {
  n_sites <- nrow(space); n_times <- nrow(time)
  # Reconstruct X (sites × times) from vec(t(X)) = v
  X <- t(matrix(v, nrow = n_times, ncol = n_sites))
  # Apply K to X: vec(t(space %*% X %*% t(time))) = (space ⊗ time) vec(t(X))
  Y <- space %*% X %*% t(time)
  as.vector(t(Y))
}

#' Fill observed values into a full vector
#'
#' Put the observed entries back into their full-length vector,
#' filling missing positions with zeros.
#'
#' @param x_obs Numeric vector of observed values (length \eqn{m}).
#' @param obs_idx Integer indices (length \eqn{m}) of observed positions in the
#'   full vector.
#' @param N Total length of the full vector.
#'
#' @return A numeric vector of length `N` with `x_obs` scattered at `obs_idx`.
fill_vector <- function(x_obs, obs_idx, N){
  v <- numeric(N)
  v[obs_idx] <- x_obs
  v
}

#' Observed-system matvec: (S K S^T + diag(noise)) v
#'
#' Takes an observed-length vector and applies the GP covariance
#' plus a per-observation noise (nugget), all without building any big matrices.
#'
#' Technically: for \eqn{K = \mathrm{space}\,\otimes\,\mathrm{time}}, returns
#' \deqn{S\,K\,S^{\mathsf T}\,v \;+\; \operatorname{diag}(\sigma^2)\,v,}
#' i.e., the observed block of the GP plus a diagonal nugget. Implemented
#' matrix-free as \code{kron_mv(with_nas(v, obs_idx, N), space_mat, time_mat)[obs_idx] + noise_var * v},
#' where \eqn{S^{\mathsf T}} “scatters’’ into the full vector and \eqn{\sigma^2}
#' denotes the per-observation noise.
#'
#' @param v Numeric vector of length \eqn{m} (observed entries).
#' @param obs_idx Integer indices of observed entries in the full vector.
#' @param N Total length of the full vector.
#' @param space_mat Spatial kernel matrix.
#' @param time_mat Temporal kernel matrix.
#' @param noise_var Scalar or length-\eqn{m} numeric nugget on the observed scale.
#'
#' @return A numeric vector of length \eqn{m}, equal to \eqn{(S K S^\top + D)v}.
Amv <- function(v, obs_idx, N, space_mat, time_mat, noise_var) {
  #browser()
  kron_mv(fill_vector(v, obs_idx, N), space_mat, time_mat)[obs_idx] + noise_var * v
}

#' Diagonal (Jacobi) preconditioner application
#'
#' Divides by an approximation to the diagonal of the system,
#' which makes the iterative solver converge faster.
#'
#' Technically: applies \eqn{M^{-1} v \approx v / \mathrm{diag}(A)}, where
#' \eqn{A = S K S^\top + \mathrm{diag}(\text{noise})} and
#' \eqn{\mathrm{diag}(K) = \mathrm{diag}(space) \otimes \mathrm{diag}(time)}.
#'
#' @param v Numeric vector to precondition (length \eqn{m}).
#' @param kdiag_full Vector \eqn{\mathrm{diag}(K)} of length \eqn{N}
#'   (typically from `as.vector(kronecker(diag(space), diag(time)))`).
#' @param obs_idx Integer indices of observed entries in the full vector.
#' @param noise_var Scalar or length-\eqn{m} numeric nugget to add to the diagonal.
#'
#' @return A numeric vector of length \eqn{m}: elementwise `v / (diagA + 1e-12)`.
M_inv <- function(v, kdiag_full, obs_idx, noise_var) {
  v / (kdiag_full[obs_idx] + noise_var + 1e-12)
}

#' Preconditioned Conjugate Gradient (PCG) solver for the observed system
#'
#' In plain terms: solves the big linear system that gives the GP weights using
#' only matrix–vector products—no huge matrices, no explicit inverse.
#'
#' Technically: solves \eqn{(S K S^\top + \mathrm{diag}(\text{noise}))\,x = b}
#' by PCG, using `Amv` for matrix–vector products and `M_inv` as a Jacobi
#' preconditioner. Stops when the relative residual falls below `tol` or after
#' `maxit` iterations (issues a warning on `maxit`).
#'
#' @param b Right-hand side vector (observed length \eqn{m}).
#' @param obs_idx Integer indices of observed entries in the full vector.
#' @param N Total length of the full vector.
#' @param space_mat Spatial kernel matrix.
#' @param time_mat Temporal kernel matrix.
#' @param noise_var Scalar or length-\eqn{m} nugget on the observed scale.
#' @param kdiag_full Vector \eqn{\mathrm{diag}(K)} of length \eqn{N}.
#' @param tol Relative residual tolerance for convergence (default `1e-8`).
#' @param maxit Maximum number of iterations (default `10000`).
#'
#' @return Numeric solution vector `x` of length \eqn{m}.
pcg <- function(b, obs_idx, N, space_mat, time_mat, noise_var, kdiag_full, tol = 1e-8, maxit = 10000) {
  x <- numeric(length(b))
  r <- b - Amv(x, obs_idx, N, space_mat, time_mat, noise_var)
  z <- M_inv(r, kdiag_full, obs_idx, noise_var)
  p <- z
  rz_old <- sum(r * z)
  for (it in seq_len(maxit)) {
    if(it == maxit){
      warning("maxit reached")
    }
    Ap <- Amv(p, obs_idx, N, space_mat, time_mat, noise_var)
    alpha <- rz_old / sum(p * Ap)
    x <- x + alpha * p
    r <- r - alpha * Ap
    if (sqrt(sum(r * r)) <= tol * sqrt(sum(b * b))) break
    z <- M_inv(r, kdiag_full, obs_idx, noise_var)
    rz_new <- sum(r * z)
    beta <- rz_new / rz_old
    p <- z + beta * p
    rz_old <- rz_new
  }
  x
}

#' Fit a spatiotemporal GP to log–counts (matrix-free PCG)
#'
#' Fits a separable spatiotemporal Gaussian process to the working
#' log–count scale and returns posterior means/uncertainties plus simple
#' predictive summaries. The solve uses a **matrix-free** preconditioned
#' conjugate gradient (PCG) on the observed system with a Kronecker
#' covariance \eqn{K = K_{\mathrm{space}} \otimes K_{\mathrm{time}}}.
#'
#' @param obs_data Tibble/data.frame with one row per site–time. Must contain
#'   at least columns `id`, `t`, `y_obs` (counts, may be `NA`),
#'   `mu_infer` (site mean on log scale), and `f_infer`
#'   (working response, e.g. `log(y_obs+1) - mu_infer`).
#' @param coordinates Tibble/data.frame with site metadata used by the spatial
#'   kernel (e.g. `id`, `lat`, `lon`).
#' @param hyperparameters Numeric length-3 vector:
#'   `c(space_length_scale, time_periodic_scale, time_long_term_scale)`.
#' @param n Integer, number of sites.
#' @param nt Integer, number of time points.
#' @return A tibble with the original rows plus columns:
#'   `z_est` (posterior mean on log scale), `tausq` (approx. variance),
#'   `z_min`, `z_max` (95% CI), lognormal summaries on the count scale
#'   `lambda_est`, `lambda_min`, `lambda_max`, and simple Negative-Binomial
#'   prediction summaries `negbin_size`, `negbin_prob`, `pred_Q2.5`,
#'   `pred_Q25`, `data_Q50`, `pred_Q75`, `pred_Q97.5`, plus `log_p`,
#'   `y_mode`, and `surprisal`.
#' @export
fit <- function(obs_data, coordinates, hyperparameters, n, nt){
  # Build kernels
  time_mat  <- time_kernel(
    times = 1:nt,
    periodic_scale = hyperparameters[2],
    long_term_scale = hyperparameters[3],
    period = 52
  )
  space_mat <- space_kernel(
    coordinates = coordinates, length_scale = hyperparameters[1]
  )

  obs_idx <- which(!is.na(obs_data$y_obs))  # ensure ordering matches vec stacking
  N <- n * nt

  # Centered observed response (log with eps, minus mu_infer)
  y_obs <- obs_data$f_infer[obs_idx]

  # Heteroscedastic nugget on the log scale:
  ## y_obs = log(Y+1) − mu_infer is noisy (Poisson + log). By the delta method,
  ## Var[log(Y+1) | λ] ≈ λ/(λ+1)^2. Using λ̂_i = exp(mu_infer_i) gives
  ## noise_var_i = λ̂_i/(λ̂_i+1)^2, which stabilises (S K S^T + diag(noise_var))
  ## and avoids overfitting log-count noise.
  lam_hat <- exp(obs_data$mu_infer[obs_idx])
  noise_var <- lam_hat / (lam_hat + 1)^2

  # Jacobi (diagonal) preconditioner:
  ## For A = S K S^T [+ diag(noise)], approximate A^{-1} with 1/diag(A).
  ## Since diag(K) = diag(space) ⊗ diag(time), build it once, then select observed entries.
  kdiag_full <- as.vector(kronecker(diag(space_mat), diag(time_mat)))

  # Solve for α using PCG:
  ## α are the GP conditioning weights solving (S K S^T + diag(noise_var)) α = y_obs.
  ## They map observations to the posterior mean via f_hat = K S^T α.
  ## PCG (Preconditioned Conjugate Gradient) is an iterative solver that needs only
  ## matrix–vector products with A (no full matrix), using M_inv as a diagonal
  ## preconditioner to speed convergence.
  alpha <- pcg(y_obs, obs_idx, N, space_mat, time_mat, noise_var, kdiag_full, tol = 1e-6)

  # Posterior mean at all entries: f_hat = K S^T alpha  (then add mu)
  f_hat <- kron_mv(fill_vector(alpha, obs_idx, N), space_mat, time_mat)

  # Write back z_est per row
  obs_data$z_est <- f_hat + obs_data$mu_infer

  # Variance (Hessian diag at ML): use your efficient diag trick, global form
  space_inv <- chol2inv(chol(regularise(space_mat)))
  time_inv  <- chol2inv(chol(regularise(time_mat)))
  kron_diag_inv <- as.vector(kronecker(diag(space_inv), diag(time_inv)))
  hess_diag <- -kron_diag_inv - exp(obs_data$z_est)  # evaluated at f_hat
  obs_data$tausq <- pmax(0, -1 / hess_diag)

  fit_data <- obs_data |>
    dplyr::mutate(
      # Overdisperson estimates, by site
      mean_lambda = mean(exp(.data$z_est)),
      var_y = stats::var(y_obs, na.rm = TRUE),
      overdispersion = (.data$mean_lambda^2) / (.data$var_y - .data$mean_lambda),
      .by = .data$id
    ) |>
    dplyr::mutate(
      # Compute the 95% confidence interval for z
      z_min = .data$z_est - 1.96 * sqrt(.data$tausq),
      z_max = .data$z_est + 1.96 * sqrt(.data$tausq),

      # Convert this uncertainty to the count scale. We cannot simply exponentiate
      # z_min and z_max because exponentiation is nonlinear, and normal confidence
      # intervals don’t transform correctly.
      # Instead, we compute quantiles of the corresponding lognormal distribution.
      lambda_est = stats::qlnorm(0.5, meanlog = .data$z_est, sdlog = sqrt(.data$tausq)),
      lambda_min = stats::qlnorm(0.025, meanlog = .data$z_est, sdlog = sqrt(.data$tausq)),
      lambda_max = stats::qlnorm(0.975, meanlog = .data$z_est, sdlog = sqrt(.data$tausq)),

      # Prediction intervals assuming Negative Binomially distributed data.
      # The Negative Binomial accounts for overdispersion beyond Poisson variability.
      negbin_size = .data$lambda_est / .data$overdispersion,  # Size parameter for NB
      negbin_prob = .data$negbin_size / (.data$negbin_size + .data$lambda_est), # NB probability

      # Compute quantiles of the Negative Binomial distribution as prediction intervals.
      pred_Q2.5 = stats::qnbinom(0.025, size = .data$negbin_size, prob = .data$negbin_prob),
      pred_Q25 = stats::qnbinom(0.25, size = .data$negbin_size, prob = .data$negbin_prob),
      data_Q50 = stats::qnbinom(0.5, size = .data$negbin_size, prob = .data$negbin_prob),
      pred_Q75 = stats::qnbinom(0.75, size = .data$negbin_size, prob = .data$negbin_prob),
      pred_Q97.5 = stats::qnbinom(0.975, size = .data$negbin_size, prob = .data$negbin_prob)
    ) |>
    dplyr::mutate(
      # Surprisal (information theory)
      log_p = stats::dnbinom(.data$y_obs, size = .data$negbin_size, prob = .data$negbin_prob, log = TRUE),
      y_mode = floor((.data$negbin_size - 1) * (1 - .data$negbin_prob) / .data$negbin_prob),
      log_p_mode = stats::dnbinom(.data$y_mode,  size = .data$negbin_size, prob = .data$negbin_prob, log = TRUE),
      surprisal = 1 - exp(.data$log_p - .data$log_p_mode) # 0: most expected, 1: highly surprising
    )

  return(fit_data)
}
