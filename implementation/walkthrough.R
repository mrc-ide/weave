# =============================================================================
# weave -- kernel-hyperparameter walkthrough
# =============================================================================
#
# WHAT THIS DOES
#   Estimates the separable-GP kernel hyperparameters from count data, then uses
#   them to predict the latent rate. The steps:
#
#     1. Build a plug-in latent field from the counts: the per-site centred and
#        scaled log1p(y).
#     2. Fit (length_scale, periodic_scale, long_term_scale) plus a noise/nugget
#        ratio by maximising the GP marginal likelihood of that field. The
#        Kronecker eigendecomposition gives the exact log-determinant and
#        quadratic form in O(n^3 + nt^3), so the full (n*nt)-square covariance is
#        never formed and the global variance is profiled out analytically.
#     3. Predict the latent rate with a closed-form separable-GP smoother and
#        turn it into a count prediction interval.
#
#   The script simulates data with a known truth and PRINTS (does not save):
#     Plot 1  simulated true mean + observations (held-out cells in red)
#     Plot 2  fitted vs true spatial and temporal kernels
#     Plot 3  Plot 1 with the predicted mean (line) + 95% prediction interval
#     A short numeric report (estimated vs true hyperparameters, coverage)
#
# CAVEATS / APPROXIMATIONS
#   - Plug-in, not fully Bayesian: it conditions on a single noisy estimate of
#     the latent field rather than integrating the field out. This attenuates
#     the length scales (periodic_scale / long_term_scale tend to come out low);
#     the nugget mitigates but does not remove it.
#   - log1p(y) is a crude stand-in for the latent log-rate (poor at low counts).
#   - Per-site scaling folds all per-site variance into one global sigma^2 and
#     amplifies noise at low-count sites.
#   - The nugget is a single homoscedastic noise term; real count noise is
#     heteroscedastic.
#   - Missing cells are mean-imputed, so they contribute no uncertainty.
#   - It returns a point (MAP) estimate -- no hyperparameter uncertainty.
#
# HOW TO RUN
#   From the project root:  source("implementation/walkthrough.R")
# =============================================================================

# -----------------------------------------------------------------------------
# 0. Setup
# -----------------------------------------------------------------------------
devtools::load_all(quiet = TRUE)

suppressPackageStartupMessages(library(ggplot2))

set.seed(20260610)

# --- simulation helpers (inlined so this script is fully self-contained) -----
# Previously in implementation/simulation.R (now archived). They draw a
# Negative-Binomial GP and apply a clustered missingness pattern.

# Ground-truth counts: y_st ~ NB(r, mean = exp(mu_s + f_st)), latent field
# f ~ N(0, space_k (x) time_k). One row per (site, time), site-major with time
# varying fastest (matching quick_mvnorm()).
simulate_data <- function(n, nt, coordinates, space_k, time_k, r = Inf) {
  f <- quick_mvnorm(space_k, time_k) # length n*nt, time fastest
  mu_full <- rep(coordinates$mu, each = nt) # site intercepts
  lambda <- exp(mu_full + f)
  y <- if (is.finite(r)) {
    stats::rnbinom(n * nt, size = r, mu = lambda)
  } else {
    stats::rpois(n * nt, lambda)
  }
  data.frame(
    id = factor(rep(seq_len(n), each = nt), levels = seq_len(n)),
    t = rep(seq_len(nt), times = n),
    lambda = lambda,
    y = y
  )
}

# A 0/1 sequence with clustered runs (whole stretches missing, not random cells).
generate_clustered_binary <- function(n, p_one, p_switch) {
  result <- numeric(n)
  result[1] <- stats::rbinom(1, 1, p_one)
  for (i in 2:n) {
    result[i] <- if (stats::runif(1) < p_switch) {
      stats::rbinom(1, 1, p_one)
    } else {
      result[i - 1]
    }
  }
  result
}

# Mask the truth with clustered missingness; observed counts (NA = missing).
observed_data <- function(data, p_one, p_switch) {
  miss <- generate_clustered_binary(nrow(data), p_one, p_switch)
  y_obs <- data$y
  y_obs[miss == 1] <- NA
  data.frame(id = data$id, t = data$t, y_obs = y_obs)
}


# -----------------------------------------------------------------------------
# 1. Controls
# -----------------------------------------------------------------------------
n <- 100 # number of sites (health facilities)
nt <- 52 * 5 # number of time points (3 yrs weekly)
period <- 52 # seasonal period (weeks/cycle)

true_length_scale <- 2 # spatial smoothness (distance units)
true_periodic_scale <- 1.1 # how "peaky" the season is
true_long_term_scale <- 150 # long-run trend smoothness
true_r <- 15 # NB dispersion (smaller = heavier tail)

show_missingness <- TRUE # draw held-out (missing) truth in red
p_one <- 0.1 # missingness controls (see observed_data)
p_switch <- 0.3
plot_sites <- 1:50 # sites shown in the per-site panels (default: all; set e.g. 1:12 to subset)


# -----------------------------------------------------------------------------
# 2. Simulate ground truth
# -----------------------------------------------------------------------------
# Model:   y_st ~ NB(r, mean = lambda_st = exp(mu_s + f_st)),
#          f ~ N(0, K_space(theta) (x) K_time(theta)).
# simulate_data() draws the latent GP and the counts; observed_data() then masks
# a clustered missingness pattern (whole stretches missing, like real dropouts)
# and returns the observed counts in `y_obs` (NA where missing).
# -----------------------------------------------------------------------------
coordinates <- data.frame(
  id = factor(1:n),
  lat = runif(n, 0, 5),
  lon = runif(n, 0, 5),
  mu = log(runif(n, 10, 80)) # site-specific mean count
)

space_k <- space_kernel(coordinates, length_scale = true_length_scale)
time_k <- time_kernel(
  1:nt,
  periodic_scale = true_periodic_scale,
  long_term_scale = true_long_term_scale,
  period = period
)

true_data <- simulate_data(n, nt, coordinates, space_k, time_k, r = true_r)
obs_data <- observed_data(true_data, p_one = p_one, p_switch = p_switch)

cat(sprintf(
  "Simulated %d cells across %d sites; %d observed (%.0f%% missing).\n",
  n * nt,
  n,
  sum(!is.na(obs_data$y_obs)),
  100 * mean(is.na(obs_data$y_obs))
))


# -----------------------------------------------------------------------------
# 3. Plot 1: simulated true mean + observations
# -----------------------------------------------------------------------------
# Black line  = true latent mean lambda (what we want to recover).
# Black points = observed counts.
# Red points   = held-out truth at the MISSING cells (only if show_missingness),
#                so you can see what the model has to reconstruct.
# -----------------------------------------------------------------------------
hf_labeller <- function(value) paste("Site", value)

# tidy frames keyed by integer site id + time, restricted to the plotted sites
truth_df <- data.frame(
  id = as.integer(true_data$id),
  t = true_data$t,
  lambda = true_data$lambda,
  y = true_data$y
)
obs_df <- data.frame(
  id = as.integer(obs_data$id),
  t = obs_data$t,
  y_obs = obs_data$y_obs
)
missing_df <- merge(truth_df, obs_df, by = c("id", "t"))
missing_df <- missing_df[is.na(missing_df$y_obs), ] # held-out truth

sub <- function(d) d[d$id %in% plot_sites, ]

base_plot <- ggplot() +
geom_line(
  data = sub(truth_df),
  aes(t, lambda),
  colour = "black",
  linewidth = 0.4
) +
geom_point(data = sub(obs_df), aes(t, y_obs), size = 0.5, colour = "grey20") +
facet_wrap(~id, scales = "free_y", labeller = labeller(id = hf_labeller)) +
labs(
  x = "Week",
  y = "Cases",
  title = "Simulated truth: black line = true mean, points = observed counts"
) +
theme_bw() +
theme(
  strip.background = element_rect(fill = "white", colour = "grey60"),
  strip.text = element_text(size = 7, face = "bold")
)

print(base_plot)


# -----------------------------------------------------------------------------
# 4. Fit the kernel hyperparameters
# -----------------------------------------------------------------------------
# Maximise the GP marginal likelihood of the plug-in field to estimate the three
# length scales and the noise/nugget ratio (profiled global variance).
# -----------------------------------------------------------------------------
fit_time <- system.time(
  est <- infer_kernel_params(obs_data, coordinates, nt = nt, period = period)
)

cat(sprintf("\nFit runtime: %.2f s\n", fit_time[3]))
report <- data.frame(
  parameter = c("length_scale", "periodic_scale", "long_term_scale"),
  truth = c(true_length_scale, true_periodic_scale, true_long_term_scale),
  estimate = round(
    c(est$length_scale, est$periodic_scale, est$long_term_scale),
    3
  )
)
print(report, row.names = FALSE)
cat(sprintf(
  "nugget_ratio (noise/signal) = %.3f   profiled sigma^2 = %.3f\n",
  est$nugget_ratio,
  est$sigma2
))


# -----------------------------------------------------------------------------
# 5. Plot 2: the fitted kernels in space and time
# -----------------------------------------------------------------------------
# The shapes implied by the estimated hyperparameters. The temporal kernel is
# the periodic component modulated by the long-term RBF decay.
# -----------------------------------------------------------------------------
# Build the spatial and temporal correlation curves for a given parameter set.
kernel_curves <- function(
  length_scale,
  periodic_scale,
  long_term_scale,
  label
) {
  sp <- data.frame(distance = seq(0, 5, length.out = 200), type = label)
  sp$correlation <- rbf_kernel(sp$distance, theta = length_scale)
  tm <- data.frame(lag = seq(0, nt - 1, length.out = 600), type = label)
  tm$correlation <- periodic_kernel(
    tm$lag,
    alpha = periodic_scale,
    period = period
  ) *
  rbf_kernel(tm$lag, theta = long_term_scale)
  list(space = sp, time = tm)
}

fit_curves <- kernel_curves(
  est$length_scale,
  est$periodic_scale,
  est$long_term_scale,
  "Estimated"
)
true_curves <- kernel_curves(
  true_length_scale,
  true_periodic_scale,
  true_long_term_scale,
  "True"
)
space_curve <- rbind(fit_curves$space, true_curves$space)
time_curve <- rbind(fit_curves$time, true_curves$time)

kernel_cols <- c(Estimated = "deeppink", True = "black")
kernel_ltys <- c(Estimated = "solid", True = "dashed")

space_kernel_plot <- ggplot(
  space_curve,
  aes(distance, correlation, colour = type, linetype = type)
) +
geom_line(linewidth = 1) +
scale_colour_manual(values = kernel_cols, name = NULL) +
scale_linetype_manual(values = kernel_ltys, name = NULL) +
labs(
  x = "Spatial distance",
  y = "Correlation",
  title = "Spatial kernel: estimated (pink) vs true (black dashed)"
) +
ylim(0, 1) +
theme_bw()

time_kernel_plot <- ggplot(
  time_curve,
  aes(lag, correlation, colour = type, linetype = type)
) +
geom_line(linewidth = 1) +
scale_colour_manual(values = kernel_cols, name = NULL) +
scale_linetype_manual(values = kernel_ltys, name = NULL) +
labs(
  x = "Temporal lag (weeks)",
  y = "Correlation",
  title = "Temporal kernel: estimated (pink) vs true (black dashed)"
) +
theme_bw()

print(space_kernel_plot)
print(time_kernel_plot)


# -----------------------------------------------------------------------------
# 6. Predict the latent rate and a count prediction interval
# -----------------------------------------------------------------------------
# Given the fitted hyperparameters, denoise the plug-in field with a closed-form
# separable-GP smoother. The noise is a scalar multiple of the identity and the
# grid is completed (missing cells mean-imputed), so everything diagonalises in
# the Kronecker eigenbasis:
#
#   posterior mean of mode (i,j):  S_ij * ghat_ij,   S_ij = lambda_ij/(lambda_ij+eta)
#   posterior var  of mode (i,j):  sigma^2 * lambda_ij * eta / (lambda_ij + eta)
#
# with lambda_ij = a_i b_j the Kronecker eigenvalues and ghat = U_s' G U_t. We
# then undo the per-site standardisation and return the posterior mean and
# variance of the log-rate.
#
# The ribbon is a PREDICTION INTERVAL for counts (not a credible interval on the
# mean): we fold observation noise into the rate posterior via the law of total
# variance and moment-match a lognormal to read off 2.5/97.5%. The NB dispersion
# r is estimated by method of moments from the observed counts (override via
# r_pred below).
#
# NB: mean-imputing missing cells with homoscedastic noise understates their
# uncertainty -- fine as a diagnostic, not a proper missing-data posterior.
# -----------------------------------------------------------------------------
gp_smoother <- function(
  obs_data,
  coordinates,
  est,
  n,
  nt,
  period,
  value = "y_obs"
) {
  ids <- sort(unique(obs_data$id))
  times <- sort(unique(obs_data$t))
  coordinates <- coordinates[match(ids, coordinates$id), , drop = FALSE]
  
  # Plug-in field + per-site centring/scaling (kept so we can undo it).
  M <- matrix(NA_real_, n, nt)
  M[cbind(
    match(obs_data$id, ids),
    match(obs_data$t, times)
  )] <- log1p(obs_data[[value]])
  row_mean <- rowMeans(M, na.rm = TRUE)
  row_mean[!is.finite(row_mean)] <- 0
  Mc <- M - row_mean
  row_sd <- apply(Mc, 1, stats::sd, na.rm = TRUE)
  row_sd[!is.finite(row_sd) | row_sd == 0] <- 1
  G <- Mc / row_sd
  G[is.na(G)] <- 0
  
  # Eigendecompositions of the fitted correlation kernels.
  eig_s <- eig_sym(space_kernel(coordinates, length_scale = est$length_scale))
  eig_t <- eig_sym(time_kernel(
    times,
    periodic_scale = est$periodic_scale,
    long_term_scale = est$long_term_scale,
    period = period
  ))
  eta <- est$nugget_ratio
  s2 <- est$sigma2
  
  lam <- outer(eig_s$values, eig_t$values) # n x nt Kronecker eigenvalues
  shrink <- lam / (lam + eta)
  pv <- s2 * lam * eta / (lam + eta) # posterior var per mode
  
  ghat <- crossprod(eig_s$vectors, G) %*% eig_t$vectors # U_s' G U_t
  Ghat <- eig_s$vectors %*% (shrink * ghat) %*% t(eig_t$vectors) # smoothed (std)
  Vstd <- (eig_s$vectors^2) %*% pv %*% t(eig_t$vectors^2) # per-cell var (std)
  
  # Undo standardisation -> log-rate scale (mu_s + f_st). Return the posterior
  # mean and variance of the log-rate so the caller can build a prediction
  # interval that also folds in observation noise.
  Z <- row_mean + row_sd * Ghat
  Vz <- (row_sd^2) * Vstd
  list(
    ids = ids,
    times = times,
    Z = Z, # posterior mean of the log-rate
    Vz = Vz, # posterior variance of the log-rate
    mean = exp(Z + Vz / 2) # posterior mean of the rate lambda
  )
}

pred <- gp_smoother(obs_data, coordinates, est, n, nt, period)

# --- estimate the NB dispersion r (method of moments on observed cells) -------
# Var(y | lambda) = lambda + lambda^2 / r, so r ~ sum(lambda^2) / sum((y-lambda)^2 - lambda).
# Set r_pred manually to override (e.g. r_pred <- true_r to use the known value).
ids_o <- match(obs_data$id, pred$ids)
tim_o <- match(obs_data$t, pred$times)
lam_o <- pred$mean[cbind(ids_o, tim_o)]
keep_o <- !is.na(obs_data$y_obs)
y_o <- obs_data$y_obs[keep_o]
lam_o <- lam_o[keep_o]
mom_den <- sum((y_o - lam_o)^2 - lam_o)
r_pred <- if (mom_den > 0) max(sum(lam_o^2) / mom_den, 0.1) else 1e6 # large r -> ~Poisson
cat(sprintf(
  "Estimated NB dispersion r (method of moments) = %.1f  (true = %.1f)\n",
  r_pred,
  true_r
))

# --- prediction interval for counts: integrate NB noise over the rate posterior
# Law of total variance with a lognormal rate posterior, then moment-match a
# lognormal to the count predictive to read off smooth 2.5/97.5% bounds.
Elam <- exp(pred$Z + pred$Vz / 2) # E[lambda]
Elam2 <- exp(2 * pred$Z + 2 * pred$Vz) # E[lambda^2]
Vlam <- exp(2 * pred$Z + pred$Vz) * (exp(pred$Vz) - 1) # Var[lambda]
pred_mean <- Elam # E[y] = E[lambda]
pred_var <- Elam + Elam2 / r_pred + Vlam # E[lambda + lambda^2/r] + Var[lambda]

mean_safe <- pmax(pred_mean, 1e-8)
ss <- log(1 + pred_var / mean_safe^2)
mln <- log(mean_safe) - ss / 2

# Long data frame keyed by integer site id + time.
pred_df <- data.frame(
  id = rep(as.integer(pred$ids), times = nt),
  t = rep(pred$times, each = n),
  mean = as.vector(pred_mean),
  lower = as.vector(stats::qlnorm(0.025, mln, sqrt(ss))),
  upper = as.vector(stats::qlnorm(0.975, mln, sqrt(ss)))
)


# -----------------------------------------------------------------------------
# 7. Plot 3: predictions on top of the truth
# -----------------------------------------------------------------------------
# Blue line  = predicted mean rate; blue ribbon = 95% credible interval for the
# latent rate. Compare against the black true-mean line and the points.
# -----------------------------------------------------------------------------
prediction_plot <- base_plot +
geom_ribbon(
  data = sub(pred_df),
  aes(t, ymin = lower, ymax = upper),
  fill = "steelblue",
  alpha = 0.25
) +
geom_line(
  data = sub(pred_df),
  aes(t, mean),
  colour = "steelblue",
  linewidth = 0.6
) +
labs(
  title = "Prediction vs truth: blue = predicted mean + 95% prediction interval (counts)",
  subtitle = if (show_missingness) {
    "black line = true mean; red = held-out truth"
  } else {
    "black line = true mean"
  }
)

print(prediction_plot)


# -----------------------------------------------------------------------------
# 8. Numeric diagnostics
# -----------------------------------------------------------------------------
# Two complementary checks at the held-out (missing) cells the model never saw:
#
#   * correlation / RMSE of the predicted mean rate vs the TRUE rate
#       -> point-prediction quality (expect the mean to track the truth well)
#   * coverage of the held-out COUNTS by the 95% prediction interval
#       -> calibration of the interval (target ~0.95)
#
# The prediction interval folds in observation noise, so it is compared against
# the true held-out COUNTS (not the latent rate). It can still miss nominal if
# the plug-in field is over/under-confident or the dispersion estimate r is off.
# -----------------------------------------------------------------------------
chk <- merge(
  missing_df[, c("id", "t", "lambda", "y")],
  pred_df,
  by = c("id", "t")
)
coverage <- mean(chk$y >= chk$lower & chk$y <= chk$upper) # held-out COUNTS in PI
rmse_held <- sqrt(mean((chk$mean - chk$lambda)^2)) # mean rate vs true rate
corr_held <- stats::cor(chk$mean, chk$lambda)
cat(sprintf("\nHeld-out cells: %d\n", nrow(chk)))
cat(sprintf(
  "corr(predicted mean rate, true rate) at held-out cells: %.3f\n",
  corr_held
))
cat(sprintf(
  "RMSE of predicted vs true rate at held-out cells:       %.2f\n",
  rmse_held
))
cat(sprintf(
  "95%% prediction-interval coverage of held-out COUNTS:    %.2f  (target ~0.95)\n",
  coverage
))

)