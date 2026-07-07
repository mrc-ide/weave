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
#     3. Predict the latent rate with gp_predict(): a matrix-free CG posterior
#        that conditions on the observed cells, then turn it into a count
#        prediction interval.
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
#   - Missing cells: with refine = TRUE (used below) hyperparameter estimation
#     fills the gaps with the GP conditional mean and refits, so missingness no
#     longer biases the length scales; the default refine = FALSE mean-imputes
#     the gaps (faster, but attenuated). Prediction (gp_predict) always
#     conditions on the observed cells, so its interval can widen over gaps.
#   - It returns a point (MAP) estimate of the hyperparameters -- their own
#     uncertainty is not propagated into the prediction.
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
nt <- 52 * 5 # number of time points (5 yrs weekly)
period <- 52 # seasonal period (weeks/cycle)

true_length_scale <- 0.5 # spatial smoothness (distance units)
true_periodic_scale <- 2 # how "peaky" the season is
true_long_term_scale <- 150 # long-run trend smoothness
true_r <- 15 # NB dispersion (smaller = heavier tail)

show_missingness <- TRUE # draw held-out (missing) truth in red
p_one <- 0.1 # missingness controls (see observed_data)
p_switch <- 0.05
plot_sites <- 1:min(n, 100) # sites shown in the per-site panels (default: all; set e.g. 1:12 to subset)

n_workers <- 1 # >1 parallelises gp_predict's posterior draws across background
# R sessions (results are identical for any value, and reproducible under the
# seed above). NOTE: workers load the INSTALLED weave, not this session's
# load_all() copy -- run devtools::install() first when setting n_workers > 1.


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

# Held-out truth at the missing cells, in red (carried into Plot 3 too, since
# the prediction plot is built on base_plot).
if (show_missingness) {
  base_plot <- base_plot +
    geom_point(
      data = sub(missing_df),
      aes(t, y),
      size = 0.6,
      colour = "red"
    ) +
    labs(subtitle = "red = held-out truth at missing cells")
}

print(base_plot)


# -----------------------------------------------------------------------------
# 4. Fit the kernel hyperparameters
# -----------------------------------------------------------------------------
# Maximise the GP marginal likelihood of the plug-in field to estimate the three
# length scales and the noise/nugget ratio (profiled global variance).
#
# refine = TRUE runs an EM-style refinement: it refits after filling the missing
# weeks with the GP conditional mean (instead of the flat per-site mean), which
# removes the downward bias the gaps would otherwise put on the length scales. It
# reuses the same fast solve as gp_predict(), so it stays cheap.
# -----------------------------------------------------------------------------
fit_time <- system.time(
  est <- infer_kernel_params(
    obs_data,
    coordinates,
    nt = nt,
    period = period,
    n_sites = 50,
    refine = TRUE
  )
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
print("")

# -----------------------------------------------------------------------------
# 6. Predict the latent rate and a count prediction interval
# -----------------------------------------------------------------------------
# gp_predict() conditions a separable GP on the OBSERVED cells only (it does not
# mean-impute the gaps), reusing the matrix-free CG machinery:
#
#   * posterior MEAN of the field -- a single CG solve, so the rate line is
#     smooth and independent of the number of draws.
#   * posterior VARIANCE          -- an exact closed-form "no gaps" part (via
#     the Kronecker eigendecomposition) plus a missing-data correction
#     estimated from `n_draws` paired perturbation draws (one CG solve each;
#     this is the expensive part). Each draw is paired with an exact
#     complete-grid twin sharing its random numbers (a control variate), so
#     most of the Monte-Carlo noise cancels and modest n_draws give tight
#     intervals -- cells far from any gap are essentially exact.
#
# The latent-rate posterior is then combined with Negative-Binomial observation
# noise (law of total variance + a lognormal moment-match) to give a 95% count
# prediction interval; the dispersion r is estimated from the data unless given.
#
# Because it conditions on the observed set, missing cells are filled by genuine
# GP interpolation and the interval can widen over gaps. (With a separable
# kernel a gap at one site is largely pinned down by other sites still reporting
# those weeks, so the widening is largest for region-wide blackouts.)
#
# Cost scales steeply with the number of sites (each draw is a full CG solve);
# reduce `n` at the top of the script or `n_draws` here to experiment quickly.
#
# The draws parallelise via future::plan(), controlled by `n_workers` in the
# Controls block: results are identical for every n_workers, and the braille
# progress bar shows only when running serially (workers can't tick it).
# -----------------------------------------------------------------------------
if (n_workers > 1) {
  future::plan(future::multisession, workers = n_workers)
} else {
  future::plan(future::sequential)
}
pred <- gp_predict(
  obs_data,
  coordinates,
  hyperparameters = est,
  nt = nt,
  period = period,
  n_draws = 100
)
future::plan(future::sequential) # back to serial for the rest of the session
pred_df <- transform(pred, id = as.integer(id))
cat(sprintf(
  "Dispersion r used for the interval (estimated) = %.1f  (true = %.1f)\n",
  attr(pred, "r"),
  true_r
))


# -----------------------------------------------------------------------------
# 7. Plot 3: predictions on top of the truth
# -----------------------------------------------------------------------------
# Blue line = posterior mean rate; blue ribbon = 95% count prediction interval.
# Compare against the black true-mean line, the grey observed points, and the
# red held-out truth.
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
    aes(t, rate),
    colour = "steelblue",
    linewidth = 0.6
  ) +
  labs(
    title = "Prediction vs truth: blue = posterior mean rate + 95% prediction interval",
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
# At the held-out (missing) cells the model never saw:
#   * correlation / RMSE of the predicted rate vs the TRUE rate (point quality)
#   * coverage of the held-out COUNTS by the 95% interval (calibration, ~0.95)
# -----------------------------------------------------------------------------
chk <- merge(
  missing_df[, c("id", "t", "lambda", "y")],
  pred_df,
  by = c("id", "t")
)
cat(sprintf("\nHeld-out cells: %d\n", nrow(chk)))
cat(sprintf(
  "corr(predicted rate, true rate) at held-out cells: %.3f\n",
  stats::cor(chk$rate, chk$lambda)
))
cat(sprintf(
  "RMSE of predicted vs true rate at held-out cells:   %.2f\n",
  sqrt(mean((chk$rate - chk$lambda)^2))
))
cat(sprintf(
  "95%% prediction-interval coverage of held-out COUNTS: %.2f  (target ~0.95)\n",
  mean(chk$y >= chk$lower & chk$y <= chk$upper)
))
