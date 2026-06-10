# =============================================================================
# weave -- PFN Phase 1 dry-run walkthrough
#
# Trains a Prior-Fitted Network end-to-end on a tiny geometry (n=20, nt=52,
# ~10k synthetic datasets, a few epochs) and compares it head-to-head with
# fit_bayes() on a single held-out truth dataset. Targets ~10-20 min total
# on a laptop CPU.
#
# Mirrors implementation/test_walkthrough.R, which drives fit_bayes(). The
# point of this script is to be re-runnable as a regression harness during
# PFN development: re-run after any change to the simulator, model, training
# loop, or inference code and inspect the loss curve, SBC, and posterior
# overlays to confirm nothing has regressed.
#
# Sections:
#   1. Setup
#   2. Truth dataset (for head-to-head)
#   3. Synthetic training set (cached on disk)
#   4. Build + train the PFN (cached checkpoint)
#   5. SBC on a fresh held-out chunk
#   6. fit_bayes vs fit_pfn on the truth dataset
#   7. Posterior overlay + predictive coverage
#   8. Optional combined diagnostic PDF
# =============================================================================

# -----------------------------------------------------------------------------
# 1. Setup
# -----------------------------------------------------------------------------
devtools::load_all()
source("implementation/simulation.R")
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(torch)
})

set.seed(20260527)
torch_manual_seed(20260527)

n <- 20 # health facilities
nt <- 52 # one year of weekly data

# Cartesian coords -- pass `distance_fn = haversine_distance` to space_kernel
# for real lat/lon. `mu` is only used by simulate_data() (the simulation
# helpers in implementation/simulation.R require it).
coords <- data.frame(
  id = factor(1:n),
  lat = runif(n, 0, 5),
  lon = runif(n, 0, 5),
  mu = log(runif(n, 10, 80))
)

cache_dir <- "implementation/_pfn_cache"
sim_path <- file.path(cache_dir, "walkthrough_train.rds")
ckpt_path <- file.path(cache_dir, "walkthrough.pt")
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)


# -----------------------------------------------------------------------------
# 2. Truth dataset for head-to-head with fit_bayes()
# -----------------------------------------------------------------------------
# We hold one dataset out from the PFN's training pool and use it as the
# "true metric". fit_bayes runs MCMC on it; fit_pfn does one forward pass;
# the two posteriors should agree (modulo Monte Carlo noise) if the PFN is
# well-trained and well-calibrated.
true_length_scale <- 0.5
true_periodic_scale <- 1
true_long_term_scale <- 30
true_r <- 10

true_data <- simulate_data(
  n,
  nt,
  coords,
  space_kernel(coords, length_scale = true_length_scale),
  time_kernel(
    1:nt,
    periodic_scale = true_periodic_scale,
    long_term_scale = true_long_term_scale,
    period = 52
  ),
  r = true_r
)
obs_data <- observed_data(true_data, p_one = 0.2, p_switch = 0.3)

cat(sprintf(
  "Truth: %d cells, %d observed (%.0f%% missing)\n",
  n * nt,
  sum(!is.na(obs_data$y_obs)),
  100 * mean(is.na(obs_data$y_obs))
))


# -----------------------------------------------------------------------------
# 3. Synthetic training set
# -----------------------------------------------------------------------------
# pfn_simulate_batch() caches the entire batch to disk. The first run takes
# ~30s (10k datasets at 2 ms each); subsequent runs reload in <1s. Delete
# sim_path to force a re-simulation (e.g. after changing the prior).
cat("\n[3] Generating PFN training data ...\n")
train_t0 <- Sys.time()
train_batch <- pfn_simulate_batch(
  coords,
  nt = nt,
  n_datasets = 50000,
  cache_path = sim_path
)
cat(sprintf(
  "    %d datasets in %.1fs\n",
  train_batch$B,
  as.numeric(Sys.time() - train_t0, units = "secs")
))


# -----------------------------------------------------------------------------
# 4. Build + train the PFN
# -----------------------------------------------------------------------------
# Checkpoint is also cached; pfn_train() resumes from it if it exists. Delete
# ckpt_path to retrain from scratch.
cat("\n[4] Training PFN ...\n")
model <- nn_pfn(coords, nt, period = 52)
cat(sprintf("    Architecture: %d parameters\n", pfn_n_params(model)))

train_res <- pfn_train(
  model,
  train_batch,
  epochs = 5,
  batch_size = 64,
  lr = 3e-3,
  checkpoint = ckpt_path,
  verbose = TRUE
)

# Loss curve.
loss_plot <- ggplot(train_res$history, aes(x = epoch)) +
  geom_line(aes(y = train_total, colour = "train"), linewidth = 0.7) +
  geom_line(aes(y = val_total, colour = "val"), linewidth = 0.7) +
  labs(y = "loss (Gaussian NLL)", colour = NULL, title = "PFN training loss") +
  theme_bw() +
  theme(legend.position = "bottom")


# -----------------------------------------------------------------------------
# 5. SBC on a fresh held-out chunk
# -----------------------------------------------------------------------------
# 200 datasets drawn from the same prior, never seen during training. If the
# rank quantiles are non-uniform, the predicted posterior is mis-calibrated
# -- read off the diagnostic from the histogram shape (see R/pfn_sbc.R).
#
# At this scale (5 epochs, 10k datasets) we expect theta and r to be roughly
# calibrated and mu/f to still be a bit off. Both should improve with more
# training; if they're consistently broken after a full run, the architecture
# needs work.
cat("\n[5] SBC on 200 fresh synthetic datasets ...\n")
sbc <- pfn_sbc(model, n_datasets = 200, seed = 1)
print(sbc)
# Uncomment for an interactive plot:
# plot(sbc)

# -----------------------------------------------------------------------------
# 6. fit_bayes vs fit_pfn on the truth dataset
# -----------------------------------------------------------------------------
# fit_bayes runs on n=20, nt=52 in ~30-60s; fit_pfn is a single forward pass
# so it's effectively instantaneous. The speedup multiplier matters at scale
# (it's what motivated the whole PFN exercise), but the more important thing
# at the dry-run scale is whether the posteriors agree.
cat("\n[6] Head-to-head: fit_bayes vs fit_pfn\n")

bayes_t <- system.time({
  bayes_fit <- fit_bayes(
    obs_data,
    n_sweeps = 1500,
    burnin = 500,
    n_chains = 2,
    period = 52,
    store_f = "thin",
    n_thin_f = 500,
    verbose = FALSE
  )
})
cat(sprintf("    fit_bayes: %.1fs\n", bayes_t[3]))

pfn_t <- system.time({
  pfn_fit <- fit_pfn(obs_data, weights = ckpt_path, n_post = 1000, seed = 1)
})
cat(sprintf(
  "    fit_pfn  : %.2fs    (%.0fx speedup)\n",
  pfn_t[3],
  bayes_t[3] / pfn_t[3]
))


# -----------------------------------------------------------------------------
# 7. Posterior overlay + predictive coverage
# -----------------------------------------------------------------------------
truth_df <- data.frame(
  parameter = c("length_scale", "periodic_scale", "long_term_scale", "r"),
  value = c(
    true_length_scale,
    true_periodic_scale,
    true_long_term_scale,
    true_r
  )
)

bayes_post_idx <- (bayes_fit$burnin + 1):bayes_fit$n_sweeps

bayes_post_df <- data.frame(
  fit = "bayes",
  length_scale = as.vector(bayes_fit$theta_trace[
    bayes_post_idx,
    "length_scale",
  ]),
  periodic_scale = as.vector(bayes_fit$theta_trace[
    bayes_post_idx,
    "periodic_scale",
  ]),
  long_term_scale = as.vector(bayes_fit$theta_trace[
    bayes_post_idx,
    "long_term_scale",
  ]),
  r = as.vector(bayes_fit$r_trace[bayes_post_idx, ])
) |>
  tidyr::pivot_longer(-fit, names_to = "parameter", values_to = "value")

pfn_post_df <- data.frame(
  fit = "pfn",
  length_scale = pfn_fit$theta_samples[, "length_scale"],
  periodic_scale = pfn_fit$theta_samples[, "periodic_scale"],
  long_term_scale = pfn_fit$theta_samples[, "long_term_scale"],
  r = pfn_fit$r_samples
) |>
  tidyr::pivot_longer(-fit, names_to = "parameter", values_to = "value")

post_df <- dplyr::bind_rows(bayes_post_df, pfn_post_df)

posterior_plot <- ggplot(post_df, aes(x = value, fill = fit, colour = fit)) +
  geom_density(alpha = 0.3) +
  geom_vline(
    data = truth_df,
    aes(xintercept = value),
    colour = "red",
    linetype = "dashed",
    linewidth = 0.5
  ) +
  facet_wrap(~parameter, scales = "free") +
  theme_bw() +
  theme(legend.position = "bottom") +
  labs(
    title = "Marginal posteriors on the truth dataset",
    subtitle = "Red dashed = truth"
  )

# 80% predictive coverage check on held-out cells.
held_out_idx <- which(is.na(obs_data$y_obs))
true_y_held <- true_data$y[order(true_data$id, true_data$t)][held_out_idx]

bayes_yrep <- posterior_predict(bayes_fit)
pfn_yrep <- posterior_predict(pfn_fit)

cov80 <- function(yrep, idx, truth) {
  q <- apply(yrep[, idx], 2, stats::quantile, probs = c(0.1, 0.9))
  mean(truth >= q[1, ] & truth <= q[2, ])
}

cat("\n[7] 80% predictive coverage on held-out cells:\n")
cat(sprintf(
  "    fit_bayes : %.2f\n",
  cov80(bayes_yrep, held_out_idx, true_y_held)
))
cat(sprintf(
  "    fit_pfn   : %.2f\n",
  cov80(pfn_yrep, held_out_idx, true_y_held)
))


# -----------------------------------------------------------------------------
# 8. Optional combined diagnostic PDF
# -----------------------------------------------------------------------------
if (requireNamespace("patchwork", quietly = TRUE)) {
  out_pdf <- file.path(cache_dir, "walkthrough_diagnostic.pdf")
  combined <- patchwork::wrap_plots(
    list(loss_plot, posterior_plot),
    ncol = 1,
    heights = c(1, 2)
  )
  ggsave(out_pdf, combined, width = 10, height = 10)
  cat(sprintf("\nDiagnostic plot saved to %s\n", out_pdf))
}

cat("\n=== Walkthrough complete ===\n")
