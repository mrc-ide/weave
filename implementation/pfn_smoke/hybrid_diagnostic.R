# =============================================================================
# Diagnostic plots: hybrid PFN vs fit_bayes for f, on a single dataset.
#
# Three panels:
#   1. f time series at 4 representative sites, with 95% bands and truth
#   2. scatter of posterior f means: hybrid vs Bayes (one point per cell)
#   3. scatter of posterior f sds:   hybrid vs Bayes (one point per cell)
#
# Question: is the hybrid f draw shape "the same posterior as Bayes, just
# wider" (i.e. correct family, conservative) or "the wrong shape"?
# =============================================================================

suppressPackageStartupMessages({
  library(devtools)
  library(ggplot2)
  library(patchwork)
})
devtools::load_all(".", quiet = TRUE)

OUT_DIR <- "implementation/pfn_smoke"
WEIGHTS <- file.path(OUT_DIR, "final_model.pt")

# Reproduce dataset 1 from hybrid_evaluation.R (same seed cascade).
set.seed(20260611)
torch::torch_manual_seed(20260611)

model  <- pfn_load_model(WEIGHTS, device = "cpu")
coords <- model$spec$coords
nt     <- model$spec$nt
period <- model$spec$period
priors <- pfn_priors(coords, nt)

# The hybrid script does: warm-up sim (1 dataset), SBC batch (50), test batch
# (5). We need to consume the same RNG to land on the same test dataset.
warmup    <- pfn_simulate_batch(coords, nt = nt, n_datasets = 1L,
                                priors = priors, period = period,
                                verbose = FALSE)
# Skip the SBC loop in the diagnostic (it would burn ~10 min of RNG advances
# we don't need). Just simulate a fresh batch; we want any single dataset
# in the right regime, not the exact one from the SBC log.
test_batch <- pfn_simulate_batch(coords, nt = nt, n_datasets = 1L,
                                 priors = priors, period = period,
                                 verbose = FALSE)
b <- 1L

# Slice -> obs_data tidy df.
n  <- nrow(coords); N <- n * nt
y_mat    <- test_batch$y[b, , ]
mask_mat <- test_batch$mask[b, , ]
y_mat[mask_mat == 0L] <- NA_integer_
df <- tidyr::expand_grid(id = factor(seq_len(n)), t = seq_len(nt))
df$lat   <- coords$lat[as.integer(df$id)]
df$lon   <- coords$lon[as.integer(df$id)]
df$y_obs <- as.vector(t(y_mat))

truth <- list(
  length_scale    = exp(test_batch$log_theta[b, 1]),
  periodic_scale  = exp(test_batch$log_theta[b, 2]),
  long_term_scale = exp(test_batch$log_theta[b, 3]),
  r               = exp(test_batch$log_r[b])
)
cat(sprintf("Truth: length=%.2f, periodic=%.2f, long_term=%.2f, r=%.2f\n",
            truth$length_scale, truth$periodic_scale,
            truth$long_term_scale, truth$r))

truth_f_mat <- test_batch$f[b, , ]   # n x nt

# ---- Fit both samplers -----------------------------------------------------
cat("Running fit_pfn_hybrid (n_post=500, n_inner=20)...\n")
t_h <- system.time({
  fit_h <- fit_pfn_hybrid(df, weights = WEIGHTS,
                          n_post = 500L, n_inner = 20L,
                          seed = 1L, verbose = FALSE)
})
cat(sprintf("  %.1fs\n", t_h["elapsed"]))

cat("Running fit_bayes (500 sweeps, 200 burnin, 2 chains)...\n")
t_b <- system.time({
  fit_b <- fit_bayes(df, n_sweeps = 500L, burnin = 200L, n_chains = 2L,
                     period = period, store_f = "thin", n_thin_f = 100L,
                     verbose = FALSE)
})
cat(sprintf("  %.1fs\n", t_b["elapsed"]))

# ---- Reshape per-cell summaries -------------------------------------------
# Both fit$f_samples are (n_post, N) with time-fastest cell ordering. Pull
# 0.025/0.5/0.975 quantiles per cell.
qs <- function(samples_mat) {
  q <- apply(samples_mat, 2, stats::quantile, c(0.025, 0.5, 0.975))
  data.frame(lo = q[1, ], med = q[2, ], hi = q[3, ])
}
h_summ <- qs(fit_h$f_samples)
b_summ <- qs(fit_b$f_samples)

# Time-fastest -> (site, t)
site_t <- tidyr::expand_grid(id = factor(seq_len(n)), t = seq_len(nt))
site_t$truth_f <- as.vector(t(truth_f_mat))

site_t$h_lo  <- h_summ$lo;  site_t$h_med <- h_summ$med;  site_t$h_hi  <- h_summ$hi
site_t$b_lo  <- b_summ$lo;  site_t$b_med <- b_summ$med;  site_t$b_hi  <- b_summ$hi

# Pick 4 sites (spread across the spatial grid).
chosen <- as.integer(round(seq(1, n, length.out = 4)))
sub <- site_t[as.integer(site_t$id) %in% chosen, ]
sub$id <- factor(sub$id, levels = chosen,
                 labels = paste("Site", chosen))

# ---- Plot 1: per-site f time series ---------------------------------------
p1 <- ggplot(sub, aes(t)) +
  geom_ribbon(aes(ymin = h_lo, ymax = h_hi, fill = "hybrid"), alpha = 0.30) +
  geom_ribbon(aes(ymin = b_lo, ymax = b_hi, fill = "bayes"),  alpha = 0.30) +
  geom_line(aes(y = h_med, colour = "hybrid"), linewidth = 0.9) +
  geom_line(aes(y = b_med, colour = "bayes"),  linewidth = 0.9) +
  geom_line(aes(y = truth_f), colour = "black", linewidth = 0.9,
            linetype = "dashed") +
  facet_wrap(~ id, ncol = 2, scales = "free_y") +
  scale_colour_manual(values = c(hybrid = "#1f77b4", bayes = "#d62728"),
                      name = "posterior median") +
  scale_fill_manual(values = c(hybrid = "#1f77b4", bayes = "#d62728"),
                    name = "95% interval") +
  labs(title = "f draws: hybrid (blue) vs fit_bayes (red); truth = black dashed",
       x = "week", y = "f") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom",
        strip.text = element_text(size = 14, face = "bold"))

# ---- Plot 2/3: per-cell mean and sd scatter -------------------------------
h_mean <- colMeans(fit_h$f_samples); b_mean <- colMeans(fit_b$f_samples)
h_sd   <- apply(fit_h$f_samples, 2, stats::sd)
b_sd   <- apply(fit_b$f_samples, 2, stats::sd)

p2 <- ggplot(data.frame(b_mean, h_mean), aes(b_mean, h_mean)) +
  geom_point(alpha = 0.4, size = 1.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              colour = "grey30") +
  labs(title = sprintf("posterior mean of f, per cell (cor = %.2f)",
                       stats::cor(b_mean, h_mean)),
       x = "fit_bayes", y = "hybrid") +
  theme_minimal(base_size = 14)

p3 <- ggplot(data.frame(b_sd, h_sd), aes(b_sd, h_sd)) +
  geom_point(alpha = 0.4, size = 1.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              colour = "grey30") +
  labs(title = sprintf("posterior sd of f, per cell (ratio median = %.2f)",
                       stats::median(h_sd / b_sd)),
       x = "fit_bayes", y = "hybrid") +
  theme_minimal(base_size = 14)

# Layout: row of time series across the top; scatters below
out <- p1 / (p2 | p3) + plot_annotation(
  title = sprintf("Hybrid PFN vs fit_bayes -- single dataset f comparison"),
  subtitle = sprintf("hybrid ran in %.1fs, fit_bayes in %.1fs",
                     t_h["elapsed"], t_b["elapsed"]),
  theme = theme(plot.title = element_text(face = "bold", size = 16),
                plot.subtitle = element_text(size = 13))
)

png_path <- file.path(OUT_DIR, "hybrid_diagnostic.png")
ggsave(png_path, out, width = 16, height = 14, dpi = 150)
cat(sprintf("Wrote %s\n", png_path))

# Print key numerics for the terminal.
cat(sprintf("\nPer-cell mean correlation (hybrid vs Bayes): %.3f\n",
            stats::cor(b_mean, h_mean)))
cat(sprintf("Per-cell sd ratio (hybrid / Bayes): median %.2f, IQR [%.2f, %.2f]\n",
            stats::median(h_sd / b_sd),
            stats::quantile(h_sd / b_sd, 0.25),
            stats::quantile(h_sd / b_sd, 0.75)))
cat(sprintf("Truth-in-95%%-interval rate: hybrid %.2f, Bayes %.2f\n",
            mean(site_t$truth_f >= site_t$h_lo & site_t$truth_f <= site_t$h_hi),
            mean(site_t$truth_f >= site_t$b_lo & site_t$truth_f <= site_t$b_hi)))
