# =============================================================================
# Final PFN evaluation: SBC on the trained moderate-scale model + a
# coverage / runtime comparison against fit_bayes() on 5 fresh test datasets.
#
# Inputs:
#   - implementation/pfn_smoke/final_model.pt  (50k sims x 100 epochs, weights = (1,1,1))
#
# Outputs (in implementation/pfn_smoke/):
#   - final_sbc.rds       : pfn_sbc object (200 fresh datasets)
#   - final_comparison.rds: list with the per-dataset results
#   - findings.md         : human-readable summary
#
# Run from the package root:
#   R_LIBS_USER="C:/Users/pwinskil/Documents/r_packages_arm64" \
#     Rscript implementation/pfn_smoke/final_comparison.R
# =============================================================================

suppressPackageStartupMessages({
  library(devtools)
})
devtools::load_all(".", quiet = TRUE)

OUT_DIR <- "implementation/pfn_smoke"
WEIGHTS <- file.path(OUT_DIR, "final_model.pt")
stopifnot(file.exists(WEIGHTS))

set.seed(20260610)
torch::torch_manual_seed(20260610)

# ---- Load model + pull geometry --------------------------------------------
model <- pfn_load_model(WEIGHTS, device = "cpu")
coords <- model$spec$coords
nt     <- model$spec$nt
period <- model$spec$period
cat(sprintf("Model geometry: n = %d sites, nt = %d, period = %d\n",
            nrow(coords), nt, period))

# ---- Step 1: SBC on the trained model --------------------------------------
cat("\n=== SBC (200 fresh datasets) ===\n")
t_sbc <- system.time({
  sbc <- pfn_sbc(model, n_datasets = 200L, seed = 1L)
})
cat(sprintf("SBC took %.1fs\n", t_sbc["elapsed"]))
print(sbc)
saveRDS(sbc, file.path(OUT_DIR, "final_sbc.rds"))

# ---- Step 2: simulate 5 fresh test datasets, fit both samplers -------------
cat("\n=== Comparison: PFN vs fit_bayes() on 5 fresh datasets ===\n")
N_TEST <- 5L
priors_train <- pfn_priors(coords, nt)
test_batch   <- pfn_simulate_batch(coords, nt = nt, n_datasets = N_TEST,
                                   priors = priors_train, period = period,
                                   verbose = FALSE)

# Build a tidy obs_data frame from one slice of the batch.
slice_to_obs_data <- function(batch, b) {
  n  <- nrow(batch$coords)
  nt <- batch$nt
  y    <- batch$y[b, , ]
  mask <- batch$mask[b, , ]
  # NA at unobserved cells.
  y_obs_mat <- y
  y_obs_mat[mask == 0L] <- NA_integer_
  # tidyr::expand_grid then attach: id varies slowest, t fastest -- match
  # the convention build_design() expects.
  df <- tidyr::expand_grid(id = factor(seq_len(n)), t = seq_len(nt))
  df$lat   <- batch$coords$lat[as.integer(df$id)]
  df$lon   <- batch$coords$lon[as.integer(df$id)]
  # t(matrix(...)) is unnecessary: y_obs_mat is already n x nt so unrolling
  # row-wise with as.vector(t(.)) gives id-then-t order.
  df$y_obs <- as.vector(t(y_obs_mat))
  df
}

# Truth vectors stored alongside each result, for coverage scoring.
truth_for <- function(batch, b) {
  list(
    length_scale    = exp(batch$log_theta[b, 1]),
    periodic_scale  = exp(batch$log_theta[b, 2]),
    long_term_scale = exp(batch$log_theta[b, 3]),
    r               = exp(batch$log_r[b])
  )
}

ci_covers <- function(samples, truth, lvl = 0.95) {
  q <- stats::quantile(samples, c((1 - lvl) / 2, 1 - (1 - lvl) / 2),
                       na.rm = TRUE)
  unname(truth >= q[1] && truth <= q[2])
}

# Coverage rows are filled as: (dataset, param, pfn_cover, bayes_cover,
# pfn_width, bayes_width).
results <- list()

for (b in seq_len(N_TEST)) {
  cat(sprintf("\n-- Dataset %d/%d --\n", b, N_TEST))
  obs_data <- slice_to_obs_data(test_batch, b)
  truth    <- truth_for(test_batch, b)
  cat(sprintf("  truth: length=%.2f, periodic=%.2f, long_term=%.2f, r=%.2f\n",
              truth$length_scale, truth$periodic_scale,
              truth$long_term_scale, truth$r))

  # PFN
  t_pfn <- system.time({
    fit_p <- fit_pfn(obs_data, weights = WEIGHTS, n_post = 1000L, seed = b)
  })
  cat(sprintf("  fit_pfn: %.2fs\n", t_pfn["elapsed"]))

  # fit_bayes() at the smoke-test scale.
  t_b <- system.time({
    fit_b <- fit_bayes(obs_data,
                       n_sweeps = 500L, burnin = 200L,
                       n_chains = 2L, period = period,
                       store_f  = "thin", n_thin_f = 100L,
                       verbose  = FALSE)
  })
  cat(sprintf("  fit_bayes: %.2fs\n", t_b["elapsed"]))

  # ---- Per-parameter coverage ---------------------------------------------
  # fit_bayes exposes the full sweep-by-chain trace, not a flat samples
  # matrix; pool post-burnin draws across chains so the coverage and width
  # estimates use every retained sweep.
  keep    <- (fit_b$burnin + 1L):fit_b$n_sweeps
  b_theta <- do.call(rbind, lapply(seq_len(fit_b$n_chains), function(ch) {
    fit_b$theta_trace[keep, , ch]
  }))
  b_r <- as.vector(fit_b$r_trace[keep, , drop = FALSE])

  cover_one <- function(name, p_samples, b_samples) {
    pw <- diff(stats::quantile(p_samples, c(0.025, 0.975)))
    bw <- diff(stats::quantile(b_samples, c(0.025, 0.975)))
    data.frame(
      dataset       = b,
      param         = name,
      truth         = truth[[name]],
      pfn_cover     = ci_covers(p_samples, truth[[name]]),
      bayes_cover   = ci_covers(b_samples, truth[[name]]),
      pfn_width     = unname(pw),
      bayes_width   = unname(bw),
      pfn_median    = stats::median(p_samples),
      bayes_median  = stats::median(b_samples),
      row.names     = NULL
    )
  }
  cov_df <- rbind(
    cover_one("length_scale",
              fit_p$theta_samples[, "length_scale"],
              b_theta[, "length_scale"]),
    cover_one("periodic_scale",
              fit_p$theta_samples[, "periodic_scale"],
              b_theta[, "periodic_scale"]),
    cover_one("long_term_scale",
              fit_p$theta_samples[, "long_term_scale"],
              b_theta[, "long_term_scale"]),
    cover_one("r", fit_p$r_samples, b_r)
  )

  # ---- Posterior-predictive coverage on observed cells --------------------
  yrep_p <- posterior_predict(fit_p)
  yrep_b <- posterior_predict(fit_b)
  obs_idx <- which(!is.na(obs_data$y_obs))
  y_observed <- obs_data$y_obs[obs_idx]
  ppc_cov <- function(yrep) {
    q <- t(apply(yrep[, obs_idx, drop = FALSE], 2, stats::quantile,
                 c(0.025, 0.975)))
    mean(y_observed >= q[, 1] & y_observed <= q[, 2])
  }

  results[[b]] <- list(
    dataset    = b,
    truth      = truth,
    cov_df     = cov_df,
    runtime_pfn   = unname(t_pfn["elapsed"]),
    runtime_bayes = unname(t_b["elapsed"]),
    ppc_pfn       = ppc_cov(yrep_p),
    ppc_bayes     = ppc_cov(yrep_b)
  )
}

saveRDS(results, file.path(OUT_DIR, "final_comparison.rds"))

# ---- Step 3: aggregate + write findings ------------------------------------
all_cov <- do.call(rbind, lapply(results, `[[`, "cov_df"))
agg <- aggregate(cbind(pfn_cover, bayes_cover, pfn_width, bayes_width) ~ param,
                 data = all_cov, FUN = mean)
agg <- agg[match(c("length_scale", "periodic_scale", "long_term_scale", "r"),
                 agg$param), ]
rt <- data.frame(
  dataset   = seq_len(N_TEST),
  pfn_sec   = vapply(results, `[[`, numeric(1), "runtime_pfn"),
  bayes_sec = vapply(results, `[[`, numeric(1), "runtime_bayes"),
  pfn_ppc   = vapply(results, `[[`, numeric(1), "ppc_pfn"),
  bayes_ppc = vapply(results, `[[`, numeric(1), "ppc_bayes")
)

cat("\n=== Aggregate coverage (5 datasets) ===\n")
print(agg, row.names = FALSE)
cat("\n=== Runtime + PPC coverage ===\n")
print(rt, row.names = FALSE)

# ---- Findings markdown -----------------------------------------------------
findings <- c(
  "# PFN proof-of-concept: final findings",
  "",
  sprintf("Trained on **50,000 simulated datasets** for **100 epochs** at"),
  sprintf("(n = %d sites, nt = %d weeks, period = %d). Loss weights = (1, 1, 1).",
          nrow(coords), nt, period),
  sprintf("Model checkpoint: `%s`.", WEIGHTS),
  "",
  "## SBC on 200 fresh datasets",
  "",
  "KS p-value vs Uniform(0,1) per parameter (lower = miscalibrated):",
  "",
  "| parameter | KS p |",
  "|---|---|",
  sprintf("| length_scale    | %.3f |", sbc$ks$length_scale),
  sprintf("| periodic_scale  | %.3f |", sbc$ks$periodic_scale),
  sprintf("| long_term_scale | %.3f |", sbc$ks$long_term_scale),
  sprintf("| log_r           | %.3f |", sbc$ks$log_r),
  sprintf("| mu_s (pooled)   | %.3g |", sbc$ks$mu_s),
  sprintf("| f (pooled)      | %.3g |", sbc$ks$f),
  "",
  "## Comparison on 5 fresh held-out datasets",
  "",
  "95% credible-interval coverage and mean width (native scale), per parameter:",
  "",
  "| parameter | PFN cov | Bayes cov | PFN width | Bayes width |",
  "|---|---|---|---|---|",
  sprintf("| %-15s | %.2f | %.2f | %.2f | %.2f |",
          agg$param, agg$pfn_cover, agg$bayes_cover,
          agg$pfn_width, agg$bayes_width),
  "",
  "Per-dataset runtime and posterior-predictive coverage on observed cells:",
  "",
  "| dataset | PFN sec | Bayes sec | PFN PPC | Bayes PPC |",
  "|---|---|---|---|---|",
  sprintf("| %d | %.2f | %.2f | %.2f | %.2f |",
          rt$dataset, rt$pfn_sec, rt$bayes_sec, rt$pfn_ppc, rt$bayes_ppc),
  "",
  "## Interpretation",
  "",
  sprintf("- **Speed**: PFN inference is %.0fx faster than `fit_bayes()` (mean %.2fs vs %.2fs).",
          mean(rt$bayes_sec) / mean(rt$pfn_sec),
          mean(rt$pfn_sec), mean(rt$bayes_sec)),
  "- **Theta calibration**: the four scalar hyperparameters",
  "  show whether scaling training data + epochs tightened SBC. Compare KS",
  "  values against the small-scale baseline (8k sims, 40 epochs).",
  "- **mu_s / f calibration**: these remained miscalibrated in the small-",
  "  scale runs because the per-cell diagonal Gaussian head cannot express",
  "  the 200-dim correlated posterior. If the KS p-values here are still",
  "  ~0, that's confirmation that scale alone won't fix it -- the head",
  "  architecture needs to be richer (low-rank+diagonal, or normalising",
  "  flow) for follow-up work.",
  ""
)
writeLines(findings, file.path(OUT_DIR, "findings.md"))
cat(sprintf("\nWrote %s\n", file.path(OUT_DIR, "findings.md")))
