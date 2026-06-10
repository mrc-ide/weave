# =============================================================================
# Hybrid PFN evaluation: SBC + PFN/hybrid/fit_bayes comparison
#
# Builds on the previous evaluation (implementation/pfn_smoke/final_comparison.R)
# adding the hybrid path. Reads final_model.pt; writes:
#   - hybrid_sbc.rds        : rank-based SBC for the hybrid (50 datasets)
#   - hybrid_comparison.rds : 5-dataset PFN / hybrid / fit_bayes results
#   - hybrid_findings.md    : human summary
#
# Run from package root:
#   R_LIBS_USER="C:/Users/pwinskil/Documents/r_packages_arm64" \
#     Rscript implementation/pfn_smoke/hybrid_evaluation.R
# =============================================================================

suppressPackageStartupMessages({
  library(devtools)
})
devtools::load_all(".", quiet = TRUE)

OUT_DIR <- "implementation/pfn_smoke"
WEIGHTS <- file.path(OUT_DIR, "final_model.pt")
stopifnot(file.exists(WEIGHTS))

set.seed(20260611)
torch::torch_manual_seed(20260611)

model  <- pfn_load_model(WEIGHTS, device = "cpu")
coords <- model$spec$coords
nt     <- model$spec$nt
period <- model$spec$period
cat(sprintf("Geometry: n = %d, nt = %d, period = %d\n",
            nrow(coords), nt, period))

# ---------------------------------------------------------------------------
# Step 0: warm-up timing -- one hybrid run on a synthetic dataset so we know
# how long SBC will take.
# ---------------------------------------------------------------------------
cat("\n=== Warm-up timing for fit_pfn_hybrid ===\n")
priors  <- pfn_priors(coords, nt)
warmup  <- pfn_simulate_batch(coords, nt = nt, n_datasets = 1L,
                              priors = priors, period = period,
                              verbose = FALSE)
slice_to_obs_data <- function(batch, b) {
  n  <- nrow(batch$coords)
  nt <- batch$nt
  y    <- batch$y[b, , ]
  mask <- batch$mask[b, , ]
  y_obs_mat <- y
  y_obs_mat[mask == 0L] <- NA_integer_
  df <- tidyr::expand_grid(id = factor(seq_len(n)), t = seq_len(nt))
  df$lat   <- batch$coords$lat[as.integer(df$id)]
  df$lon   <- batch$coords$lon[as.integer(df$id)]
  df$y_obs <- as.vector(t(y_obs_mat))
  df
}
obs_w <- slice_to_obs_data(warmup, 1L)
t_warm <- system.time({
  fit_w <- fit_pfn_hybrid(obs_w, weights = WEIGHTS,
                          n_post = 200L, n_inner = 10L,
                          seed = 0L, verbose = FALSE)
})
cat(sprintf("Hybrid (n_post=200, n_inner=10) per dataset: %.1fs\n",
            t_warm["elapsed"]))
cat(sprintf("  inner PCG mean iters: %.1f\n", mean(fit_w$pcg_iters_mean)))

# Decide SBC budget from the warm-up
SBC_N        <- 50L
SBC_N_POST   <- 200L
SBC_N_INNER  <- 10L
estimated_min <- SBC_N * t_warm["elapsed"] / 60
cat(sprintf("Planned SBC: %d datasets => ~%.1f min\n",
            SBC_N, estimated_min))

# ---------------------------------------------------------------------------
# Step 1: rank-based SBC for the hybrid
#
# For each of SBC_N synthetic datasets:
#   - simulate (theta, mu_s, r, f, y) from the priors
#   - run hybrid inference to get N_POST posterior samples of each scalar
#   - compute the rank of the truth among the samples (in [0, 1])
# Under perfect calibration, ranks should be Uniform(0,1) per parameter.
# ---------------------------------------------------------------------------
cat("\n=== SBC: hybrid path ===\n")
sbc_batch <- pfn_simulate_batch(coords, nt = nt, n_datasets = SBC_N,
                                priors = priors, period = period,
                                verbose = FALSE)

# Rank quantile of `truth` among `samples` (with continuity correction so
# the values land in (0, 1) rather than {0, 1} at the extremes -- KS hates
# atoms there).
rank_of <- function(samples, truth) {
  (sum(samples < truth) + stats::runif(1)) / (length(samples) + 1)
}

rank_theta  <- matrix(NA_real_, SBC_N, 4L,
                      dimnames = list(NULL, c("length_scale", "periodic_scale",
                                              "long_term_scale", "log_r")))
# For mu_s and f, store per-cell ranks pooled across datasets (vector grows).
rank_mu_all <- numeric(0)
rank_f_all  <- numeric(0)

pb <- progress::progress_bar$new(
  format = "  SBC [:bar] :percent  dataset :current/:total  eta :eta",
  total = SBC_N, clear = FALSE, width = 70)

for (b in seq_len(SBC_N)) {
  obs <- slice_to_obs_data(sbc_batch, b)
  fit <- fit_pfn_hybrid(obs, weights = WEIGHTS,
                        n_post = SBC_N_POST, n_inner = SBC_N_INNER,
                        seed = b, verbose = FALSE)

  truth_theta <- exp(sbc_batch$log_theta[b, ])         # native scale
  truth_logr  <- sbc_batch$log_r[b]
  truth_mu    <- sbc_batch$mu_s[b, ]                   # length n
  truth_f     <- as.vector(t(sbc_batch$f[b, , ]))      # length N

  rank_theta[b, 1] <- rank_of(fit$theta_samples[, "length_scale"],
                              truth_theta[1])
  rank_theta[b, 2] <- rank_of(fit$theta_samples[, "periodic_scale"],
                              truth_theta[2])
  rank_theta[b, 3] <- rank_of(fit$theta_samples[, "long_term_scale"],
                              truth_theta[3])
  rank_theta[b, 4] <- rank_of(log(fit$r_samples), truth_logr)

  rank_mu_all <- c(rank_mu_all,
                   vapply(seq_len(design_n <- nrow(coords)),
                          function(i) rank_of(fit$mu_samples[, i], truth_mu[i]),
                          numeric(1)))
  rank_f_all <- c(rank_f_all,
                  vapply(seq_along(truth_f),
                         function(j) rank_of(fit$f_samples[, j], truth_f[j]),
                         numeric(1)))
  pb$tick()
}

ks_p <- function(v) suppressWarnings(stats::ks.test(v, "punif")$p.value)
sbc_ks <- list(
  length_scale    = ks_p(rank_theta[, "length_scale"]),
  periodic_scale  = ks_p(rank_theta[, "periodic_scale"]),
  long_term_scale = ks_p(rank_theta[, "long_term_scale"]),
  log_r           = ks_p(rank_theta[, "log_r"]),
  mu_s            = ks_p(rank_mu_all),
  f               = ks_p(rank_f_all)
)
sbc_result <- list(
  rank_theta = rank_theta,
  rank_mu    = rank_mu_all,
  rank_f     = rank_f_all,
  ks         = sbc_ks,
  n_datasets = SBC_N,
  n_post     = SBC_N_POST,
  n_inner    = SBC_N_INNER
)
saveRDS(sbc_result, file.path(OUT_DIR, "hybrid_sbc.rds"))

cat("\nHybrid SBC KS p-values vs Uniform(0,1) (lower => miscalibrated):\n")
for (nm in names(sbc_ks)) {
  p   <- sbc_ks[[nm]]
  tag <- if (is.na(p)) "" else if (p < 0.01) "  *** miscalibrated"
         else if (p < 0.05) "  *  borderline" else ""
  cat(sprintf("  %-18s  %.3f%s\n", nm, p, tag))
}

# ---------------------------------------------------------------------------
# Step 2: comparison on 5 fresh held-out datasets
# ---------------------------------------------------------------------------
cat("\n=== Comparison: PFN / hybrid / fit_bayes on 5 fresh datasets ===\n")
N_TEST <- 5L
test_batch <- pfn_simulate_batch(coords, nt = nt, n_datasets = N_TEST,
                                 priors = priors, period = period,
                                 verbose = FALSE)

ci_covers <- function(samples, truth, lvl = 0.95) {
  q <- stats::quantile(samples, c((1 - lvl) / 2, 1 - (1 - lvl) / 2),
                       na.rm = TRUE)
  unname(truth >= q[1] && truth <= q[2])
}

cover_row <- function(b, param, truth, samples_list) {
  do.call(data.frame, c(list(dataset = b, param = param, truth = truth),
    setNames(lapply(samples_list, function(s) ci_covers(s, truth)),
             paste0(names(samples_list), "_cover")),
    setNames(lapply(samples_list, function(s)
      unname(diff(stats::quantile(s, c(0.025, 0.975))))),
             paste0(names(samples_list), "_width")),
    setNames(lapply(samples_list, stats::median),
             paste0(names(samples_list), "_median"))
  ))
}

results <- list()
for (b in seq_len(N_TEST)) {
  cat(sprintf("\n-- Dataset %d/%d --\n", b, N_TEST))
  obs   <- slice_to_obs_data(test_batch, b)
  truth <- list(
    length_scale    = exp(test_batch$log_theta[b, 1]),
    periodic_scale  = exp(test_batch$log_theta[b, 2]),
    long_term_scale = exp(test_batch$log_theta[b, 3]),
    r               = exp(test_batch$log_r[b])
  )
  cat(sprintf("  truth: length=%.2f, periodic=%.2f, long_term=%.2f, r=%.2f\n",
              truth$length_scale, truth$periodic_scale,
              truth$long_term_scale, truth$r))

  t_p <- system.time({
    fit_p <- fit_pfn(obs, weights = WEIGHTS, n_post = 1000L, seed = b)
  })
  cat(sprintf("  fit_pfn         : %.2fs\n", t_p["elapsed"]))

  t_h <- system.time({
    fit_h <- fit_pfn_hybrid(obs, weights = WEIGHTS,
                            n_post = 500L, n_inner = 20L,
                            seed = b, verbose = FALSE)
  })
  cat(sprintf("  fit_pfn_hybrid  : %.2fs\n", t_h["elapsed"]))

  t_b <- system.time({
    fit_b <- fit_bayes(obs,
                       n_sweeps = 500L, burnin = 200L,
                       n_chains = 2L, period = period,
                       store_f  = "thin", n_thin_f = 100L,
                       verbose  = FALSE)
  })
  cat(sprintf("  fit_bayes       : %.2fs\n", t_b["elapsed"]))

  # Pool fit_bayes post-burnin theta + r across chains.
  keep <- (fit_b$burnin + 1L):fit_b$n_sweeps
  b_theta <- do.call(rbind, lapply(seq_len(fit_b$n_chains), function(ch) {
    fit_b$theta_trace[keep, , ch]
  }))
  b_r <- as.vector(fit_b$r_trace[keep, , drop = FALSE])

  # Per-parameter coverage rows.
  cov_df <- do.call(rbind, list(
    cover_row(b, "length_scale", truth$length_scale, list(
      pfn    = fit_p$theta_samples[, "length_scale"],
      hybrid = fit_h$theta_samples[, "length_scale"],
      bayes  = b_theta[, "length_scale"])),
    cover_row(b, "periodic_scale", truth$periodic_scale, list(
      pfn    = fit_p$theta_samples[, "periodic_scale"],
      hybrid = fit_h$theta_samples[, "periodic_scale"],
      bayes  = b_theta[, "periodic_scale"])),
    cover_row(b, "long_term_scale", truth$long_term_scale, list(
      pfn    = fit_p$theta_samples[, "long_term_scale"],
      hybrid = fit_h$theta_samples[, "long_term_scale"],
      bayes  = b_theta[, "long_term_scale"])),
    cover_row(b, "r", truth$r, list(
      pfn    = fit_p$r_samples,
      hybrid = fit_h$r_samples,
      bayes  = b_r))
  ))

  # Posterior-predictive coverage on observed cells.
  yrep_p <- posterior_predict(fit_p)
  yrep_h <- posterior_predict(fit_h)
  yrep_b <- posterior_predict(fit_b)
  obs_cells <- which(!is.na(obs$y_obs))
  y_obs     <- obs$y_obs[obs_cells]
  ppc <- function(yrep) {
    q <- t(apply(yrep[, obs_cells, drop = FALSE], 2, stats::quantile,
                 c(0.025, 0.975)))
    mean(y_obs >= q[, 1] & y_obs <= q[, 2])
  }

  results[[b]] <- list(
    dataset   = b,
    truth     = truth,
    cov_df    = cov_df,
    runtime   = c(pfn = unname(t_p["elapsed"]),
                  hybrid = unname(t_h["elapsed"]),
                  bayes  = unname(t_b["elapsed"])),
    ppc       = c(pfn = ppc(yrep_p), hybrid = ppc(yrep_h),
                  bayes = ppc(yrep_b))
  )
}
saveRDS(results, file.path(OUT_DIR, "hybrid_comparison.rds"))

# ---------------------------------------------------------------------------
# Step 3: aggregate + write findings
# ---------------------------------------------------------------------------
all_cov <- do.call(rbind, lapply(results, `[[`, "cov_df"))
agg <- aggregate(cbind(pfn_cover, hybrid_cover, bayes_cover,
                       pfn_width, hybrid_width, bayes_width) ~ param,
                 data = all_cov, FUN = mean)
agg <- agg[match(c("length_scale", "periodic_scale", "long_term_scale", "r"),
                 agg$param), ]
rt <- data.frame(
  dataset    = seq_len(N_TEST),
  pfn_sec    = vapply(results, function(r) r$runtime["pfn"],    numeric(1)),
  hybrid_sec = vapply(results, function(r) r$runtime["hybrid"], numeric(1)),
  bayes_sec  = vapply(results, function(r) r$runtime["bayes"],  numeric(1)),
  pfn_ppc    = vapply(results, function(r) r$ppc["pfn"],    numeric(1)),
  hybrid_ppc = vapply(results, function(r) r$ppc["hybrid"], numeric(1)),
  bayes_ppc  = vapply(results, function(r) r$ppc["bayes"],  numeric(1))
)
cat("\n=== Aggregate coverage (5 datasets) ===\n")
print(agg, row.names = FALSE)
cat("\n=== Runtime + PPC coverage ===\n")
print(rt, row.names = FALSE)

findings <- c(
  "# Hybrid PFN evaluation",
  "",
  sprintf("Same trained model as before (`final_model.pt`, 50k sims x 100 epochs,"),
  sprintf("n = %d, nt = %d, period = %d).", nrow(coords), nt, period),
  "",
  "## Hybrid SBC (rank-based, 50 fresh datasets, 200 posterior draws each)",
  "",
  "| parameter | KS p vs Uniform |",
  "|---|---|",
  sprintf("| length_scale    | %.3f |", sbc_ks$length_scale),
  sprintf("| periodic_scale  | %.3f |", sbc_ks$periodic_scale),
  sprintf("| long_term_scale | %.3f |", sbc_ks$long_term_scale),
  sprintf("| log_r           | %.3f |", sbc_ks$log_r),
  sprintf("| mu_s            | %.3g |", sbc_ks$mu_s),
  sprintf("| f               | %.3g |", sbc_ks$f),
  "",
  "## Comparison: 95% credible interval coverage and width (5 fresh datasets)",
  "",
  "| parameter | PFN cov | Hybrid cov | Bayes cov | PFN width | Hybrid width | Bayes width |",
  "|---|---|---|---|---|---|---|",
  sprintf("| %-15s | %.2f | %.2f | %.2f | %.2f | %.2f | %.2f |",
          agg$param,
          agg$pfn_cover, agg$hybrid_cover, agg$bayes_cover,
          agg$pfn_width, agg$hybrid_width, agg$bayes_width),
  "",
  "## Runtime + posterior-predictive coverage",
  "",
  "| dataset | PFN s | Hybrid s | Bayes s | PFN PPC | Hybrid PPC | Bayes PPC |",
  "|---|---|---|---|---|---|---|",
  sprintf("| %d | %.2f | %.2f | %.2f | %.2f | %.2f | %.2f |",
          rt$dataset, rt$pfn_sec, rt$hybrid_sec, rt$bayes_sec,
          rt$pfn_ppc, rt$hybrid_ppc, rt$bayes_ppc),
  "",
  sprintf("Hybrid mean speedup vs fit_bayes: **%.1fx**",
          mean(rt$bayes_sec) / mean(rt$hybrid_sec)),
  ""
)
writeLines(findings, file.path(OUT_DIR, "hybrid_findings.md"))
cat(sprintf("\nWrote %s\n", file.path(OUT_DIR, "hybrid_findings.md")))
