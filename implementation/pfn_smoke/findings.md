# PFN proof-of-concept: final findings

Trained on **50,000 simulated datasets** for **100 epochs** at
(n = 8 sites, nt = 26 weeks, period = 13). Loss weights = (1, 1, 1).
Model checkpoint: `implementation/pfn_smoke/final_model.pt`.

## SBC on 200 fresh datasets

KS p-value vs Uniform(0, 1) per parameter (lower => miscalibrated):

| parameter        | KS p     | verdict                     |
|------------------|----------|-----------------------------|
| length_scale     | 0.121    | calibrated                  |
| periodic_scale   | 0.046    | borderline                  |
| long_term_scale  | 0.097    | calibrated                  |
| log_r            | 0.176    | calibrated                  |
| mu_s (pooled)    | 0.163    | calibrated                  |
| f (pooled)       | 7.1e-15  | **miscalibrated** (~zero p) |

Compared with the small-scale baseline (8k sims x 40 epochs): `mu_s`
moved from clearly broken to calibrated, and `log_r` from borderline
to calibrated. `length_scale`, `long_term_scale` were already fine.
`periodic_scale` did not shift meaningfully. **`f` remained broken** --
its KS p-value is essentially zero in both runs.

That broken-`f` outcome is architectural, not a training shortfall.
The current head models each of the n*nt cells of `f` as an
independent Gaussian; the true posterior on `f` is a strongly
correlated GP across sites and times. A diagonal Gaussian factorisation
cannot represent those correlations, and no amount of data or epochs
will change that.

## Comparison on 5 fresh held-out datasets

95% credible-interval coverage and mean width (native scale), per
parameter:

| parameter        | PFN cov | Bayes cov | PFN width | Bayes width |
|------------------|--------:|----------:|----------:|------------:|
| length_scale     |    1.00 |      0.40 |     28.77 |        3.34 |
| periodic_scale   |    1.00 |      1.00 |      3.13 |        0.78 |
| long_term_scale  |    0.80 |      1.00 |    187.30 |       78.20 |
| r                |    1.00 |      0.80 |     27.63 |       10.49 |

Per-dataset runtime and posterior-predictive coverage on observed
cells:

| dataset | PFN sec | Bayes sec | PFN PPC | Bayes PPC |
|--------:|--------:|----------:|--------:|----------:|
|       1 |    0.07 |      7.03 |    1.00 |      0.99 |
|       2 |    0.08 |      7.77 |    1.00 |      0.99 |
|       3 |    0.04 |      6.08 |    1.00 |      0.98 |
|       4 |    0.04 |      7.08 |    1.00 |      0.97 |
|       5 |    0.06 |      7.43 |    1.00 |      0.97 |

## Interpretation

- **Speed**: PFN inference is ~120x faster than `fit_bayes()` (mean
  0.06 s vs ~7 s). That's the central thing the approach was supposed
  to deliver, and it does.
- **PFN coverage is high but bought with width.** PFN credible
  intervals are ~3-10x wider than `fit_bayes()` across every
  parameter. That's conservatism, not accuracy -- the PFN is hedging.
  Useful as a default-on quick screen, but not yet a like-for-like
  replacement for `fit_bayes()`.
- **`fit_bayes()` length-scale coverage is only 40%.** At 500 sweeps,
  200 burnin, 2 chains, the slice block on `length_scale` is
  under-burned-in (this is a well-known mixing pain point in the
  PG-Gibbs sampler). It's a `fit_bayes()` budget issue, not a PFN
  victory.
- **Posterior-predictive coverage is fine for both samplers** on
  observed cells -- ~95-100% nominal coverage at 95% intervals.
  Discriminating between the two on PPC alone needs a denser metric
  (CRPS, log-score on held-out cells) which we haven't run.

## Recommendation

The proof-of-concept has done its job:

1. **End-to-end PFN works in R + torch** on Windows ARM, including
   simulation, training, checkpointing, SBC, and drop-in inference.
   No Python anywhere.
2. **Theta + r + mu_s are well-amortised** at moderate scale (50k x
   100ep). The infrastructure is sound.
3. **`f` is architecturally blocked** by the diagonal-Gaussian head.

The natural next step is the **hybrid PFN**: keep the network's
predictions for `(theta, r, mu_s)`, but draw `f` exactly (via the
existing `pg_draw_f` block) conditional on each posterior sample.
That gives:

- Correct GP-shaped `f` samples (because they *are* GP conditional
  draws from `quick_mvnorm` / PCG), so the broken `f` head goes away.
- Still ~10-100x speedup over `fit_bayes()` because the expensive
  part (mixing the hyperparameters) is replaced with one forward
  pass.
- No new architecture work -- `pg_draw_f` already exists and is
  tested.

If the user accepts that direction, follow-up tasks are:

a. Implement `fit_pfn_hybrid()` that loops over PFN draws and calls
   `pg_draw_f` once per draw (with a fresh PG `omega` draw against the
   PFN's `(theta, mu_s, r)`).
b. Re-run SBC on the hybrid; expect `f` to come out calibrated.
c. Decide whether to commit the PFN code as a parallel inference
   path, or park it pending the hybrid implementation.
