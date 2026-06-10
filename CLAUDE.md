# weave — project context

R package for full Bayesian inference on health-facility count data. The
generative model is

```
y_{i,t} ~ NB(r, mu_{i,t}),   log(mu_{i,t}) = f_{i,t} + mu_s_i
f       ~ N(0, K_space(theta_s) ⊗ K_time(theta_t))
```

— a separable spatio-temporal Gaussian process under a Negative-Binomial
observation model, with clustered missingness in practice.

## How to run things

```bash
# Tests
R_LIBS_USER="C:/Users/pwinskil/Documents/r_packages_arm64" \
  Rscript -e "devtools::test()"

# R CMD check
R_LIBS_USER="C:/Users/pwinskil/Documents/r_packages_arm64" \
  Rscript -e "devtools::check()"

# Full walkthrough (≈2 minutes on a laptop)
R_LIBS_USER="C:/Users/pwinskil/Documents/r_packages_arm64" \
  Rscript implementation/test_walkthrough.R
```

`Rscript` invoked from a non-interactive shell doesn't pick up
`R_LIBS_USER` from `~/.Renviron`, so always set it explicitly — see also
`~/.claude/CLAUDE.md` (global) for the same convention.

## Architecture

Single public inference path: **`fit_bayes()`**, a multi-chain Pólya-Gamma
augmented Gibbs sampler.

```
R/fit_bayes.R       Public entry: validate inputs, run chains
       │                 (serial or future.apply parallel),
       │                 pool draws, compute Rhat / PCG health.
       │
       ▼
R/pg_sampler.R      One sweep = five block updates:
                       1. omega | f, mu, r, y     (Pólya-Gamma)
                       2. f     | omega, mu, theta (exact Gaussian
                                                    via perturbation
                                                    sampler + PCG)
                       3. mu_s  | omega, f, y, r   (closed-form Normal)
                       4. theta | f                (three univariate
                                                    slice samplers, with
                                                    off-axis kernel cache
                                                    and Cholesky-in-slice)
                       5. r     | y, psi           (adaptive RW-MH on log r)
       │
       ▼
R/kron.R            Kronecker linear algebra. safe_chol / safe_eigen
R/pcg.R             with auto-jitter retry. Generic PCG used by f-block;
R/sample.R          quick_mvnorm uses Cholesky of the two small kernels.
R/kernel.R          RBF + periodic time kernel; haversine + Euclidean
                    spatial distance.
R/data.R            build_design() packs tidy obs_data into the design list
                    the sampler consumes.

R/posterior.R       Welford pooling across chains; posterior_predict().
R/diagnose.R        diagnose_bayes() — Rhat, ESS, KS-vs-prior per
                    parameter; flags prior-driven or poorly-mixed.

implementation/test_walkthrough.R   Heavily commented end-to-end demo
implementation/simulation.R         Ground-truth generators (NOT in pkg
                                      namespace — sourced by demos/tests)
implementation/bench_theta_block.R  Benchmark for the θ-slice block
```

## Things explored but not taken further

Three exploratory inference paths were tried and rejected on this branch.
Each is preserved somewhere — branch name / commit hash listed — so we
can pick them back up if a future use case justifies the work.

### 1. Elliptical Slice Sampler (ESS) — alternative f-block

**What:** Replace the PG-augmented exact Gaussian f-draw with Murray-
Adams-MacKay (2010) ESS, which only needs a prior draw and 2–5 likelihood
evaluations per step. Avoids PG augmentation, perturbation sampler, and
PCG entirely.

**Why dropped:** ESS mode-trapped catastrophically on θ. Rhat on
`length_scale` stayed at ~90 across 4 chains regardless of sweep count
(tried up to 7× more iterations). Within-chain ESS was healthy (~1100)
but chains found different posterior basins and never visited each
other. The fundamental problem is that ESS proposes f rotations on an
ellipse anchored at the current f, and the θ-slice's log-target is
driven by `f' K(θ)^-1 f`; the two updates are too coupled to escape
local basins. PG-Gibbs's exact f-draw decouples f from θ, which is what
makes the θ chain mix globally.

**Where:** Removed in commit "Remove ESS sampler and associated files"
(part of the now-orphaned exploration). The code (R/ess_sampler.R,
R/fit_bayes_ess.R, tests, parallel walkthrough, A/B benchmark) is no
longer on any branch; see the historical plan section in
`~/.claude/plans/ok-can-we-create-wiggly-token.md` for the full design.

### 2. Pure PFN (Prior-Fitted Network) — amortised one-shot inference

**What:** A neural network (1D conv temporal encoder + multi-head spatial
attention) trained on 50k simulated `(y, theta, r, mu_s, f)` draws to
predict diagonal-Gaussian posteriors for each component in one forward
pass. ~1400 lines of R + torch. End-to-end working, including SBC.

**Why dropped:** Three issues, only one fatal.
1. **Speed win is real** (~120× faster than fit_bayes at smoke scale)
   *but* only useful if you have many datasets to infer at the same
   geometry. The package's primary workflow is fitting a single
   posterior, where pure PFN's training-cost upfront doesn't amortise.
2. **Intervals are 8–10× wider than `fit_bayes`** on every parameter.
   Calibrated by SBC (truth in interval at the right rate) but
   conservative to the point of being uninformative for parameter
   inference on θ.
3. **f draws are jagged white noise.** The per-cell diagonal-Gaussian
   head cannot represent a strongly correlated 200-dim GP-shaped
   posterior. SBC on f returns KS p ≈ 0. This is architectural; no
   amount of training fixes it. Would require a richer head (low-rank +
   diagonal Gaussian, or normalising flow) to fix.

**Where:** Archived on `origin/pfn-poc-archive`. Includes:
- All `R/pfn_*.R` and `R/fit_pfn.R`
- The trained `final_model.pt` checkpoint (50k sims × 100 epochs)
- `implementation/pfn_smoke/findings.md` for the full write-up
- `hybrid_diagnostic.png` for the visual comparison

To revisit: `git checkout pfn-poc-archive`. To rebuild from scratch
without resuming the archive: re-run `implementation/pfn_smoke/
final_comparison.R` from that branch.

### 3. Hybrid PFN — amortised θ/r/mu_s + exact f-draw

**What:** Sample `(theta, r, mu_s)` from the PFN's diagonal-Gaussian
heads (which SBC says are calibrated), then draw f exactly via the
existing `pg_draw_f()` conditional on those. ~50 new lines composing
existing tested primitives.

**Why dropped:** Built to fix the jagged-f problem in path 2. It does
fix it — diagnostic plots show hybrid `f` draws indistinguishable from
`fit_bayes` (per-cell mean cor = 0.96, sd ratio median = 0.98,
coverage 96% vs nominal 95%). But the runtime story is bad: at smoke
scale (n = 8, nt = 26) the hybrid is **6× slower** than `fit_bayes`,
because each of the 500 posterior draws needs a fresh inner Gibbs
trajectory and the inner Gibbs's per-iteration cost is comparable to
`fit_bayes`'s per-sweep cost. We never tested whether the comparison
flips at production scale (n ≈ 1000) — the user's actual use case
hasn't been articulated clearly enough to justify a second retraining
run at that scale.

Also: the θ marginals from the PFN are still 8× wider than
`fit_bayes`'s, even though hybrid `f` matches. So the hybrid is correct
for predictions, but its θ posterior is still uninformative compared
to `fit_bayes`'s.

**Where:** Same archive — `origin/pfn-poc-archive`. Code in
`R/fit_pfn_hybrid.R`; comparison in `implementation/pfn_smoke/
hybrid_findings.md`.

## What's a sensible "next try" if we want to come back

Ranked by promise:

1. **Low-rank + diagonal Gaussian head on f**, retrained. The current
   diagonal head failed because it can't represent correlations. A head
   that predicts `Sigma = D + V V^T` with `V` of rank 16–32 captures
   the dominant correlation modes of a GP and is in the right
   distributional family. ~50–80 lines of new torch code + retraining
   (~30 min). Likely fixes pure-PFN's `f` problem outright while
   keeping the speed.
2. **Multivariate Gaussian head on `(θ, log r)`**, retrained. The
   current diagonal head treats the four scalar hyperparameters as
   independent, but they are correlated in the posterior. A 4×4
   Cholesky-factor head would tighten the marginals materially.
   ~30 lines + retraining.
3. **Production-scale benchmark** of the hybrid vs `fit_bayes` at
   n ≈ 200–500. The hybrid's edge gets bigger as n grows because
   `fit_bayes`'s θ-slice cost scales as n³. At smoke scale we saw the
   worst case for the hybrid.

## R environment

User runs Windows ARM. Local R toolchain:
- R binary: `C:/Program Files/R-aarch64/R-4.4.2/bin/Rscript`
- Library path: `C:/Users/pwinskil/Documents/r_packages_arm64`

torch 0.17.0 is confirmed to work on this platform (verified during the
PFN exploration). If you ever need the PFN branch again, `library(torch);
torch::torch_tensor(1)` will pull libtorch binaries on first call.
