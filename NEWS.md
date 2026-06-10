# weave 0.3.0

## Breaking changes

* **Removed `fit()`**, the legacy fast point estimator under a linearised-
  Gaussian surrogate likelihood. `fit_bayes()` is now the only inference
  entry point. The Kronecker linear-algebra primitives (`kron_*` in
  `R/kron.R`) and the PCG solver (`R/pcg.R`) that `fit()` shared with
  `fit_bayes()` are retained; only the public estimator is gone.
* Removed the `weave:::hutchinson_diag()` internal helper that powered the
  stochastic posterior-variance estimate inside `fit()`.

# weave 0.2.0

## New features

* **`fit_bayes()`** — full Bayesian inference for the NB-GP model via a
  Pólya–Gamma augmented Gibbs sampler. No Laplace approximation: the
  latent-field block is an exact Gaussian draw via matrix-free PCG, the
  hyperparameters are updated by univariate slice samplers on the log scale,
  and the NB dispersion `r` is updated by adaptive random-walk MH on `log r`.
* **`build_design()`** — single utility that packs a tidy data frame into
  the design list both `fit()` and `fit_bayes()` consume (flat count
  vector, observed-index vector, per-site coordinates, initial site
  intercepts).
* **`haversine_distance()`** — great-circle distance helper for use with
  `space_kernel(distance_fn = haversine_distance)` when working with real
  geographic coordinates.
* **`posterior_predict()`** — draws of `y_rep` from a `weave_bayes` fit for
  posterior-predictive checks.
* `simulate_data()` (in `implementation/simulation.R`) now accepts an `r`
  argument, defaulting to `Inf` for the Poisson limit.

## Breaking changes

* `fit()` no longer reads `mu_infer` / `f_infer` from the input data; it
  derives the per-site offset internally from `y_obs` (or `n`). Old
  pipelines that hand-constructed `obs_data` with these columns need to be
  re-plumbed; new code should pass a tidy data frame and let `build_design()`
  do the bookkeeping.
* `fit()`'s `hyperparameters` argument is now a named list
  `list(length_scale, periodic_scale, long_term_scale)`. Length-3 numeric
  vectors are still accepted with a deprecation warning.
* Removed: `infer_space_kernel_params()`, `infer_time_kernel_params()` —
  empirical-correlation hyperparameter inference. Replaced by the slice
  updates inside `fit_bayes()`.
* Removed: the dead `llh()` Hessian-builder and the Jacobi `M_inv()`
  preconditioner from `fit.R`.
* `data_process()` no longer computes `observed_sigmasq` or `start_par`;
  it now produces a single `mu_init` column (renamed from `observed_mu`).

## Bug fixes

* **Posterior variance in `fit()`** was incorrect: the previous expression
  conflated `diag(Σ⁻¹)` with `diag(Σ⁻¹)⊗diag(Σ⁻¹)` and then took `−1/H_ii`
  rather than `(−H⁻¹)_ii`. Replaced with a stochastic Hutchinson diagonal
  of the linearised posterior covariance, computed via PCG against
  Rademacher probes.
* PCG's `maxit` warning now fires only on actual non-convergence (it used
  to fire when `it == maxit` even if that iteration converged).
* `regularise()` now uses a relative jitter (`lambda * mean(diag(x))`)
  rather than an absolute one — scale-safe.

## Internals

* Linear-algebra primitives split out of `fit.R` into `R/kron.R`
  (Kronecker matvec, eigen cache, Kronecker-eigen solver, log-det,
  quadratic form) and `R/pcg.R` (generic PCG, observed-system matvec,
  Kronecker-eigen and Jacobi preconditioner closures).
* PCG in both `fit()` and the f-block of `fit_bayes()` now uses the
  Kronecker-eigen preconditioner by default.

## Dependencies

* New `Imports`: `BayesLogit` (Pólya–Gamma sampler — supports real-valued
  shape, required for NB with non-integer `r`), `progress`, `dplyr`,
  `tidyr`, `stats`.
* New `Suggests`: `coda`, `lifecycle`, `ggplot2`.
