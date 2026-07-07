# Posterior variance of the latent field (exact part + control-variate draws)

Estimates the per-cell posterior variance of the GP field conditioned on
the observed cells, splitting it into an exact term and a small
Monte-Carlo correction: \$\$\operatorname{Var}(f) =
V\_{\mathrm{complete}} + \mathbb{E}\\\left\[d\_{\mathrm{obs}}^2 -
d\_{\mathrm{complete}}^2\right\].\$\$

## Usage

``` r
gp_posterior_var(
  obs_idx,
  N,
  space_mat,
  time_mat,
  noise_var,
  Rs_chol,
  Rt_chol,
  n_draws,
  tol = 0.001,
  progress_bar = NULL
)
```

## Arguments

- obs_idx:

  Integer indices of observed cells in the full vector.

- N:

  Total number of cells (`n * nt`).

- space_mat:

  Spatial kernel matrix (variance-scaled).

- time_mat:

  Temporal kernel matrix.

- noise_var:

  Scalar observation-noise variance \\\nu\\.

- Rs_chol, Rt_chol:

  Upper Cholesky factors of `space_mat` / `time_mat`.

- n_draws:

  Number of paired draws for the missing-data correction.

- tol:

  CG tolerance for the draw solves.

- progress_bar:

  Optional progress bar (from `make_curve_bar()`); ticked once per draw.
  Only effective under a single-worker plan (parallel workers cannot
  tick a bar in the calling session).

## Value

Numeric vector of length `N`: the posterior variance of the field.

## Details

\\V\_{\mathrm{complete}} = \operatorname{diag}(K - K(K+\nu I)^{-1}K)\\
is the posterior variance had *every* cell been observed – available in
closed form through the Kronecker eigendecomposition, no draws needed.
Missing cells change the variance only locally, so the draws are spent
purely on that correction: each zero-mean perturbation draw
\\d\_{\mathrm{obs}}\\ (Papandreou & Yuille 2010, one CG solve) is paired
with an exact complete-grid twin \\d\_{\mathrm{complete}}\\ built from
the *same* random numbers \\u, e\\, so their difference is nearly
noise-free away from gaps (a control variate). Modest draw counts
therefore give variances that plain Monte-Carlo would need hundreds of
draws to match, and cells far from any gap are essentially exact.

The draws are independent, so they run through
[`future.apply::future_lapply()`](https://future.apply.futureverse.org/reference/future_lapply.html):
serial under the default
[`future::plan()`](https://future.futureverse.org/reference/plan.html),
parallel when the caller selects a multi-worker plan.
`future.seed = TRUE` gives every draw its own pre-generated
L'Ecuyer-CMRG stream, so results are reproducible under
[`set.seed()`](https://rdrr.io/r/base/Random.html) and identical for
every backend and worker count.
