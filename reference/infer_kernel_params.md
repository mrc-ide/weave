# Quick exact-marginal-likelihood estimate of the kernel hyperparameters

Estimates `(length_scale, periodic_scale, long_term_scale)` plus a
noise-to-signal ratio by maximising the exact GP marginal likelihood of
a plug-in latent field, using the Kronecker eigendecomposition so the
full `(n * nt)`-square covariance is never formed. Fast and
deterministic – no MCMC, no iterative solver.

## Usage

``` r
infer_kernel_params(
  obs_data,
  coordinates,
  nt,
  period,
  value = "y_obs",
  standardise = TRUE,
  priors = default_kernel_priors(),
  start = c(length_scale = 1, periodic_scale = 1, long_term_scale = 100, nugget_ratio =
    0.1),
  n_sites = NULL,
  refine = FALSE,
  refine_iter = 3L
)
```

## Arguments

- obs_data:

  Data frame with `id` (site), `t` (time) and the count column named by
  `value`. `t` is a numeric time index whose *differences* encode real
  elapsed time, so gaps and uneven spacing between time points are
  modelled as genuine time distances (use e.g. weeks or days since a
  reference).
  [`gp_predict()`](https://mrc-ide.github.io/weave/reference/gp_predict.md)
  must be given the same `t` encoding.

- coordinates:

  Site coordinates: a data frame with `id`, `lon`, `lat`, one row per
  site in `obs_data` (rows are matched by `id`, so order does not
  matter).

- nt:

  Number of time points.

- period:

  Period of the seasonal cycle, in the same units as `t`.

- value:

  Name of the count column (default `"y_obs"`).

- standardise:

  Logical; standardise the plug-in field per site (default `TRUE`).

- priors:

  Log-normal priors, see
  [`default_kernel_priors()`](https://mrc-ide.github.io/weave/reference/default_kernel_priors.md).

- start:

  Named/length-4 starting values on the natural scale (`length_scale`,
  `periodic_scale`, `long_term_scale`, `nugget_ratio`).

- n_sites:

  Optional integer. If supplied and smaller than the number of sites,
  the hyperparameters are estimated from a random subsample of this many
  sites. The kernel hyperparameters are shared, population-level
  quantities, so a representative site subsample estimates the same
  length-scales at a fraction of the \\O(n^3)\\ cost – useful for very
  large site counts. Default `NULL` uses all sites. The subsample is
  drawn from the current RNG state, so set a seed beforehand (e.g.
  [`set.seed()`](https://rdrr.io/r/base/Random.html)) for a reproducible
  estimate. Note: this subsamples *sites* only, not time points (the
  temporal kernel needs the full series to resolve the periodic and
  long-term scales).

- refine:

  Logical; if `TRUE`, run `refine_iter` EM-style refinement passes that
  re-fit after filling the gaps with the GP conditional mean (see
  Details). Default `FALSE` (the fast single-pass estimate). Recommended
  when missingness is non-trivial.

- refine_iter:

  Number of refinement passes when `refine = TRUE` (default `3`).
  Ignored when `refine = FALSE`.

## Value

A list with `length_scale`, `periodic_scale`, `long_term_scale`,
`nugget_ratio`, the profiled `sigma2`, the maximised `log_posterior`,
and `convergence` (the `optim` code; 0 = success).

## Details

This is the recommended quick estimator when a fast hyperparameter
estimate is wanted (e.g. as a starting point for a downstream sampler,
or as a standalone summary).

Strictly, the score maximised is the marginal likelihood *plus* weakly
informative log-normal priors on the four parameters (a MAP estimate;
see
[`default_kernel_priors()`](https://mrc-ide.github.io/weave/reference/default_kernel_priors.md)).
The priors act as mild regularisation that keeps weakly identified
parameters – notably `long_term_scale` – away from the boundary; pass
`priors` to change them.

Set `refine` to enable an EM-style refinement that removes the bias
missing cells introduce. Each pass refits after replacing the gaps with
the GP posterior (conditional) mean under the current estimate – a
correlation-aware fill, not the flat mean-imputation – using the same
matrix-free CG solve as
[`gp_predict()`](https://mrc-ide.github.io/weave/reference/gp_predict.md).
The expensive observed-cell solve runs only once per pass (not inside
the optimiser), so it stays cheap, and it typically converges in 2-3
passes to the estimate you would get with no missing data at all. It
does not remove the intrinsic plug-in attenuation (conditioning on a
noisy field rather than integrating the latent field out), only the part
caused by the gaps.
