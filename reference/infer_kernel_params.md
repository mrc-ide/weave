# Quick exact-marginal-likelihood estimate of the kernel hyperparameters

Estimates \`(length_scale, periodic_scale, long_term_scale)\` plus a
noise-to-signal ratio by maximising the exact GP marginal likelihood of
a plug-in latent field, using the Kronecker eigendecomposition so the
full \`(n \* nt)\`-square covariance is never formed. Fast and
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
    0.1)
)
```

## Arguments

- obs_data:

  Data frame with \`id\` (site), \`t\` (time) and the count column named
  by \`value\`.

- coordinates:

  Site coordinates (data frame with \`lon\`, \`lat\`), ordered to match
  \`sort(unique(obs_data\$id))\`.

- nt:

  Number of time points.

- period:

  Period of the seasonal cycle.

- value:

  Name of the count column (default \`"y_obs"\`).

- standardise:

  Logical; standardise the plug-in field per site (default \`TRUE\`).

- priors:

  Log-normal priors, see \[default_kernel_priors()\].

- start:

  Named/length-4 starting values on the natural scale (\`length_scale\`,
  \`periodic_scale\`, \`long_term_scale\`, \`nugget_ratio\`).

## Value

A list with \`length_scale\`, \`periodic_scale\`, \`long_term_scale\`,
\`nugget_ratio\`, the profiled \`sigma2\`, the maximised
\`log_posterior\`, and \`convergence\` (the \`optim\` code; 0 =
success).

## Details

This is the recommended quick estimator when a fast hyperparameter
estimate is wanted (e.g. as a starting point for a downstream sampler,
or as a standalone summary).
