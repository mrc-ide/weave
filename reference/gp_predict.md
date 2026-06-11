# Predict the latent rate and a count prediction interval (PCG)

Given kernel hyperparameters (e.g. from \[infer_kernel_params()\]),
predicts the latent rate \\\lambda = e^{\mu_s + f\_{st}}\\ at every
site-by-time cell by conditioning a separable Gaussian process on the
observed counts only. The posterior mean is obtained from a single
matrix-free PCG solve (so it is smooth and deterministic); the posterior
variance is estimated from \`n_draws\` perturbation draws. The
latent-rate posterior is then combined with Negative-Binomial
observation noise (law of total variance, lognormal moment-match) to
give a 95

## Usage

``` r
gp_predict(
  obs_data,
  coordinates,
  hyperparameters,
  nt,
  period,
  n_draws = 100,
  r = NULL,
  value = "y_obs",
  standardise = TRUE,
  pcg_tol = 1e-06
)
```

## Arguments

- obs_data:

  Data frame with \`id\` (site), \`t\` (time) and a count column named
  by \`value\` (\`NA\` where missing).

- coordinates:

  Site coordinates (data frame with \`id\`, \`lon\`, \`lat\`).

- hyperparameters:

  A list with elements \`length_scale\`, \`periodic_scale\`,
  \`long_term_scale\`, \`nugget_ratio\` and \`sigma2\` – the value
  returned by \[infer_kernel_params()\].

- nt:

  Number of time points.

- period:

  Period of the seasonal cycle.

- n_draws:

  Number of posterior draws used to estimate the variance (the
  prediction interval). Controls only the interval, not the mean. Use
  \`0\` to return the smooth posterior-mean rate only (one solve, no
  interval).

- r:

  Negative-Binomial dispersion for the count interval. If \`NULL\`
  (default) it is estimated by method of moments from the observed
  counts.

- value:

  Name of the count column (default \`"y_obs"\`).

- standardise:

  Logical; standardise the plug-in field per site (default \`TRUE\`),
  matching \[infer_kernel_params()\].

- pcg_tol:

  Convergence tolerance for the PCG solves.

## Value

A data frame with one row per cell and columns \`id\`, \`t\`, \`rate\`
(posterior point estimate of \\\lambda\\), and – when \`n_draws \>= 1\`
– \`lower\` and \`upper\` (the 95 used and \`n_draws\` are attached as
attributes.

## Details

Because the fit conditions on the observed set, missing cells are filled
by genuine GP interpolation and their prediction interval can widen over
gaps – unlike a completed-grid smoother that mean-imputes the gaps.

## See also

\[infer_kernel_params()\]
