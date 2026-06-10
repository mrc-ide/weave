# Fit hyperparmeters

Fit hyperparmeters

## Usage

``` r
fit(
  obs_data,
  coordinates,
  nt,
  period,
  n_sites,
  mask_prop = 0.2,
  verbose = FALSE,
  par0 = c(1, 5, 100),
  lower = c(1e-04, 0.8, 52 * 1.5),
  upper = c(2, 10, 500)
)
```

## Arguments

- obs_data:

  Data frame of observations with at least \`id\`, \`y_obs\`,
  \`mu_infer\`, and \`f_infer\` columns.

- coordinates:

  Data frame of site coordinates containing an \`id\` column that
  matches the site identifiers in \`obs_data\`.

- nt:

  Number of time points per site.

- period:

  Period used for the temporal kernel.

- n_sites:

  Number of sites to sample for optimisation.

- mask_prop:

  Proportion of observed points per site to hold out.

- verbose:

  Logical; whether to print objective values during optimisation.

- par0:

  Initial hyperparameter values passed to \`optim()\`.

- lower:

  Lower bounds for the hyperparameters passed to \`optim()\`.

- upper:

  Upper bounds for the hyperparameters passed to \`optim()\`.

## Value

A list as returned by \`stats::optim()\`.
