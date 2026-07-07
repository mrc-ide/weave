# Process raw epidemiological data for the GP model

Validates and prepares input data into the shape consumed by
[`infer_kernel_params()`](https://mrc-ide.github.io/weave/reference/infer_kernel_params.md)
and
[`gp_predict()`](https://mrc-ide.github.io/weave/reference/gp_predict.md):
it completes the site-by-time grid, drops sites that cannot be modelled,
assigns a factor site `id`, and returns the observations and coordinates
as separate, ready-to-use frames.

## Usage

``` r
data_process(data, ..., drop_zero = FALSE)
```

## Arguments

- data:

  A data frame containing site identifiers, time `t`, counts `n`, and
  coordinates `lat` and `lon`.

- ...:

  Bare (unquoted) column names that jointly identify a site, e.g.
  `data_process(raw, region, facility_name)`.

- drop_zero:

  Logical; passed to
  [`data_missing()`](https://mrc-ide.github.io/weave/reference/data_missing.md)
  – also drop sites whose observed counts sum to zero (default `FALSE`).

## Value

A list with three elements ready to pass to the model functions:

- `obs_data`:

  Observations with the site-identifier columns, the factor `id`, time
  `t`, and the count column `y_obs` (`NA` where missing).

- `coordinates`:

  One row per site with `id`, `lon` and `lat`.

- `nt`:

  The number of time points.
