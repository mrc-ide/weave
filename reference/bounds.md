# Compute quantile bounds for GP-based count draws

This function generates draws from a Gaussian process state, applies
Poisson observation noise, and computes quantile bounds of the resulting
counts for each observation.

## Usage

``` r
bounds(
  state,
  n_lambda = 30,
  n_draw = 100,
  quantiles = c(0.025, 0.25, 0.75, 0.975)
)
```

## Arguments

- state:

  A fitted GP state object used as input to
  [`gp_draw()`](https://mrc-ide.github.io/weave/reference/gp_draw.md).

- n_lambda:

  Integer. Number of latent Gaussian process draws.

- n_draw:

  Integer. Number of replicate draws with Poisson noise applied.

- quantiles:

  Numeric vector. Quantiles to compute for the count draws (values
  between 0 and 1).

## Value

A data frame with one row per observation. Contains the columns:

- `id` Observation identifier.

- `t` Observation time index.

- Additional columns for each quantile requested, named `q{quantile}`
  (e.g. `q0.025`, `q0.5`).

## See also

[`gp_draw`](https://mrc-ide.github.io/weave/reference/gp_draw.md),
[`quantile`](https://rdrr.io/r/stats/quantile.html)
