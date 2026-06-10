# Estimate the temporal kernel

Builds a temporal covariance matrix by combining periodic and long-term
RBF components with a nugget term for numerical stability.

## Usage

``` r
time_kernel(
  times,
  periodic_scale,
  long_term_scale,
  nugget = 1e-09,
  period = 52
)
```

## Arguments

- times:

  A numeric vector of time indices.

- periodic_scale:

  A positive numeric scalar controlling the periodic variation.

- long_term_scale:

  A positive numeric scalar for the long-term length scale.

- nugget:

  A non-negative numeric scalar added to the diagonal for numerical
  stability.

- period:

  A positive numeric scalar giving the period of the seasonal component.

## Value

A positive-definite matrix representing temporal covariance.
