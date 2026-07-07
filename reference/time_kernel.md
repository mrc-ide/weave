# Build the temporal correlation matrix

Builds a temporal correlation matrix by combining periodic and long-term
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

  A positive numeric scalar controlling how sharply correlation falls
  within each seasonal cycle: smaller values allow sharp seasonal peaks;
  larger values give a gentler, smoother cycle.

- long_term_scale:

  A positive numeric scalar for the long-term length-scale.

- nugget:

  A non-negative numeric scalar added to the diagonal for numerical
  stability.

- period:

  A positive numeric scalar giving the period of the seasonal component,
  in the same units as `times`.

## Value

A positive-definite matrix representing temporal correlation.
