# Estimate the spatial kernel

Builds a spatial covariance matrix using an RBF kernel with a nugget
term for numerical stability.

## Usage

``` r
space_kernel(coordinates, length_scale, nugget = 1e-09)
```

## Arguments

- coordinates:

  A data frame with columns \`lon\` and \`lat\` in degrees.

- length_scale:

  A positive numeric scalar for the spatial length scale.

- nugget:

  A non-negative numeric scalar added to the diagonal for numerical
  stability.

## Value

A positive-definite matrix representing spatial covariance.
