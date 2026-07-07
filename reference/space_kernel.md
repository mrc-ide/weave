# Build the spatial correlation matrix

Builds a spatial correlation matrix using an RBF kernel with a nugget
term for numerical stability. Distances are Euclidean in the coordinate
units (see
[`get_spatial_distance()`](https://mrc-ide.github.io/weave/reference/get_spatial_distance.md)),
so `length_scale` is in those same units.

## Usage

``` r
space_kernel(coordinates, length_scale, nugget = 1e-09)
```

## Arguments

- coordinates:

  A data frame with columns `lon` and `lat`.

- length_scale:

  A positive numeric scalar for the spatial length-scale.

- nugget:

  A non-negative numeric scalar added to the diagonal for numerical
  stability.

## Value

A positive-definite matrix representing spatial correlation.
