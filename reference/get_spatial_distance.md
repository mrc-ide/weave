# Pairwise spatial distances

Computes pairwise Euclidean distances between locations, in the units of
the coordinates. No great-circle correction is applied: with raw
longitude/latitude degrees, one degree of longitude shrinks with
latitude, so for large or high-latitude extents project the coordinates
first (e.g. to km) and interpret `length_scale` in those units.

## Usage

``` r
get_spatial_distance(coordinates)
```

## Arguments

- coordinates:

  A data frame with columns `lon` and `lat`.

## Value

A symmetric matrix of pairwise spatial distances.
