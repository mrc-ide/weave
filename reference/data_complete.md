# Complete site-time combinations

Adds rows for all combinations of site identifiers and time `t`.

## Usage

``` r
data_complete(data, ...)
```

## Arguments

- data:

  A data frame containing site identifiers, time `t`, counts `n`, and
  coordinates `lat` and `lon`.

- ...:

  Bare (unquoted) column names that jointly identify a site, e.g.
  `region, facility_name`.

## Value

A data frame with missing site-time combinations filled in and `n` set
to `NA`.
