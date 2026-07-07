# Order data and assign identifiers

Arranges data by site and time and creates a factor `id` per site.

## Usage

``` r
data_order_index(data, ...)
```

## Arguments

- data:

  A data frame containing site identifiers and time `t`.

- ...:

  Bare (unquoted) column names that jointly identify a site, e.g.
  `region, facility_name`.

## Value

A data frame ordered by site and time with an `id` column.
