# Order data and assign identifiers

Arranges data by site and time and creates a factor \`id\` per site.

## Usage

``` r
data_order_index(data, ...)
```

## Arguments

- data:

  A data frame containing site identifiers and time \`t\`.

- ...:

  Columns identifying sites passed to \[dplyr::group_by()\] (unquoted).

## Value

A data frame ordered by site and time with an \`id\` column.
