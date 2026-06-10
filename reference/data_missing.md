# Drop sites with missing data

Removes sites lacking counts or coordinates and reports them.

## Usage

``` r
data_missing(data, ...)
```

## Arguments

- data:

  A data frame containing site identifiers, time \`t\`, counts \`n\`,
  and coordinates \`lat\` and \`lon\`.

- ...:

  Columns identifying sites passed to \[dplyr::group_by()\] (unquoted).

## Value

A data frame with problem sites removed.
