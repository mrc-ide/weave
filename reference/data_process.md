# Process raw epidemiological data

Validates and prepares input data for modelling.

## Usage

``` r
data_process(data, ...)
```

## Arguments

- data:

  A data frame containing site identifiers, time \`t\`, counts \`n\`,
  and coordinates \`lat\` and \`lon\`.

- ...:

  Columns identifying sites passed to \[dplyr::group_by()\] (unquoted).

## Value

A processed data frame ready for model fitting.
