# Calculate observed summary statistics

Computes log-scale mean and variance summaries for each site.

## Usage

``` r
data_observed_summary(data, ...)
```

## Arguments

- data:

  A data frame containing site identifiers, time \`t\`, counts \`n\`,
  and coordinates \`lat\` and \`lon\`.

- ...:

  Columns identifying sites passed to \[dplyr::group_by()\] (unquoted).

## Value

A data frame with \`observed_mu\` and \`observed_sigmasq\` columns.
