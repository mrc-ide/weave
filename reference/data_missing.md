# Drop sites with missing data

Removes sites that cannot be modelled and reports them. Sites with no
observed counts (all \`n\` missing) or without coordinates (all
\`lat\`/\`lon\` missing) are always dropped. Sites whose observed counts
are all zero are dropped only when \`drop_zero = TRUE\`.

Under the \`log1p\` per-site-centred rate model an all-zero site is
valid low-rate data, not a defect, so it is retained by default; set
\`drop_zero = TRUE\` if all-zero facilities should be treated as
reporting artefacts and removed.

## Usage

``` r
data_missing(data, ..., drop_zero = FALSE)
```

## Arguments

- data:

  A data frame containing site identifiers, time \`t\`, counts \`n\`,
  and coordinates \`lat\` and \`lon\`.

- ...:

  Columns identifying sites passed to \[dplyr::group_by()\] (unquoted).

- drop_zero:

  Logical; also drop sites whose observed counts sum to zero (default
  \`FALSE\`).

## Value

A data frame with problem sites removed.
