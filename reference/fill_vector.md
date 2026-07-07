# Fill observed values into a full vector

Put the observed entries back into their full-length vector, filling
missing positions with zeros.

## Usage

``` r
fill_vector(x_obs, obs_idx, N)
```

## Arguments

- x_obs:

  Numeric vector of observed values (length \\m\\).

- obs_idx:

  Integer indices (length \\m\\) of observed positions in the full
  vector.

- N:

  Total length of the full vector.

## Value

A numeric vector of length `N` with `x_obs` scattered at `obs_idx`.
