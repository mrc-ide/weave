# Observed-system matvec: (S K S^T + diag(noise)) v

Takes an observed-length vector and applies the GP covariance plus a
per-observation noise (nugget), all without building any big matrices.

## Usage

``` r
Amv(v, obs_idx, N, space_mat, time_mat, noise_var)
```

## Arguments

- v:

  Numeric vector of length \\m\\ (observed entries).

- obs_idx:

  Integer indices of observed entries in the full vector.

- N:

  Total length of the full vector.

- space_mat:

  Spatial kernel matrix.

- time_mat:

  Temporal kernel matrix.

- noise_var:

  Scalar or length-\\m\\ numeric nugget on the observed scale.

## Value

A numeric vector of length \\m\\, equal to \\(S K S^\top + D)v\\.

## Details

Technically: for \\K = \mathrm{space}\\\otimes\\\mathrm{time}\\, returns
\$\$S\\K\\S^{\mathsf T}\\v \\+\\ \operatorname{diag}(\sigma^2)\\v,\$\$
i.e., the observed block of the GP plus a diagonal nugget. Implemented
matrix-free as
`kron_mv(with_nas(v, obs_idx, N), space_mat, time_mat)[obs_idx] + noise_var * v`,
where \\S^{\mathsf T}\\ “scatters’’ into the full vector and
\\\sigma^2\\ denotes the per-observation noise.
