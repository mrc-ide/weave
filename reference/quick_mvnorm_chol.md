# Quick multivariate normal draw over two dimensions (Cholesky precomputed)

As
[`quick_mvnorm()`](https://mrc-ide.github.io/weave/reference/quick_mvnorm.md),
but taking precomputed Cholesky factors so repeated draws (e.g. the
perturbation draws in
[`gp_predict()`](https://mrc-ide.github.io/weave/reference/gp_predict.md))
skip the factorisation cost.

## Usage

``` r
quick_mvnorm_chol(space_chol, time_chol)
```

## Arguments

- space_chol:

  Upper-triangular Cholesky factor of the space kernel matrix, as
  returned by [`chol()`](https://rdrr.io/r/base/chol.html). Passing the
  lower-triangular factor gives silently wrong draws.

- time_chol:

  Upper-triangular Cholesky factor of the time kernel matrix, as
  returned by [`chol()`](https://rdrr.io/r/base/chol.html).

## Value

A numeric vector of length `nrow(space_chol) * nrow(time_chol)`, ordered
site by site with time varying fastest (matching
`kronecker(space, time)`).
