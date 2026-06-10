# Quick multivariate normal samples over two dimensions (cholesky precomputed)

This is equivalent to estimating the full spatio-temporal covariance
matrix and sampling from the multivariate normal distribution: full_k
\<- kronecker(dist_k, time_k) f \<- mvrnorm(1, rep(0, n \* nt), full_k)

## Usage

``` r
quick_mvnorm_chol(space_chol, time_chol)
```

## Arguments

- space_chol:

  Cholesky decomposition of sapace kernel matrix

- time_chol:

  Cholesky decomposition of time kernel matrix
