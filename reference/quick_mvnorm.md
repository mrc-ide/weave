# Quick multivariate normal samples over two dimensions

This is equivalent to estimating the full spatio-temporal covariance
matrix and sampling from the multivariate normal distribution: full_k
\<- kronecker(dist_k, time_k) f \<- mvrnorm(1, rep(0, n \* nt), full_k)

## Usage

``` r
quick_mvnorm(space, time)
```

## Arguments

- space:

  Space kernel matrix

- time:

  Time kernel matrix
