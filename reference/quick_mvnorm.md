# Quick multivariate normal draw over two dimensions

Draws one sample from a zero-mean Gaussian with separable space-time
covariance, without ever forming the full matrix. This is equivalent to
forming the full spatio-temporal covariance and drawing from the
multivariate normal distribution:

## Usage

``` r
quick_mvnorm(space, time)
```

## Arguments

- space:

  Space kernel matrix.

- time:

  Time kernel matrix.

## Value

A numeric vector of length `nrow(space) * nrow(time)`, ordered site by
site with time varying fastest (matching `kronecker(space, time)`).

## Details

    full_k <- kronecker(space, time)
    f <- MASS::mvrnorm(1, rep(0, nrow(space) * nrow(time)), full_k)
