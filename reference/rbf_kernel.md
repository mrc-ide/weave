# Radial basis function kernel

Computes the radial basis function (RBF) kernel for a distance vector or
matrix: \$\$k(d) = \exp\left(-\frac{d^2}{2\theta^2}\right).\$\$ This is
the correlation form of the kernel (\\k(0) = 1\\); any global variance
is applied separately.

## Usage

``` r
rbf_kernel(x, theta)
```

## Arguments

- x:

  A numeric vector or matrix of distances.

- theta:

  A positive numeric scalar giving the length-scale parameter (the
  \\\ell\\ of textbook presentations): correlation decays with distance
  over this scale, so larger values give smoother functions.

## Value

A numeric vector or matrix with RBF kernel values.
