# Periodic kernel

Computes a periodic kernel for a distance vector or matrix: \$\$k(d) =
\exp\left(-\frac{2\sin^2(\pi d / p)}{\alpha^2}\right).\$\$

## Usage

``` r
periodic_kernel(x, alpha, period)
```

## Arguments

- x:

  A numeric vector or matrix of distances.

- alpha:

  A positive numeric scalar controlling how sharply correlation falls
  within each cycle: smaller values allow sharp seasonal peaks; larger
  values give a gentler, smoother cycle.

- period:

  A positive numeric scalar giving the period \\p\\.

## Value

A numeric vector or matrix of periodic kernel values.
