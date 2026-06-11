# Kronecker diagonal of a separable kernel

A helper: the diagonal of \`space ⊗ time\` without allocating the big
dense \`kronecker()\` product.

## Usage

``` r
kdiag_from_factors(space_diag, time_diag, n, nt)
```

## Arguments

- space_diag:

  Space matrix diagonal

- time_diag:

  Time matrix diagonal

- n:

  Number of sites

- nt:

  Number of times

## Value

Kronecker diagonal
