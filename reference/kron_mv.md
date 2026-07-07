# Fast Kronecker–product matrix–vector multiply (times vary fastest)

In plain terms: multiplies a big covariance \`K = space ⊗ time\` by a
vector without ever forming \`K\`, using a reshape–multiply–reshape
trick.

## Usage

``` r
kron_mv(v, space, time)
```

## Arguments

- v:

  Numeric vector of length \`nrow(space) \* nrow(time)\`, ordered with
  times varying fastest within site.

- space:

  Spatial kernel matrix (size \\n \times n\\). Must be symmetric
  (kernels are, by construction) – the implementation relies on
  \`t(space) == space\`.

- time:

  Temporal kernel matrix (size \\nt \times nt\\).

## Value

A numeric vector the same length as \`v\`.

## Details

Technically: for \\v = \mathrm{vec}(X^\top)\\ with \`times\` varying
fastest, computes \\(space \otimes time)\\v =
\mathrm{vec}\\\big((space\\X\\time^\top)^\top\big)\\.
