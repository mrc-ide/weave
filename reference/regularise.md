# Add a small ridge to a square matrix

In plain terms: this adds a tiny value to the diagonal so the matrix is
better-behaved numerically (e.g., invertible and Cholesky-able).

## Usage

``` r
regularise(x, lambda = 1e-05)
```

## Arguments

- x:

  A square numeric matrix.

- lambda:

  Non-negative ridge (diagonal) value to add. Default \`1e-5\`.

## Value

A matrix the same size as \`x\` with \`lambda\` added to the diagonal.

## Details

Technically: returns \\X + \lambda I\\, which improves condition number
and ensures positive definiteness when \\\lambda \> 0\\.
