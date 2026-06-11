# Preconditioned Conjugate Gradient (PCG) solver for the observed system

In plain terms: solves the big linear system that gives the GP weights
using only matrix–vector products—no huge matrices, no explicit inverse.

## Usage

``` r
pcg(
  b,
  obs_idx,
  N,
  space_mat,
  time_mat,
  noise_var,
  kdiag_full,
  tol = 1e-08,
  maxit = 10000
)
```

## Arguments

- b:

  Right-hand side vector (observed length \\m\\).

- obs_idx:

  Integer indices of observed entries in the full vector.

- N:

  Total length of the full vector.

- space_mat:

  Spatial kernel matrix.

- time_mat:

  Temporal kernel matrix.

- noise_var:

  Scalar or length-\\m\\ nugget on the observed scale.

- kdiag_full:

  Vector \\\mathrm{diag}(K)\\ of length \\N\\.

- tol:

  Relative residual tolerance for convergence (default \`1e-8\`).

- maxit:

  Maximum number of iterations (default \`10000\`).

## Value

Numeric solution vector \`x\` of length \\m\\.

## Details

Technically: solves \\(S K S^\top + \mathrm{diag}(\text{noise}))\\x =
b\\ by PCG, using \`Amv\` for matrix–vector products and a Jacobi
(diagonal) preconditioner \\M^{-1} v \approx v / \mathrm{diag}(A)\\,
gathered once up front. Stops when the relative residual falls below
\`tol\` or after \`maxit\` iterations (issues a warning on \`maxit\`).
