# Conjugate Gradient (CG) solver for the observed system

In plain terms: solves the big linear system that gives the GP weights
using only matrix-vector products – no huge matrices, no explicit
inverse.

## Usage

``` r
cg(b, obs_idx, N, space_mat, time_mat, noise_var, tol = 1e-08, maxit = 10000)
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

- tol:

  Relative residual tolerance for convergence (default `1e-8`).

- maxit:

  Maximum number of iterations (default `10000`).

## Value

Numeric solution vector `x` of length \\m\\.

## Details

Technically: solves \\(S K S^\top + \mathrm{diag}(\text{noise}))\\x =
b\\ by plain CG, using `Amv` for matrix-vector products. Stops when the
relative residual falls below `tol` or after `maxit` iterations (issues
a warning on `maxit`).

Deliberately unpreconditioned: the separable kernel is built from
*correlation* matrices (unit diagonal plus a constant nugget) and the
noise is a scalar, so \\\mathrm{diag}(A)\\ is exactly constant. A Jacobi
(diagonal) preconditioner therefore only rescales the residual and
leaves the CG iterates unchanged – it is an exact no-op here, so don't
add one.
