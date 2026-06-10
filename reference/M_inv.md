# Diagonal (Jacobi) preconditioner application

Divides by an approximation to the diagonal of the system, which makes
the iterative solver converge faster.

## Usage

``` r
M_inv(v, kdiag_full, obs_idx, noise_var)
```

## Arguments

- v:

  Numeric vector to precondition (length \\m\\).

- kdiag_full:

  Vector \\\mathrm{diag}(K)\\ of length \\N\\ (typically from
  \`as.vector(kronecker(diag(space), diag(time)))\`).

- obs_idx:

  Integer indices of observed entries in the full vector.

- noise_var:

  Scalar or length-\\m\\ numeric nugget to add to the diagonal.

## Value

A numeric vector of length \\m\\: elementwise \`v / (diagA + 1e-12)\`.

## Details

Technically: applies \\M^{-1} v \approx v / \mathrm{diag}(A)\\, where
\\A = S K S^\top + \mathrm{diag}(\text{noise})\\ and \\\mathrm{diag}(K)
= \mathrm{diag}(space) \otimes \mathrm{diag}(time)\\.
