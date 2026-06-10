# Exact separable-GP marginal log-likelihood of a field, with a nugget

Evaluates the exact log-density of \`g ~ N(0, sigma^2 ((R_s (x) R_t) +
eta I))\` given the eigendecompositions of the spatial and temporal
correlation kernels, via the Kronecker log-determinant and
quadratic-form identities. The global variance \`sigma^2\` is profiled
out (concentrated log-likelihood) and the profiled value is attached as
\`attr(., "sigma2")\`.

## Usage

``` r
gp_marginal_loglik(g, n, nt, eig_s, eig_t, eta)
```

## Arguments

- g:

  Plug-in field, length \`n \* nt\`, ordered sites x times (time
  fastest).

- n:

  Number of sites.

- nt:

  Number of time points.

- eig_s:

  Eigendecomposition (\`eigen\` object) of the spatial correlation
  kernel.

- eig_t:

  Eigendecomposition (\`eigen\` object) of the temporal correlation
  kernel.

- eta:

  Noise-to-signal ratio (nugget), a non-negative scalar.

## Value

The concentrated log-likelihood (numeric scalar), with the profiled
\`sigma2\` attached as an attribute.
