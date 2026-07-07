# Exact separable-GP marginal log-likelihood of a field, with a nugget

Evaluates the exact log-density of
`g ~ N(0, sigma^2 ((R_s (x) R_t) + eta I))` given the
eigendecompositions of the spatial and temporal correlation kernels, via
the Kronecker log-determinant and quadratic-form identities. The global
variance `sigma^2` is profiled out (concentrated log-likelihood) and the
profiled value is attached as `attr(., "sigma2")`.

## Usage

``` r
gp_marginal_loglik(g, n, nt, eig_s, eig_t, eta)
```

## Arguments

- g:

  Plug-in field, length `n * nt`, ordered sites x times (time fastest),
  as returned by
  [`build_plugin_field()`](https://mrc-ide.github.io/weave/reference/build_plugin_field.md).

- n:

  Number of sites.

- nt:

  Number of time points.

- eig_s:

  Eigendecomposition (`eigen` object) of the spatial correlation kernel.

- eig_t:

  Eigendecomposition (`eigen` object) of the temporal correlation
  kernel.

- eta:

  Noise-to-signal ratio (nugget), a non-negative scalar.

## Value

The concentrated log-likelihood (numeric scalar), with the profiled
`sigma2` attached as an attribute.

## Details

[`infer_kernel_params()`](https://mrc-ide.github.io/weave/reference/infer_kernel_params.md)
maximises this quantity for you; call it directly when you want to score
a candidate set of hyperparameters yourself (e.g. profiling a likelihood
surface). The kernels should be correlation matrices (unit diagonal) as
built by
[`space_kernel()`](https://mrc-ide.github.io/weave/reference/space_kernel.md)
/
[`time_kernel()`](https://mrc-ide.github.io/weave/reference/time_kernel.md),
decomposed with `eigen(., symmetric = TRUE)`.

## Examples

``` r
n <- 5; nt <- 20
coords <- data.frame(id = 1:n, lon = runif(n), lat = runif(n))
g <- stats::rnorm(n * nt)
eig_s <- eigen(space_kernel(coords, length_scale = 1), symmetric = TRUE)
eig_t <- eigen(time_kernel(1:nt, periodic_scale = 1, long_term_scale = 50,
                           period = 52), symmetric = TRUE)
ll <- gp_marginal_loglik(g, n, nt, eig_s, eig_t, eta = 0.1)
ll
#> [1] -158.7387
#> attr(,"sigma2")
#> [1] 9.878955
attr(ll, "sigma2")  # the profiled global variance
#> [1] 9.878955
```
