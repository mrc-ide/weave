# Posterior mean of the intensity surface

Computes the posterior mean of \\\lambda = \exp(z)\\ where \\z = f +
\mu\_{\mathrm{infer}}\\. The latent field \\f\\ follows a zero–mean GP
with separable covariance \\K = K\_{\mathrm{space}} \otimes
K\_{\mathrm{time}}\\. A Gaussian working likelihood on the log scale
with heteroscedastic diagonal variance \\D =
\mathrm{diag}(\text{noise\\var})\\ is used. This evaluates
\$\$f\_{\text{hat}} = K S^\top (S K S^\top + D)^{-1}
y\_{\text{work}},\$\$ then returns \\\exp(f\_{\text{hat}} +
\mu\_{\mathrm{infer}})\\.

## Usage

``` r
gp_posterior_mean(state, tol = 1e-06)
```

## Arguments

- state:

  A sampler state created by
  [`gp_build_state()`](https://mrc-ide.github.io/weave/reference/gp_build_state.md),
  containing at least `space_mat`, `time_mat`, `obs_idx`, `N`, `y_work`,
  `noise_var`, `kdiag_full`, `A_solve`, and `mu_infer`. The vector
  layout is sites × times with time varying fastest.

- tol:

  Convergence tolerance passed to the inner PCG solve.

## Value

A numeric vector of length `state$N` giving the posterior mean intensity
\\\lambda\\ in the same ordering as `state$mu_infer`.
