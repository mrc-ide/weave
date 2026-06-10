# One posterior draw of the intensity surface

Generates a single posterior sample of \\\lambda = \exp(z)\\ under the
Gaussian working model on the log scale with heteroscedastic variance
\\D\\. It draws a GP prior realisation \\\eta \sim \mathcal N(0, K)\\,
adds working-likelihood noise \\\varepsilon \sim \mathcal N(0, D)\\ on
the observed log scale, solves \\(S K S^\top + D)\alpha\_{\sim} =
(S\eta + \varepsilon) - y\_{\text{work}}\\, and returns \\\exp(\eta - K
S^\top \alpha\_{\sim} + \mu\_{\mathrm{infer}})\\.

## Usage

``` r
gp_draw(state, tol = 1e-06)
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

A numeric vector of length `state$N` containing one posterior draw of
\\\lambda\\ in the same ordering as `state$mu_infer`.
