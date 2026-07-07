# Predict the latent rate and a count prediction interval (CG)

Given kernel hyperparameters (e.g. from
[`infer_kernel_params()`](https://mrc-ide.github.io/weave/reference/infer_kernel_params.md)),
predicts the latent rate \\\lambda = e^{\mu_s + f\_{st}}\\ at every
site-by-time cell by conditioning a separable Gaussian process on the
observed counts only. The posterior mean is obtained from a single
matrix-free CG solve (so it is smooth and deterministic). The posterior
variance splits into an exact closed-form "no gaps" part (via the
Kronecker eigendecomposition) plus a missing-data correction estimated
from `n_draws` paired perturbation draws (a control variate; see
[`gp_posterior_var()`](https://mrc-ide.github.io/weave/reference/gp_posterior_var.md)),
so modest draw counts give tight intervals. The latent-rate posterior is
then combined with Negative-Binomial observation noise (law of total
variance, lognormal moment-match) to give a 95% count prediction
interval.

## Usage

``` r
gp_predict(
  obs_data,
  coordinates,
  hyperparameters,
  nt,
  period,
  n_draws = 100,
  r = NULL,
  value = "y_obs",
  standardise = TRUE,
  cg_tol = 1e-06,
  cg_draw_tol = 0.001,
  progress = TRUE
)
```

## Arguments

- obs_data:

  Data frame with `id` (site), `t` (time) and a count column named by
  `value` (`NA` where missing). `t` is a numeric time index whose
  *differences* encode real elapsed time, so gaps and uneven spacing
  between time points are modelled as genuine time distances (use e.g.
  weeks or days since a reference). Must use the same `t` encoding as
  [`infer_kernel_params()`](https://mrc-ide.github.io/weave/reference/infer_kernel_params.md).

- coordinates:

  Site coordinates (data frame with `id`, `lon`, `lat`).

- hyperparameters:

  A list with elements `length_scale`, `periodic_scale`,
  `long_term_scale`, `nugget_ratio` and `sigma2` – the value returned by
  [`infer_kernel_params()`](https://mrc-ide.github.io/weave/reference/infer_kernel_params.md).

- nt:

  Number of time points.

- period:

  Period of the seasonal cycle, in the same units as `t`.

- n_draws:

  Number of paired posterior draws used to estimate the missing-data
  correction to the variance (the prediction interval). Controls only
  the interval, not the mean. Use `0` to return the smooth
  posterior-mean rate only (one solve, no interval). Because most of the
  variance is computed exactly and the draws only estimate the gap
  correction, modest values (25–100) already give tight intervals.

- r:

  Negative-Binomial dispersion for the count interval. If `NULL`
  (default) it is estimated by method of moments from the observed
  counts. `Inf` is valid and means Poisson observation noise (no
  overdispersion); the method-of-moments estimate returns `Inf` itself
  when the data show no excess variance, so `attr(., "r")` on the result
  can be `Inf`. The distribution affects only the interval width – the
  `rate` column does not depend on it.

- value:

  Name of the count column (default `"y_obs"`).

- standardise:

  Logical; standardise the plug-in field per site (default `TRUE`),
  matching
  [`infer_kernel_params()`](https://mrc-ide.github.io/weave/reference/infer_kernel_params.md).

- cg_tol:

  Convergence tolerance for the single posterior-mean CG solve (the
  deterministic part of the prediction).

- cg_draw_tol:

  Convergence tolerance for the `n_draws` perturbation-draw CG solves.
  Deliberately looser than `cg_tol`: the draws only feed a Monte-Carlo
  variance whose own relative error is \\\approx
  1/\sqrt{2\\(n\_{draws}-1)}\\ (about 5% at 200 draws), so solver error
  below that is wasted work. At the default `1e-3` the posterior sd
  typically changes by well under 1% relative to a tight solve, while
  the draw loop needs roughly half the CG iterations.

- progress:

  Logical; show a progress bar over the posterior-draw loop (the
  expensive part). Defaults to `TRUE`, but the bar is drawn only in an
  interactive UTF-8 / truecolor terminal – it stays silent in scripts,
  knitr, logs, and CI. Under a multi-worker
  [`future::plan()`](https://future.futureverse.org/reference/plan.html)
  the bar is suppressed (workers cannot tick it). Set `FALSE` to disable
  it entirely.

## Value

A data frame with one row per site-by-time cell (site-week) and columns
`id`, `t`, `rate` (posterior point estimate of \\\lambda\\), and – when
`n_draws >= 1` – `lower` and `upper` (the 95% count prediction
interval). The dispersion `r` used and `n_draws` are attached as
attributes.

## Details

Because the fit conditions on the observed set, missing cells are filled
by genuine GP interpolation and their prediction interval can widen over
gaps – unlike a completed-grid smoother that mean-imputes the gaps.

Predictions are made at every `(id, t)` cell present in `obs_data`. To
predict at time points with no data anywhere (e.g. a future week),
append rows with `NA` counts for those times (every site) and increase
`nt` accordingly; the interval widens with distance from the data. Treat
extrapolation beyond the observed range with the usual caution.

## Parallel execution

The `n_draws` perturbation draws are independent CG solves and run
through the future framework. By default (no
[`future::plan()`](https://future.futureverse.org/reference/plan.html)
set) they run serially, exactly as before. To spread them across CPU
cores, set a plan before calling and reset it after:

    future::plan(future::multisession, workers = 4)
    pred <- gp_predict(...)
    future::plan(future::sequential)

Results are identical for every backend and worker count, and
reproducible under [`set.seed()`](https://rdrr.io/r/base/Random.html)
(each draw gets its own pre-generated L'Ecuyer-CMRG stream). Parallelism
pays off when `n * n_draws` is large: each worker costs about a second
to start and receives the kernel matrices once. For very large site
counts you may need to raise `options(future.globals.maxSize = ...)`
(the matrices shipped to workers are ~40 MB at 1000 sites x 260 weeks;
the default cap is 500 MiB). If R uses a multithreaded BLAS (e.g.
OpenBLAS/MKL), cap its threads inside workers to avoid oversubscription;
R's shipped BLAS is single-threaded, so by default there is nothing to
do. See
[`vignette("parallel", package = "weave")`](https://mrc-ide.github.io/weave/articles/parallel.md)
for a walkthrough.

## See also

[`infer_kernel_params()`](https://mrc-ide.github.io/weave/reference/infer_kernel_params.md)
