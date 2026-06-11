# weave ![weave package logo](reference/figures/Weave.png)

`weave` interpolates noisy, gappy spatio-temporal count data — such as
malaria cases reported by health facilities — with a separable Gaussian
process. It estimates how quickly counts become uncorrelated across
space and time, then predicts the underlying rate at every site and
week, fills in missing observations, and attaches an honest prediction
interval.

## Installation

You can install the development version of weave from
[GitHub](https://github.com/) with:

``` r

# install.packages("pak")
pak::pak("mrc-ide/weave")
```

## Quick start

Estimate the kernel hyperparameters from observed counts, then predict
the latent rate (with a 95% prediction interval) at every site-week:

``` r

library(weave)

# estimate the spatial + temporal kernel length-scales
est <- infer_kernel_params(obs_data, coordinates, nt = nt, period = 52)

# predict the rate, fill the gaps, and attach a prediction interval
pred <- gp_predict(obs_data, coordinates, hyperparameters = est,
                   nt = nt, period = 52, n_draws = 100)
```

Here `obs_data` has columns `id` (site), `t` (week) and `y_obs` (count,
`NA` where missing), and `coordinates` gives each site’s `lon`/`lat`.

## Learn more

- [Gaussian processes: a gentle
  introduction](https://mrc-ide.github.io/weave/articles/gaussian-processes.html)
  — what a GP is, the role of the kernel, and conditioning on data.
- [The weave
  walkthrough](https://mrc-ide.github.io/weave/articles/walkthrough.html)
  — the full method, from estimating the kernels to predicting rates
  with honest uncertainty.
