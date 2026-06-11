
<!-- README.md is generated from README.Rmd. Please edit that file -->

# weave <img src="man/figures/Weave.png" align="right" width=30% height=30% alt="weave package logo" />

<!-- badges: start -->

[![R-CMD-check](https://github.com/mrc-ide/weave/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/mrc-ide/weave/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/mrc-ide/weave/graph/badge.svg)](https://app.codecov.io/gh/mrc-ide/weave)
<!-- badges: end -->

`weave` turns noisy, gappy spatio-temporal count data — such as malaria
cases reported by health facilities — into smooth estimates of the
underlying rate. It fills in missing weeks and attaches an honest
prediction interval to every estimate.

## How it works

`weave` treats the counts as arising from a smooth latent surface that
varies across space and time, modelled as a separable Gaussian process —
a spatial kernel combined with a temporal kernel (periodic season ×
long-term drift). In three steps it:

1.  **builds a quick plug-in field** from the counts (per-site
    standardised `log(1 + y)`);
2.  **estimates the kernel hyperparameters** — how fast counts
    decorrelate across space, within the season, and over the long term
    — by maximising the Gaussian-process marginal likelihood, kept fast
    by the kernel’s Kronecker structure; and
3.  **predicts the latent rate** at every site and week, conditioning on
    the observed cells (so gaps are genuinely interpolated rather than
    guessed) and combining rate uncertainty with Negative-Binomial count
    noise to form a 95% prediction interval.

## Installation

You can install the development version of weave from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("mrc-ide/weave")
```

## Getting started

`weave` needs two inputs: `obs_data`, a data frame with columns `id`
(site), `t` (week) and `y_obs` (the count, `NA` where missing); and
`coordinates`, giving each site’s `lon`/`lat`. The workflow is two calls
— estimate the kernels, then predict:

``` r
library(weave)

# 1. estimate the spatial + temporal kernel length-scales
est <- infer_kernel_params(obs_data, coordinates, nt = nt, period = 52)

# 2. predict the rate, fill the gaps, and attach a 95% prediction interval
pred <- gp_predict(obs_data, coordinates, hyperparameters = est,
                   nt = nt, period = 52, n_draws = 100)
```

`pred` has one row per site-week, with the posterior `rate` and a
`lower`/`upper` prediction interval.

## Learn more

- [Gaussian processes: a gentle
  introduction](https://mrc-ide.github.io/weave/articles/gaussian-processes.html)
  — what a Gaussian process is, the role of the kernel, and conditioning
  on data.
- [The weave
  walkthrough](https://mrc-ide.github.io/weave/articles/walkthrough.html)
  — the full method worked end to end, from estimating the kernels to
  predicting rates with honest uncertainty.
