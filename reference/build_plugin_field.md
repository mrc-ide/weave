# Build a plug-in latent field from observed counts

Forms a cheap estimate of the latent log-intensity field as the per-site
centred (and optionally scaled) \`log1p\` of the observed counts.
Missing cells are mean-imputed (0 after centring) and therefore
contribute nothing to the marginal likelihood.

## Usage

``` r
build_plugin_field(obs_data, n, nt, value = "y_obs", standardise = TRUE)
```

## Arguments

- obs_data:

  Data frame with \`id\` (site), \`t\` (time) and the count column named
  by \`value\`.

- n:

  Number of sites.

- nt:

  Number of time points.

- value:

  Name of the count column (default \`"y_obs"\`).

- standardise:

  Logical; scale each site to unit variance after centring (default
  \`TRUE\`).

## Value

A numeric vector of length \`n \* nt\`, ordered sites x times (time
fastest).

## Details

Per-site centring removes the site intercept \`mu_s\`; per-site scaling
homogenises per-site variances so a single global \`sigma^2\` and the
\*correlation\* kernels apply.
