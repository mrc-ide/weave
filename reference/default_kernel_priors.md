# Default priors for the kernel hyperparameters

Weakly-informative log-normal priors (i.e. Normal priors on the log
scale) on \`length_scale\`, \`periodic_scale\`, \`long_term_scale\` and
the noise-to-signal ratio \`nugget_ratio\`. They act as mild
regularisation on an otherwise maximum-likelihood fit, keeping
weakly-identified parameters (notably \`long_term_scale\`) away from the
boundary.

## Usage

``` r
default_kernel_priors()
```

## Value

A named list of \`list(meanlog, sdlog)\` priors.
