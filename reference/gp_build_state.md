# Build state object

This convenience wrapper assembles the ingredients the Gaussian process
(GP) model needs so downstream functions can work with a single, tidy
bundle.

## Usage

``` r
gp_build_state(obs_data, coordinates, hyperparameters, n, nt, period = 52)
```

## Arguments

- obs_data:

  Observation data

- coordinates:

  Site coordinates

- hyperparameters:

  Vector of hyperparameters

- n:

  N sites

- nt:

  N times

- period:

  Period

## Value

List of pre-computaed GP elements

## Details

It builds the spatial and temporal kernels, indexes the observed counts,
prepares working responses and variances for the Poisson log-scale
approximation, and returns the matrices and closures required for GP
inference.
