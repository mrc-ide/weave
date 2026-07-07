# Package index

## Core workflow

Prepare the data, estimate the kernel hyperparameters, predict.

- [`data_process()`](https://mrc-ide.github.io/weave/reference/data_process.md)
  : Process raw epidemiological data for the GP model
- [`infer_kernel_params()`](https://mrc-ide.github.io/weave/reference/infer_kernel_params.md)
  : Quick exact-marginal-likelihood estimate of the kernel
  hyperparameters
- [`gp_predict()`](https://mrc-ide.github.io/weave/reference/gp_predict.md)
  : Predict the latent rate and a count prediction interval (CG)

## Kernels and distances

- [`space_kernel()`](https://mrc-ide.github.io/weave/reference/space_kernel.md)
  : Build the spatial correlation matrix
- [`time_kernel()`](https://mrc-ide.github.io/weave/reference/time_kernel.md)
  : Build the temporal correlation matrix
- [`rbf_kernel()`](https://mrc-ide.github.io/weave/reference/rbf_kernel.md)
  : Radial basis function kernel
- [`periodic_kernel()`](https://mrc-ide.github.io/weave/reference/periodic_kernel.md)
  : Periodic kernel
- [`get_spatial_distance()`](https://mrc-ide.github.io/weave/reference/get_spatial_distance.md)
  : Pairwise spatial distances
- [`get_temporal_distance()`](https://mrc-ide.github.io/weave/reference/get_temporal_distance.md)
  : Pairwise temporal distances

## Building blocks

The pieces the workflow is assembled from, exported for direct use.

- [`build_plugin_field()`](https://mrc-ide.github.io/weave/reference/build_plugin_field.md)
  : Build a plug-in latent field from observed counts
- [`gp_marginal_loglik()`](https://mrc-ide.github.io/weave/reference/gp_marginal_loglik.md)
  : Exact separable-GP marginal log-likelihood of a field, with a nugget
- [`default_kernel_priors()`](https://mrc-ide.github.io/weave/reference/default_kernel_priors.md)
  : Default priors for the kernel hyperparameters
- [`quick_mvnorm()`](https://mrc-ide.github.io/weave/reference/quick_mvnorm.md)
  : Quick multivariate normal draw over two dimensions
- [`quick_mvnorm_chol()`](https://mrc-ide.github.io/weave/reference/quick_mvnorm_chol.md)
  : Quick multivariate normal draw over two dimensions (Cholesky
  precomputed)
