# =============================================================================
# PFN torch model
#
# Architecture (CPU-friendly, fixed-geometry):
#
#   y, mask: (B, n, nt)
#       |
#       v
#   [Temporal encoder]      -- 1D conv over time, weight-shared across sites
#       |                      Input channels: [log1p(y)*mask, mask]
#       v
#   per-cell features (B, n, d, nt)
#   per-site features (B, n, d)       (= mean over time of per-cell)
#       |
#       v
#   [Spatial block]         -- single multi-head attention layer with
#       |                      additive bias derived from the (fixed)
#       v                      spatial distance matrix.
#   contextualised site embeddings (B, n, d)
#       |
#       +-> [theta+r head] (mean over sites -> Gaussian over log-params, R^4)
#       +-> [mu_s head]    (per site -> Gaussian, R^n)
#       +-> [f head]       (per cell -> Gaussian, R^{n*nt})
#
# Two structural choices vs a vanilla transformer:
#
#   * Tokens are SITES (n ~ 100-200), not cells (N ~ 30k). Attention is
#     O(n^2 d), which is what makes CPU training tractable.
#
#   * The fixed spatial geometry is baked in as an attention bias buffer
#     (-d_ij^2 normalised so mean off-diagonal is 1), scaled by a learnable
#     per-head scalar. This replaces positional encoding and gives the model
#     a clean prior on which sites should attend to which (cf. DeepRV
#     §4.2.3 kernel-attention bias).
#
# Outputs are diagonal Gaussians for simplicity. If calibration ends up
# weak we'll swap the theta+r head for a normalizing flow without touching
# the rest of the pipeline.
# =============================================================================


# Internal: shared 1D conv stack over time, applied site-wise.
.nn_pfn_temporal_encoder <- function(d_out, n_in_channels = 2,
                                     hidden_channels = 32,
                                     kernel_size = 7) {
  # Lazy nn_module() construction; caller wraps with torch::nn_module().
  list(
    n_in_channels   = n_in_channels,
    hidden_channels = hidden_channels,
    d_out           = d_out,
    kernel_size     = kernel_size
  )
}


#' Build the PFN torch model for a fixed facility geometry
#'
#' Returns an `nn_module` whose `$forward(y, mask)` consumes a batched
#' `(B, n, nt)` count tensor and mask, and returns a list of posterior
#' parameters:
#'
#'   - `theta_mean`, `theta_logvar` -- (B, 4) Gaussian over
#'     `log(length_scale, periodic_scale, long_term_scale, r)`
#'   - `mu_mean`, `mu_logvar`       -- (B, n)
#'   - `f_mean`, `f_logvar`         -- (B, n, nt)
#'
#' All Gaussian, all diagonal. Replace the head later for richer families
#' without touching the encoder/spatial block.
#'
#' @param coords Data frame with `lat`, `lon`; defines the fixed geometry.
#' @param nt Number of timepoints.
#' @param period Periodic kernel period (default 52). Used to inform the
#'   temporal kernel size; not strictly required by the architecture but
#'   convenient to keep on the module for inference.
#' @param d_model Embedding dimension (default 64).
#' @param n_heads Attention heads (default 4); `d_model` must be divisible
#'   by `n_heads`.
#' @param n_spatial_blocks Number of stacked spatial attention blocks
#'   (default 1).
#' @param kernel_size Time-conv kernel size (default 7).
#' @param hidden_channels Conv hidden channel count (default 32).
#'
#' @return An `nn_module` instance. Move to GPU with `model$to(device = "cuda")`.
#' @export
nn_pfn <- function(coords, nt, period = 52,
                   d_model          = 64,
                   n_heads          = 4,
                   n_spatial_blocks = 1L,
                   kernel_size      = 7,
                   hidden_channels  = 32) {
  if (!requireNamespace("torch", quietly = TRUE)) {
    stop("`torch` is required for nn_pfn(). ",
         "Install with install.packages('torch') and torch::install_torch().")
  }
  if (d_model %% n_heads != 0L) {
    stop("`d_model` (", d_model, ") must be divisible by `n_heads` (",
         n_heads, ").")
  }

  n <- nrow(coords)

  # Spatial distance bias buffer: -d^2, normalised so the off-diagonal mean
  # is 1, then sign-flipped (negative => closer sites get higher attention).
  dist_mat <- get_spatial_distance(coords[, c("lon", "lat")])
  d2       <- dist_mat^2
  off_mean <- mean(d2[upper.tri(d2)])
  if (!is.finite(off_mean) || off_mean <= 0) off_mean <- 1
  d_bias <- -d2 / off_mean

  module <- torch::nn_module(
    classname = "nn_pfn",

    initialize = function() {
      # Stash a few sizes on `self` so we can use them in forward().
      self$n <- n
      self$nt <- nt
      self$period <- period
      self$d_model <- d_model
      self$n_heads <- n_heads
      self$d_head <- d_model %/% n_heads

      # Full construction args captured so checkpoints can round-trip the
      # architecture (see pfn_train()/fit_pfn()). Plain R list, not a buffer
      # -- doesn't move with $to(device), which is fine since these are
      # CPU-side metadata.
      self$spec <- list(
        coords           = coords,
        nt               = nt,
        period           = period,
        d_model          = d_model,
        n_heads          = n_heads,
        n_spatial_blocks = n_spatial_blocks,
        kernel_size      = kernel_size,
        hidden_channels  = hidden_channels
      )

      # -- Temporal encoder ---------------------------------------------------
      pad <- kernel_size %/% 2L
      self$tconv1 <- torch::nn_conv1d(2L, hidden_channels,
                                      kernel_size = kernel_size,
                                      padding = pad)
      self$tconv2 <- torch::nn_conv1d(hidden_channels, d_model,
                                      kernel_size = kernel_size,
                                      padding = pad)
      self$tact   <- torch::nn_gelu()

      # -- Spatial block(s) ---------------------------------------------------
      # Fixed distance bias buffer. Registered properly so it moves with
      # $to(device) and persists in $state_dict() for checkpointing.
      self$register_buffer(
        "d_bias",
        torch::torch_tensor(d_bias, dtype = torch::torch_float())
      )

      # Stack of spatial blocks. Each block consumes the shared d_bias buffer.
      self$spatial_blocks <- torch::nn_module_list(
        lapply(seq_len(n_spatial_blocks), function(.) {
          .pfn_make_spatial_block(d_model, n_heads)
        })
      )

      # -- Heads --------------------------------------------------------------
      # theta+r: pooled site embedding -> 4 means + 4 log-vars
      self$theta_head <- torch::nn_sequential(
        torch::nn_linear(d_model, d_model),
        torch::nn_gelu(),
        torch::nn_linear(d_model, 8L)
      )
      # mu_s: per-site -> 2
      self$mu_head <- torch::nn_linear(d_model, 2L)
      # f: per-cell, combines site context (broadcast over time) with
      # per-cell temporal features (output of conv stack).
      self$f_head <- torch::nn_sequential(
        torch::nn_linear(2L * d_model, d_model),
        torch::nn_gelu(),
        torch::nn_linear(d_model, 2L)
      )
    },

    forward = function(y, mask) {
      # y, mask: (B, n, nt) float tensors. mask is 1 = observed, 0 = missing.

      sizes <- y$size()
      B  <- sizes[1L]; n_ <- sizes[2L]; nt_ <- sizes[3L]

      # ---- Temporal encoder -------------------------------------------------
      # log1p(y) gated by mask, plus mask as a separate channel so the
      # network knows what's missing vs. genuine zeros.
      y_log <- torch::torch_log1p(y) * mask
      x     <- torch::torch_stack(list(y_log, mask), dim = 3L) # (B, n, 2, nt)
      x     <- x$reshape(c(B * n_, 2L, nt_))                   # (B*n, 2, nt)

      x <- self$tact(self$tconv1(x))
      x <- self$tact(self$tconv2(x))                            # (B*n, d, nt)

      per_cell <- x$reshape(c(B, n_, self$d_model, nt_))        # (B,n,d,nt)
      per_site <- per_cell$mean(dim = 4L)                        # (B,n,d)

      # ---- Spatial block(s) -------------------------------------------------
      h <- per_site
      for (i in seq_len(length(self$spatial_blocks))) {
        h <- self$spatial_blocks[[i]](h, self$d_bias)
      }                                                          # (B,n,d)

      # ---- theta+r head -----------------------------------------------------
      pooled       <- h$mean(dim = 2L)                           # (B, d)
      theta_out    <- self$theta_head(pooled)                    # (B, 8)
      theta_mean   <- theta_out[, 1:4]
      theta_logvar <- theta_out[, 5:8]

      # ---- mu_s head --------------------------------------------------------
      mu_out    <- self$mu_head(h)                               # (B, n, 2)
      mu_mean   <- mu_out[, , 1]
      mu_logvar <- mu_out[, , 2]

      # ---- f head -----------------------------------------------------------
      # Broadcast site context (B, n, d) -> (B, n, d, nt), concat with
      # per-cell features (B, n, d, nt), then permute to (B, n, nt, 2d) for
      # the Linear stack.
      h_bcast <- h$unsqueeze(4L)$expand(c(B, n_, self$d_model, nt_))
      f_in    <- torch::torch_cat(list(h_bcast, per_cell), dim = 3L)
      f_in    <- f_in$permute(c(1L, 2L, 4L, 3L))                 # (B,n,nt,2d)
      f_out   <- self$f_head(f_in)                                # (B,n,nt,2)
      f_mean    <- f_out[, , , 1]
      f_logvar  <- f_out[, , , 2]

      list(
        theta_mean   = theta_mean,
        theta_logvar = theta_logvar,
        mu_mean      = mu_mean,
        mu_logvar    = mu_logvar,
        f_mean       = f_mean,
        f_logvar     = f_logvar
      )
    }
  )

  module()
}


# Internal: build a single spatial attention block with externally-supplied
# distance-bias buffer. Pre-norm transformer block layout (LN -> attn -> +
# -> LN -> ffn -> +).
.pfn_make_spatial_block <- function(d_model, n_heads) {
  d_head <- d_model %/% n_heads

  torch::nn_module(
    classname = "pfn_spatial_block",
    initialize = function() {
      self$d_model <- d_model
      self$n_heads <- n_heads
      self$d_head  <- d_head

      self$qkv      <- torch::nn_linear(d_model, 3L * d_model)
      self$out_proj <- torch::nn_linear(d_model, d_model)

      # Per-head amplitude on the spatial bias. Initialise small so the
      # network starts close to plain self-attention and learns to use the
      # geometric prior as needed.
      self$bias_scale <- torch::nn_parameter(
        torch::torch_full(c(n_heads), 0.1)
      )

      self$ln1 <- torch::nn_layer_norm(d_model)
      self$ln2 <- torch::nn_layer_norm(d_model)
      self$ffn <- torch::nn_sequential(
        torch::nn_linear(d_model, 4L * d_model),
        torch::nn_gelu(),
        torch::nn_linear(4L * d_model, d_model)
      )
    },
    forward = function(x, d_bias) {
      # x: (B, n, d_model); d_bias: (n, n) on same device.
      sizes <- x$size()
      B <- sizes[1L]; n_ <- sizes[2L]; d <- sizes[3L]

      h <- self$ln1(x)
      qkv <- self$qkv(h)$reshape(c(B, n_, 3L, self$n_heads, self$d_head))
      qkv <- qkv$permute(c(3L, 1L, 4L, 2L, 5L))   # (3, B, H, n, d_head)
      q <- qkv[1]; k <- qkv[2]; v <- qkv[3]

      scale  <- 1 / sqrt(self$d_head)
      scores <- torch::torch_matmul(q, k$transpose(3L, 4L)) * scale
      # scores: (B, H, n, n)

      # Spatial bias: (H, n, n) = bias_scale[h] * d_bias  ->  broadcast over B.
      bias <- self$bias_scale$reshape(c(self$n_heads, 1L, 1L)) *
              d_bias$unsqueeze(1L)
      scores <- scores + bias$unsqueeze(1L)

      attn <- torch::nnf_softmax(scores, dim = 4L)
      out  <- torch::torch_matmul(attn, v)             # (B, H, n, d_head)
      out  <- out$permute(c(1L, 3L, 2L, 4L))$reshape(c(B, n_, d))
      out  <- self$out_proj(out)

      x <- x + out
      x <- x + self$ffn(self$ln2(x))
      x
    }
  )()
}


#' Count trainable parameters in a PFN module
#'
#' Sanity-check helper for the CPU budget. Use during architecture changes
#' to ensure the parameter count stays in the target range.
#'
#' @param model An `nn_module` returned by [nn_pfn()].
#' @return Integer scalar.
#' @export
pfn_n_params <- function(model) {
  sum(vapply(model$parameters, function(p) p$numel(), integer(1)))
}
