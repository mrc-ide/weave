# =============================================================================
# PFN training loop
#
# Loss:
#
#   L = -log q_phi(log_theta, log_r | y, mask)        Gaussian NLL on R^4
#       -log q_phi(mu_s          | y, mask)            Gaussian NLL on R^n
#       -log q_phi(f             | y, mask)            Gaussian NLL on R^{n*nt}
#
# Per-cell f and per-site mu_s NLLs are MEAN-pooled within each dataset
# (instead of summed) so the three loss terms have comparable magnitudes
# regardless of geometry. This keeps the optimiser well-conditioned across
# different (n, nt).
#
# We drop the log(2*pi) constant -- it doesn't affect the gradient and the
# remaining terms are easier to read in logs.
# =============================================================================


# Internal: Gaussian NLL with a soft floor on the predicted log-variance.
#
# The log_var floor is the standard trick to keep early-training instability
# from blowing up: when the network predicts log_var = -10 and the residual
# is large, the loss explodes and the gradient with it. Clamping log_var at
# the lower end (default -7 => stdev ~ 1e-3.5) caps the worst case without
# distorting the well-conditioned regime.
.gaussian_nll <- function(target, mean, log_var, log_var_floor = -7) {
  log_var <- log_var$clamp(min = log_var_floor)
  0.5 * (log_var + (target - mean)^2 / log_var$exp())
}


#' Compute the PFN loss for one minibatch
#'
#' Public so [fit_pfn()] and SBC code can reuse the exact same loss
#' decomposition during validation. Returns the scalar total loss plus
#' per-component scalar terms for logging.
#'
#' `target_log4` is the (B, 4) concatenation of `log_theta` and `log_r`,
#' which is what the model's `theta_*` head predicts. Pass it through
#' explicitly so callers control the ordering.
#'
#' @param out List of forward outputs from [nn_pfn()].
#' @param target_log4 (B, 4) tensor: `cbind(log_theta, log_r)`.
#' @param target_mu_s (B, n) tensor.
#' @param target_f (B, n, nt) tensor.
#' @param weights Optional named list of scalar weights on each loss
#'   component (`theta`, `mu`, `f`). Defaults to all 1.
#' @return A list with `loss` (scalar tensor with grad) and `parts` (named
#'   numeric vector for logging).
#' @export
pfn_loss <- function(out, target_log4, target_mu_s, target_f,
                     weights = list(theta = 1, mu = 1, f = 1)) {
  l_theta <- .gaussian_nll(target_log4, out$theta_mean,   out$theta_logvar)$mean()
  l_mu    <- .gaussian_nll(target_mu_s, out$mu_mean,      out$mu_logvar)$mean()
  l_f     <- .gaussian_nll(target_f,    out$f_mean,       out$f_logvar)$mean()

  total <- weights$theta * l_theta + weights$mu * l_mu + weights$f * l_f

  list(
    loss  = total,
    parts = c(theta = as.numeric(l_theta),
              mu    = as.numeric(l_mu),
              f     = as.numeric(l_f))
  )
}


#' Train a PFN model on a simulated batch
#'
#' Splits the input batch into a train/validation slice, runs `epochs`
#' passes over the train slice with minibatches of `batch_size`, and
#' tracks per-epoch component losses on both slices. Checkpoints to
#' `checkpoint` after each epoch (atomic write via a sidecar `.tmp` rename
#' so an interrupted run never leaves a half-written file).
#'
#' Resume: if `checkpoint` exists at start and `resume = TRUE`, the model
#' state and optimiser state are loaded from it and training continues
#' from the saved epoch.
#'
#' @param model An `nn_module` from [nn_pfn()].
#' @param batch A list from [pfn_simulate_batch()] (R arrays). Converted to
#'   tensors on `device` once at the start.
#' @param epochs Number of passes over the train slice.
#' @param batch_size Minibatch size.
#' @param lr Initial learning rate.
#' @param val_frac Fraction of `batch` held out for validation (default 0.1).
#' @param grad_clip Max gradient norm before clipping (default 1.0). Pass
#'   `NULL` to disable.
#' @param weight_decay AdamW weight decay (default 1e-4).
#' @param loss_weights Named list of component weights; see [pfn_loss()].
#' @param device Torch device string (`"cpu"` by default; `"cuda"` once GPU
#'   becomes available).
#' @param checkpoint Path to a `.pt` file. If non-NULL, model + optimiser
#'   state are saved after every epoch.
#' @param resume If TRUE and `checkpoint` exists, load it and continue from
#'   the saved epoch.
#' @param verbose Show a per-epoch progress message.
#'
#' @return A list with:
#'   - `history` : data frame of per-epoch loss components (train + val)
#'   - `model`   : the (in-place trained) model
#'   - `epoch`   : final epoch index reached
#' @export
pfn_train <- function(model, batch,
                      epochs       = 2L,
                      batch_size   = 32L,
                      lr           = 3e-3,
                      val_frac     = 0.1,
                      grad_clip    = 1.0,
                      weight_decay = 1e-4,
                      loss_weights = list(theta = 1, mu = 1, f = 1),
                      device       = "cpu",
                      checkpoint   = NULL,
                      resume       = TRUE,
                      verbose      = TRUE) {
  if (!requireNamespace("torch", quietly = TRUE)) {
    stop("`torch` is required for pfn_train().")
  }

  # Move tensors + model to device once.
  model$to(device = device)
  tt <- pfn_to_tensors(batch, device = device)

  B <- tt$B
  n_val <- max(1L, as.integer(round(val_frac * B)))
  n_train <- B - n_val
  if (n_train < batch_size) {
    warning("`batch_size` (", batch_size, ") larger than train slice (",
            n_train, "); using batch_size = ", n_train, ".")
    batch_size <- n_train
  }

  # Deterministic split: last n_val datasets are validation. The simulator
  # output has no ordering structure (each draw is independent), so this is
  # equivalent to a random split without needing seed bookkeeping.
  train_idx <- seq_len(n_train)
  val_idx   <- (n_train + 1L):B

  target_log4_all <- torch::torch_cat(
    list(tt$log_theta, tt$log_r$unsqueeze(2L)), dim = 2L)

  optimizer <- torch::optim_adamw(model$parameters, lr = lr,
                                  weight_decay = weight_decay)

  start_epoch <- 1L
  loaded_ckpt <- NULL
  if (!is.null(checkpoint) && resume && file.exists(checkpoint)) {
    if (verbose) message("[pfn_train] resuming from ", checkpoint)
    loaded_ckpt <- torch::torch_load(checkpoint)
    # torch_save() strips the data.frame class from history on round-trip;
    # restore it so callers can ggplot/print() the result.
    if (!is.data.frame(loaded_ckpt$history)) {
      loaded_ckpt$history <- as.data.frame(loaded_ckpt$history)
    }
    model$load_state_dict(loaded_ckpt$model)
    optimizer$load_state_dict(loaded_ckpt$optim)
    start_epoch <- loaded_ckpt$epoch + 1L
    if (start_epoch > epochs) {
      message("[pfn_train] checkpoint already at epoch ", loaded_ckpt$epoch,
              " >= requested epochs ", epochs, "; nothing to do.")
      return(list(history = loaded_ckpt$history, model = model,
                  epoch = loaded_ckpt$epoch))
    }
  }

  history <- if (!is.null(loaded_ckpt)) {
    loaded_ckpt$history
  } else {
    data.frame(
      epoch       = integer(0),
      train_total = numeric(0), train_theta = numeric(0),
      train_mu    = numeric(0), train_f     = numeric(0),
      val_total   = numeric(0), val_theta   = numeric(0),
      val_mu      = numeric(0), val_f       = numeric(0),
      sec         = numeric(0)
    )
  }

  # Helper: average loss over a slice of indices (no grad).
  eval_slice <- function(idx) {
    model$eval()
    on.exit(model$train())
    n_eval <- length(idx)
    if (n_eval == 0) return(c(total = NA, theta = NA, mu = NA, f = NA))
    parts_acc <- c(theta = 0, mu = 0, f = 0)
    n_chunks  <- 0L
    chunk_sz  <- min(batch_size, n_eval)
    chunks <- split(idx, ceiling(seq_along(idx) / chunk_sz))
    with_no_grad <- torch::with_no_grad
    with_no_grad({
      for (chunk in chunks) {
        b <- list(
          y         = tt$y[chunk, , ],
          mask      = tt$mask[chunk, , ],
          log_theta = tt$log_theta[chunk, ],
          log_r     = tt$log_r[chunk],
          mu_s      = tt$mu_s[chunk, ],
          f         = tt$f[chunk, , ]
        )
        target_log4 <- target_log4_all[chunk, ]
        out <- model(b$y, b$mask)
        lo  <- pfn_loss(out, target_log4, b$mu_s, b$f,
                        weights = loss_weights)
        parts_acc <- parts_acc + lo$parts
        n_chunks  <- n_chunks + 1L
      }
    })
    parts <- parts_acc / n_chunks
    c(total = sum(loss_weights[c("theta","mu","f")] |> unlist() * parts),
      parts)
  }

  for (epoch in seq.int(start_epoch, epochs)) {
    t0 <- Sys.time()
    model$train()

    perm <- sample(train_idx)
    chunks <- split(perm, ceiling(seq_along(perm) / batch_size))

    train_acc <- c(theta = 0, mu = 0, f = 0)
    n_chunks  <- 0L

    pb <- if (verbose) progress::progress_bar$new(
      format = sprintf("  epoch %d/%d [:bar] :percent loss::loss",
                       epoch, epochs),
      total  = length(chunks), clear = FALSE, width = 70) else NULL

    for (chunk in chunks) {
      target_log4 <- target_log4_all[chunk, ]
      out <- model(tt$y[chunk, , ], tt$mask[chunk, , ])
      lo  <- pfn_loss(out,
                      target_log4,
                      tt$mu_s[chunk, ],
                      tt$f[chunk, , ],
                      weights = loss_weights)

      optimizer$zero_grad()
      lo$loss$backward()
      if (!is.null(grad_clip)) {
        torch::nn_utils_clip_grad_norm_(model$parameters, max_norm = grad_clip)
      }
      optimizer$step()

      train_acc <- train_acc + lo$parts
      n_chunks  <- n_chunks + 1L
      if (!is.null(pb)) pb$tick(tokens = list(
        loss = sprintf("%.3f", as.numeric(lo$loss))
      ))
    }

    train_parts <- train_acc / n_chunks
    train_total <- sum(unlist(loss_weights[c("theta","mu","f")]) * train_parts)

    val_metrics <- eval_slice(val_idx)
    sec <- as.numeric(Sys.time() - t0, units = "secs")

    history <- rbind(history, data.frame(
      epoch       = epoch,
      train_total = train_total,
      train_theta = train_parts["theta"],
      train_mu    = train_parts["mu"],
      train_f     = train_parts["f"],
      val_total   = val_metrics["total"],
      val_theta   = val_metrics["theta"],
      val_mu      = val_metrics["mu"],
      val_f       = val_metrics["f"],
      sec         = sec,
      row.names   = NULL
    ))

    if (verbose) {
      message(sprintf(
        "  epoch %d: train %.3f (theta %.3f / mu %.3f / f %.3f) | val %.3f | %.1fs",
        epoch, train_total,
        train_parts["theta"], train_parts["mu"], train_parts["f"],
        val_metrics["total"], sec))
    }

    if (!is.null(checkpoint)) {
      dir.create(dirname(checkpoint), recursive = TRUE,
                 showWarnings = FALSE)
      tmp <- paste0(checkpoint, ".tmp")
      torch::torch_save(
        list(model   = model$state_dict(),
             optim   = optimizer$state_dict(),
             epoch   = epoch,
             history = history,
             spec    = model$spec),
        tmp
      )
      file.rename(tmp, checkpoint)   # atomic on Windows + POSIX
    }
  }

  list(history = history, model = model, epoch = epochs)
}
