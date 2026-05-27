test_that("diagnose_bayes returns expected columns and types", {
  skip_if_not_installed("BayesLogit")
  skip_if_not(file.exists("../../implementation/simulation.R"),
              "simulation helpers not found")
  source("../../implementation/simulation.R", local = TRUE)

  set.seed(2026)
  n <- 10; nt <- 26
  coords <- data.frame(id = factor(seq_len(n)),
                       lat = runif(n, 0, 5), lon = runif(n, 0, 5),
                       mu = log(runif(n, 5, 50)))
  sk <- space_kernel(coords, length_scale = 1.5)
  tk <- time_kernel(seq_len(nt), periodic_scale = 1,
                    long_term_scale = 50, period = 13)
  sim <- simulate_data(n, nt, coords, sk, tk, r = 20)
  obs <- observed_data(sim, p_one = 0.2, p_switch = 0.3)
  fb <- fit_bayes(obs, n_sweeps = 100, burnin = 30, n_chains = 2,
                  period = 13, verbose = FALSE)

  d <- diagnose_bayes(fb)
  expect_s3_class(d, "weave_diagnosis")
  expect_setequal(d$param,
                  c("length_scale", "periodic_scale", "long_term_scale", "r"))
  expect_true(all(c("median", "q025", "q975", "rhat", "ess",
                    "prior_overlap_p", "is_poorly_mixed",
                    "is_prior_driven", "status") %in% names(d)))
  expect_true(is.logical(d$is_poorly_mixed))
  expect_true(is.logical(d$is_prior_driven))
  # All quantiles ordered: q025 <= median <= q975.
  expect_true(all(d$q025 <= d$median))
  expect_true(all(d$median <= d$q975))
})

test_that("diagnose_bayes flags poor mixing when Rhat is high", {
  # Build a synthetic weave_bayes-like object with deliberately stuck
  # chains so we don't have to wait for a real one.
  skip_if_not_installed("coda")
  fb <- structure(
    list(
      n_sweeps = 100, burnin = 20, n_chains = 2,
      theta_trace = array(NA_real_, c(100, 3, 2),
                          dimnames = list(NULL,
                                          c("length_scale", "periodic_scale",
                                            "long_term_scale"),
                                          c("chain1", "chain2"))),
      r_trace     = matrix(NA_real_, 100, 2),
      priors = list(
        theta = list(
          length_scale    = list(sample = function() stats::rlnorm(1, 0, 1)),
          periodic_scale  = list(sample = function() stats::rlnorm(1, 0, 1)),
          long_term_scale = list(sample = function() stats::rlnorm(1, 0, 1))
        ),
        r = list(sample = function() stats::rgamma(1, 2, 0.1))
      )
    ),
    class = c("weave_bayes", "list")
  )
  # Chain 1 oscillates around 1, chain 2 around 5 -- maximally stuck.
  set.seed(11)
  fb$theta_trace[, "length_scale", 1] <- 1 + rnorm(100, 0, 0.05)
  fb$theta_trace[, "length_scale", 2] <- 5 + rnorm(100, 0, 0.05)
  # The other params can be unmoving (Rhat undefined but won't crash).
  fb$theta_trace[, "periodic_scale", ]  <- 1 + rnorm(200, 0, 0.05)
  fb$theta_trace[, "long_term_scale", ] <- 50 + rnorm(200, 0, 0.05)
  fb$r_trace[, ] <- 10 + rnorm(200, 0, 0.05)

  d <- diagnose_bayes(fb)
  expect_true(d$is_poorly_mixed[d$param == "length_scale"])
})
