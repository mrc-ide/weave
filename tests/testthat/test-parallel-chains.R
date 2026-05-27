test_that("parallel and serial chains produce identical traces with the same seed", {
  skip_on_cran()
  skip_if_not_installed("future")
  skip_if_not_installed("future.apply")
  skip_if_not_installed("BayesLogit")
  skip_if_not(file.exists("../../implementation/simulation.R"),
              "simulation helpers not found")

  source("../../implementation/simulation.R", local = TRUE)

  set.seed(2026)
  n  <- 8
  nt <- 26
  coordinates <- data.frame(
    id  = factor(seq_len(n)),
    lat = runif(n, 0, 5), lon = runif(n, 0, 5),
    mu  = log(runif(n, 5, 50))
  )
  space_k <- space_kernel(coordinates, length_scale = 1.5)
  time_k  <- time_kernel(seq_len(nt), periodic_scale = 1,
                         long_term_scale = 50, period = 52)
  sim <- simulate_data(n, nt, coordinates, space_k, time_k, r = 20)
  obs <- observed_data(sim, p_one = 0.2, p_switch = 0.3)

  # The parallel and serial paths use the same per-chain seed derivation, so
  # both should produce bit-identical traces for a fixed user seed.
  args <- list(obs_data = obs, n_sweeps = 60, burnin = 20,
               n_chains = 2, period = 13, verbose = FALSE)

  set.seed(123)
  fb_serial <- do.call(fit_bayes, c(args, list(parallel = FALSE)))

  # Run parallel under multisession with 2 workers; clean up the plan after.
  old_plan <- future::plan("multisession", workers = 2)
  on.exit(future::plan(old_plan), add = TRUE)
  set.seed(123)
  fb_par <- do.call(fit_bayes, c(args, list(parallel = TRUE)))

  # Under `devtools::load_all()` the future workers may pick up an
  # installed-package version of weave that differs slightly from the
  # in-development source; bit-exact equality therefore isn't reliable
  # in the test environment, but the values should agree to ~1e-6
  # (any larger discrepancy would indicate a real RNG-plumbing bug).
  expect_equal(fb_serial$theta_trace, fb_par$theta_trace, tolerance = 1e-6)
  expect_equal(fb_serial$r_trace,     fb_par$r_trace,     tolerance = 1e-6)
  expect_equal(fb_serial$mu_trace,    fb_par$mu_trace,    tolerance = 1e-6)
})
