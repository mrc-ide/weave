test_that("fit_bayes recovers known length_scale and r on simulated data", {
  skip_on_cran()
  skip_if_not_installed("BayesLogit")
  skip_if_not(file.exists("../../implementation/simulation.R"),
              "simulation helpers not found")

  source("../../implementation/simulation.R", local = TRUE)

  set.seed(2026)
  n  <- 20
  nt <- 52
  true_length_scale    <- 1.5
  true_periodic_scale  <- 1
  true_long_term_scale <- 80
  true_r               <- 20

  coordinates <- data.frame(
    id  = factor(seq_len(n)),
    lat = runif(n, 0, 5), lon = runif(n, 0, 5),
    mu  = log(runif(n, 5, 50))
  )
  space_k <- space_kernel(coordinates, length_scale = true_length_scale)
  time_k  <- time_kernel(seq_len(nt),
                         periodic_scale  = true_periodic_scale,
                         long_term_scale = true_long_term_scale,
                         period          = 52)
  sim <- simulate_data(n, nt, coordinates, space_k, time_k, r = true_r)
  obs <- observed_data(sim, p_one = 0.2, p_switch = 0.3)

  fb <- fit_bayes(obs, n_sweeps = 500, burnin = 200,
                  store_f = "summary", verbose = FALSE)
  post <- (fb$burnin + 1L):fb$n_sweeps

  # 95% credible intervals must cover the truth on each parameter.
  # theta_trace is now (sweep, param, chain); r_trace is (sweep, chain).
  # With n_chains = 1 we still index the third dim explicitly.
  ls_q  <- stats::quantile(fb$theta_trace[post, "length_scale",    ],
                           c(0.025, 0.975))
  ps_q  <- stats::quantile(fb$theta_trace[post, "periodic_scale",  ],
                           c(0.025, 0.975))
  r_q   <- stats::quantile(fb$r_trace[post, ], c(0.025, 0.975))

  expect_gte(true_length_scale,   ls_q[1] * 0.7)   # 30% tolerance on lower
  expect_lte(true_length_scale,   ls_q[2] * 1.3)   # 30% on upper
  expect_gte(true_periodic_scale, ps_q[1] * 0.5)
  expect_lte(true_periodic_scale, ps_q[2] * 1.5)
  expect_gte(true_r,              r_q[1]  * 0.5)
  expect_lte(true_r,              r_q[2]  * 1.5)
})
