# Hybrid PFN evaluation

Same trained model as before (`final_model.pt`, 50k sims x 100 epochs,
n = 8, nt = 26, period = 13).

## Hybrid SBC (rank-based, 50 fresh datasets, 200 posterior draws each)

| parameter | KS p vs Uniform |
|---|---|
| length_scale    | 0.570 |
| periodic_scale  | 0.061 |
| long_term_scale | 0.123 |
| log_r           | 0.016 |
| mu_s            | 0.037 |
| f               | 4.04e-196 |

## Comparison: 95% credible interval coverage and width (5 fresh datasets)

| parameter | PFN cov | Hybrid cov | Bayes cov | PFN width | Hybrid width | Bayes width |
|---|---|---|---|---|---|---|
| length_scale    | 1.00 | 1.00 | 1.00 | 23.34 | 21.47 | 2.54 |
| periodic_scale  | 1.00 | 1.00 | 0.80 | 3.80 | 3.56 | 0.95 |
| long_term_scale | 1.00 | 1.00 | 1.00 | 180.62 | 181.60 | 67.25 |
| r               | 1.00 | 1.00 | 1.00 | 30.62 | 30.29 | 19.13 |

## Runtime + posterior-predictive coverage

| dataset | PFN s | Hybrid s | Bayes s | PFN PPC | Hybrid PPC | Bayes PPC |
|---|---|---|---|---|---|---|
| 1 | 0.07 | 107.82 | 11.60 | 1.00 | 0.99 | 0.96 |
| 2 | 0.06 | 57.93 | 12.42 | 1.00 | 1.00 | 0.99 |
| 3 | 0.09 | 79.36 | 15.09 | 1.00 | 0.99 | 0.98 |
| 4 | 0.11 | 66.27 | 13.11 | 1.00 | 1.00 | 0.99 |
| 5 | 0.11 | 52.47 | 10.45 | 1.00 | 1.00 | 1.00 |

Hybrid mean speedup vs fit_bayes: **0.2x**

