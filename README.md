# adaptivejump

`adaptivejump` implements the methods in *Adaptive Test for Jump* by Huifang
Ma and Long Feng.

## Methods

- `aj_test()` implements the frictionless Ait-Sahalia--Jacod (AJ)
  power-variation ratio test.
- `lm_test()` implements the frictionless Lee--Mykland (LM) max test with
  Gumbel or pathwise-bootstrap calibration.
- `adaptive_jump_test()` returns AJ, LM, and their Cauchy combination (CC).
- `ajj_test()` implements the pre-averaged noise-robust ratio statistic (AJJ).
- `lm_noise_test()` implements the noisy local-average LM statistic with
  LM (2012) Table 5 block-size interpolation, TSRSV spot-volatility
  standardization, and a fully recursive parametric bootstrap.
- `adaptive_jump_test_noise()` returns AJJ, noisy LM, and CC.

The package also exposes `estimate_noise_sd()`, `estimate_spot_volatility()`,
`simulate_jump_path()`, `run_jump_simulation()`, `clean_intraday_data()`,
`resample_last_tick()`, `read_tos_archive()`, `daily_jump_tests()`, and
`bh_select()`.

## Installation

```r
remotes::install_github("flnankai/adaptivejump")
```

## Basic usage

```r
library(adaptivejump)

sim <- simulate_jump_path(model = "dense", theta = 2, mark_sd = 1)
adaptive_jump_test(exp(sim$latent_log_price), bootstrap_rep = 199)

noisy <- simulate_jump_path(
  model = "sparse", lambda = 1.5, mark_sd = sqrt(0.05), noise_sd = 0.005
)
adaptive_jump_test_noise(
  exp(noisy$observed_log_price), bootstrap_rep = 199, critical_value_c = 4
)
```

## Reproduction

The package contains no empirical raw data.  The public ThinkOrSwim archive
used for the six-asset application is available from Kaggle:
[Intraday market data](https://www.kaggle.com/datasets/brtnsmth/intraday-market-data).
Download it manually, then provide its local path to `read_tos_archive()` or
`inst/reproduce/07_empirical_3s.R`.

| Script | Main result |
| --- | --- |
| `01_noiseless_size.R` | Frictionless null size |
| `02_noiseless_dense.R` | Dense-local frictionless power |
| `02b_noiseless_dense_one_second_var1.R` | One-second dense local power with `Y ~ N(0, 1)` |
| `03_noiseless_sparse.R` | Sparse finite-activity frictionless power |
| `04_noisy_size.R` | Noisy null size |
| `05_noisy_dense.R` | Four-panel noisy dense-local power curves |
| `06_noisy_sparse.R` | Four-panel noisy sparse power curves |
| `07_empirical_3s.R` | Three-second noise-robust daily tests for six assets |
| `08_empirical_coarse.R` | One-, two-, and five-minute frictionless daily tests |

The simulation defaults match the paper: 5,000 size replications and 1,000
power replications.  Each driver accepts a smaller count for a smoke run and a
worker count for Windows parallel execution.  Bootstrap replications are never
split into batches.

See [REPRODUCTION.md](REPRODUCTION.md) and `inst/reproduce/README.md` for
commands and input/output conventions.

## References

- Ait-Sahalia, Y. and Jacod, J. (2009). Testing for jumps in a discretely
  observed process. *Annals of Statistics*.
- Ait-Sahalia, Y., Jacod, J. and Li, J. (2012). Testing for jumps in noisy
  high frequency data. *Journal of Econometrics*.
- Lee, S. S. and Mykland, P. A. (2008). Jumps in financial markets: A new
  nonparametric test and jump dynamics. *Review of Financial Studies*.
- Lee, S. S. and Mykland, P. A. (2012). Jumps in equilibrium prices and market
  microstructure noise. *Journal of Econometrics*.
