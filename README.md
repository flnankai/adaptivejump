# adaptivejump

`adaptivejump` implements the methods in *Adaptive Test for Jump* by Huifang
Ma and Long Feng:

- `aj_test()` implements the frictionless Ait-Sahalia--Jacod (AJ) power-
  variation ratio test.
- `lm_test()` implements the frictionless Lee--Mykland (LM) max test with
  Gumbel or pathwise-bootstrap calibration.
- `adaptive_jump_test()` returns AJ, LM, and their Cauchy combination (CC).
- `ajj_test()` implements the pre-averaged Ait-Sahalia--Jacod--Li-style
  noise-robust ratio statistic (AJJ).
- `lm_noise_test()` implements the noisy local-average LM statistic with the
  LM (2012) Table 5 block-size interpolation, TSRSV spot-volatility
  standardization, and a fully recursive parametric bootstrap.
- `adaptive_jump_test_noise()` returns AJJ, noisy LM, and CC.

The package contains no proprietary or large empirical data.  It includes
small synthetic examples only.  The empirical study uses the public ThinkOrSwim
three-second data published as Kaggle's
[Intraday market data](https://www.kaggle.com/datasets/brtnsmth/intraday-market-data).
Download the archive manually, then pass its local path to `read_tos_archive()`
or to the scripts in `inst/reproduce/`.

## Installation

```r
remotes::install_github("flnankai/adaptivejump")
```

For a local checkout:

```r
install.packages("adaptivejump_0.1.0.tar.gz", repos = NULL, type = "source")
```

## Basic usage

```r
library(adaptivejump)

sim <- simulate_jump_path(model = "dense", theta = 2, mark_sd = 1)
out <- adaptive_jump_test(exp(sim$latent_log_price), bootstrap_rep = 199)
out

noisy <- simulate_jump_path(model = "sparse", lambda = 1.5, mark_sd = sqrt(0.05),
                            noise_sd = 0.005)
out_noise <- adaptive_jump_test_noise(
  exp(noisy$observed_log_price),
  bootstrap_rep = 199,
  critical_value_c = 4
)
out_noise
```

## Reproducing the paper

The scripts in `inst/reproduce/` are standalone drivers after the package is
installed.  They generate output files in a user-selected directory and do
not read or modify the original working-paper directory.

| Script | Main result |
| --- | --- |
| `01_noiseless_size.R` | Frictionless null size |
| `02_noiseless_power.R` | Dense and sparse frictionless power |
| `03_noisy_size.R` | Noisy null size under Gaussian and standardized `t8` noise |
| `04_noisy_dense_power.R` | Four-panel noisy dense-local power curves |
| `05_noisy_sparse_power.R` | Four-panel noisy sparse power curves |
| `06_empirical_3s.R` | Three-second noise-robust daily tests for six assets |
| `07_empirical_coarse.R` | One-, two-, and five-minute frictionless daily tests |

The full paper settings use 5,000 size replications and 1,000 power
replications.  Each script accepts a smaller replication count for a quick
smoke run.  Parallel workers are controlled with a `workers` command-line
argument; no bootstrap batching is used.

## Data format

Daily data passed to `daily_jump_tests()` must have columns `date`, `time`, and
`price`.  The data utilities retain the 09:30:00--16:00:00 U.S. regular session,
discard days with missing or constant prices, and use the last available price
in each resampling interval.

## References

- Ait-Sahalia, Y. and Jacod, J. (2009). Testing for jumps in a discretely
  observed process. *Annals of Statistics*.
- Ait-Sahalia, Y., Jacod, J. and Li, J. (2012). Testing for jumps in noisy
  high frequency data. *Journal of Econometrics*.
- Lee, S. S. and Mykland, P. A. (2008). Jumps in financial markets: A new
  nonparametric test and jump dynamics. *Review of Financial Studies*.
- Lee, S. S. and Mykland, P. A. (2012). Jumps in equilibrium prices and market
  microstructure noise. *Journal of Econometrics*.
