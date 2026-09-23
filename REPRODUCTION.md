# Reproduction Guide

The current reproduction drivers are in `inst/reproduce/`:

| Driver | Result |
| --- | --- |
| `01_noiseless_size.R` | Frictionless null size |
| `02_noiseless_dense.R` | Frictionless dense-local power |
| `03_noiseless_sparse.R` | Frictionless sparse finite-activity power |
| `04_noisy_size.R` | Noisy null size |
| `05_noisy_dense.R` | Noisy dense-local power panels |
| `06_noisy_sparse.R` | Noisy sparse power panels |
| `07_empirical_3s.R` | Three-second noise-robust six-asset study |
| `08_empirical_coarse.R` | One-, two-, and five-minute frictionless six-asset study |

The empirical archive is deliberately not in the repository.  Download it from
<https://www.kaggle.com/datasets/brtnsmth/intraday-market-data>, then call the
three-second driver with the local archive path.  The script creates cleaned
regular-session data and all subsequent drivers consume those local CSV files.

All scripts use the full bootstrap replication count passed to the test
functions.  They intentionally have no bootstrap batch-size argument.
