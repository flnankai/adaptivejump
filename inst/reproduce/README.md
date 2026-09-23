# Reproduction scripts

Install `adaptivejump` first, then run these scripts with `Rscript`.  The
simulation defaults are the full paper replication counts: 5,000 for size and
1,000 for power.  Replace the second command-line argument by `100` for a
short smoke run.

The public empirical archive is not included in this repository.  Download it
from [Kaggle: Intraday market data](https://www.kaggle.com/datasets/brtnsmth/intraday-market-data)
and supply the local zip-file path to `07_empirical_3s.R`.  It contains weekly
three-second ThinkOrSwim snapshots.  The scripts retain only 09:30--16:00 U.S.
regular-session observations, remove non-weekdays and invalid daily paths, and
then use the last observation in each coarser interval.

```sh
Rscript 01_noiseless_size.R output/noiseless-size 100 30 199
Rscript 02_noiseless_dense.R output/noiseless-dense 100 30 199
Rscript 03_noiseless_sparse.R output/noiseless-sparse 100 30 199
Rscript 04_noisy_size.R output/noisy-size 100 30 199
Rscript 05_noisy_dense.R output/noisy-dense 100 30 199 "5,10"
Rscript 06_noisy_sparse.R output/noisy-sparse 100 30 199 "5,10"
Rscript 07_empirical_3s.R /path/to/archive.zip output/empirical-3s 30 199
Rscript 08_empirical_coarse.R output/empirical-3s/clean_data output/empirical-coarse 30 199
```

All bootstrap replications are evaluated together within each test call.  The
scripts deliberately do not expose or use a bootstrap batch-size parameter.
