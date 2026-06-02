#! /bin/bash

# By default, the p=25, q=50, n=200 settings are generated.
# For additional results, uncomment the desired setting.

# p=25, q=50
Rscript generate-data-simple.R --p=25 --q=50 --nobs=200 --delta=0
# Rscript generate-data-simple.R --p=25 --q=50 --nobs=200 --delta=0.25
# Rscript generate-data-simple.R --p=25 --q=50 --nobs=200 --delta=0.5
# Rscript generate-data-simple.R --p=25 --q=50 --nobs=200 --delta=0.75
Rscript generate-data-simple.R --p=25 --q=50 --nobs=200 --delta=1
# Rscript generate-data-simple.R --p=25 --q=50 --nobs=400 --delta=0
# Rscript generate-data-simple.R --p=25 --q=50 --nobs=400 --delta=1

# p=25, q=100
# Rscript generate-data-simple.R --p=25 --q=100 --nobs=200 --delta=0
# Rscript generate-data-simple.R --p=25 --q=100 --nobs=200 --delta=1
# Rscript generate-data-simple.R --p=25 --q=100 --nobs=400 --delta=0
# Rscript generate-data-simple.R --p=25 --q=100 --nobs=400 --delta=1

# varying SNR
# Rscript generate-data-varying-snr.R