#! /bin/bash

# By default, cspine and RegGMM are run on the p=25, q=50, n=200 settings.
# For additional results, uncomment the desired setting.

# p=25, q=50, n=200
Rscript run-sim.R --p=25 --q=50 --nobs=200 --delta=0
# Rscript run-sim.R --p=25 --q=50 --nobs=200 --delta=0.25
# Rscript run-sim.R --p=25 --q=50 --nobs=200 --delta=0.5
# Rscript run-sim.R --p=25 --q=50 --nobs=200 --delta=0.75
Rscript run-sim.R --p=25 --q=50 --nobs=200 --delta=1
# Rscript run-sim-oracle.R --p=25 --q=50 --nobs=200
# Rscript run-sim-mtgmmreg.R --p=25 --q=50 --nobs=200 --delta=0
# Rscript run-sim-mtgmmreg.R --p=25 --q=50 --nobs=200 --delta=0.25
# Rscript run-sim-mtgmmreg.R --p=25 --q=50 --nobs=200 --delta=0.5
# Rscript run-sim-mtgmmreg.R --p=25 --q=50 --nobs=200 --delta=0.75
# Rscript run-sim-mtgmmreg.R --p=25 --q=50 --nobs=200 --delta=1

# p=25, q=50, n=400
# Rscript run-sim.R --p=25 --q=50 --nobs=400 --delta=0
# Rscript run-sim.R --p=25 --q=50 --nobs=400 --delta=1
# Rscript run-sim-oracle.R --p=25 --q=50 --nobs=400
# Rscript run-sim-mtgmmreg.R --p=25 --q=50 --nobs=400 --delta=0
# Rscript run-sim-mtgmmreg.R --p=25 --q=50 --nobs=400 --delta=1

# p=25, q=100, n=200
# Rscript run-sim.R --p=25 --q=100 --nobs=200 --delta=0
# Rscript run-sim.R --p=25 --q=100 --nobs=200 --delta=1
# Rscript run-sim-oracle.R --p=25 --q=100 --nobs=200
# Rscript run-sim-mtgmmreg.R --p=25 --q=100 --nobs=200 --delta=0
# Rscript run-sim-mtgmmreg.R --p=25 --q=100 --nobs=200 --delta=1

# p=25, q=100, n=400
# Rscript run-sim.R --p=25 --q=100 --nobs=400 --delta=0
# Rscript run-sim.R --p=25 --q=100 --nobs=400 --delta=1
# Rscript run-sim-oracle.R --p=25 --q=100 --nobs=400
# Rscript run-sim-mtgmmreg.R --p=25 --q=100 --nobs=400 --delta=0
# Rscript run-sim-mtgmmreg.R --p=25 --q=100 --nobs=400 --delta=1

# varying SNR (p=25, q=50, n=200)
# Rscript run-sim-varying-snr.R --p=25 --q=50 --nobs=200 --c=0.29
# Rscript run-sim-varying-snr.R --p=25 --q=50 --nobs=200 --c=0.41
# Rscript run-sim-varying-snr.R --p=25 --q=50 --nobs=200 --c=0.59
# Rscript run-sim-varying-snr.R --p=25 --q=50 --nobs=200 --c=0.87
# Rscript run-sim-varying-snr.R --p=25 --q=50 --nobs=200 --c=1.35
# Rscript run-sim-oracle-varying-snr.R --p=25 --q=50 --nobs=200 --c=0.29
# Rscript run-sim-oracle-varying-snr.R --p=25 --q=50 --nobs=200 --c=0.41
# Rscript run-sim-oracle-varying-snr.R --p=25 --q=50 --nobs=200 --c=0.59
# Rscript run-sim-oracle-varying-snr.R --p=25 --q=50 --nobs=200 --c=0.87
# Rscript run-sim-oracle-varying-snr.R --p=25 --q=50 --nobs=200 --c=1.35