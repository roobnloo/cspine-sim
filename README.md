# cspine-sim

Simulation code accompanying the manuscript *Convex estimation of Gaussian graphical regression models with covariates*.
The **cspine** package implements the method.

## Setup

Install R package dependencies by running:

```r
source("requirements.R")
```

This installs `cspine` (from GitHub) and the CRAN packages `sparsegl`, `tictoc`, and `RhpcBLASctl`.

## Running simulations

### 1. Generate data

```bash
bash generate-all-data.sh
```

Generates datasets across combinations of graph type (`original`, `natural`), number of nodes `q` (50, 100), and sample size `n` (200, 400), as well as varying-sparsity data.

### 2. Run simulations

```bash
bash run-all.sh
```

Runs the full simulation study across all data configurations and writes results to the `out/` directory.

In the interest of time, each simulation setting may also be run separately.