# Install required packages for cspine-sim
# Run this script once before running any simulations.

if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes", repos = "https://cloud.r-project.org")

# --- GitHub packages ---

if (!requireNamespace("cspine", quietly = TRUE))
  remotes::install_github("roobnloo/sglssnal")

# --- CRAN packages ---

cran_pkgs <- c(
  "Rcpp",
  "tictoc",
  "Matrix",
  "RSpectra",
  "sparsegl"
)

to_install <- cran_pkgs[!sapply(cran_pkgs, requireNamespace, quietly = TRUE)]
if (length(to_install) > 0)
  install.packages(to_install, repos = "https://cloud.r-project.org")
