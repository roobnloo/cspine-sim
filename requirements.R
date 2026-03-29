# Install required packages for cspine-sim
# Run this script once before running any simulations.

if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")

# --- GitHub packages ---

if (!requireNamespace("cspine", quietly = TRUE))
  remotes::install_github("cspine/roobnloo")

# --- CRAN packages ---

cran_pkgs <- c(
  "sparsegl",
  "tictoc",
  "RhpcBLASctl"
)

to_install <- cran_pkgs[!sapply(cran_pkgs, requireNamespace, quietly = TRUE)]
if (length(to_install) > 0)
  install.packages(to_install)
