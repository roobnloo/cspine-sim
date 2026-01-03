library(cspine)
source("performance.R")
source("gmmreg.R")

args <- commandArgs(trailingOnly = TRUE)
sparsity_setting <- args[1]
stopifnot(sparsity_setting %in% c("beta", "gamma"))

p <- 25
q <- 10
n <- 300
model <- "natural"

all_data <- readRDS(sprintf("data/p%dq%d-n%d-%s-varying-sparsity.rds", p, q, n, model))
if (sparsity_setting == "gamma") {
  all_data <- readRDS(sprintf("data/p%dq%d-n%d-%s-varying-sparsity-gamma.rds", p, q, n, model))
}

nrep <- length(all_data)
metrics <- c(
  "tpr", "fpr", "tpr_pop", "fpr_pop", "tpr_cov", "fpr_cov",
  "beta_err", "omega_err", "gamma_err", "mean_err", "omega_tpr", "omega_fpr", paste0("s_", sparsity_setting)
)
c_perf <- matrix(nrow = nrep, ncol = length(metrics), dimnames = list(NULL, metrics))
g_perf <- matrix(nrow = nrep, ncol = length(metrics), dimnames = list(NULL, metrics))
dir.create("./out", showWarnings = FALSE)
outpath <- file.path("out", sprintf("p%dq%d-n%d-%s-varying-sparsity-result", p, q, n, model))
if (sparsity_setting == "gamma") {
  outpath <- file.path("out", sprintf("p%dq%d-n%d-%s-varying-sparsity-result-gamma", p, q, n, model))
}

for (i in seq_len(nrep)) {
  message("Rep ", i)
  s <- all_data[[i]]
  sparsity <- sum(abs(s$tb) > 0)
  if (sparsity_setting == "gamma") {
    sparsity <- sum(abs(s$mg) > 0)
  }
  tictoc::tic()
  g_result <- gmmreg(s$X, s$U, ncores = 13)
  tictoc::toc()
  pgs <- performance(g_result, s, s$tb, s$mg)
  g_perf[i, ] <- c(pgs, sparsity)
  saveRDS(g_perf[1:i, ], paste0(outpath, "-RegGMM.rds"))

  tictoc::tic()
  c_result <- cspine(s$X, s$U, ncores = 13)
  tictoc::toc()
  pcs <- performance(c_result, s, s$tb, s$mg)
  c_perf[i, ] <- c(pcs, sparsity)
  saveRDS(c_perf[1:i, ], paste0(outpath, "-cspine.rds"))

  print(rbind(round(g_perf[i, ], 3), round(c_perf[i, ], 3)))
}
