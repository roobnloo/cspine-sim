library(cspine)
source("performance.R")
source("gmmreg.R")

set.seed(145461)

p <- 25
q <- 50
n <- 200
model <- "original"

all_data <- readRDS(sprintf("data/p%dq%d-n%d-%s-varying-snr.rds", p, q, n, model))

nrep <- length(all_data)
metrics <- c(
  "tpr", "fpr", "tpr_pop", "fpr_pop", "tpr_cov", "fpr_cov",
  "beta_err", "omega_err", "gamma_err", "mean_err", "omega_tpr", "omega_fpr", "snr"
)
c_perf <- matrix(nrow = nrep, ncol = length(metrics), dimnames = list(NULL, metrics))
g_perf <- matrix(nrow = nrep, ncol = length(metrics), dimnames = list(NULL, metrics))
dir.create("./out", showWarnings = FALSE)
outpath <- file.path("out", sprintf("p%dq%d-n%d-%s-varying-snr-result", p, q, n, model))

for (i in seq_len(nrep)) {
  message("Rep ", i)
  s <- all_data[[i]]
  snr <- mean(s$snr)
  tictoc::tic()
  g_result <- gmmreg(s$X, s$U, ncores = 13)
  tictoc::toc()
  pgs <- performance(g_result, s, s$tb, s$mg)
  g_perf[i, ] <- c(pgs, snr)
  saveRDS(g_perf[1:i, ], paste0(outpath, "-RegGMM.rds"))

  tictoc::tic()
  c_result <- cspine(s$X, s$U, ncores = 13)
  tictoc::toc()
  pcs <- performance(c_result, s, s$tb, s$mg)
  c_perf[i, ] <- c(pcs, snr)
  saveRDS(c_perf[1:i, ], paste0(outpath, "-cspine.rds"))

  print(rbind(round(g_perf[i, ], 3), round(c_perf[i, ], 3)))
}
