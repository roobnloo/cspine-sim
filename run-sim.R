library(cspine)
source("performance.R")
source("gmmreg.R")

args <- commandArgs(trailingOnly = TRUE)
p <- as.integer(args[1])
q <- as.integer(args[2])
n <- as.integer(args[3])
setting <- args[4]
stopifnot(setting %in% c("natural", "original"))
# p <- 25; q<-100; n<-200; setting <- "original"

setting_str <- sprintf("p%dq%d-n%d-%s", p, q, n, setting)
data_file <- file.path("data", paste0(setting_str, ".rds"))
generated <- readRDS(data_file)
tb_true <- generated[[length(generated)]]$tb
mg_true <- generated[[length(generated)]]$mg
nrep <- length(generated) - 1

metrics <- c("tpr", "fpr", "tpr_pop", "fpr_pop", "tpr_cov", "fpr_cov", "beta_err", "omega_err", "mean_err")
reggmm_results <- matrix(nrow = nrep, ncol = 9)
colnames(reggmm_results) <- metrics
cspine_results <- matrix(nrow = nrep, ncol = 9)
colnames(cspine_results) <- metrics

dir.create("./out", showWarnings = FALSE)
for (i in seq_len(nrep)) {
  message("Rep ", i)
  s <- generated[[i]]
  tictoc::tic()
  g_result <- gmmreg(s$X, s$U, ncores = 13)
  tictoc::toc()
  pgs <- performance(g_result, s, tb_true, mg_true)
  reggmm_results[i, ] <- pgs
  tictoc::tic()
  c_result <- cspine(s$X, s$U, ncores = 13)
  tictoc::toc()
  pcs <- performance(c_result, s, tb_true, mg_true)
  cspine_results[i, ] <- pcs
  print(rbind(round(pgs, 3), round(pcs, 3)))
  saveRDS(reggmm_results[1:i, ], file.path("out", paste0(setting_str, "-result-RegGMM.rds")))
  saveRDS(cspine_results[1:i, ], file.path("out", paste0(setting_str, "-result-cspine.rds")))
}
