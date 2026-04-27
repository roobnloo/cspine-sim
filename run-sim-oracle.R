# Usage: Rscript run-sim-oracle.R --p=25 --q=50 --nobs=200 [--nrep=100]
# Runs oracle RegGMM only (delta=1, true mg supplied via mg_oracle)
source("impl/gmmreg_ssnal.R")
source("performance.R")
RhpcBLASctl::omp_set_num_threads(1)
RhpcBLASctl::blas_set_num_threads(1)

args <- commandArgs(trailingOnly = TRUE)

parse_arg <- function(args, key, default = NULL) {
  pat <- paste0("^--", key, "=(.+)$")
  m <- regmatches(args, regexpr(pat, args, perl = TRUE))
  if (length(m) == 0L) default else sub(pat, "\\1", m)
}

p <- as.integer(parse_arg(args, "p"))
q <- as.integer(parse_arg(args, "q"))
nobs <- as.integer(parse_arg(args, "nobs"))
nrep <- as.integer(parse_arg(args, "nrep", default = 100))
for (req in c("p", "q", "nobs")) {
  if (is.na(get(req))) stop("Required argument missing: --", req)
}

true_param_path <- file.path("data", sprintf("coef_p%dq%d.rds", p, q))
if (!file.exists(true_param_path)) stop("Generate data first.")
true_param <- tryCatch(
  readRDS(true_param_path),
  error = function(e) stop("Generate data first.")
)
setting_str <- sprintf("p%dq%d-n%d-d1.00", p, q, nobs)
metrics <- c(
  "tpr", "fpr", "tpr_pop", "fpr_pop", "tpr_cov", "fpr_cov", "beta_err", "gamma_err",
  "mu_err", "omega_err", "omega_tpr", "omega_fpr"
)

out_dir <- file.path("out", setting_str)
dir.create(out_dir, showWarnings = FALSE)

message(sprintf("Settings: p=%d, q=%d, nobs=%d, nrep=%d (oracle RegGMM, delta=1)", p, q, nobs, nrep))
message("Output directory: ", out_dir)

datasets <- readRDS(file.path("data", sprintf("%s.rds", setting_str)))

oracle_csv <- file.path(out_dir, "result-RegGMM-oracle.csv")
write(paste(metrics, collapse = ","), oracle_csv)

for (i in seq_len(nrep)) {
  message("Rep ", i)
  x_mat <- datasets[[i]]$X
  u_mat <- datasets[[i]]$U
  mu_true <- datasets[[i]]$mu
  omega_true <- datasets[[i]]$omega
  i_u <- cbind(1, u_mat)
  tictoc::tic()
  g_result <- gmmreg_ssnal(
    x_mat, u_mat,
    alpha = 0.75, nl1 = 100, lambda_factor = 0.1, num_cores = 25,
    mg_oracle = true_param$mg
  )
  tictoc::toc()
  g_est <- est_omega_mu(g_result$beta, g_result$gamma, i_u, u_mat, p, nobs, "original")
  pgs <- c(
    performance(g_result$beta, g_result$gamma, true_param$tb, true_param$mg),
    performance_supp(g_est$mu, g_est$omega, mu_true, omega_true)
  )
  print(round(pgs, 3))
  write(paste(pgs, collapse = ","), oracle_csv, append = TRUE)
}
