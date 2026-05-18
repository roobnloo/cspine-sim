# Usage: Rscript run-sim.R --p=25 --q=50 --nobs=200 --delta=1 [--nrep=100] [--start_id=1]
source("impl/cspine_ssnal.R")
source("impl/gmmreg_ssnal.R")
source("performance.R")
if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
  RhpcBLASctl::omp_set_num_threads(1)
  RhpcBLASctl::blas_set_num_threads(1)
}

args <- commandArgs(trailingOnly = TRUE)

parse_arg <- function(args, key, default = NULL) {
  pat <- paste0("^--", key, "=(.+)$")
  m <- regmatches(args, regexpr(pat, args, perl = TRUE))
  if (length(m) == 0L) default else sub(pat, "\\1", m)
}

p <- as.integer(parse_arg(args, "p"))
q <- as.integer(parse_arg(args, "q"))
nobs <- as.integer(parse_arg(args, "nobs"))
delta <- as.numeric(parse_arg(args, "delta"))
nrep <- as.integer(parse_arg(args, "nrep", default = 100))
start_id <- as.integer(parse_arg(args, "start_id", default = 1))
for (req in c("p", "q", "nobs", "delta")) {
  if (is.na(get(req))) stop("Required argument missing: --", req)
}
if (delta < 0 || delta > 1) stop("--delta must be between 0 and 1 inclusive")

true_param_path <- file.path("data", sprintf("coef_p%dq%d.rds", p, q))
if (!file.exists(true_param_path)) stop("Generate data first.")
true_param <- tryCatch(
  readRDS(true_param_path),
  error = function(e) stop("Generate data first.")
)
setting_str <- sprintf("p%dq%d-n%d-d%.2f", p, q, nobs, delta)
metrics <- c(
  "tpr", "fpr", "tpr_pop", "fpr_pop", "tpr_cov", "fpr_cov", "beta_err", "rel_beta_err", "gamma_err",
  "mu_err", "omega_err", "omega_tpr", "omega_fpr"
)

out_dir <- file.path("out", setting_str)
dir.create(out_dir, showWarnings = FALSE)

message(sprintf("Settings: p=%d, q=%d, nobs=%d, delta=%.2f, nrep=%d, start_id=%d", p, q, nobs, delta, nrep, start_id))
message("Output directory: ", out_dir)

datasets <- readRDS(file.path("data", sprintf("%s.rds", setting_str)))

reggmm_csv <- file.path(out_dir, "result-RegGMM.csv")
cspine_csv <- file.path(out_dir, "result-cspine.csv")
if (start_id == 1) {
  write(paste(metrics, collapse = ","), reggmm_csv)
  write(paste(metrics, collapse = ","), cspine_csv)
} else {
  if (!file.exists(reggmm_csv)) stop("--start_id != 1 but output CSV does not exist: ", reggmm_csv)
  if (!file.exists(cspine_csv)) stop("--start_id != 1 but output CSV does not exist: ", cspine_csv)
}

for (i in seq(start_id, nrep)) {
  message("Rep ", i)
  x_mat <- datasets[[i]]$X
  u_mat <- datasets[[i]]$U
  mu_true <- datasets[[i]]$mu
  omega_true <- datasets[[i]]$omega
  i_u <- cbind(1, u_mat)
  tictoc::tic()
  g_result <- gmmreg_ssnal(x_mat, u_mat, alpha = 0.75, nl1 = 100, lambda_factor = 0.1, num_cores = 25)
  tictoc::toc()
  g_est <- est_omega_mu(g_result$beta, g_result$gamma, i_u, u_mat, p, nobs, "original")
  pgs <- c(
    performance(g_result$beta, g_result$gamma, true_param$tb, true_param$mg),
    performance_supp(g_est$mu, g_est$omega, mu_true, omega_true)
  )
  tictoc::tic()
  c_result <- cspine_ssnal(x_mat, u_mat, alpha = 0.75, nl1 = 100, lambda_factor = 0.1, num_cores = 25)
  tictoc::toc()
  c_est <- est_omega_mu(c_result$beta_raw, c_result$gamma, i_u, u_mat, p, nobs, "natural")
  pcs <- c(
    performance(c_result$beta_raw, c_result$gamma, true_param$tb, true_param$mg),
    performance_supp(c_est$mu, c_est$omega, mu_true, omega_true)
  )
  metric_mat <- rbind(round(pgs, 3), round(pcs, 3))
  rownames(metric_mat) <- c("RegGMM", "cspine")
  print(metric_mat)
  write(paste(pgs, collapse = ","), reggmm_csv, append = TRUE)
  write(paste(pcs, collapse = ","), cspine_csv, append = TRUE)
}
