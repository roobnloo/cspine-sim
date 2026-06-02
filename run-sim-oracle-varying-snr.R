# Usage: Rscript run-sim-oracle-varying-snr.R --p=25 --q=50 --nobs=200 --c=0.29 [--nrep=50] [--start_id=1]
# Runs oracle RegGMM only (delta=1, true mg supplied via mg_oracle)
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
c_val <- as.numeric(parse_arg(args, "c"))
nrep <- as.integer(parse_arg(args, "nrep", default = 50))
start_id <- as.integer(parse_arg(args, "start_id", default = 1))
for (req in c("p", "q", "nobs", "c_val")) {
  if (is.na(get(req))) stop("Required argument missing: --", sub("_val", "", req))
}

# ---- Load true parameters (fixed Gamma and B_h) ----
true_param_path <- file.path("data", sprintf("coef_p%dq%d.rds", p, q))
if (!file.exists(true_param_path)) stop("True parameter file not found: ", true_param_path)
true_param <- readRDS(true_param_path)
tb <- true_param$tb # p × p × (q+1)
mg <- true_param$mg # p × q

# Pre-compute vectorized Omega pieces for computing omega_true on-the-fly
Bmat <- matrix(tb[, , 2:(q + 1)], p * p, q)
B0_vec <- as.vector(tb[, , 1])
diag_idx <- seq(1L, p * p, by = p + 1L)

# ---- Load data ----
c_str <- gsub("\\.", "p", as.character(c_val))
data_file <- file.path("data", sprintf("p%dq%d-n%d-varying-snr-c%s.rds", p, q, nobs, c_str))
if (!file.exists(data_file)) stop("Data file not found: ", data_file)
all_data <- readRDS(data_file)
nrep <- min(nrep, length(all_data))

# ---- Output setup ----
metrics <- c(
  "tpr", "fpr", "tpr_pop", "fpr_pop", "tpr_cov", "fpr_cov", "beta_err", "rel_beta_err",
  "gamma_err", "mu_err", "omega_err", "omega_tpr", "omega_fpr", "snr"
)
setting_str <- sprintf("p%dq%d-n%d-varying-snr-c%s", p, q, nobs, c_str)
out_dir <- file.path("out", setting_str)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

message(sprintf(
  "Settings: p=%d, q=%d, nobs=%d, c=%.2f, nrep=%d, start_id=%d (oracle RegGMM)",
  p, q, nobs, c_val, nrep, start_id
))
message("Output directory: ", out_dir)

oracle_csv <- file.path(out_dir, "result-RegGMM-oracle.csv")
if (start_id == 1L) {
  write(paste(metrics, collapse = ","), oracle_csv)
} else {
  if (!file.exists(oracle_csv)) stop("--start_id != 1 but output CSV does not exist: ", oracle_csv)
}

# ---- Main loop ----
for (i in seq(start_id, nrep)) {
  message("Rep ", i)
  s <- all_data[[i]]
  x_mat <- s$X
  u_mat <- s$U # already scaled by c
  i_u <- cbind(1, u_mat)
  snr <- s$snr

  # Compute true mu and omega from fixed params and scaled U
  mu_true <- t(mg %*% t(u_mat)) # n × p
  ov <- B0_vec + Bmat %*% t(u_mat) # (p^2) × n
  ov[diag_idx, ] <- 1
  omega_true <- array(ov, dim = c(p, p, nobs)) # p × p × n

  tictoc::tic()
  g_result <- gmmreg_ssnal(
    x_mat, u_mat,
    alpha = 0.75, nl1 = 100, lambda_factor = 0.1, num_cores = 25,
    mg_oracle = true_param$mg
  )
  tictoc::toc()
  g_est <- est_omega_mu(g_result$beta, g_result$gamma, i_u, u_mat, p, nobs, "original")
  pgs <- c(
    performance(g_result$beta, g_result$gamma, tb, mg),
    performance_supp(g_est$mu, g_est$omega, mu_true, omega_true),
    snr
  )
  print(round(pgs, 3))
  write(paste(pgs, collapse = ","), oracle_csv, append = TRUE)
}
