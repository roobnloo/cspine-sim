# Usage: Rscript run-sim-varying-c.R --nobs=200 [--nrep=30]
source("impl/cspine_ssnal.R")
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

nobs <- as.integer(parse_arg(args, "nobs"))
nrep <- as.integer(parse_arg(args, "nrep", default = 30L))
if (is.na(nobs)) stop("Required argument missing: --nobs")

coef_param <- readRDS(file.path("data", "coef_p25q10_vary_c.rds"))
tB_all <- coef_param$tB_all
mg <- coef_param$mg
c_scale <- coef_param$c_scale

p <- dim(tB_all)[1]
q <- dim(tB_all)[3] - 1L

metrics <- c(
  "tpr", "fpr", "tpr_pop", "fpr_pop", "tpr_cov", "fpr_cov", "beta_err", "rel_beta_err", "gamma_err",
  "mu_err", "omega_err", "omega_tpr", "omega_fpr"
)

for (c_idx in seq_along(c_scale)) {
  c_val <- c_scale[c_idx]
  setting_str <- sprintf("vary-c-c%.2f-n%d", c_val, nobs)
  data_file <- file.path("data", sprintf("%s.rds", setting_str))
  if (!file.exists(data_file)) stop("Data file not found: ", data_file, ". Run generate-data-varying-c.R first.")

  datasets <- readRDS(data_file)
  tb <- tB_all[, , , c_idx]

  out_dir <- file.path("out", setting_str)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  message(sprintf(
    "\n=== c_scale=%.2f (p=%d, q=%d, nobs=%d, nrep=%d) ===",
    c_val, p, q, nobs, nrep
  ))
  message("Output directory: ", out_dir)

  reggmm_csv <- file.path(out_dir, "result-RegGMM.csv")
  cspine_csv <- file.path(out_dir, "result-cspine.csv")
  oracle_csv <- file.path(out_dir, "result-RegGMM-oracle.csv")
  write(paste(metrics, collapse = ","), reggmm_csv)
  write(paste(metrics, collapse = ","), cspine_csv)
  write(paste(metrics, collapse = ","), oracle_csv)

  for (i in seq_len(nrep)) {
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
      performance(g_result$beta, g_result$gamma, tb, mg),
      performance_supp(g_est$mu, g_est$omega, mu_true, omega_true)
    )

    tictoc::tic()
    c_result <- cspine_ssnal(x_mat, u_mat, alpha = 0.75, nl1 = 100, lambda_factor = 0.1, num_cores = 25)
    tictoc::toc()
    c_est <- est_omega_mu(c_result$beta_raw, c_result$gamma, i_u, u_mat, p, nobs, "natural")
    pcs <- c(
      performance(c_result$beta_raw, c_result$gamma, tb, mg),
      performance_supp(c_est$mu, c_est$omega, mu_true, omega_true)
    )

    tictoc::tic()
    o_result <- gmmreg_ssnal(
      x_mat, u_mat,
      alpha = 0.75, nl1 = 100, lambda_factor = 0.1, num_cores = 25,
      mg_oracle = mg
    )
    tictoc::toc()
    o_est <- est_omega_mu(o_result$beta, o_result$gamma, i_u, u_mat, p, nobs, "original")
    pos <- c(
      performance(o_result$beta, o_result$gamma, tb, mg),
      performance_supp(o_est$mu, o_est$omega, mu_true, omega_true)
    )

    metric_mat <- rbind(round(pgs, 3), round(pcs, 3), round(pos, 3))
    rownames(metric_mat) <- c("RegGMM", "cspine", "oracle")
    print(metric_mat)

    write(paste(pgs, collapse = ","), reggmm_csv, append = TRUE)
    write(paste(pcs, collapse = ","), cspine_csv, append = TRUE)
    write(paste(pos, collapse = ","), oracle_csv, append = TRUE)
  }
}
