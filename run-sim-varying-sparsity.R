# Usage: Rscript run-sim-varying-sparsity.R --p=25 --q=10 --nobs=300 --sparsity_setting=beta [--nrep=100]
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

p <- as.integer(parse_arg(args, "p"))
q <- as.integer(parse_arg(args, "q"))
nobs <- as.integer(parse_arg(args, "nobs"))
sparsity_setting <- parse_arg(args, "sparsity_setting")
nrep <- as.integer(parse_arg(args, "nrep", default = 100))
for (req in c("p", "q", "nobs")) {
  if (is.na(get(req))) stop("Required argument missing: --", req)
}
if (is.null(sparsity_setting)) stop("Required argument missing: --sparsity_setting")
stopifnot(sparsity_setting %in% c("beta", "gamma"))

model <- "natural"

data_file <- if (sparsity_setting == "gamma") {
  sprintf("data/p%dq%d-n%d-%s-varying-sparsity-gamma.rds", p, q, nobs, model)
} else {
  sprintf("data/p%dq%d-n%d-%s-varying-sparsity.rds", p, q, nobs, model)
}
all_data <- readRDS(data_file)

metrics <- c(
  "tpr", "fpr", "tpr_pop", "fpr_pop", "tpr_cov", "fpr_cov", "beta_err", "rel_beta_err", "gamma_err",
  "mu_err", "omega_err", "omega_tpr", "omega_fpr", paste0("s_", sparsity_setting)
)

setting_str <- sprintf("p%dq%d-n%d-%s-varying-sparsity-%s", p, q, nobs, model, sparsity_setting)
out_dir <- file.path("out", setting_str)
dir.create(out_dir, showWarnings = FALSE)

message(sprintf("Settings: p=%d, q=%d, nobs=%d, sparsity_setting=%s, nrep=%d", p, q, nobs, sparsity_setting, nrep))
message("Output directory: ", out_dir)

reggmm_csv <- file.path(out_dir, "result-RegGMM.csv")
cspine_csv <- file.path(out_dir, "result-cspine.csv")
write(paste(metrics, collapse = ","), reggmm_csv)
write(paste(metrics, collapse = ","), cspine_csv)

for (i in seq_len(nrep)) {
  message("Rep ", i)
  s <- all_data[[i]]
  x_mat <- s$X
  u_mat <- s$U
  mu_true <- s$mu
  omega_true <- s$omega
  i_u <- cbind(1, u_mat)
  sparsity <- if (sparsity_setting == "gamma") sum(abs(s$mg) > 0) else sum(abs(s$tb) > 0)

  tictoc::tic()
  g_result <- gmmreg_ssnal(x_mat, u_mat, alpha = 0.75, nl1 = 100, lambda_factor = 0.1, num_cores = 25)
  tictoc::toc()
  g_est <- est_omega_mu(g_result$beta, g_result$gamma, i_u, u_mat, p, nobs, "original")
  pgs <- c(
    performance(g_result$beta, g_result$gamma, s$tb, s$mg),
    performance_supp(g_est$mu, g_est$omega, mu_true, omega_true),
    sparsity
  )

  tictoc::tic()
  c_result <- cspine_ssnal(x_mat, u_mat, alpha = 0.75, nl1 = 100, lambda_factor = 0.1, num_cores = 25)
  tictoc::toc()
  c_est <- est_omega_mu(c_result$beta_raw, c_result$gamma, i_u, u_mat, p, nobs, "natural")
  pcs <- c(
    performance(c_result$beta_raw, c_result$gamma, s$tb, s$mg),
    performance_supp(c_est$mu, c_est$omega, mu_true, omega_true),
    sparsity
  )

  metric_mat <- rbind(round(pgs, 3), round(pcs, 3))
  rownames(metric_mat) <- c("RegGMM", "cspine")
  print(metric_mat)
  write(paste(pgs, collapse = ","), reggmm_csv, append = TRUE)
  write(paste(pcs, collapse = ","), cspine_csv, append = TRUE)
}
