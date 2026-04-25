# Usage: Rscript generate-data-simple.R --p=25 --q=50 --nobs=200 --delta=1 [--nrep=100]
suppressPackageStartupMessages(library(MASS))
RhpcBLASctl::blas_set_num_threads(1)
RhpcBLASctl::omp_set_num_threads(1)
source("impl/cspine_ssnal.R")
source("impl/gmmreg_ssnal.R")
source("performance.R")

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
for (req in c("p", "q", "nobs", "delta")) {
  if (is.na(get(req))) stop("Required argument missing: --", req)
}
if (delta < 0 || delta > 1) stop("--delta must be between 0 and 1 inclusive")

true_param <- readRDS(file.path("data", sprintf("coef_p%dq%d.rds", p, q)))
tb <- true_param$tb
mg <- true_param$mg

data_dir <- file.path("data", sprintf("p%dq%d-n%d-d%.2f", p, q, nobs, delta))
dir.create(data_dir, showWarnings = FALSE)

message(sprintf("Generating %d datasets with p=%d, q=%d, nobs=%d...", nrep, p, q, nobs))

repi <- 0L
attempt <- 0L
while (repi < nrep) {
  attempt <- attempt + 1L
  set.seed(42 + attempt)

  u_mat <- matrix(sample(c(0L, 1L), nobs * q, replace = TRUE), nobs, q)
  cont_idx <- sample(seq_len(q), q / 2, replace = FALSE)
  u_mat[, cont_idx] <- matrix(runif(nobs * q / 2), nobs, q / 2)
  u_mat <- apply(u_mat, 2, scale)
  i_u <- cbind(1, u_mat)

  # X^(i) ~ N(Gamma U^(i), Omega(U^(i))^{-1})
  # tb stores regression betas; Omega = -tB %*% iU
  x_mat <- matrix(0, nobs, p)
  valid <- TRUE
  for (i in seq_len(nobs)) {
    omega <- -apply(tb, c(1, 2), function(b) b %*% i_u[i, ])
    diag(omega) <- 1
    if (any(eigen(omega, symmetric = TRUE, only.values = TRUE)$values <= 0)) {
      valid <- FALSE
      break
    }
    sigma <- solve(omega)
    mu <- delta * mg %*% u_mat[i, ] + (1 - delta) * sigma %*% mg %*% u_mat[i, ]
    x_mat[i, ] <- mvrnorm(1, mu, sigma)
  }
  if (!valid) next

  repi <- repi + 1L
  write.table(as.data.frame(x_mat),
    sep = ",",
    file.path(data_dir, sprintf("X_%d.csv", repi)),
    row.names = FALSE, col.names = FALSE
  )
  write.table(as.data.frame(u_mat),
    sep = ",",
    file.path(data_dir, sprintf("U_%d.csv", repi)),
    row.names = FALSE, col.names = FALSE
  )
  message(repi, " ", appendLF = FALSE)
}
message("\nFinished.")
