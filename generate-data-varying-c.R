# Usage: Rscript generate-data-varying-c.R --c_scale=1 --nobs=200 [--nrep=30]
suppressPackageStartupMessages(library(MASS))

args <- commandArgs(trailingOnly = TRUE)

parse_arg <- function(args, key, default = NULL) {
  pat <- paste0("^--", key, "=(.+)$")
  m <- regmatches(args, regexpr(pat, args, perl = TRUE))
  if (length(m) == 0L) default else sub(pat, "\\1", m)
}

c_val <- as.numeric(parse_arg(args, "c_scale"))
nobs <- as.integer(parse_arg(args, "nobs"))
nrep <- as.integer(parse_arg(args, "nrep", default = 30L))
for (req in c("c_val", "nobs")) {
  if (is.na(get(req))) stop("Required argument missing: --", req)
}

true_param <- readRDS(file.path("data", "coef_p25q10_vary_c.rds"))
tB_all <- true_param$tB_all
mg <- true_param$mg
c_scale <- true_param$c_scale

c_idx <- which(abs(c_scale - c_val) < 1e-9)
if (length(c_idx) == 0L) {
  stop("--c_scale=", c_val, " not found in precomputed levels: ", paste(c_scale, collapse = ", "))
}
c_idx <- c_idx[1L]

tb <- tB_all[, , , c_idx]
p <- dim(tb)[1]
q <- dim(tb)[3] - 1L

out_file <- file.path("data", sprintf("vary-c-c%.2f-n%d.rds", c_val, nobs))
message(sprintf(
  "Generating %d datasets (p=%d, q=%d, nobs=%d, c_scale=%.2f)...",
  nrep, p, q, nobs, c_val
))

datasets <- vector("list", nrep)
for (repi in seq_len(nrep)) {
  success <- FALSE
  for (attempt in seq_len(50L)) {
    set.seed(42 + (repi - 1L) * 10L + attempt)

    u_mat <- matrix(sample(c(0L, 1L), nobs * q, replace = TRUE), nobs, q)
    cont_idx <- sample(seq_len(q), q / 2, replace = FALSE)
    u_mat[, cont_idx] <- matrix(runif(nobs * q / 2), nobs, q / 2)
    u_mat <- apply(u_mat, 2, scale)
    i_u <- cbind(1, u_mat)

    x_mat <- matrix(0, nobs, p)
    mu_mat <- matrix(0, nobs, p)
    omega_arr <- array(0, dim = c(p, p, nobs))
    valid <- TRUE
    for (i in seq_len(nobs)) {
      omega <- -apply(tb, c(1, 2), function(b) b %*% i_u[i, ])
      diag(omega) <- 1
      if (any(eigen(omega, symmetric = TRUE, only.values = TRUE)$values <= 0)) {
        valid <- FALSE
        break
      }
      sigma <- solve(omega)
      mu <- mg %*% u_mat[i, ]
      x_mat[i, ] <- mvrnorm(1, mu, sigma)
      mu_mat[i, ] <- mu
      omega_arr[, , i] <- omega
    }
    if (!valid) next

    datasets[[repi]] <- list(X = x_mat, U = u_mat, mu = mu_mat, omega = omega_arr)
    success <- TRUE
    break
  }
  if (!success) stop("Failed to generate valid dataset for rep ", repi, " after 50 attempts.")
  message(repi, " ", appendLF = FALSE)
}
saveRDS(datasets, out_file)
message("\nFinished.")
