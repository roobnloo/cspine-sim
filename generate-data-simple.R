# Usage: Rscript generate-data-simple.R --p=25 --q=50 --nobs=200 --delta=1 [--nrep=100]
suppressPackageStartupMessages(library(MASS))

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

out_file <- file.path("data", sprintf("p%dq%d-n%d-d%.2f.rds", p, q, nobs, delta))

message(sprintf("Generating %d datasets with p=%d, q=%d, nobs=%d, delta=%.2f...", nrep, p, q, nobs, delta))

datasets <- vector("list", nrep)
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
    mu <- (1 - delta) * mg %*% u_mat[i, ] + delta * sigma %*% mg %*% u_mat[i, ]
    x_mat[i, ] <- mvrnorm(1, mu, sigma)
    mu_mat[i, ] <- mu
    omega_arr[, , i] <- omega
  }
  if (!valid) next

  repi <- repi + 1L
  datasets[[repi]] <- list(X = x_mat, U = u_mat, mu = mu_mat, omega = omega_arr)
  message(repi, " ", appendLF = FALSE)
}
saveRDS(datasets, out_file)
message("\nFinished.")
