suppressPackageStartupMessages(library(MASS))

# ---- Parameters ----
p <- 25L
q <- 50L
n <- 200L
n_rep <- 100L
c_vals <- round(sqrt(c(0.01, seq(0.025, 0.125, 0.025)) / 0.083), 2)
c_vals <- round(sqrt(0.01 / 0.083), 2)
pd_tol <- 1e-6
pd_fail_thresh <- 0.01
master_seed <- 8472L

# ---- Load fixed true parameters ----
true_param <- readRDS(file.path("data", sprintf("coef_p%dq%d.rds", p, q)))
tb <- true_param$tb
mg <- true_param$mg

# ---- Pre-compute vectorized Omega pieces (reused across all c and reps) ----
# Bmat[,h] = vec(B_h), so contrib = Bmat %*% t(U_scaled) gives vec(Omega_i) contributions
Bmat <- matrix(tb[, , 2:(q + 1)], p * p, q) # (p^2) x q
B0_vec <- as.vector(tb[, , 1]) # p^2, zero diagonal
diag_idx <- seq(1L, p * p, by = p + 1L) # indices of diagonal in vec

# ---- Seed matrix: one seed per (replicate, c-level), all drawn from master seed ----
set.seed(master_seed)
all_seeds <- matrix(
  sample.int(1e8, length(c_vals) * n_rep),
  nrow = n_rep, ncol = length(c_vals)
)

# ---- Helper: generate one replicate ----
gen_replicate <- function(c_val, seed) {
  set.seed(seed)

  U <- matrix(0, n, q)
  U[, seq_len(q / 2)] <- matrix(sample(0:1, n * (q / 2), replace = TRUE), n, q / 2)
  U[, seq(q / 2 + 1L, q)] <- apply(matrix(runif(n * (q / 2)), n, q / 2), 2, scale)
  U_scaled <- c_val * U

  # Vectorized Omega construction
  contrib <- Bmat %*% t(U_scaled) # (p^2) x n
  omega_vecs <- -(B0_vec + contrib) # (p^2) x n
  omega_vecs[diag_idx, ] <- 1 # set diagonal = 1

  pd_fail <- logical(n)
  X <- matrix(NA_real_, n, p)
  snr_num <- numeric(n)
  snr_den <- numeric(n)

  for (i in seq_len(n)) {
    omega_i <- matrix(omega_vecs[, i], p, p)
    ev_min <- min(eigen(omega_i, symmetric = TRUE, only.values = TRUE)$values)
    if (ev_min <= pd_tol) {
      pd_fail[i] <- TRUE
      next
    }
    sigma_i <- solve(omega_i)
    mu_i <- drop(mg %*% U_scaled[i, ])
    X[i, ] <- mvrnorm(1, mu_i, sigma_i)
    snr_num[i] <- sum(mu_i^2)
    snr_den[i] <- sum(diag(sigma_i))
  }

  pd_fail_rate <- mean(pd_fail)
  if (pd_fail_rate > pd_fail_thresh) {
    return(list(failed = TRUE, pd_fail_rate = pd_fail_rate))
  }

  valid <- !pd_fail
  snr <- mean(snr_num[valid]) / mean(snr_den[valid])

  list(
    U            = U_scaled,
    X            = X,
    snr          = snr,
    c            = c_val,
    seed         = seed,
    pd_fail_rate = pd_fail_rate
  )
}

# ---- Main loop ----
dir.create("data", showWarnings = FALSE)

avg_snrs <- numeric(length(c_vals))
for (ci in seq_along(c_vals)) {
  c_val <- c_vals[ci]
  message(sprintf("\n--- c = %.2f ---", c_val))
  reps <- vector("list", n_rep)
  skipped <- FALSE

  for (r in seq_len(n_rep)) {
    res <- gen_replicate(c_val, all_seeds[r, ci])
    if (isTRUE(res$failed)) {
      warning(sprintf(
        "c=%.2f rep=%d: PD failure rate %.1f%% > 1%%, skipping this c level",
        c_val, r, 100 * res$pd_fail_rate
      ))
      skipped <- TRUE
      break
    }
    reps[[r]] <- res
    cat(sprintf("  rep %2d  snr=%.4f  pd_fail=%.4f\n", r, res$snr, res$pd_fail_rate))
  }

  if (!skipped) {
    avg_snr <- mean(sapply(reps, `[[`, "snr"))
    avg_snrs[ci] <- avg_snr
    message(sprintf("  avg snr = %.4f", avg_snr))
    c_str <- gsub("\\.", "p", as.character(c_val))
    outpath <- file.path("data", sprintf("p%dq%d-n%d-varying-snr-c%s.rds", p, q, n, c_str))
    # saveRDS(reps, outpath)
    message("Saved ", outpath)
  }
}
