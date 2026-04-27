suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(tibble))

performance <- function(tb_hat, mg_hat, tb_true, mg_true) {
  stopifnot(all(dim(tb_hat) == dim(tb_true)))
  stopifnot(all(dim(mg_hat) == dim(mg_true)))

  metrics <- c(
    "tpr", "fpr", "tpr_pop", "fpr_pop",
    "tpr_cov", "fpr_cov", "beta_err", "gamma_err"
  )
  stats <- numeric(length(metrics))
  names(stats) <- metrics

  stats["tpr"] <- sum(tb_hat != 0 & tb_true != 0) / sum(tb_true != 0)
  stats["fpr"] <- sum(tb_hat != 0 & tb_true == 0) / sum(tb_true == 0)

  stats["tpr_pop"] <- sum(tb_hat[, , 1] != 0 & tb_true[, , 1] != 0) /
    sum(tb_true[, , 1] != 0)
  stats["fpr_pop"] <- sum(tb_hat[, , 1] != 0 & tb_true[, , 1] == 0) /
    sum(tb_true[, , 1] == 0)

  stats["tpr_cov"] <- sum(tb_hat[, , -1] != 0 & tb_true[, , -1] != 0) /
    sum(tb_true[, , -1] != 0)
  stats["fpr_cov"] <- sum(tb_hat[, , -1] != 0 & tb_true[, , -1] == 0) /
    sum(tb_true[, , -1] == 0)

  stats["beta_err"] <- sqrt(sum((tb_hat - tb_true)^2))
  stats["gamma_err"] <- sqrt(sum((mg_hat - mg_true)^2))

  stats
}

performance_supp <- function(mu_hat, omega_hat, mu_true, omega_true) {
  mu_err <- sqrt(mean((mu_hat - mu_true)^2))
  omega_err <- sqrt(mean((omega_hat - omega_true)^2))

  p <- dim(omega_hat)[1]
  nobs <- dim(omega_hat)[3]
  off_diag <- !diag(p)

  tpr_fpr <- vapply(seq_len(nobs), function(k) {
    true_pos <- omega_true[, , k][off_diag] != 0
    hat_pos <- omega_hat[, , k][off_diag] != 0
    c(
      sum(hat_pos & true_pos) / max(sum(true_pos), 1L),
      sum(hat_pos & !true_pos) / max(sum(!true_pos), 1L)
    )
  }, numeric(2))

  c(
    mu_err = mu_err, omega_err = omega_err,
    tpr = mean(tpr_fpr[1, ]), fpr = mean(tpr_fpr[2, ])
  )
}

est_omega_mu <- function(beta, gamma, i_u, u_mat, p, nobs, method = c("original", "natural")) {
  method <- match.arg(method)
  beta_mat <- matrix(beta, nrow = p * p, ncol = ncol(i_u))
  omega_hat <- array(-(beta_mat %*% t(i_u)), dim = c(p, p, nobs))

  diag_idx <- cbind(rep(seq_len(p), nobs), rep(seq_len(p), nobs), rep(seq_len(nobs), each = p))
  omega_hat[diag_idx] <- 1

  min_eigs <- apply(omega_hat, 3, function(om) {
    min(eigen(om, symmetric = TRUE, only.values = TRUE)$values)
  })
  for (k in which(min_eigs <= 0)) {
    om <- omega_hat[, , k]
    inflate <- 1e-6
    while (min(eigen(om, symmetric = TRUE, only.values = TRUE)$values) <= 0) {
      diag(om) <- diag(om) + inflate
      inflate <- inflate * 10
    }
    omega_hat[, , k] <- om
  }

  mu_hat <- gamma %*% t(u_mat)
  if (method == "natural") {
    mu_hat <- vapply(seq_len(nobs), function(k) solve(omega_hat[, , k], mu_hat[, k]), numeric(p))
  }

  list(omega = omega_hat, mu = t(mu_hat))
}
