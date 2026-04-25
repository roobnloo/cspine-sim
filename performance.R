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
