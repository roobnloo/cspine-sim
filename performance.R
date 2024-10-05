suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(tibble))

performance <- function(fit, s, tb_true, mg_true, simple = FALSE) {
  stopifnot(all(dim(fit$beta) == dim(tb_true)))
  p <- dim(tb_true)[1]
  n <- nrow(s$X)
  stats <- vector(length = 9)
  names(stats) <- c(
    "tpr", "fpr", "tpr_pop", "fpr_pop", "tpr_cov", "fpr_cov",
    "beta_err", "omega_err", "mean_err"
  )
  stats["tpr"] <- sum(fit$beta != 0 & tb_true != 0) / sum(tb_true != 0)
  stats["fpr"] <- sum(fit$beta != 0 & tb_true == 0) / sum(tb_true == 0)

  stats["tpr_pop"] <- sum(fit$beta[, , 1] != 0 & tb_true[, , 1] != 0) / sum(tb_true[, , 1] != 0)
  stats["fpr_pop"] <- sum(fit$beta[, , 1] != 0 & tb_true[, , 1] == 0) / sum(tb_true[, , 1] == 0)

  stats["tpr_cov"] <- sum(fit$beta[, , -1] != 0 & tb_true[, , -1] != 0) / sum(tb_true[, , -1] != 0)
  stats["fpr_cov"] <- sum(fit$beta[, , -1] != 0 & tb_true[, , -1] == 0) / sum(tb_true[, , -1] == 0)

  beta_err <- 0
  for (i in 1:p) {
    beta_err <- beta_err + sqrt(sum((tb_true[i, -i, ] + fit$beta_raw[i, -i, ])^2))
  }
  stats["beta_err"] <- beta_err

  if (simple) {
    # Only return TPR, FPR, and beta_err
    return(stats[c("tpr", "fpr", "tpr_pop", "fpr_pop", "tpr_cov", "fpr_cov", "beta_err")])
  }

  omega_err <- 0
  mean_err <- 0
  iu <- cbind(1, s$U)

  for (i in 1:n) {
    omega <- apply(tb_true, c(1, 2), \(b) b %*% iu[i, ])
    diag(omega) <- 0

    p <- predict(fit, s$U[i, ])
    omhat <- p$precision
    diag(omhat) <- 0
    omega_err <- omega_err + sum((omega - omhat)^2) / n

    mean_err <- mean_err + sum((s$mumx[i, ] - p$mean)^2) / n
    # print(p$mean / n)
  }
  stats["omega_err"] <- omega_err
  stats["mean_err"] <- mean_err

  return(stats)
}

beta_viz <- function(beta_mx, title = "", limits = NULL, guides = T,
                     tileborder = T, fill_legend = T) {
  d <- nrow(beta_mx)
  beta0 <- as_tibble(cbind(
    expand.grid(rev(seq_len(d)), seq_len(d)),
    c(beta_mx)
  )) |>
    setNames(c("row", "col", "value"))

  tilecolor <- ifelse(tileborder, "gray30", "white")
  p <- ggplot(beta0, mapping = aes(x = col, y = row, fill = value)) +
    geom_tile() +
    scale_fill_gradient2(limits = limits) +
    coord_fixed() +
    labs(title = title) +
    theme_minimal()
  if (!fill_legend) {
    p <- p + guides(fill = "none")
  } else if (!guides) {
    p <- p + guides(x = "none", y = "none", fill = "none") +
      labs(x = NULL, y = NULL)
  }

  return(p)
}

beta_viz_compare <- function(computed, actual, cov_lbl, guides = T,
                             tileborder = T) {
  lim <- max(abs(c(as.numeric(computed), as.numeric(actual))))
  comp <- beta_viz(computed,
    title = paste("Estimated", cov_lbl),
    limits = c(-lim, lim),
    guides = guides,
    tileborder = tileborder,
    fill_legend = F
  )
  act <- beta_viz(actual,
    title = paste("True", cov_lbl),
    limits = c(-lim, lim),
    guides = guides,
    tileborder = tileborder
  )
  comp + act
}

beta_viz_list <- function(beta_list, cov_lbl, guides = T, tileborder = T) {
  lim <- max(abs(unlist(beta_list)))

  plots <- vector(mode = "list", length = length(beta_list))
  for (i in seq_along(beta_list)) {
    plots[[i]] <- beta_viz(beta_list[[i]],
      title = paste(names(beta_list)[i], cov_lbl),
      limits = c(-lim, lim),
      tileborder = tileborder
    ) +
      guides(x = "none", y = "none") +
      labs(x = NULL, y = NULL)
    if (i != length(beta_list)) {
      plots[[i]] <- plots[[i]] + guides(fill = "none")
    }
  }

  Reduce(`+`, plots)
}

cv_result_plot <- function(cv_result, node) {
  df <- expand.grid(
    lambda = seq_along(cv_result$lambdapath[, node]),
    sglmix = seq_along(cv_result$sglmixpath),
    gmix = seq_along(cv_result$gmixpath)
  )
  df$error <- as.numeric(cv_result$cv_mse[node, , , ])
  opt_lambda_idx <- cv_result$cv_lambda_idx[node]
  opt_gmix_idx <- cv_result$cv_gmix_idx[node]
  p <- ggplot(df) +
    geom_tile(
      color = "gray80",
      mapping = aes(x = gmix, y = lambda, fill = error)
    ) +
    geom_point(
      mapping = aes(x = x, y = y),
      data = data.frame(x = opt_gmix_idx, y = opt_lambda_idx),
      color = "tomato"
    ) +
    scale_x_continuous(
      breaks = seq_along(cv_result$gmixpath),
      labels = cv_result$gmixpath
    ) +
    scale_y_reverse(
      breaks = seq_along(cv_result$lambdapath[, node]),
      labels = round(cv_result$lambdapath[, node], 3)
    ) +
    scale_fill_gradient(low = "white", high = "black") +
    facet_wrap(~sglmix, nrow = 1) +
    guides(x = "none", y = "none", fill = "none") +
    coord_fixed() +
    theme_classic()
  return(p)
}

gamma_viz <- function(gamma_mx, title = "", limits = NULL) {
  d <- nrow(gamma_mx)
  p <- ncol(gamma_mx)
  gamma_tbl <- as_tibble(cbind(
    expand.grid(rev(seq_len(d)), seq_len(p)),
    c(gamma_mx)
  )) |>
    setNames(c("row", "col", "value"))

  ggplot(gamma_tbl, mapping = aes(x = col, y = row, fill = value)) +
    # geom_tile(color = "gray30") +
    geom_tile() +
    scale_fill_gradient2(limits = limits) +
    coord_fixed() +
    labs(title = title) +
    theme_minimal()
}

gamma_viz_compare <- function(computed, actual) {
  lim <- max(abs(c(as.numeric(computed), as.numeric(actual))))
  computed_plot <- gamma_viz(computed, "Computed gamma", c(-lim, lim)) +
    guides(fill = "none")
  actual_plot <- gamma_viz(actual, "Actual gamma", c(-lim, lim))
  computed_plot + actual_plot
}

gamma_viz_list <- function(gamma_list) {
  lim <- max(abs(unlist(gamma_list)))
  plots <- vector(mode = "list", length = length(gamma_list))
  for (i in seq_along(gamma_list)) {
    plots[[i]] <- gamma_viz(
      gamma_list[[i]],
      paste(names(gamma_list)[i], "gamma"), c(-lim, lim)
    ) +
      labs(x = NULL, y = NULL) +
      guides(x = "none", y = "none")
    if (i != length(gamma_list)) {
      plots[[i]] <- plots[[i]] + guides(fill = "none")
    }
  }
  Reduce(`+`, plots) + plot_layout(ncol = 1)
}
