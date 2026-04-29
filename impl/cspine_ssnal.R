library(sglssnal)

source("impl/helpers.R")

#' Nodewise sparse-group lasso regression for the CSPINE model
#'
#' @param responses n x p matrix of responses
#' @param u_cov n x q design matrix (without intercept column)
#' @param alpha sparse-group lasso mixing parameter in [0,1]; alpha*lambda is
#'   the L1 penalty and (1-alpha)*lambda is the group L2 penalty
#' @param nl1 number of lambda values in the cross-validation grid
#' @param lambda_factor ratio of smallest to largest lambda in the grid
#' @param nfolds number of cross-validation folds
#' @param stoptol_cv solver convergence tolerance used during cross-validation
#' @param stoptol_final solver convergence tolerance for the final fit
#' @param deterministic_folds if TRUE use MATLAB-style mod-based fold assignment
#' @param verbose print solver progress for the final fit of each node
#' @param maxit maximum solver iterations
#' @param standardize if TRUE (default) standardize each design-matrix column
#'   to mean zero and unit variance before the nodewise regression, ensuring
#'   gamma and beta are penalized on the same scale; returned coefficients are
#'   adjusted back to the original scale (divided by each column's sd)
#' @return list with:
#'   beta: p x p x (q+1) symmetrized precision-matrix coefficients (-beta_hat/sigma2)
#'   beta_raw: p x p x (q+1) symmetrized raw regression coefficients
#'   gamma: p x q mean regression coefficients
#'   sigma2: length-p vector of estimated residual variances
cspine_ssnal <- function(
    responses,
    u_cov,
    alpha,
    nl1,
    lambda_factor,
    nfolds = 5L,
    stoptol_cv = 1e-4,
    stoptol_final = 1e-6,
    deterministic_folds = TRUE,
    num_cores = 1L,
    verbose = FALSE,
    maxit = 5000L,
    standardize = TRUE) {
  n <- nrow(responses)
  p <- ncol(responses)
  q <- ncol(u_cov)
  dim_beta <- (p - 1L) * (q + 1L)
  dim_coef <- q + dim_beta

  # Expanded design: n x p*(q+1); uses original responses (no centering)
  inter_m <- intxmx(responses, u_cov)

  # Group structure: 2q+1 groups
  #   Groups 1..q:     gamma variables, size 1 each,  pfgroup=0 (pure L1)
  #   Group  q+1:      beta_0 variables, size p-1,    pfgroup=0 (pure L1)
  #   Groups q+2..2q+1: beta_1..beta_q, size p-1 each, pfgroup=1 (L1 + group L2)
  grp_vec <- seq_len(dim_coef)
  grp_starts <- c(seq_len(q), q + 1L + (p - 1L) * seq.int(0L, q))
  grp_ends <- c(seq_len(q), q + (p - 1L) * seq.int(1L, q + 1L))
  grp_idx <- rbind(grp_starts, grp_ends)
  pfgroup <- c(rep(0L, q + 1L), rep(1L, q))
  # pfgroup <- c(rep(1L, q), 0L, rep(1L, q))

  if (deterministic_folds) {
    fold_id <- (seq.int(0L, n - 1L) %% nfolds) + 1L
  } else {
    fold_id <- cut(sample(seq_len(n)), nfolds, labels = FALSE)
  }

  fit_node <- function(j) {
    # Design: [u_cov | X_{-j} | u1*X_{-j} | ... | uq*X_{-j}]
    remove_cols <- j + (0:q) * p
    a <- cbind(u_cov, inter_m[, -remove_cols, drop = FALSE])
    b_j <- responses[, j]

    standardize_cols <- function(mat, means, sds) {
      sweep(sweep(mat, 2L, means), 2L, sds, "/")
    }

    if (standardize) {
      a_col_means <- colMeans(a)
      a_col_sds <- pmax(apply(a, 2L, sd), .Machine$double.eps)
      a_fit <- standardize_cols(a, a_col_means, a_col_sds)
    } else {
      a_fit <- a
    }

    lam1_max <- max(abs(crossprod(a_fit, b_j)))
    lambda1_seq <- lam1_max *
      exp(seq(log(1), log(lambda_factor), length.out = nl1))

    cv_error <- matrix(0, nrow = nl1, ncol = length(alpha))
    for (fold in seq_len(nfolds)) {
      train_idx <- which(fold_id != fold)
      test_idx <- which(fold_id == fold)
      if (standardize) {
        fold_means <- colMeans(a[train_idx, , drop = FALSE])
        fold_sds <- pmax(
          apply(a[train_idx, , drop = FALSE], 2L, sd),
          .Machine$double.eps
        )
        a_train <- standardize_cols(a[train_idx, , drop = FALSE], fold_means, fold_sds)
        a_test <- standardize_cols(a[test_idx, , drop = FALSE], fold_means, fold_sds)
      } else {
        a_train <- a[train_idx, , drop = FALSE]
        a_test <- a[test_idx, , drop = FALSE]
      }
      b_train <- b_j[train_idx]
      b_test <- b_j[test_idx]

      for (k in seq_along(alpha)) {
        alpha_k <- alpha[k]
        warm_x <- NULL
        warm_y <- NULL
        warm_z <- NULL
        for (i in seq_len(nl1)) {
          lambda_r <- lambda1_seq[i] / alpha_k
          result <- sglssnal::sglssnal(
            A        = a_train,
            b        = b_train,
            grp_vec  = grp_vec,
            grp_idx  = grp_idx,
            lambda   = lambda_r,
            alpha    = alpha_k,
            pfgroup  = pfgroup,
            stoptol  = stoptol_cv,
            printyes = FALSE,
            maxit    = maxit,
            x0       = warm_x,
            y0       = warm_y,
            z0       = warm_z
          )
          warm_x <- result$x
          warm_y <- result$y
          warm_z <- result$z
          res <- a_test %*% result$x - b_test
          cv_error[i, k] <- cv_error[i, k] + sum(res^2)
        }
      }
    }

    best_ij <- which(cv_error == min(cv_error), arr.ind = TRUE)[1L, , drop = FALSE]
    best_i <- best_ij[1L]
    best_k <- best_ij[2L]
    alpha_best <- alpha[best_k]
    lambda_r_best <- lambda1_seq[best_i] / alpha_best

    result_final <- sglssnal::sglssnal(
      A        = a_fit,
      b        = b_j,
      grp_vec  = grp_vec,
      grp_idx  = grp_idx,
      lambda   = lambda_r_best,
      alpha    = alpha_best,
      pfgroup  = pfgroup,
      stoptol  = stoptol_final,
      printyes = verbose,
      maxit    = maxit
    )

    coef_std <- result_final$x
    # Residuals from the standardized fitted values (avoids the centering
    # offset that would arise from using a %*% coef_j directly).
    resid_j <- b_j - as.numeric(a_fit %*% coef_std)
    # Adjust back: divide by sd so coef_j applies to the original unscaled a.
    coef_j <- if (standardize) coef_std / a_col_sds else coef_std
    df_j <- result_final$info$nnz
    sigma2_j <- sum(resid_j^2) / max(n - df_j, 1L)

    message(j, " ", appendLF = FALSE)
    list(
      gamma         = coef_j[seq_len(q)],
      beta          = coef_j[seq.int(q + 1L, dim_coef)],
      sigma2        = sigma2_j,
      cv_lambda_idx = best_i,
      cv_lambda     = lambda1_seq[best_i],
      cv_alpha_idx  = best_k,
      cv_alpha      = alpha_best,
      lambda        = lambda1_seq,
      cv_error      = cv_error
    )
  }

  message("Running CSPINE nodewise regressions...")
  node_fits <- parallel::mclapply(seq_len(p), fit_node, mc.cores = num_cores)
  message("\nFinished regressions.")

  nl2 <- length(alpha)
  ghat_mx <- matrix(0, nrow = p, ncol = q)
  sigma2 <- numeric(p)
  cv_lambda_idx <- integer(p)
  cv_lambda <- numeric(p)
  cv_alpha_idx <- integer(p)
  cv_alpha <- numeric(p)
  lambda_path <- matrix(0, nrow = p, ncol = nl1)
  cvm <- array(0, dim = c(nl1, nl2, p))
  BB_raw <- array(0, dim = c(p, p, q + 1L))
  BB_scaled <- array(0, dim = c(p, p, q + 1L))

  for (j in seq_len(p)) {
    fit_j <- node_fits[[j]]
    ghat_mx[j, ] <- fit_j$gamma
    sigma2[j] <- fit_j$sigma2
    cv_lambda_idx[j] <- fit_j$cv_lambda_idx
    cv_lambda[j] <- fit_j$cv_lambda
    cv_alpha_idx[j] <- fit_j$cv_alpha_idx
    cv_alpha[j] <- fit_j$cv_alpha
    lambda_path[j, ] <- fit_j$lambda
    cvm[, , j] <- fit_j$cv_error
    beta_j_mx <- matrix(fit_j$beta, nrow = p - 1L, ncol = q + 1L)
    rows_j <- setdiff(seq_len(p), j)
    BB_raw[rows_j, j, ] <- beta_j_mx
    BB_scaled[rows_j, j, ] <- -beta_j_mx / fit_j$sigma2
  }

  Bhat <- array(0, dim = c(p, p, q + 1L))
  Bhat_raw <- array(0, dim = c(p, p, q + 1L))
  for (k in seq_len(q + 1L)) {
    Bhat[, , k] <- symmetrize(BB_scaled[, , k])
    Bhat_raw[, , k] <- symmetrize(BB_raw[, , k])
  }

  list(
    beta = Bhat,
    beta_raw = Bhat_raw,
    gamma = ghat_mx,
    sigma2 = sigma2,
    cv_lambda_idx = cv_lambda_idx,
    cv_lambda = cv_lambda,
    cv_alpha_idx = cv_alpha_idx,
    cv_alpha = cv_alpha,
    lambda = lambda_path,
    alpha = alpha,
    cvm = cvm
  )
}
