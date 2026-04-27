library(sglssnal)
library(sparsegl)
library(Matrix)

source("impl/helpers.R")

#' Multi-task sparse-group lasso regression
#'
#' Stacks all p nodewise regression problems into one block-diagonal system
#' and tunes a single shared lambda via one joint 5-fold CV. This differs from
#' gmmreg_ssnal, which runs p independent CVs and selects p separate lambdas.
#'
#' @param responses n x p matrix of responses
#' @param u_cov n x q design matrix (without intercept)
#' @param alpha sparse-group lasso mixing parameter in [0,1]
#' @param nl1 number of lambda values in the cross-validation grid
#' @param lambda_factor ratio of smallest to largest lambda in the grid
#' @param nfolds number of cross-validation folds
#' @param stoptol_cv solver tolerance during cross-validation
#' @param stoptol_final solver tolerance for the final fit
#' @param deterministic_folds if TRUE use MATLAB-style mod-based fold assignment
#' @param verbose print solver progress for the final fit
#' @param maxit maximum solver iterations
#' @param mg_oracle optional p x q matrix of true means; skips stage 1 if given
#' @param standardize if TRUE standardize each design-matrix column per node
#' @param num_cores number of parallel cores for CV folds (via mclapply)
#' @return list with beta (p x p x (q+1)), gamma (p x q), and CV diagnostics
mtgmmreg_ssnal <- function(
    responses,
    u_cov,
    alpha,
    nl1,
    lambda_factor,
    nfolds = 5L,
    stoptol_cv = 1e-4,
    stoptol_final = 1e-6,
    deterministic_folds = TRUE,
    verbose = FALSE,
    maxit = 5000L,
    mg_oracle = NULL,
    standardize = TRUE,
    num_cores = 1L) {
  n <- nrow(responses)
  p <- ncol(responses)
  q <- ncol(u_cov)
  node_dim <- (p - 1L) * (q + 1L)

  # Stage 1: estimate mean matrix (identical to gmmreg_ssnal)
  if (!is.null(mg_oracle)) {
    ghat_mx <- mg_oracle
  } else {
    ghat_mx <- matrix(0, nrow = p, ncol = q)
    message("Stage 1: estimating mean matrix...")
    for (node in seq_len(p)) {
      result <- cv.sparsegl(
        u_cov, responses[, node], seq_len(q),
        asparse = 1, intercept = FALSE, standardize = TRUE
      )
      ghat_mx[node, ] <- as.numeric(coef(result, s = "lambda.min")[-1])
      message(node, " ", appendLF = FALSE)
    }
    message("\nFinished stage 1.")
  }

  z0 <- responses - u_cov %*% t(ghat_mx)

  i_u <- cbind(1, u_cov)
  inter_m <- Reduce(cbind, lapply(seq_len(q + 1L), \(k) z0 * i_u[, k]))

  # Build per-node design matrices and responses
  a_raw_list <- vector("list", p)
  a_std_list <- vector("list", p)
  b_list <- vector("list", p)
  col_sds_list <- vector("list", p)

  standardize_cols <- function(mat, means, sds) {
    sweep(sweep(mat, 2L, means), 2L, sds, "/")
  }

  for (j in seq_len(p)) {
    remove_cols <- j + (0:q) * p
    a_j <- inter_m[, -remove_cols, drop = FALSE]
    a_raw_list[[j]] <- a_j
    b_list[[j]] <- z0[, j] - mean(z0[, j])
    if (standardize) {
      col_means <- colMeans(a_j)
      col_sds <- pmax(apply(a_j, 2L, sd), .Machine$double.eps)
      a_std_list[[j]] <- standardize_cols(a_j, col_means, col_sds)
      col_sds_list[[j]] <- col_sds
    } else {
      a_std_list[[j]] <- a_j
    }
  }

  # Block-diagonal sparse design matrix and stacked response
  A_full <- Matrix::bdiag(a_std_list)
  b_full <- unlist(b_list)

  # Group structure: q+1 groups of size p*(p-1), coupling all nodes per covariate.
  # Block-diagonal has node-first ordering; grp_vec permutes to covariate-first
  # so that covariate k's columns from all p nodes form one contiguous group.
  grp_vec <- unlist(lapply(seq_len(q + 1L), function(k) {
    unlist(lapply(seq_len(p), function(j) {
      (j - 1L) * node_dim + (k - 1L) * (p - 1L) + seq_len(p - 1L)
    }))
  }))
  mt_grp_size <- p * (p - 1L)
  grp_starts <- 1L + mt_grp_size * seq.int(0L, q)
  grp_ends <- mt_grp_size * seq.int(1L, q + 1L)
  grp_idx <- rbind(grp_starts, grp_ends)
  pfgroup <- c(0, rep(1, q))

  # Lambda sequence derived from the combined system
  lam1_max <- max(abs(as.numeric(Matrix::crossprod(A_full, b_full))))
  lambda1_seq <- lam1_max * exp(seq(log(1), log(lambda_factor), length.out = nl1))

  # CV fold assignment: same obs-level folds applied to all p nodes
  if (deterministic_folds) {
    cut <- (seq.int(0L, n - 1L) %% nfolds) + 1L
  } else {
    cut <- cut(sample(seq_len(n)), nfolds, labels = FALSE)
  }

  message("Running MtGMMReg joint CV (", nfolds, " folds, ", num_cores, " cores)...")

  run_fold <- function(fold) {
    train_obs <- which(cut != fold)
    test_obs <- which(cut == fold)

    a_tr_list <- a_te_list <- vector("list", p)
    for (j in seq_len(p)) {
      a_j_raw <- a_raw_list[[j]]
      if (standardize) {
        fm <- colMeans(a_j_raw[train_obs, , drop = FALSE])
        fsd <- pmax(
          apply(a_j_raw[train_obs, , drop = FALSE], 2L, sd),
          .Machine$double.eps
        )
        a_tr_list[[j]] <- standardize_cols(a_j_raw[train_obs, ], fm, fsd)
        a_te_list[[j]] <- standardize_cols(a_j_raw[test_obs, ], fm, fsd)
      } else {
        a_tr_list[[j]] <- a_j_raw[train_obs, ]
        a_te_list[[j]] <- a_j_raw[test_obs, ]
      }
    }
    A_tr <- Matrix::bdiag(a_tr_list)
    A_te <- Matrix::bdiag(a_te_list)
    b_tr <- unlist(lapply(b_list, `[`, train_obs))
    b_te <- unlist(lapply(b_list, `[`, test_obs))

    fold_error <- matrix(0, nrow = nl1, ncol = length(alpha))
    for (k in seq_along(alpha)) {
      alpha_k <- alpha[k]
      warm_x <- warm_y <- warm_z <- NULL
      for (i in seq_len(nl1)) {
        result <- sglssnal::sglssnal(
          A        = A_tr,
          b        = b_tr,
          grp_vec  = grp_vec,
          grp_idx  = grp_idx,
          lambda   = lambda1_seq[i] / alpha_k,
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
        res <- A_te %*% result$x - b_te
        fold_error[i, k] <- fold_error[i, k] + sum(res^2)
      }
    }
    message("fold ", fold, " done")
    fold_error
  }

  fold_errors <- parallel::mclapply(seq_len(nfolds), run_fold, mc.cores = num_cores)
  cv_error <- Reduce(`+`, fold_errors)

  best_ij <- which(cv_error == min(cv_error), arr.ind = TRUE)[1L, , drop = FALSE]
  best_i <- best_ij[1L]
  best_k <- best_ij[2L]
  alpha_best <- alpha[best_k]
  lambda_r_best <- lambda1_seq[best_i] / alpha_best

  message(
    "Final fit (lambda=", signif(lambda1_seq[best_i], 3),
    ", alpha=", alpha_best, ")..."
  )
  result_final <- sglssnal::sglssnal(
    A        = A_full,
    b        = b_full,
    grp_vec  = grp_vec,
    grp_idx  = grp_idx,
    lambda   = lambda_r_best,
    alpha    = alpha_best,
    pfgroup  = pfgroup,
    stoptol  = stoptol_final,
    printyes = verbose,
    maxit    = maxit
  )
  x_coef <- result_final$x

  # Extract per-node coefficients and un-standardize
  beta_raw <- matrix(0, nrow = p, ncol = node_dim)
  for (j in seq_len(p)) {
    coef_j <- x_coef[((j - 1L) * node_dim + 1L):(j * node_dim)]
    beta_raw[j, ] <- if (standardize) coef_j / col_sds_list[[j]] else coef_j
  }

  # Organize into p x p x (q+1) array and symmetrize (same as gmmreg_ssnal)
  BB <- array(0, dim = c(p, p, q + 1L))
  for (j in seq_len(p)) {
    BB[setdiff(seq_len(p), j), j, ] <- matrix(beta_raw[j, ], nrow = p - 1L, ncol = q + 1L)
  }
  Bhat <- array(0, dim = c(p, p, q + 1L))
  for (k in seq_len(q + 1L)) {
    Bhat[, , k] <- symmetrize(BB[, , k])
  }

  list(
    beta          = Bhat,
    gamma         = ghat_mx,
    cv_lambda_idx = best_i,
    cv_lambda     = lambda1_seq[best_i],
    cv_alpha_idx  = best_k,
    cv_alpha      = alpha_best,
    lambda        = lambda1_seq,
    alpha         = alpha
  )
}
