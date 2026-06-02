library(sglssnal)
library(sparsegl)

source("impl/helpers.R")

predict.gmmreg <- function(fit, newcovar) {
  q <- dim(fit$gamma)[2]
  dim(newcovar) <- NULL
  if (length(newcovar) != q) {
    stop("Expected covariate vector of length ", q, ".")
  }
  omega <- apply(fit$beta, c(1, 2), \(b) b %*% c(1, newcovar))
  diag(omega) <- 1 / fit$sigma2
  mu <- fit$gamma %*% newcovar
  list(precision = omega, mean = mu)
}

#' Nodewise sparse-group lasso regression (faithful port of MATLAB GMMReg.m)
#'
#' @param z0 n x p matrix of centered responses
#' @param u_cov n x q design matrix (without intercept column)
#' @param alpha sparse-group lasso mixing parameter in [0,1]; alpha*lambda is
#'   the L1 penalty and (1-alpha)*lambda is the group L2 penalty
#' @param nl1 number of lambda values in the cross-validation grid
#' @param lambda_factor ratio of smallest to largest lambda in the grid
#' @param nfolds number of cross-validation folds
#' @param stoptol_cv solver convergence tolerance used during cross-validation
#' @param stoptol_final solver convergence tolerance for the final fit
#' @param deterministic_folds if TRUE use MATLAB-style mod-based fold
#'   assignment; if FALSE use random assignment
#' @param verbose print solver progress for the final fit of each node
#' @param maxit maximum solver iterations
#' @param mg_oracle optional p x q matrix of true means; if non-NULL, responses
#'   are centered using this matrix and stage 1 is skipped (oracle estimator)
#' @param standardize if TRUE (default) standardize each design-matrix column
#'   to mean zero and unit variance before the nodewise regression, ensuring
#'   beta coefficients are penalized on the same scale; returned coefficients
#'   are adjusted back to the original scale (divided by each column's sd)
#' @return p x (p-1)*(q+1) matrix; row j holds regression coefficients for
#'   predicting z0[,j] from all other variables interacted with u_cov
gmmreg_ssnal <- function(
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
    mg_oracle = NULL,
    standardize = TRUE) {
  n <- nrow(responses)
  p <- ncol(responses)
  q <- ncol(u_cov)
  dim_coef <- (p - 1L) * (q + 1L)

  # Stage 1: estimate mean matrix via nodewise lasso (cv.sparsegl)
  # If mg_oracle is provided, use it directly (oracle estimator where mg is known)
  if (!is.null(mg_oracle)) {
    ghat_mx <- mg_oracle
  } else {
    ghat_mx <- matrix(0, nrow = p, ncol = q)
    message("Stage 1: estimating mean matrix...")
    nodewise_gamma <- function(node) {
      result <- cv.sparsegl(
        u_cov, responses[, node], seq_len(q),
        asparse = 1, intercept = FALSE, standardize = TRUE
      )
      message(node, " ", appendLF = FALSE)
      as.numeric(coef(result, s = "lambda.min")[-1])
    }

    if (num_cores > 1L) {
      step1 <- parallel::mclapply(seq_len(p), nodewise_gamma, mc.cores = num_cores)
    } else {
      step1 <- lapply(seq_len(p), nodewise_gamma)
    }
    message("\nFinished stage 1.")

    for (node in seq_len(p)) {
      ghat_mx[node, ] <- step1[[node]]
    }
    rm(step1)
    gc()
  }

  # Center responses
  z0 <- responses - u_cov %*% t(ghat_mx)

  # Expanded design matrix: n x p*(q+1)
  # Block k: z0 * i_u[,k] elementwise, where i_u = [1, u_cov]
  i_u <- cbind(1, u_cov)
  inter_m <- Reduce(cbind, lapply(seq_len(q + 1L), \(k) z0 * i_u[, k]))

  # Group structure: q+1 consecutive groups of size p-1
  # Group 1 is the intercept block — no group penalty (pfgroup = 0)
  grp_vec <- seq_len(dim_coef)
  grp_starts <- 1L + (p - 1L) * seq.int(0L, q)
  grp_ends <- (p - 1L) * seq.int(1L, q + 1L)
  grp_idx <- rbind(grp_starts, grp_ends)
  pfgroup <- c(0, rep(1, q))

  # CV fold assignment: deterministic mode matches MATLAB mod(0:n-1, 5)+1
  if (deterministic_folds) {
    cut <- (seq.int(0L, n - 1L) %% nfolds) + 1L
  } else {
    cut <- cut(sample(seq_len(n)), nfolds, labels = FALSE)
  }

  fit_node <- function(j) {
    # Remove columns that predict variable j: j, p+j, 2p+j, ..., q*p+j
    remove_cols <- j + (0:q) * p
    a <- inter_m[, -remove_cols, drop = FALSE]
    b_j <- z0[, j] - mean(z0[, j])

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

    # Lambda sequence: lam1_max * exp(linspace(log(1), log(lambda_factor), nl1))
    lam1_max <- max(abs(crossprod(a_fit, b_j)))
    lambda1_seq <- lam1_max *
      exp(seq(log(1), log(lambda_factor), length.out = nl1))

    # Lipschitz constant: largest eigenvalue of a*a', tol=1e-3 matches MATLAB
    # aat_fun <- function(x) as.numeric(a %*% crossprod(a, x))
    # lip <- RSpectra::eigs_sym(
    #   aat_fun,
    #   n = n, k = 1L, which = "LA",
    #   opts = list(retvec = FALSE, tol = 1e-3)
    # )$values

    # CV: outer loop over folds, middle over alpha, inner over lambda.
    # Warm starts reset per alpha path within each fold (matches MATLAB Cerror(i,k)).
    cv_error <- matrix(0, nrow = nl1, ncol = length(alpha))

    for (fold in seq_len(nfolds)) {
      train_idx <- which(cut != fold)
      test_idx <- which(cut == fold)
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
            # Lip      = lip,
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

    # Select first minimiser (matches MATLAB [id1,id2]=find(Cerror==min(min(Cerror))))
    best_ij <- which(cv_error == min(cv_error), arr.ind = TRUE)[1L, , drop = FALSE]
    best_i <- best_ij[1L]
    best_k <- best_ij[2L]
    alpha_best <- alpha[best_k]
    lambda_r_best <- lambda1_seq[best_i] / alpha_best

    # Final estimation on full data with tighter tolerance
    result_final <- sglssnal::sglssnal(
      A        = a_fit,
      b        = b_j,
      grp_vec  = grp_vec,
      grp_idx  = grp_idx,
      lambda   = lambda_r_best,
      alpha    = alpha_best,
      pfgroup  = pfgroup,
      stoptol  = stoptol_final,
      # Lip      = lip,
      printyes = verbose,
      maxit    = maxit
    )

    coef_std <- result_final$x
    # Adjust back: divide by sd so coefficients apply to the original unscaled a.
    coef_j <- if (standardize) coef_std / a_col_sds else coef_std

    message(j, " ", appendLF = FALSE)
    list(
      x             = coef_j,
      cv_lambda_idx = best_i,
      cv_lambda     = lambda1_seq[best_i],
      cv_alpha_idx  = best_k,
      cv_alpha      = alpha_best,
      lambda        = lambda1_seq,
      cv_error      = cv_error
    )
  }

  message("Running GMMReg nodewise regressions...")
  rows <- parallel::mclapply(seq_len(p), fit_node, mc.cores = num_cores)
  message("\nFinished regressions.")

  nl2 <- length(alpha)
  cv_lambda_idx <- integer(p)
  cv_lambda <- numeric(p)
  cv_alpha_idx <- integer(p)
  cv_alpha <- numeric(p)
  lambda_path <- matrix(0, nrow = p, ncol = nl1)
  cvm <- array(0, dim = c(nl1, nl2, p))
  for (j in seq_len(p)) {
    cv_lambda_idx[j] <- rows[[j]]$cv_lambda_idx
    cv_lambda[j] <- rows[[j]]$cv_lambda
    cv_alpha_idx[j] <- rows[[j]]$cv_alpha_idx
    cv_alpha[j] <- rows[[j]]$cv_alpha
    lambda_path[j, ] <- rows[[j]]$lambda
    cvm[, , j] <- rows[[j]]$cv_error
  }
  beta <- t(do.call(cbind, lapply(rows, `[[`, "x")))

  # Organize into p x p x (q+1) array: column j gets the (p-1) x (q+1) block
  # of regression coefficients for node j, placed at rows != j
  BB <- array(0, dim = c(p, p, q + 1L))
  for (j in seq_len(p)) {
    BB[setdiff(seq_len(p), j), j, ] <- matrix(beta[j, ], nrow = p - 1L, ncol = q + 1L)
  }

  # Symmetrize each slice with the "and" (min-magnitude) rule
  Bhat <- array(0, dim = c(p, p, q + 1L))
  for (k in seq_len(q + 1L)) {
    Bhat[, , k] <- symmetrize(BB[, , k])
  }
  list(
    beta          = Bhat,
    gamma         = ghat_mx,
    cv_lambda_idx = cv_lambda_idx,
    cv_lambda     = cv_lambda,
    cv_alpha_idx  = cv_alpha_idx,
    cv_alpha      = cv_alpha,
    lambda        = lambda_path,
    alpha         = alpha,
    cvm           = cvm
  )
}
