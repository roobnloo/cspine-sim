library(sparsegl)

intxmx <- function(X, U) {
  q <- ncol(U)
  iU <- cbind(1, U)
  result <- lapply(seq_len(q + 1), \(j) {
    X * iU[, j]
  })
  result <- Reduce(cbind, result)
  result
}

symmetrize <- function(mx, rule = "and") {
  if (rule == "and") {
    result <- mx * (abs(mx) < t(abs(mx))) + t(mx) * (t(abs(mx)) < abs(mx))
  } else {
    result <- mx * (abs(mx) >= t(abs(mx))) + t(mx) * (t(abs(mx)) >= abs(mx))
  }
  return(result)
}

predict.gmmreg <- function(fit, newcovar) {
  q <- dim(fit$gamma)[2]
  dim(newcovar) <- NULL
  if (length(newcovar) != q) {
    stop("Expected covariate vector of length ", q, ".")
  }
  omega <- apply(fit$beta, c(1, 2), \(b) b %*% c(1, newcovar))
  diag(omega) <- 1 / fit$sigma2

  mu <- fit$gamma %*% newcovar
  return(list(precision = omega, mean = mu))
}


gmmreg <- function(
    responses, covariates, asparse = seq(0.1, 1, by = 0.1),
    nlambda = 100, lam_max = NULL, lambda_factor = 1e-4,
    nfolds = 5, verbose = FALSE, ncores = 1, skip_stage1 = FALSE) {
  stopifnot(
    is.matrix(responses), is.matrix(covariates),
    nrow(responses) == nrow(covariates),
    all(asparse >= 0), all(asparse <= 1)
  )
  p <- ncol(responses)
  q <- ncol(covariates)
  n <- nrow(responses)
  nvars <- (p - 1) * (q + 1)
  nasparse <- length(asparse)

  # Estimate mean matrix
  ghat_mx <- matrix(0, nrow = p, ncol = q)
  g0 <- numeric(p)

  if (!skip_stage1) {
    message("Stage 1")

    nodewise_gamma <- function(node) {
      result <- cv.sparsegl(covariates, responses[, node], seq_len(q),
        asparse = 1, intercept = TRUE, standardize = TRUE
      )
      message(paste(node, " "), appendLF = FALSE)
      gamma <- as.numeric(coef(result, s = "lambda.min"))
      return(gamma)
    }

    if (ncores > 1) {
      step1_result <- parallel::mclapply(seq_len(p), nodewise_gamma, mc.cores = ncores)
    } else {
      step1_result <- lapply(seq_len(p), nodewise_gamma)
    }

    message("\nFinished stage 1")

    for (node in seq_len(p)) {
      ghat_mx[node, ] <- step1_result[[node]][-1]
      g0[node] <- step1_result[[node]][1] # intercept
    }

    rm(step1_result)
    gc()
  }

  # Initialize covariate array
  # Includes the population matrix, hence +1
  bhat_tens <- array(0, dim = c(p, p, q + 1))

  # Estimated variances
  varhat <- vector(length = p)

  Z <- responses - (g0 + covariates %*% t(ghat_mx))
  intmx <- intxmx(Z, covariates)

  # foldid <- sample(cut(seq_len(n), nfolds, labels = FALSE))
  foldid <- rep(1:5, each = n / nfolds)
  nodewise_beta <- function(node) {
    y <- Z[, node]
    y <- y - mean(y)
    mx <- intmx[, -(seq(0, q) * p + node)]
    mxs <- as.matrix(mx %*% Matrix::Diagonal(x = 1 / sqrt(Matrix::colSums(mx^2))))

    if (is.null(lam_max)) {
      lam1_max <- max(abs(crossprod(mxs, y)))
    }
    lambda1 <- lam1_max * exp(seq(log(1), log(lambda_factor), length = nlambda))

    # There are (q + 1) groups and the size of each group is p-1
    grp_idx <- rep(1:(q + 1), each = p - 1)

    cvm_mx <- matrix(nrow = nlambda, ncol = nasparse)
    betas <- matrix(nrow = nvars, ncol = nasparse)
    varhat <- matrix(nrow = nlambda, ncol = nasparse)
    mse <- numeric(nasparse)
    pf_group <- c(0, rep(1, q))

    for (asid in seq_len(nasparse)) {
      cv_result <- cv.sparsegl(
        mx, y, grp_idx,
        asparse = asparse[asid],
        pf_group = pf_group,
        foldid = foldid,
        lambda = lambda1
      )
      cvm_mx[, asid] <- cv_result$cvm
      lam_min_ind <- which.min(cv_result$cvm)
      fit <- cv_result$sparsegl.fit
      betas[, asid] <- as.numeric(fit$beta[, lam_min_ind])
      mse[asid] <- fit$mse[lam_min_ind]
    }
    cvind <- arrayInd(which.min(cvm_mx), dim(cvm_mx))
    alpha_min_ind <- cvind[2]
    beta_cv <- betas[, alpha_min_ind]

    # Compute nnz based on largest magnitude coefficients
    bcs <- cumsum(sort(abs(beta_cv), decreasing = TRUE))
    nnz <- which(bcs >= 0.999 * sum(abs(beta_cv)))[1]
    if (nnz >= n) {
      sigma2 <- 1
    } else {
      sigma2 <- mse[alpha_min_ind] * n / abs(n - nnz)
    }

    message(paste(node, " "), appendLF = FALSE)
    return(list(
      beta = beta_cv,
      varhat = sigma2,
      cvind = cvind,
      mse = cvm_mx,
      y = y,
      mx = mx
    ))
  }

  message("Stage 2")

  if (ncores > 1) {
    result <- parallel::mclapply(seq_len(p), nodewise_beta, mc.cores = ncores)
  } else {
    result <- lapply(seq_len(p), nodewise_beta)
  }

  message("\nFinished stage 2")

  reg <- list(
    y = matrix(nrow = p, ncol = n),
    mx = array(dim = c(p, n, (p - 1) * (q + 1)))
  )

  cv_mse <- array(dim = c(p, nlambda, nasparse))
  for (node in seq_len(p)) {
    varhat[node] <- result[[node]]$varhat
    bhat_tens[node, -node, ] <- result[[node]]$beta
    cv_mse[node, , ] <- result[[node]]$mse
    reg$y[node, ] <- result[[node]]$y
    reg$mx[node, , ] <- result[[node]]$mx
  }

  bhat_symm <- array(0, dim = c(p, p, q + 1))
  for (h in seq_len(q + 1)) {
    bhat_symm[, , h] <- symmetrize(-diag(1 / varhat) %*% bhat_tens[, , h])
  }

  result <- list(
    gamma = ghat_mx,
    beta = bhat_symm,
    beta_raw = bhat_tens,
    sigma2 = varhat,
    cv_mse = cv_mse,
    reg = reg,
    intmx = intmx
  )
  class(result) <- "gmmreg"
  result
}
