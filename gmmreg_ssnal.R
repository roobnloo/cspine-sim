library(sglssnal)

# [X, u1 X, ..., uq X]
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


gmmreg_ssnal <- function(
    responses, covariates, asparse = seq(0.1, 1, by = 0.1),
    nlambda = 100, lam_max = NULL, lambda_factor = 0.01,
    nfolds = 5, verbose = FALSE, ncores = 1, skip_stage1 = FALSE) {
  stopifnot(
    is.matrix(responses), is.matrix(covariates),
    nrow(responses) == nrow(covariates),
    all(asparse > 0), all(asparse <= 1)
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
      result <- glmnet::cv.glmnet(covariates, responses[, node], alpha = 1)
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
  sigma2 <- vector(length = p)

  Z <- responses - (g0 + covariates %*% t(ghat_mx))
  zsd <- sqrt(colSums(Z^2) / n)
  zsd[abs(zsd) < 1e-9] <- 1
  umean <- colMeans(covariates)
  u <- sweep(covariates, 2, umean, "-")
  usd <- sqrt(colSums(u^2) / n)
  usd[abs(usd) < 1e-9] <- 1
  u <- sweep(u, 2, usd, "/")

  # foldid <- rep(1:5, each = n / nfolds)
  nodewise_beta <- function(node) {
    y <- Z[, node]
    y <- y - mean(y)
    z_scale <- sweep(Z, 2, zsd, "/")[, -node]
    mx <- intxmx(z_scale, u)

    if (is.null(lam_max)) {
      mina <- min(asparse)
      if (mina == 0) {
        mina <- 1
      }
      lam1_max <- max(abs(crossprod(mx, y))) / mina
    }
    lambda1 <- lam1_max * exp(seq(log(1), log(lambda_factor), length = nlambda))

    # There are (q + 1) groups and the size of each group is p-1
    grp_vec <- seq(1, nvars)
    start_ids <- 1 + (p - 1) * seq(0, q)
    end_ids <- start_ids + p - 2
    grp_idx <- rbind(start_ids, end_ids)
    pfgroup <- c(0, rep(1, q))

    foldid <- cut(sample(seq_len(n)), nfolds, labels = FALSE)
    cv_result <- sglssnal::cv.sglssnal(
      mx, y, grp_vec, grp_idx, asparse, lambda1,
      pfgroup = pfgroup, foldid = foldid, quietall = TRUE
    )

    nnz <- cv_result$info$nnz
    sigma2 <- 1
    if (n > nnz) {
      sigma2 <- cv_result$info$mse * n / abs(n - nnz)
    }

    x_unstd <- cv_result$x
    for (h in seq_len(q)) {
      x_unstd[1:(p - 1)] <- x_unstd[1:(p - 1)] -
        x_unstd[(p - 1) * h + 1:(p - 1)] * umean[h] / usd[h]
    }
    x_unstd[1:(p - 1)] <- x_unstd[1:(p - 1)] / zsd[-node]

    for (h in seq_len(q)) {
      x_unstd[(p - 1) * h + 1:(p - 1)] <- x_unstd[(p - 1) * h + 1:(p - 1)] / (usd[h] * zsd[-node])
    }
    x_unstd[abs(x_unstd) < 1e-9] <- 0

    message(paste(node, " "), appendLF = FALSE)
    return(list(
      beta = x_unstd,
      sigma2 = sigma2,
      cvm = cv_result$cv_info$cvm,
      cv_idx = cv_result$cv_info$cv_idx,
      mse = cv_result$info$mse,
      y = y,
      mx = mx,
      lambda = lambda1
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
  cv_idx <- matrix(nrow = p, ncol = 2)
  lambdas <- matrix(nrow = p, ncol = nlambda)
  for (node in seq_len(p)) {
    sigma2[node] <- result[[node]]$sigma2
    bhat_tens[node, -node, ] <- result[[node]]$beta
    cv_mse[node, , ] <- result[[node]]$cvm
    cv_idx[node, ] <- result[[node]]$cv_idx
    reg$y[node, ] <- result[[node]]$y
    reg$mx[node, , ] <- result[[node]]$mx
    lambdas[node, ] <- result[[node]]$lambda
  }

  bhat_symm <- array(0, dim = c(p, p, q + 1))
  for (h in seq_len(q + 1)) {
    bhat_symm[, , h] <- symmetrize(-diag(1 / sigma2) %*% bhat_tens[, , h])
  }

  result <- list(
    gamma = ghat_mx,
    beta = bhat_symm,
    beta_raw = bhat_tens,
    sigma2 = sigma2,
    cv_mse = cv_mse,
    cv_idx = cv_idx,
    reg = reg,
    lambdas = lambdas
  )
  class(result) <- "gmmreg"
  result
}
