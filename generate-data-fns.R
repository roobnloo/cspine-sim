suppressPackageStartupMessages(library(MASS))
suppressPackageStartupMessages(library(igraph))

x_options <- function(
    model = c("original", "natural"),
    cov_type = c("binary", "mixed"),
    sigma2 = 1) {
  list(
    model = match.arg(model), cov_type = match.arg(cov_type), sigma2 = sigma2
  )
}

generate_x <- function(n, tb, mg, opts = x_options()) {
  p <- nrow(mg)
  q <- ncol(mg)
  U <- matrix(sample(c(0, 1), n * q, replace = TRUE), n, q)
  if (opts$cov_type == "mixed") {
    U2 <- apply(matrix(runif(n * q), n, q), 2, scale)
    cind <- sample(1:q, q / 2, replace = FALSE)
    U[, cind] <- U2[, cind] # select continuous covariates
  }
  iU <- cbind(rep(1, n), U)
  X <- matrix(0, n, p)
  mumx <- matrix(0, n, p) # mean vector for each observation
  snr <- numeric(n)

  for (i in 1:n) {
    omega <- apply(tb, c(1, 2), \(b) b %*% iU[i, ])
    diag(omega) <- opts$sigma2
    sigma <- solve(omega)
    mu <- mg %*% U[i, ]
    if (opts$model == "natural") {
      mu <- sigma %*% mu
      snr[i] <- sum(diag(crossprod(sigma %*% mg))) / (sum(diag(sigma)))
    } else {
      snr[i] <- sum(diag(crossprod(mg))) / (sum(diag(sigma)))
    }

    mumx[i, ] <- mu
    X[i, ] <- mvrnorm(1, mu, sigma)
  }

  return(list(X = X, U = U, mumx = mumx, snr = snr))
}

gamma_options <- function(p, q, prob = 0.2, frob = NULL) {
  list(p = p, q = q, prob = prob, frob = frob)
}

generate_mg <- function(opts = gamma_options()) {
  p <- opts$p
  q <- opts$q
  mg_bern <- matrix(rbinom(p * q, 1, opts$prob), p, q)
  mg <- rnorm(p * q) * mg_bern
  mg <- mg * mg_bern
  if (!is.null(opts$frob)) {
    mg <- mg / sqrt(sum(mg^2))
    mg <- mg * opts$frob
  }
  return(mg)
}

beta_options <- function(
    p = 25, q = 50, ve = 0.01, nz_cov = 5, lim = c(0.35, 0.5), pwr0 = 1, sigma2 = 1) {
  list(
    p = p, q = q, ve = ve, nz_cov = nz_cov, lim = lim, pwr0 = pwr0, sigma2 = sigma2
  )
}

generate_tb <- function(opts = beta_options()) {
  p <- opts$p
  q <- opts$q

  tB <- array(0, dim = c(p, p, q + 1)) #+1 for intercept
  l <- opts$lim[1]
  u <- opts$lim[2]

  # degs<-c(2,2,1,1,1,1,rep(0,p-6))
  # g <- sample_degseq(degs, method="simple.no.multiple") #generate with given degrees

  g1 <- sample_pa(p, power = opts$pwr0, directed = FALSE) # scale-free network
  A <- as_adjacency_matrix(g1, sparse = FALSE)
  rind <- sample(1:p, p, replace = FALSE)
  A <- A[rind, rind]
  tb <- matrix(0, p, p)
  tb[lower.tri(A) & A > 0] <- sample(
    c(
      runif(sum(A), -u, -l),
      runif(sum(A), l, u)
    ),
    sum(A) / 2,
    replace = FALSE
  )
  tB[, , 1] <- tb + t(tb)

  for (j in seq(2, opts$nz_cov + 1)) {
    all_zero <- TRUE
    while (all_zero) {
      g2 <- sample_gnp(p, opts$ve, directed = FALSE, loops = FALSE) # random network
      A <- as_adjacency_matrix(g2, sparse = FALSE)
      rind <- sample(1:p, p, replace = FALSE)
      A <- A[rind, rind]
      tb <- matrix(0, p, p)
      tb[lower.tri(A) & A > 0] <- sample(
        c(
          runif(sum(A), -u, -l),
          runif(sum(A), l, u)
        ),
        sum(A) / 2,
        replace = FALSE
      )
      tB[, , j] <- tb + t(tb)
      if (any(tB[, , j] != 0)) {
        all_zero <- FALSE
      }
    }
  }

  tB_temp <- array(0, dim = c(p, p, q + 1))
  for (j in 1:p) {
    # ensure diagonal dominance
    tB_temp[j, , ] <- tB[j, , ] / (sum(abs(tB[j, , ])) * 1.5)
  }

  for (j in 1:(q + 1)) {
    tB[, , j] <- (tB_temp[, , j] + t(tB_temp[, , j])) / 2
  }

  return(tB)
}
