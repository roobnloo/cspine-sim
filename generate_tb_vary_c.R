library(MASS)
library(igraph)

p <- 25
q <- 10
c_scale <- seq(0.5, 1.3, by = 0.2)
nc <- length(c_scale)

l <- 0.35
u <- 0.5

generate_tB <- function(p, q, c_val) {
  tB <- array(0, dim = c(p, p, q + 1))

  g1 <- sample_pa(p, power = 2.5, directed = FALSE)
  A <- as_adjacency_matrix(g1, sparse = FALSE)
  rind <- sample(1:p, p, replace = FALSE)
  A <- A[rind, rind]
  tb <- matrix(0, p, p)
  tb[lower.tri(A) & A > 0] <- sample(
    c(runif(sum(A), -u, -l), runif(sum(A), l, u)),
    sum(A) / 2,
    replace = FALSE
  )
  tB[, , 1] <- tb + t(tb)

  for (j in 2:6) {
    g2 <- sample_gnp(p, 0.01, directed = FALSE, loops = FALSE)
    A <- as_adjacency_matrix(g2, sparse = FALSE)
    rind <- sample(1:p, p, replace = FALSE)
    A <- A[rind, rind]
    tb <- matrix(0, p, p)
    tb[lower.tri(A) & A > 0] <- sample(
      c(runif(sum(A), -u, -l), runif(sum(A), l, u)),
      sum(A) / 2,
      replace = FALSE
    )
    tB[, , j] <- tb + t(tb)
  }

  tB_temp <- array(0, dim = c(p, p, q + 1))
  for (j in 1:p) {
    s <- sum(abs(tB[, j, ]))
    if (s > 0) tB_temp[, j, ] <- tB[, j, ] / s / 1.5
  }
  for (j in 1:(q + 1)) {
    tB[, , j] <- (tB_temp[, , j] + t(tB_temp[, , j])) / 2
  }

  tB * c_val
}

set.seed(42)

tB_all <- array(0, dim = c(p, p, q + 1, nc))
for (i in seq_along(c_scale)) {
  tB_all[, , , i] <- generate_tB(p, q, c_scale[i])
}

mg <- matrix(0, p, q)
mg[sample(1:(p * q), p * q * 0.1)] <- 0.15

saveRDS(
  list(tB_all = tB_all, mg = mg, c_scale = c_scale),
  file.path("data", sprintf("coef_p%dq%d_vary_c.rds", p, q))
)
