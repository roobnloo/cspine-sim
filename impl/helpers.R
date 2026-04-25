# [X, u1 X, ..., uq X]
intxmx <- function(x_mat, u_mat) {
  i_u <- cbind(1, u_mat)
  result <- lapply(seq_len(ncol(i_u)), \(j) {
    x_mat * i_u[, j]
  })
  Reduce(cbind, result)
}

symmetrize <- function(mx, rule = "and") {
  if (rule == "and") {
    result <- mx * (abs(mx) < t(abs(mx))) + t(mx) * (t(abs(mx)) < abs(mx))
  } else {
    result <- mx * (abs(mx) >= t(abs(mx))) + t(mx) * (t(abs(mx)) >= abs(mx))
  }
  result
}
