source("generate-data-fns.R")

p <- 25
n <- 200
setting <- "original"
set.seed(3242113)

nrep <- 1
x_opts <- x_options(setting, "mixed")
q_seq <- seq(50, 100, by = 5)
generated <- vector(mode = "list", length = nrep * length(q_seq))
gid <- 1
message("Generating data with varying sparsities...")
for (rep in seq_len(nrep)) {
  for (q in q_seq) {
    mg_opts <- gamma_options(p, q, prob = 0.3, frob = 6)
    mg <- generate_mg(mg_opts)
    tb <- generate_tb(beta_options(p, q))
    gx <- generate_x(n, tb, mg, x_opts)
    gx$tb <- tb
    gx$mg <- mg
    generated[[gid]] <- gx
    gid <- gid + 1
    cat(q, " ")
  }
}

dir.create("./data", showWarnings = FALSE)
outpath <- file.path("data", sprintf("p%dq%d-n%d-%s-varying-q.rds", p, q, n, setting))
saveRDS(generated, outpath)
message("\nGenerated data saved to ", outpath)
