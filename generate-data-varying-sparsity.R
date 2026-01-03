source("generate-data-fns.R")

p <- 25
q <- 10
n <- 300
setting <- "natural"
set.seed(725726)

nrep <- 50
x_opts <- x_options(setting, "mixed")
mg_opts <- gamma_options(p, q, prob = 0.1, frob = 4)
mg <- generate_mg(mg_opts)
s_gamma <- sum(abs(mg) > 0)
ves <- c(0.005, 0.01, 0.02, 0.03, 0.04, 0.05)
generated <- vector(mode = "list", length = nrep * length(ves))
gid <- 1
message("Generating data with varying sparsities...")
for (rep in seq_len(nrep)) {
  for (ve in ves) {
    tb <- generate_tb(beta_options(p, q, ve))
    gx <- generate_x(n, tb, mg, x_opts)
    gx$tb <- tb
    gx$mg <- mg
    generated[[gid]] <- gx
    gid <- gid + 1
    cat(sum(abs(tb) > 0), " ")
  }
}

dir.create("./data", showWarnings = FALSE)
outpath <- file.path("data", sprintf("p%dq%d-n%d-%s-varying-sparsity.rds", p, q, n, setting))
saveRDS(generated, outpath)
message("\nGenerated data saved to ", outpath)

set.seed(13151)
nrep <- 30
p_gamma <- seq(0.05, 0.8, length = 10)
generated <- vector(mode = "list", length = nrep * length(p_gamma))
gid <- 1
tb <- generate_tb(beta_options(p, q, ve = 0.01))

message("Generating data with varying gamma sparsities...")
for (rep in seq_len(nrep)) {
  for (pg in p_gamma) {
    mg <- generate_mg(gamma_options(p, q, prob = pg, frob = 4))
    gx <- generate_x(n, tb, mg, x_opts)
    gx$tb <- tb
    gx$mg <- mg
    generated[[gid]] <- gx
    gid <- gid + 1
    cat(sum(abs(mg) > 0), " ")
  }
}

dir.create("./data", showWarnings = FALSE)
outpath <- file.path("data", sprintf("p%dq%d-n%d-%s-varying-sparsity-gamma.rds", p, q, n, setting))
saveRDS(generated, outpath)
message("\nGenerated data saved to ", outpath)

set.seed(991145)
nrep <- 30
q <- 50
setting <- "original"
x_opts <- x_options(setting, "mixed")
p_gamma <- seq(0.01, 0.8, length = 10)
generated <- vector(mode = "list", length = nrep * length(p_gamma))
gid <- 1
tb <- generate_tb(beta_options(p, q, ve = 0.01))
message("Generating data with varying gamma sparsities original model...")
for (rep in seq_len(nrep)) {
  for (pg in p_gamma) {
    mg <- generate_mg(gamma_options(p, q, prob = pg, frob = 6))
    gx <- generate_x(n, tb, mg, x_opts)
    gx$tb <- tb
    gx$mg <- mg
    generated[[gid]] <- gx
    gid <- gid + 1
    cat(sum(abs(mg) > 0), " ")
  }
}

dir.create("./data", showWarnings = FALSE)
outpath <- file.path("data", sprintf("p%dq%d-n%d-%s-varying-sparsity-gamma.rds", p, q, n, setting))
saveRDS(generated, outpath)
message("\nGenerated data saved to ", outpath)
