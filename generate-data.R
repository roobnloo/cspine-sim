source("generate-data-fns.R")

args <- commandArgs(trailingOnly = TRUE)
p <- as.integer(args[1])
q <- as.integer(args[2])
n <- as.integer(args[3])
setting <- args[4]
stopifnot(setting %in% c("natural", "original"))
seed <- as.numeric(args[5])

message("Generating with ", sprintf("(%d, %d, %d, %s)...", p, q, n, setting))
# p <- 25; q<-100; n<-200; setting <- "original"

n_rep <- 100
x_opts <- x_options(setting, "mixed")
tb_opts <- beta_options(p, q, ve = 0.01)

# Scale Frobenius norm of mg to control SNR.
mg_scale <- 6
if (setting == "natural") {
  mg_scale <- 4
}
mg_opts <- gamma_options(p, q, prob = 0.3, frob = mg_scale)

set.seed(seed)
tb <- generate_tb(tb_opts)
mg <- generate_mg(mg_opts)

generated <- vector(mode = "list", length = n_rep + 1)
for (i in seq_len(n_rep)) {
  generated[[i]] <- generate_x(n, tb, mg, x_opts)
  cat(i, " ")
}
generated[[n_rep + 1]] <- list(tb = tb, mg = mg)

dir.create("./data", showWarnings = FALSE)
outpath <- file.path("data", sprintf("p%dq%d-n%d-%s.rds", p, q, n, setting))
saveRDS(generated, outpath)
message("\nGenerated data saved to ", outpath)
