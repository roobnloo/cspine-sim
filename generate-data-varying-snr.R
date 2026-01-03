source("generate-data-fns.R")

p <- 25
q <- 50
n <- 200
setting <- "original"
set.seed(38951)

nrep <- 10
x_opts <- x_options(setting, "mixed")
frobs <- runif(nrep, 3, 10)
generated <- vector(mode = "list", length = nrep)
gid <- 1
message("Generating data with SNR...")
while (gid <= nrep) {
  mg_opts <- gamma_options(p, q, prob = 0.3, frob = frobs[gid])
  mg <- generate_mg(mg_opts)
  tb <- generate_tb(beta_options(p, q, 0.01))
  gx <- tryCatch(
    {
      generate_x(n, tb, mg, x_opts)
    },
    error = function(e) {
      # message("Error in generating data: ", e)
      return(NULL)
    }
  )
  if (is.null(gx)) {
    next
  }
  gx$tb <- tb
  gx$mg <- mg
  generated[[gid]] <- gx
  gid <- gid + 1
  cat(mean(gx$snr), " ")
}

dir.create("./data", showWarnings = FALSE)
outpath <- file.path("data", sprintf("p%dq%d-n%d-%s-varying-snr.rds", p, q, n, setting))
saveRDS(generated, outpath)
message("\nGenerated data saved to ", outpath)
