suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggthemr))

c_results <- readRDS("out/p25q10-n300-natural-varying-sparsity-result-cspine.rds")
g_results <- readRDS("out/p25q10-n300-natural-varying-sparsity-result-RegGMM.rds")

ggthemr::ggthemr("fresh")
tibble(
  s_beta = c_results[, "s_beta"], err =
    c_results[, "beta_err"] + c_results[, "gamma_err"]
) |>
  ggplot(mapping = aes(x = s_beta, y = err)) +
  geom_point(alpha = 0.7, shape = 4, color = swatch()[3]) +
  # stat_function(fun = \(x) val(theta, x), linetype = "dashed") +
  ylab(expression(beta[~ ~err] + gamma[~err])) +
  xlab(expression(s[B]))
ggsave("s_beta_err.png")

c_results_g <- readRDS("out/p25q10-n300-natural-varying-sparsity-result-gamma-cspine.rds")
g_results_g <- readRDS("out/p25q10-n300-natural-varying-sparsity-result-gamma-RegGMM.rds")
tibble(
  s_gamma = c_results_g[, "s_gamma"],
  err = c_results_g[, "beta_err"] + c_results_g[, "gamma_err"]
) |>
  ggplot(mapping = aes(x = s_gamma, y = err)) +
  geom_point(alpha = 0.7, shape = 4, color = swatch()[3]) +
  # stat_function(fun = \(x) val(theta, x), linetype = "dashed") +
  ylab(expression(beta[~ ~err] + gamma[~err])) +
  xlab(expression(s[Gamma]))
ggsave("s_gamma_err.png")


# val <- function(theta, x) {
#   m1 <- theta[1]
#   m2 <- theta[2]
#   b1 <- theta[3]
#   b2 <- theta[4]
#   sqrt(x * m1 + b1) * m2 + b2
# }

# fn <- function(theta, x, y) {
#   sum((y - (val(theta, x)))^2)
# }

# inds <- which(df$s_beta < 80)
# opt <- optim(rep(1, 4), fn, NULL, df$s_beta[inds], df$err[inds])
# (theta <- opt$par)
# axis <- seq(50, 200)

c_results_g <- readRDS("out/p25q50-n300-original-varying-sparsity-result-gamma-cspine.rds")
g_results_g <- readRDS("out/p25q50-n300-original-varying-sparsity-result-gamma-RegGMM.rds")
ggthemr::ggthemr("fresh")
tibble(
  q = c(c_results_g[, "s_gamma"], g_results_g[, "s_gamma"]),
  err = c(c_results_g[, "beta_err"], g_results_g[, "beta_err"]),
  tpr = c(c_results_g[, "tpr"], g_results_g[, "tpr"]),
  method = rep(c("cspine", "RegGMM"), each = nrow(c_results_g))
) |>
  ggplot(mapping = aes(x = q, y = err, color = method, shape = method)) +
  geom_jitter(alpha = 0.8, width = 20) +
  geom_smooth(se = FALSE) +
  scale_shape_manual(values = c(4, 21)) +
  scale_color_manual(values = swatch()[3:2]) +
  ylab(expression(beta[~ ~err])) +
  xlab(expression(s[Gamma])) +
  theme(plot.margin = margin(l = 10, r = 5))
# ggsave("s_gamma_err.png")
