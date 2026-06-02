suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(tidyr))
ggthemr::ggthemr("fresh")

dirs <- list.dirs("out", recursive = FALSE)

# ---- Varying delta ----

parse_dir_name <- function(name) {
  m <- regmatches(name, regexec("p(\\d+)q(\\d+)-n(\\d+)-d([0-9.]+)", name))[[1]]
  if (length(m) == 0) {
    return(NULL)
  }
  list(p = as.integer(m[2]), q = as.integer(m[3]), n = as.integer(m[4]), d = as.numeric(m[5]))
}

rows <- list()
for (d in dirs) {
  params <- tryCatch(parse_dir_name(basename(d)), error = function(e) NULL)
  if (is.null(params)) next
  if (params$p != 25 || params$q != 50 || params$n != 200) next
  method_files <- list.files(d, pattern = "^result-.*\\.csv$", full.names = FALSE)
  methods <- sub("^result-(.+)\\.csv$", "\\1", method_files)
  for (method in methods) {
    f <- file.path(d, paste0("result-", method, ".csv"))
    if (!file.exists(f)) next
    df <- read.csv(f)
    df$p <- params$p
    df$q <- params$q
    df$n <- params$n
    df$d <- params$d
    df$method <- method
    rows <- c(rows, list(df))
  }
}

if (length(rows) > 0) {
  method_labels <- c(cspine = "cspine", RegGMM = "RegGMM", mtRegGMM = "mt-RegGMM")
  data <- bind_rows(rows) |>
    mutate(d = factor(d)) |>
    filter(method %in% names(method_labels))

  if (nrow(data) > 0) {
    present_methods <- intersect(names(method_labels), unique(data$method))
    data$method <- factor(data$method, levels = present_methods, labels = method_labels[present_methods])

    p_tpr <- ggplot(data, aes(x = d, y = tpr, fill = method)) +
      geom_boxplot(outlier.size = 0.8, position = position_dodge(0.8)) +
      labs(x = expression(delta), y = "TPR", fill = "Method") +
      scale_fill_brewer(palette = "Set2") +
      theme(legend.position = "bottom", legend.direction = "horizontal")

    p_beta <- ggplot(data, aes(x = d, y = beta_err, fill = method)) +
      geom_boxplot(outlier.size = 0.8, position = position_dodge(0.8)) +
      labs(x = expression(delta), y = expression(beta ~ error), fill = "Method") +
      scale_fill_brewer(palette = "Set2") +
      theme(legend.position = "bottom", legend.direction = "horizontal")

    # ggsave("out/plot_tpr.pdf", p_tpr, width = 7, height = 4)
    ggsave("out/plot_beta_err.pdf", p_beta, width = 7, height = 4)
    message("Saved out/plot_beta_err.pdf")
  }
}

# ---- Varying SNR ----

parse_snr_dir <- function(name) {
  m <- regmatches(name, regexec("p(\\d+)q(\\d+)-n(\\d+)-varying-snr-c(\\d+)p(\\d+)", name))[[1]]
  if (length(m) == 0) {
    return(NULL)
  }
  list(
    p = as.integer(m[2]),
    q = as.integer(m[3]),
    n = as.integer(m[4]),
    c = as.numeric(paste0(m[5], ".", m[6]))
  )
}

rows <- list()
for (d in dirs) {
  params <- tryCatch(parse_snr_dir(basename(d)), error = function(e) NULL)
  if (is.null(params)) next
  method_files <- list.files(d, pattern = "^result-.*\\.csv$", full.names = FALSE)
  methods <- sub("^result-(.+)\\.csv$", "\\1", method_files)
  for (method in methods) {
    f <- file.path(d, paste0("result-", method, ".csv"))
    if (!file.exists(f)) next
    df <- read.csv(f)
    df$p <- params$p
    df$q <- params$q
    df$n <- params$n
    df$c <- params$c
    df$method <- method
    rows <- c(rows, list(df))
  }
}

if (length(rows) > 0) {
  data <- bind_rows(rows) |>
    filter(method %in% c("cspine", "RegGMM"))

  if (nrow(data) > 0) {
    data$method <- as.factor(data$method)

    snrs <- c(0.01, seq(0.025, 0.1, by = 0.025))
    snr_means <- data.frame(c = sort(unique(data$c)), snr_mean = snrs)
    data <- left_join(data, snr_means, by = "c") |>
      mutate(log_snr = log10(snr_mean), snr_mean = factor(snr_mean, levels = snrs))

    p_beta <- ggplot(data, aes(x = snr_mean, y = beta_err, fill = method, group = interaction(c, method))) +
      geom_boxplot() +
      labs(x = "SNR", y = expression(beta ~ error), fill = "Method") +
      scale_fill_brewer(palette = "Set2") +
      theme(legend.position = "bottom", legend.direction = "horizontal")

    ggsave("out/plot_snr_beta_err.pdf", p_beta, width = 7, height = 4)
    message("Saved out/plot_snr_beta_err.pdf")
  }
}
