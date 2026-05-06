library(dplyr)
library(stringr)

out_dir <- "out"
target_d <- 0.00

parse_dir_name <- function(name) {
  m <- regmatches(name, regexec("p(\\d+)q(\\d+)-n(\\d+)-d([0-9.]+)", name))[[1]]
  list(p = as.integer(m[2]), q = as.integer(m[3]), n = as.integer(m[4]), d = as.numeric(m[5]))
}

dirs <- list.dirs(out_dir, recursive = FALSE)
dirs <- dirs[grepl("^p", basename(dirs))]
rows_mean <- list()
rows_sd <- list()

for (d in dirs) {
  params <- tryCatch(parse_dir_name(basename(d)), error = function(e) NULL)
  if (is.null(params) || is.na(params$d)) next
  if (abs(params$d - target_d) > 1e-9) next

  method_files <- list.files(d, pattern = "^result-.*\\.csv$", full.names = FALSE)
  methods <- sub("^result-(.+)\\.csv$", "\\1", method_files)

  for (method in methods) {
    f <- file.path(d, paste0("result-", method, ".csv"))
    df <- read.csv(f)
    means <- as.list(colMeans(df, na.rm = TRUE))
    sds <- as.list(apply(df, 2, sd, na.rm = TRUE))
    base <- c(params, method = method)
    rows_mean <- c(rows_mean, list(as.data.frame(c(base, means))))
    rows_sd <- c(rows_sd, list(as.data.frame(c(base, sds))))
  }
}

mean_df <- bind_rows(rows_mean) |>
  # filter(n == 400) |>
  select(p, q, n, d, method, everything()) |>
  arrange(p, q, n, d, method)

sd_df <- bind_rows(rows_sd) |>
  # filter(n == 400) |>
  select(p, q, n, d, method, everything()) |>
  arrange(p, q, n, d, method)

fmt_val <- function(x) {
  if (is.na(x)) "--" else sprintf("$%1.3f$", x)
}
fmt_sd <- function(x) {
  if (is.na(x)) "" else sprintf("$(%1.3f)$", x)
}

# bold_dir: "G" = bold max, "L" = bold min, "-" = no bold
# embolden_methods: subset of method names to compete for bold (others are never bolded)
embolden <- function(method_lines, metric_vals, bold_dirs,
                     embolden_methods = names(method_lines)) {
  stopifnot(length(method_lines) == length(metric_vals))
  stopifnot(length(bold_dirs) == length(metric_vals[[1]]))

  compete_ids <- which(names(method_lines) %in% embolden_methods)

  for (j in seq_along(bold_dirs)) {
    bd <- bold_dirs[j]
    if (bd == "-") next
    col_vals <- sapply(metric_vals[compete_ids], function(v) v[j])
    if (any(is.na(col_vals))) next # skip bolding if any competing method is NA

    if (bd == "G") {
      local_ids <- which(col_vals == max(col_vals))
    } else {
      local_ids <- which(col_vals == min(col_vals))
    }
    best_ids <- compete_ids[local_ids]
    formatted <- sprintf("%1.3f", col_vals[local_ids])
    for (k in seq_along(best_ids)) {
      method_lines[best_ids[k]] <- sub(
        paste0("$", formatted[k], "$"),
        sprintf("$\\mathbf{%s}$", formatted[k]),
        method_lines[best_ids[k]],
        fixed = TRUE
      )
    }
  }
  method_lines
}

make_table_d1 <- function() {
  compared_methods <- c("cspine", "RegGMM", "mtRegGMM", "RegGMM-oracle")
  metrics <- c("tpr", "fpr", "beta_err", "omega_err")
  bold_dirs <- c("G", "L", "L", "L")

  settings <- mean_df |>
    select(p, q, n) |>
    distinct() |>
    arrange(n, p, q)

  latex_code <- "\\begin{table}[ht]\n"
  latex_code <- paste0(latex_code, "\\centering\n")
  latex_code <- paste0(latex_code, "\\small\n")
  latex_code <- paste0(latex_code, "\\begin{tabular}{r|r|r|cccc}\n")
  latex_code <- paste0(latex_code, "  \\hline\n")
  latex_code <- paste0(
    latex_code,
    r"($n$ & $(p, q)$ & Method & $\text{TPR}$ & $\text{FPR}$ & $\vec\beta_\text{err}$ & $\mat\Omega_\text{err}$ \\)", "\n"
  )
  latex_code <- paste0(latex_code, "  \\hline\n")

  for (s in seq_len(nrow(settings))) {
    sp <- settings$p[s]
    sq <- settings$q[s]
    sn <- settings$n[s]
    next_n <- if (s < nrow(settings)) settings$n[s + 1] else NA

    setting_methods <- mean_df |>
      filter(p == sp, q == sq, n == sn) |>
      pull(method)
    active_methods <- intersect(compared_methods, setting_methods)
    num_active <- length(active_methods)

    method_lines <- character(num_active)
    names(method_lines) <- active_methods
    metric_vals <- vector("list", num_active)
    names(metric_vals) <- active_methods

    for (mi in seq_along(active_methods)) {
      m <- active_methods[mi]
      mr <- mean_df |> filter(p == sp, q == sq, n == sn, method == m)
      sr <- sd_df |> filter(p == sp, q == sq, n == sn, method == m)

      prefix <- if (mi == 1) {
        n_str <- if (is.na(next_n) || next_n != sn) as.character(sn) else ""
        paste0(n_str, " & \\multirow{", num_active, "}{*}{$(", sp, ", ", sq, ")$}")
      } else {
        " & "
      }

      mvals <- setNames(sapply(metrics, function(col) mr[[col]]), metrics)
      metric_vals[[mi]] <- mvals

      cells <- paste(
        mapply(function(col) {
          paste(fmt_val(mr[[col]]), fmt_sd(sr[[col]]))
        }, metrics),
        collapse = " & "
      )

      method_lines[mi] <- paste0(
        prefix, " & \\texttt{", m, "} & ", cells, " \\\\\n"
      )
    }

    # method_lines <- embolden(method_lines, metric_vals, bold_dirs,
    #   embolden_methods = setdiff(active_methods, "RegGMM-oracle")
    # )
    latex_code <- paste0(latex_code, paste0(method_lines, collapse = ""))
    hline_str <- if (!is.na(next_n) && next_n == sn) "  \\cline{2-7}\n" else "  \\hline\n"
    latex_code <- paste0(latex_code, hline_str)
  }

  latex_code <- paste0(latex_code, "\\end{tabular}\n")
  latex_code <- paste0(
    latex_code,
    r"(\caption{Mean and standard deviation of performance metrics over simulated data sets, $\delta = 1$.})", "\n"
  )
  latex_code <- paste0(latex_code, "\\label{tbl:sim-d1}\n")
  latex_code <- paste0(latex_code, "\\end{table}\n")

  cat(latex_code)
}

make_table_d1()
