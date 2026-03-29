library(tidyr)
library(dplyr)
library(stringr)

prefix <- "set0"
out_dir <- "out"

# methods <- factor(c("cspine", "RegGMM", "ANTAC", "glasso", "MB"),
#                   levels = c("cspine", "RegGMM", "ANTAC", "glasso", "MB"))
methods <- factor(c("cspine", "RegGMM"), levels = c("cspine", "RegGMM"))

num_pq <- 2 # (25, 50), (50, 50), (25, 100)
num_n <- 2 # 200, 400

meandf <- expand.grid(
  method = methods,
  p = c(25), q = c(50, 100),
  n = c(200, 400), model = c("natural", "original"),
  betaTPR = 0, betaFPR = 0, popTPR = 0, popFPR = 0,
  covTPR = 0, covFPR = 0, beta_err = 0, omega_err = 0,
  mean_err = 0, omega_tpr = 0, omega_fpr = 0
)

sddf <- expand.grid(
  method = methods,
  p = c(25), q = c(50, 100),
  n = c(200, 400), model = c("natural", "original"),
  betaTPR = 0, betaFPR = 0, popTPR = 0, popFPR = 0,
  covTPR = 0, covFPR = 0, beta_err = 0, omega_err = 0,
  mean_err = 0, omega_tpr = 0, omega_fpr = 0
)

metrics <- c(
  "betaTPR", "betaFPR", "popTPR", "popFPR", "covTPR", "covFPR",
  "beta_err", "omega_err", "mean_err", "omega_tpr", "omega_fpr"
)

for (rowid in seq_len(nrow(meandf))) {
  row <- meandf[rowid, ]
  path <- ""
  switch(as.character(row$method),
    cspine = {
      path <- sprintf(
        "%s/p%dq%d-n%d-%s-result-%s.rds", out_dir,
        row[[2]], row[[3]], row[[4]], row[[5]], row$method
      )
    },
    RegGMM = {
      path <- sprintf(
        "%s/p%dq%d-n%d-%s-result-%s.rds", out_dir,
        row[[2]], row[[3]], row[[4]], row[[5]], row$method
      )
    },
    glasso = {
      path <- sprintf(
        "%s/p%dq%d/glasso_n%dreparam%s-result.rds", prefix,
        row[[2]], row[[3]], row[[4]], row[[5]]
      )
    },
    MB = {
      path <- sprintf(
        "%s/p%dq%d/mb_n%dreparam%s-result.rds", prefix,
        row[[2]], row[[3]], row[[4]], row[[5]]
      )
    },
    ANTAC = {
      path <- sprintf(
        "%s/p%dq%d/antac_n%dreparam%s-result.rds", prefix,
        row[[2]], row[[3]], row[[4]], row[[5]]
      )
    },
    {
      stop("Invalid method!")
    }
  )

  result <- readRDS(path)

  if (row$method %in% c("cspine", "RegGMM")) {
    colnames(result)[1:6] <- metrics[1:6]
  } else {
  }

  res_mean <- apply(result, 2, mean, na.rm = TRUE) |> round(3)
  res_sd <- apply(result, 2, sd, na.rm = TRUE) |> round(3)
  meandf[rowid, ][metrics] <- res_mean[metrics]
  sddf[rowid, ][metrics] <- res_sd[metrics]
}

embolden <- function(method_line, bold_pattern) {
  # Pattern is a string of L (lesser) and G (greater), indicating which value should be bold.
  # E.g. "GLLL"
  # line1 <-  "natural & 200 & (25, 50) & \\texttt{cspine} & $0.947$ ($0.037$) & $0.002$ ($0.000$) & $4.249$ ($0.357$) "
  # line2 <- "natural & 200 & (25, 50) & \\texttt{RegGMM} & $0.822$ ($0.063$) & $0.003$ ($0.001$) & $4.767$ ($0.422$) "
  # method_line <- c(line1, line2)

  bold_pattern <- strsplit(bold_pattern, "")[[1]]
  # pattern <- paste(rep(r"(\$([0-9]+\.[0-9]+)\$)", length(bold_pattern)), collapse = " .* & ")
  select_means <- seq_along(bold_pattern) * 2 - 1

  vals_str <-
    lapply(method_line, \(line) str_match_all(line, r"(\$([0-9]+\.[0-9]+)\$)")[[1]][, -1])
  vals_str <- Reduce(rbind, vals_str)
  vals <- as.numeric(vals_str)
  dim(vals) <- dim(vals_str)

  if (ncol(vals) != length(bold_pattern)) {
    stop(sprintf(
      "Length of bold pattern (%d) does not match the number of values in the line (%d)",
      length(bold_pattern), ncol(vals)
    ))
  }

  for (i in seq_along(bold_pattern)) {
    bp <- bold_pattern[i]
    if (bp == "G") {
      maxids <- which(vals[, i] == max(vals[, i]))
      for (maxid in maxids) {
        method_line[maxid] <-
          sub(paste0("$", vals_str[maxid, i], "$"),
            sprintf("$\\mathbf{%s}$", vals_str[maxid, i]),
            method_line[maxid],
            fixed = TRUE
          )
      }
    } else if (bp == "L") {
      minids <- which(vals[, i] == min(vals[, i]))
      for (minid in minids) {
        method_line[minid] <-
          sub(vals_str[minid, i],
            sprintf("\\mathbf{%s}", vals_str[minid, i]),
            method_line[minid],
            fixed = TRUE
          )
      }
    } else {
      stop("Invalid bold pattern")
    }
  }

  return(method_line)
}

supp_table <- function(model = c("natural", "original")) {
  model <- match.arg(model)
  # compared_methods <- c("cspine", "RegGMM", "ANTAC", "glasso")
  compared_methods <- c("cspine", "RegGMM")
  num_methods <- length(compared_methods)

  # Begin LaTeX table
  latex_code <- "\\begin{table}[ht]\n"
  latex_code <- paste0(latex_code, "\\centering\n")
  latex_code <- paste0(latex_code, "\\small\n")
  latex_code <- paste0(latex_code, "\\begin{tabular}{r|r|r|rrrr}\n")
  latex_code <- paste0(latex_code, "  \\hline\n")
  latex_code <- paste0(latex_code, r"($n$ & $(p, q)$ & Method & $\mat\Omega_\text{TPR}$ & $\mat\Omega_\text{FPR}$ & $\mat\Omega_\text{err}$ & $\vec\mu_\text{err}$ \\)", "\n")
  latex_code <- paste0(latex_code, "  \\hline\n")

  # Loop through the dataframe to fill in values
  line_counter <- 1
  for (i in seq_len(nrow(meandf))) {
    line_code <- ""
    if (model != meandf$model[i]) {
      next
    }
    # if (line_counter %% (num_methods * num_pq * num_n) == 1) {
    #   if (meandf$reparam[i]) {
    #     line_code <- "natural"
    #   } else {
    #     line_code <- "original"
    #   }
    # }
    # line_code <- paste0(line_code, " & ")


    if (line_counter %% (num_methods * num_pq) == 1) {
      line_code <- paste0(line_code, meandf$n[i])
    }
    line_code <- paste0(line_code, " & ")

    if (line_counter %% num_methods == 1) {
      line_code <- paste0(line_code, "(", meandf$p[i], ", ", meandf$q[i], ")")
    }
    line_code <- paste0(line_code, " & ", sprintf("\\texttt{%s}", meandf$method[i]), " & ")

    if (!(meandf$method[i] %in% compared_methods)) {
      next
    }
    if (meandf$method[i] == "cspine") {
      method_line <- vector(length = length(compared_methods), mode = "character")
    }
    names(method_line) <- compared_methods
    latex_line <- paste0(line_code, sprintf("$%1.2f$", meandf$omega_tpr[i]), " ", sprintf("$(%1.2f)$", sddf$omega_tpr[i]), " & ")
    latex_line <- paste0(latex_line, sprintf("$%1.2f$", meandf$omega_fpr[i]), " ", sprintf("$(%1.2f)$", sddf$omega_fpr[i]), " & ")
    latex_line <- paste0(latex_line, sprintf("$%1.2f$", meandf$omega_err[i]), " ", sprintf("$(%1.2f)$", sddf$omega_err[i]), " & ")
    latex_line <- paste0(latex_line, sprintf("$%1.2f$", meandf$mean_err[i]), " ", sprintf("$(%1.2f)$", sddf$mean_err[i]), " \\\\\n")

    method_line[meandf$method[i]] <- latex_line

    if (meandf$method[i] == compared_methods[length(compared_methods)]) {
      line_res <- embolden(method_line, "GLLL")
      latex_code <- paste0(latex_code, paste0(line_res, collapse = ""))
    }
    line_counter <- line_counter + 1
  }

  # End LaTeX table
  latex_code <- paste0(latex_code, "  \\hline\n")
  latex_code <- paste0(latex_code, "\\end{tabular}\n")
  latex_code <- paste0(latex_code, "\\caption{my caption}\n")
  latex_code <- paste0(latex_code, "\\end{table}\n")

  # Print LaTeX code
  cat(latex_code)
}


main_table <- function(model = "natural") {
  compared_methods <- c("cspine", "RegGMM")
  num_methods <- length(compared_methods)

  # Begin LaTeX table
  latex_code <- "\\begin{table}[ht]\n"
  latex_code <- paste0(latex_code, "\\centering\n")
  latex_code <- paste0(latex_code, "\\small\n")
  latex_code <- paste0(latex_code, "\\begin{tabular}{r|r|r|rrrrrr}\n")
  latex_code <- paste0(latex_code, "  \\hline\n")
  latex_code <- paste0(latex_code, r"($n$ & $(p, q)$ & Method & $\text{TPR}$ & $\text{TPR}_\text{pop}$ & $\text{FPR}_\text{pop}$ & $\text{TPR}_\text{cov}$ & $\vec\beta_\text{err}$ & $\mat\Omega_\text{err}$\\)", "\n")
  latex_code <- paste0(latex_code, "  \\hline\n")

  # Loop through the dataframe to fill in values
  line_counter <- 1
  for (i in seq_len(nrow(meandf))) {
    line_code <- ""
    if (model != meandf$model[i]) {
      next
    }

    if (line_counter %% (num_methods * num_pq) == 1) {
      line_code <- paste0(line_code, meandf$n[i])
    }
    line_code <- paste0(line_code, " & ")

    if (line_counter %% num_methods == 1) {
      line_code <- paste0(line_code, "(", meandf$p[i], ", ", meandf$q[i], ")")
    }
    line_code <- paste0(line_code, " & ", sprintf("\\texttt{%s}", meandf$method[i]), " & ")

    if (!(meandf$method[i] %in% compared_methods)) {
      next
    }
    if (meandf$method[i] == "cspine") {
      method_line <- vector(length = length(compared_methods), mode = "character")
    }
    names(method_line) <- compared_methods
    latex_line <- paste0(line_code, sprintf("$%1.2f$", meandf$betaTPR[i]), " ", sprintf("$(%1.2f)$", sddf$betaTPR[i]), " & ")
    latex_line <- paste0(latex_line, sprintf("$%1.2f$", meandf$popTPR[i]), " ", sprintf("$(%1.2f)$", sddf$popTPR[i]), " &")
    latex_line <- paste0(latex_line, sprintf("$%1.2f$", meandf$popFPR[i]), " ", sprintf("$(%1.2f)$", sddf$popFPR[i]), " &")
    latex_line <- paste0(latex_line, sprintf("$%1.2f$", meandf$covTPR[i]), " ", sprintf("$(%1.2f)$", sddf$covTPR[i]), " &")
    latex_line <- paste0(latex_line, sprintf("$%1.2f$", meandf$beta_err[i]), " ", sprintf("$(%1.2f)$", sddf$beta_err[i]), " &")
    latex_line <- paste0(latex_line, sprintf("$%1.2f$", meandf$omega_err[i]), " ", sprintf("$(%1.2f)$", sddf$omega_err[i]), " \\\\\n")

    method_line[meandf$method[i]] <- latex_line

    if (meandf$method[i] == compared_methods[length(compared_methods)]) {
      line_res <- embolden(method_line, "GGLGLL")
      latex_code <- paste0(latex_code, paste0(line_res, collapse = ""))
    }
    line_counter <- line_counter + 1
  }

  # End LaTeX table
  latex_code <- paste0(latex_code, "  \\hline\n")
  latex_code <- paste0(latex_code, "\\end{tabular}\n")
  main_caption <- r"(Mean and standard error of performance metrics over $100$ data sets.)"
  latex_code <- paste0(latex_code, sprintf("\\caption{%s}\n", main_caption))
  latex_code <- paste0(latex_code, "\\label{tbl:sim}\n")
  latex_code <- paste0(latex_code, "\\end{table}\n")

  # Print LaTeX code
  cat(latex_code)
}

# main_table("natural")
# main_table("original")
supp_table("natural")
supp_table("original")
