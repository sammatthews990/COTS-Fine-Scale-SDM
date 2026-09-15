args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
v5_dir <- normalizePath(file.path(dirname(script_file), ".."), winslash = "/")
source_file <- file.path(v5_dir, "scripts", "02_fit_shared_state_120m.R")
code <- readLines(source_file, warn = FALSE)
code <- sub(
  "psi <- as.numeric\\(report\\$psi\\)",
  "state_par <- fit$par[grepl(\"^beta_state\", names(fit$par))]; psi <- stats::plogis(as.numeric(artifact$data$X %*% state_par))",
  code
)
eval(parse(text = code, keep.source = TRUE), envir = globalenv())
