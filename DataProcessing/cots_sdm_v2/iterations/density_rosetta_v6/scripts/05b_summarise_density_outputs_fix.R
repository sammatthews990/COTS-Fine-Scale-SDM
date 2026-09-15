args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
v6_dir <- normalizePath(file.path(dirname(script_file), ".."), winslash = "/")
source_file <- file.path(v6_dir, "scripts", "05_summarise_density_outputs.R")
code <- readLines(source_file, warn = FALSE)
code <- sub(
  "coefficients[[length(coefficients) + 1L]] <- tibble::tibble(",
  "beta_values <- as.numeric(model$beta); coefficients[[length(coefficients) + 1L]] <- tibble::tibble(",
  code, fixed = TRUE
)
code <- sub(
  "coefficient = as.numeric(model$beta)", "coefficient = beta_values",
  code, fixed = TRUE
)
eval(parse(text = code, keep.source = TRUE), envir = globalenv())
