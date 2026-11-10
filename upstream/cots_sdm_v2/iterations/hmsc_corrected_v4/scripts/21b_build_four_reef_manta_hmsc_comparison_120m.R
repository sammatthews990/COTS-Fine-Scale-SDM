args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
v4_dir <- normalizePath(file.path(dirname(script_file), ".."), winslash = "/")
source_file <- file.path(v4_dir, "scripts", "21_build_four_reef_manta_hmsc_comparison_120m.R")
code <- readLines(source_file, warn = FALSE)
code <- sub(
  "dimnames = list\\(NULL, missing\\)\\)\\)\\)\\)\\)",
  "dimnames = list(NULL, missing)))",
  code
)
eval(parse(text = code, keep.source = TRUE), envir = globalenv())
