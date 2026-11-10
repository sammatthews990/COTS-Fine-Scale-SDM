args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
v4_dir <- normalizePath(file.path(dirname(script_file), ".."), winslash = "/")
source_file <- file.path(v4_dir, "scripts", "13_assemble_support_predictors_200m.R")
code <- readLines(source_file, warn = FALSE)
code <- sub(
  'gaussian_cache_v2_union\\.rds',
  'gaussian_cache_v2_bound.rds',
  code
)
eval(parse(text = code, keep.source = TRUE), envir = globalenv())
