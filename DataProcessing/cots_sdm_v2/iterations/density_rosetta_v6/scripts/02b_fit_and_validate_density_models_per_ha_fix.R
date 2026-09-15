args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
v6_dir <- normalizePath(file.path(dirname(script_file), ".."), winslash = "/")
source_file <- file.path(v6_dir, "scripts", "02_fit_and_validate_density_models.R")
code <- readLines(source_file, warn = FALSE)
code <- sub(
  "ctx <- initialise_density_v6\\(script_file\\)",
  "ctx <- initialise_density_v6(script_file); source(file.path(v6_dir, 'R', 'rosetta_density_unit_fix.R'))",
  code
)
eval(parse(text = code, keep.source = TRUE), envir = globalenv())
