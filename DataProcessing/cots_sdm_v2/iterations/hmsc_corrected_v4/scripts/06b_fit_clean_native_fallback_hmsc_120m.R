args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
v4_dir <- normalizePath(file.path(dirname(script_file), ".."), winslash = "/")
fit_file <- file.path(v4_dir, "scripts", "06_fit_clean_native_fallback_hmsc_120m.R")
code <- readLines(fit_file, warn = FALSE)
code <- gsub(
  '  "BPI", "EAST", "HCU", "NORTH", "SLO", "SVF", "VCU", "VRM", "DEM",',
  '  "EAST", "HCU", "NORTH", "SLO", "SVF", "VCU", "VRM", "DEM",',
  code, fixed = TRUE
)
eval(parse(text = code, keep.source = TRUE), envir = globalenv())
