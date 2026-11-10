args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
v4_dir <- normalizePath(file.path(dirname(script_file), ".."), winslash = "/")
source_file <- file.path(v4_dir, "scripts", "08_validate_preliminary_hmsc_120m.R")
code <- readLines(source_file, warn = FALSE)
code <- gsub("120m", "200m", code, fixed = TRUE)
code <- sub(
  'c\\("enhanced_cots", "enhanced_joint", "clean_joint"\\)',
  'c("enhanced_cots_200support_pca2")',
  code
)
eval(parse(text = code, keep.source = TRUE), envir = globalenv())
