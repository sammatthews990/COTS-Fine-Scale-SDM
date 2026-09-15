args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
v5_dir <- normalizePath(file.path(dirname(script_file), ".."), winslash = "/")
source_file <- file.path(v5_dir, "scripts", "01_prepare_shared_state_data_120m.R")
code <- readLines(source_file, warn = FALSE)
code <- sub(
  '    \\.data\\$ecological_eligible,',
  '    (.data$ecological_eligible | .data$method == "cull"),',
  code
)
code <- sub(
  '    direct_detection = observed_any_v4\\(\\.data\\$direct_detection\\),',
  '    direct_detection = observed_any_v4(ifelse(.data$method == "cull", .data$high_density, .data$direct_detection)),',
  code
)
eval(parse(text = code, keep.source = TRUE), envir = globalenv())
