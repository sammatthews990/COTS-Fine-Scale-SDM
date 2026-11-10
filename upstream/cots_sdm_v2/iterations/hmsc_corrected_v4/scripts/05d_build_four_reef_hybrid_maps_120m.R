args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
v4_dir <- normalizePath(file.path(dirname(script_file), ".."), winslash = "/")
v2_dir <- normalizePath(file.path(v4_dir, "..", ".."), winslash = "/")
source(file.path(v2_dir, "R", "gaussian_terrain.R"))
source(file.path(v2_dir, "R", "mapping.R"))
map_file <- file.path(v4_dir, "scripts", "05_build_four_reef_hybrid_maps_120m.R")
code <- readLines(map_file, warn = FALSE)
replacement <- c(
  "fallback_calibration_data <- data.frame(",
  "  y = stats::qnorm(clip(enhanced_overlap[calibration_keep])),",
  "  x = stats::qnorm(clip(clean_overlap[calibration_keep]))",
  ")",
  "fallback_calibration <- stats::lm(y ~ x, data = fallback_calibration_data)",
  "calibrate_clean <- function(probability) {",
  "  result <- rep(NA_real_, length(probability))",
  "  keep <- is.finite(probability)",
  "  result[keep] <- stats::pnorm(stats::predict(",
  "    fallback_calibration,",
  "    newdata = data.frame(x = stats::qnorm(clip(probability[keep])))",
  "  ))",
  "  result",
  "}"
)
code <- c(code[seq_len(61L)], replacement, code[78L:length(code)])
eval(parse(text = code, keep.source = TRUE), envir = globalenv())
