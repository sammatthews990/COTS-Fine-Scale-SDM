args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
v6_dir <- normalizePath(file.path(dirname(script_file), ".."), winslash = "/")
source_file <- file.path(v6_dir, "scripts", "04_build_four_reef_density_maps.R")
code <- readLines(source_file, warn = FALSE)
code <- gsub(
  "align_raster_to_template(legacy, template)",
  "align_raster_to_template(legacy, rasters$salad_native)",
  code, fixed = TRUE
)
code <- gsub(
  "align_raster_to_template(terra::rast(v2_file), template)",
  "align_raster_to_template(terra::rast(v2_file), rasters$salad_native)",
  code, fixed = TRUE
)
code <- gsub(
  "align_raster_to_template(baseline, template)",
  "align_raster_to_template(baseline, rasters$salad_native)",
  code, fixed = TRUE
)
eval(parse(text = code, keep.source = TRUE), envir = globalenv())
