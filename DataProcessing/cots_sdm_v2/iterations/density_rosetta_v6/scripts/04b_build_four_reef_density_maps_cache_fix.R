args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
v6_dir <- normalizePath(file.path(dirname(script_file), ".."), winslash = "/")
source_file <- file.path(v6_dir, "scripts", "04_build_four_reef_density_maps.R")
code <- readLines(source_file, warn = FALSE)
code <- sub(
  'cache_dir = file.path\\(v4_dir, "tmp", "map_terrain_60m", slug\\)',
  'cache_dir = file.path(v6_dir, "tmp", "map_terrain_60m", slug)',
  code
)
eval(parse(text = code, keep.source = TRUE), envir = globalenv())
