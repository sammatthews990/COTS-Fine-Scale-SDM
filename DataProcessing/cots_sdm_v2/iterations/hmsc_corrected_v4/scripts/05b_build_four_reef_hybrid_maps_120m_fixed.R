args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
v4_dir <- normalizePath(file.path(dirname(script_file), ".."), winslash = "/")
v2_dir <- normalizePath(file.path(v4_dir, "..", ".."), winslash = "/")

# The original v4 map entry point omitted this helper from its isolated
# bootstrap. Load it explicitly, then execute the unchanged mapping workflow.
source(file.path(v2_dir, "R", "gaussian_terrain.R"))
source(file.path(v4_dir, "scripts", "05_build_four_reef_hybrid_maps_120m.R"))
