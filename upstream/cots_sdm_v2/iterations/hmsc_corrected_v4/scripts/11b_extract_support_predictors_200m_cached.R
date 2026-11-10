args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
v4_dir <- normalizePath(file.path(dirname(script_file), ".."), winslash = "/")
source_file <- file.path(v4_dir, "scripts", "11_extract_support_predictors_200m.R")
code <- readLines(source_file, warn = FALSE)
replacement <- c(
  'base_cache_dir <- file.path(cfg$paths$tmp_dir, "base_exact_200m")',
  'dir.create(base_cache_dir, recursive = TRUE, showWarnings = FALSE)',
  'base_rows <- lapply(seq_along(chunks), function(i) {',
  '  cache_file <- file.path(base_cache_dir, sprintf("chunk_%03d.rds", i))',
  '  if (file.exists(cache_file)) return(readRDS(cache_file))',
  '  message("Base polygon chunk ", i, " of ", length(chunks))',
  '  rows <- chunks[[i]]',
  '  polygons <- terra::vect(cell_polygons(cells[rows, , drop = FALSE], terra::crs(base)))',
  '  value <- terra::extract(base, polygons, fun = mean, na.rm = TRUE)',
  '  value <- value[, setdiff(names(value), "ID"), drop = FALSE]',
  '  result <- dplyr::bind_cols(cells[rows, , drop = FALSE], tibble::as_tibble(value))',
  '  saveRDS(result, cache_file)',
  '  result',
  '})',
  'base_values <- dplyr::bind_rows(base_rows)',
  'saveRDS(base_values, file.path(cfg$paths$artifacts_dir, "base_support_predictors_200m.rds"))'
)
code <- c(code[seq_len(32L)], replacement, code[42L:length(code)])
eval(parse(text = code, keep.source = TRUE), envir = globalenv())
