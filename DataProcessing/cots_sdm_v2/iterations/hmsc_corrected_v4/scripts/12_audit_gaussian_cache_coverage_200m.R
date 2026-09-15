args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
v4_dir <- normalizePath(file.path(dirname(script_file), ".."), winslash = "/")
v2_dir <- normalizePath(file.path(v4_dir, "..", ".."), winslash = "/")

base <- readRDS(file.path(v4_dir, "artifacts", "base_support_predictors_200m.rds"))
base$cache_cell_id <- sub("^200m\\|\\|", "", base$cell_id)

cache_dirs <- c(
  point = file.path(v2_dir, "tmp", "hmsc_v2_gaussian_point_tiles"),
  polygon = file.path(v2_dir, "tmp", "hmsc_v2_gaussian_tiles"),
  map = file.path(v2_dir, "tmp", "hmsc_v2_map_gaussian_tiles")
)

read_cache <- function(path) {
  files <- list.files(path, pattern = "\\.rds$", full.names = TRUE, recursive = TRUE)
  rows <- lapply(files, function(file) {
    value <- tryCatch(readRDS(file), error = function(e) NULL)
    if (is.null(value) || !"cell_id" %in% names(value)) return(NULL)
    value
  })
  rows <- dplyr::bind_rows(rows)
  rows[!duplicated(rows$cell_id), , drop = FALSE]
}

pieces <- lapply(names(cache_dirs), function(label) {
  value <- read_cache(cache_dirs[[label]])
  cat(label, "files/rows", length(list.files(cache_dirs[[label]], pattern = "\\.rds$", recursive = TRUE)), nrow(value), "\n")
  value$cache_source <- label
  value
})
names(pieces) <- names(cache_dirs)

coverage <- lapply(names(pieces), function(label) {
  tibble::tibble(
    cache = label,
    cache_rows = nrow(pieces[[label]]),
    matched = sum(base$cache_cell_id %in% pieces[[label]]$cell_id),
    unmatched = sum(!base$cache_cell_id %in% pieces[[label]]$cell_id)
  )
}) |>
  dplyr::bind_rows()

all_cache <- dplyr::bind_rows(pieces) |>
  dplyr::arrange(factor(.data$cache_source, levels = names(cache_dirs))) |>
  dplyr::distinct(.data$cell_id, .keep_all = TRUE)
coverage <- dplyr::bind_rows(
  coverage,
  tibble::tibble(
    cache = "union",
    cache_rows = nrow(all_cache),
    matched = sum(base$cache_cell_id %in% all_cache$cell_id),
    unmatched = sum(!base$cache_cell_id %in% all_cache$cell_id)
  )
)

saveRDS(all_cache, file.path(v4_dir, "artifacts", "gaussian_cache_v2_union.rds"))
utils::write.csv(
  coverage,
  file.path(v4_dir, "outputs", "gaussian_cache_coverage_200m.csv"),
  row.names = FALSE
)
print(coverage)
