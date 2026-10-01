# Shared helpers for the final COTS fine-scale SDM pipeline.

find_repo_root <- function(start = getwd()) {
  path <- normalizePath(start, winslash = "/", mustWork = TRUE)
  repeat {
    if (dir.exists(file.path(path, ".git"))) return(path)
    parent <- dirname(path)
    if (identical(parent, path)) stop("Could not find the repository root above ", start)
    path <- parent
  }
}

current_script_file <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) != 1L) stop("Run final pipeline scripts with Rscript.")
  normalizePath(sub("^--file=", "", file_arg), winslash = "/", mustWork = TRUE)
}

# Parse --key=value command-line arguments into a named list.
script_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  args <- args[grepl("^--[A-Za-z0-9_]+=", args)]
  keys <- sub("^--([A-Za-z0-9_]+)=.*$", "\\1", args)
  stats::setNames(as.list(sub("^--[A-Za-z0-9_]+=", "", args)), keys)
}

initialise_final <- function() {
  script_file <- current_script_file()
  final_dir <- normalizePath(file.path(dirname(script_file), ".."), winslash = "/")
  source(file.path(final_dir, "config.R"), local = FALSE)
  cfg <- final_config(final_dir)
  for (f in c("hmsc.R", "support_grid.R", "rasters.R", "figures.R", "context.R")) {
    source(file.path(final_dir, "R", f), local = FALSE)
  }
  invisible(lapply(cfg$dirs, dir.create, recursive = TRUE, showWarnings = FALSE))
  log_message("Script: ", basename(script_file))
  cfg
}

log_message <- function(...) {
  message(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), ...)
}

sha256_file <- function(path) {
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}

# Fail early if a frozen input is missing or has changed since it was frozen.
verify_frozen_input <- function(path, expected_sha256 = NA_character_) {
  if (!file.exists(path)) stop("Missing frozen input: ", path)
  if (!is.na(expected_sha256)) {
    actual <- sha256_file(path)
    if (!identical(actual, expected_sha256)) {
      stop("Checksum mismatch for ", path, "\n  expected ", expected_sha256, "\n  actual   ", actual)
    }
  }
  invisible(path)
}

write_table <- function(x, cfg, name) {
  path <- file.path(cfg$dirs$tables, name)
  utils::write.csv(x, path, row.names = FALSE)
  log_message("Wrote ", path)
  invisible(path)
}

# Upsert one row per layer into the deliverable layer manifest.
register_layer <- function(cfg, layer_id, path, description, model, units, extra = list()) {
  manifest_path <- file.path(cfg$dirs$layers, "layer_manifest.csv")
  r <- terra::rast(path)
  row <- tibble::tibble(
    layer_id = layer_id,
    file = basename(path),
    description = description,
    model = model,
    units = units,
    resolution_m = terra::res(r)[[1]],
    crs = terra::crs(r, describe = TRUE)$name,
    bytes = file.info(path)$size,
    sha256 = sha256_file(path),
    created = format(Sys.time(), "%Y-%m-%d %H:%M"),
    !!!extra
  )
  manifest <- if (file.exists(manifest_path)) {
    utils::read.csv(manifest_path, stringsAsFactors = FALSE, check.names = FALSE)
  } else NULL
  if (!is.null(manifest)) {
    manifest <- manifest[manifest$layer_id != layer_id, , drop = FALSE]
    manifest <- dplyr::bind_rows(lapply(manifest, as.character) |> tibble::as_tibble(),
                                 dplyr::mutate(row, dplyr::across(dplyr::everything(), as.character)))
  } else {
    manifest <- row
  }
  utils::write.csv(manifest[order(manifest$layer_id), ], manifest_path, row.names = FALSE)
  invisible(row)
}
