# Verify every frozen input exists and matches its pinned SHA-256 (including
# the multi-GB predictor stacks, so this takes a few minutes).

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "..", "R", "common.R"))
cfg <- initialise_final()

status <- dplyr::bind_rows(lapply(names(cfg$inputs), function(name) {
  entry <- cfg$inputs[[name]]
  exists <- file.exists(entry$path)
  actual <- if (exists) sha256_file(entry$path) else NA_character_
  tibble::tibble(
    input = name, path = entry$path, exists = exists,
    bytes = if (exists) file.info(entry$path)$size else NA_real_,
    checksum_ok = exists && identical(actual, entry$sha256)
  )
}))
write_table(status, cfg, "input_check.csv")
print(status[, c("input", "exists", "checksum_ok")], n = Inf)
if (!all(status$checksum_ok)) stop("Frozen inputs are missing or changed: ", paste(status$input[!status$checksum_ok], collapse = ", "))
log_message("All ", nrow(status), " frozen inputs verified.")
