# Stage the shared data store (one-off migration, 10 November 2026).
#
# Question: can the pipeline be rerun from a single, checksummed data store on
# the project SharePoint instead of files scattered through the repository?
# Method: copy every file the final pipeline reads (the frozen inputs pinned in
# config.R at tag v1.0-2027-products, plus the few read directly), the ADRIA
# metacommunity snapshot behind scripts 22/25, and the delivered products into
# COTS_DATA_DIR, keeping a fixed layout:
#   raw/         source data as received (cull/manta workbook, cull-site KMZ,
#                reef forecasts, reef features, geomorphic map, larval
#                concentrations, ReefScan)
#   predictors/  the two 30 m predictor stacks
#   frozen/      pinned intermediates from the exploratory work, cached
#                inputs (eReefs extract) and the ADRIA/COTSMod snapshot
#   products/    final habitat layers and 2027 cull-site rankings
# Validation: every copy is checked against its pinned SHA-256 (or ADRIA's own
# metadata checksum, or the source file's checksum when nothing is pinned).
# Files already present with the right checksum are skipped, so the script can
# be rerun. Writes DATA_MANIFEST.csv to the store and to final/.
#
# Run (from the repository root; ADRIA_DIR = the ADRIA.jl checkout holding the
# 21 October calibration runs):
#   $env:ADRIA_DIR = "C:/Users/smatthew/Desktop/GitKraken/ADRIA.jl"
#   Rscript pipeline/scripts/00a_stage_data_store.R  (run before the move: DataProcessing/cots_sdm_v2/final/...)

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "..", "R", "common.R"))
cfg <- initialise_final()
args <- script_args()
dry_run <- identical(args$dry_run, "1")

data_root <- Sys.getenv("COTS_DATA_DIR", file.path(
  Sys.getenv("USERPROFILE"), "Australian Institute of Marine Science", "AIMS Decision Support - General",
  "01_PROJECTS", "CCIP Projects", "CCIP-BR-13-Fine Scale", "Data"))
if (!dir.exists(data_root)) stop("Data store not found: ", data_root)

# Source locations: the config.R frozen at the tag, so the script still knows
# where everything came from after config.R points at the store.
tag_config <- tempfile(fileext = ".R")
system2("git", c("-C", shQuote(cfg$repo_root), "show",
                 "v1.0-2027-products:DataProcessing/cots_sdm_v2/final/config.R"), stdout = tag_config)
old <- new.env()
sys.source(tag_config, envir = old)
old_cfg <- old$final_config(cfg$final_dir)
dp <- file.path(cfg$repo_root, "DataProcessing")
adria_root <- normalizePath(Sys.getenv("ADRIA_DIR", file.path(dirname(cfg$repo_root), "ADRIA.jl")), "/")
cotsmod_root <- file.path(dirname(cfg$repo_root), "COTSMod-Calibration")

# Store path for each pinned input.
pinned_store <- c(
  predictor_stack_clean = "predictors/predictors_clean_12layer.tif",
  predictor_stack_rg = "predictors/predictors_terra_30m_full_extended.tif",
  geomorphic = "raw/geomorphic/GBR10 GBRMP Geomorphic.tif",
  reef_features = "raw/reef_features/Great_Barrier_Reef_Features.gpkg",
  manta_workbook = "raw/cull_manta/260529-COTS-Manta-Cull-RHIS-Lawrence-CSIRO.xlsx",
  cull_sites_kmz = "raw/cull_sites/Eotr_CotsCullSites_2025_11_19_1_58_PM_1.kmz",
  reef_forecast_2026 = "raw/reef_forecast/gbrPreds2026Adj_20262408.csv",
  reef_forecast_hindcast = "raw/reef_forecast/gbrPredsAdj_20262408.csv",
  larval_report_day10_shp = "raw/larval_concentration/fngbr_reef_daily_means_day10.shp",
  larval_report_day10_dbf = "raw/larval_concentration/fngbr_reef_daily_means_day10.dbf",
  larval_report_day10_shx = "raw/larval_concentration/fngbr_reef_daily_means_day10.shx",
  larval_report_day10_prj = "raw/larval_concentration/fngbr_reef_daily_means_day10.prj",
  legacy_clean_model = "frozen/legacy_sdm/fitted_model_clean_20260820_1041.rds",
  legacy_clean_raster = "frozen/legacy_sdm/COTS_prob_0.02_cpue_year2026_clean_20260820_1041.tif",
  legacy_clean_metrics = "frozen/legacy_sdm/model_performance_metrics_clean_20260820_1041.rds",
  legacy_rg_model = "frozen/legacy_sdm/fitted_model_withRG_20260820_1217.rds",
  legacy_rg_metrics = "frozen/legacy_sdm/model_performance_metrics_withRG_20260820_1217.rds",
  hmsc_enhanced_joint = "frozen/hmsc_v4/enhanced_joint_fit_data_120m.rds",
  hmsc_clean_joint = "frozen/hmsc_v4/clean_joint_fit_data_120m.rds",
  hmsc_clean_native_joint = "frozen/hmsc_v4/clean_native_joint_fit_data_120m.rds",
  hmsc_reference_cells = "frozen/hmsc_v4/corrected_cell_units_120m.rds",
  manta_tows = "frozen/hmsc_v4/manta_footprints_120m.rds",
  footprint_links = "frozen/hmsc_v3/footprint_links_120m.rds",
  assessment_dives = "frozen/prepared/cull_observations_clean.rds",
  ranking_train_site_years = "frozen/prepared/site_year_all.rds",
  assessment_sites = "frozen/final_pins/cull_sites_designated.gpkg",
  ranking_site_sample_points = "frozen/final_pins/cull_sites_points.rds",
  assessment_evidence = "frozen/final_pins/cull_site_evidence_scores.rds",
  assessment_manta_sites = "frozen/final_pins/manta_response_site_ranks.rds",
  assessment_forward = "frozen/final_pins/site_ranking_forward.rds",
  assessment_forward_coefficients = "frozen/final_pins/site_ranking_forward_coefficients.csv",
  ranking_deployed_habitat = "frozen/final_pins/cots_legacy_rg_year2026_clean_fallback_30m.tif",
  larval_report_summary = "frozen/connectivity_lizard/reef_concentration_summary.csv",
  larval_report_zones = "frozen/connectivity_lizard/reef_zones.csv",
  larval_report_releases = "frozen/connectivity_lizard/zone_by_release.csv",
  larval_report_latitude = "frozen/connectivity_lizard/latitude_trend.csv",
  report_metapop_domain = "frozen/cotsmod/lizard_cluster.gpkg",
  # The 15 October observations (22); ADRIA's sandbox copy has since been
  # replaced by the version the 21 October snapshot uses.
  mc_observations = "frozen/cotsmod/reef_cots_20261015.csv"
)
source_override <- c(
  mc_observations = file.path(cotsmod_root, "store/data/archive_ltmp/0d8871dd90a5/reef_cots.csv")
)
rel_to <- function(path, root) {
  path <- normalizePath(path, "/", mustWork = FALSE)
  root <- paste0(normalizePath(root, "/", mustWork = FALSE), "/")
  if (!startsWith(path, root)) stop(path, " is not under ", root)
  substring(path, nchar(root) + 1L)
}
# ADRIA run paths exceed Windows' 260-character limit under the SharePoint
# root, so the snapshot is one zip (members keep ADRIA-relative paths;
# initialise_final() extracts it to cache/adria). Manifest rows for members
# read "<zip>::<member>".
adria_zip <- "frozen/adria_snapshot_20261021.zip"
adria_store <- function(path) paste0(adria_zip, "::", rel_to(path, adria_root))

# ADRIA inputs: the tagged config.R resolves them under ADRIA_DIR too.
old_adria <- adria_root
adria_inputs <- setdiff(names(old_cfg$inputs), names(pinned_store))
not_adria <- adria_inputs[!startsWith(normalizePath(vapply(old_cfg$inputs[adria_inputs], `[[`, "", "path"), "/", mustWork = FALSE),
                                      normalizePath(old_adria, "/", mustWork = FALSE))]
if (length(not_adria)) stop("No store path for: ", paste(not_adria, collapse = ", "))
adria_src <- file.path(adria_root, vapply(adria_inputs, function(n) rel_to(old_cfg$inputs[[n]]$path, old_adria), ""))

files <- dplyr::bind_rows(
  tibble::tibble(input = names(pinned_store), group = "pinned",
                 source = vapply(names(pinned_store), function(n) {
                   if (n %in% names(source_override)) source_override[[n]] else old_cfg$inputs[[n]]$path
                 }, ""),
                 store = unname(pinned_store),
                 sha256 = vapply(names(pinned_store), function(n) old_cfg$inputs[[n]]$sha256, "")),
  tibble::tibble(input = adria_inputs, group = "pinned", source = adria_src,
                 store = vapply(adria_src, adria_store, ""),
                 sha256 = vapply(adria_inputs, function(n) old_cfg$inputs[[n]]$sha256, ""))
)

# Read directly by final scripts (not pinned before; pinned in config.R now).
direct <- tibble::tibble(
  input = c("site_year_first_visit", "annual_rank_predictions"),
  group = "direct",
  source = c(file.path(dp, "cots_sdm_v2/artifacts/site_year_first_visit.rds"),
             file.path(dp, "cots_sdm_v2/outputs/annual_rank_predictions.rds")),
  store = c("frozen/prepared/site_year_first_visit.rds", "frozen/prepared/annual_rank_predictions.rds"),
  sha256 = NA_character_)
folder <- function(src_dir, store_dir, group, pattern = NULL) {
  f <- list.files(src_dir, pattern = pattern, full.names = TRUE)
  tibble::tibble(input = NA_character_, group = group, source = f,
                 store = file.path(store_dir, basename(f)), sha256 = NA_character_)
}
sidecars <- dplyr::bind_rows(
  folder(file.path(dp, "data/GBR10_GBRMP_Geomorphic"), "raw/geomorphic", "sidecar"),
  folder(file.path(dp, "data/fngbr_LarvalConcentrations"), "raw/larval_concentration", "sidecar", "day10\\.cpg$"),
  folder(file.path(cfg$repo_root, "data/ReefScan"), "raw/reefscan", "direct", "\\.csv$"),
  folder(file.path(cfg$final_dir, "cache/ereefs"), "frozen/ereefs", "cache")
)

# ADRIA snapshot: the two report folders plus every file their metadata lists
# (verified against the metadata's own checksums).
cal_dir <- file.path(adria_root, rel_to(dirname(old_cfg$inputs$calibration_report_20261021$path), old_adria))
proj_dir <- file.path(dirname(cal_dir), "projections_2050")
adria_files <- dplyr::bind_rows(folder(cal_dir, adria_store(cal_dir), "adria"),
                                folder(proj_dir, adria_store(proj_dir), "adria"))
cm <- jsonlite::read_json(file.path(cal_dir, "metadata.json"))
run_dir <- function(x) file.path(adria_root, gsub("\\\\", "/", x))
cal_input_path <- function(key) {
  if (key %in% c("reef_cots.csv", "reef_manta.csv", "reef_photo_transect.csv")) return(file.path(adria_root, "sandbox/data", key))
  if (key == "targets_reef_year.csv") return(file.path(adria_root, "sandbox/data/derived/density_targets_2020s", key))
  if (startsWith(key, "cull/")) return(file.path(run_dir(cm$cull_run), sub("^cull/", "", key)))
  file.path(run_dir(cm$run), key)
}
pm <- jsonlite::read_json(file.path(proj_dir, "metadata.json"))
meta_inputs <- dplyr::bind_rows(
  tibble::tibble(source = vapply(names(cm$inputs), cal_input_path, ""), sha256 = unlist(cm$inputs)),
  tibble::tibble(source = file.path(adria_root, "sandbox/calibration/runs", names(pm$inputs)), sha256 = unlist(pm$inputs))
)
meta_inputs$input <- NA_character_
meta_inputs$group <- "adria"
meta_inputs$store <- vapply(meta_inputs$source, adria_store, "")

# Delivered products (not inputs; recorded with their checksums).
products <- dplyr::bind_rows(
  folder(cfg$dirs$layers, "products/habitat_layers", "product", "\\.(tif|csv)$"),
  folder(file.path(cfg$dirs$outputs, "site_ranking"), "products/site_ranking", "product")
)

files <- dplyr::bind_rows(files, direct, sidecars, adria_files, meta_inputs, products) |>
  dplyr::distinct(store, .keep_all = TRUE)
files[] <- lapply(files, unname)
# SharePoint rewrites the library-metadata parts of Office files after upload,
# which changes their checksums; store each one inside its own zip.
office <- grepl("\\.(xlsx|xlsm|docx|pptx)$", files$store)
files$store[office] <- paste0(files$store[office], ".zip::", basename(files$store[office]))
absent <- files$source[!file.exists(files$source)]
if (length(absent)) stop("Missing source files:\n  ", paste(absent, collapse = "\n  "))
files$bytes <- file.size(files$source)
log_message(nrow(files), " files, ", round(sum(files$bytes) / 2^30, 2), " GB to stage in ", data_root)
print(dplyr::summarise(dplyr::group_by(files, top = sub("/.*", "", store)), files = dplyr::n(), GB = round(sum(bytes) / 2^30, 2)))
if (dry_run) quit(save = "no")

# Zipped files ("<zip>::<member>"): verify every member at its source, zip
# once, then check every member again after extracting the stored zip. The
# zip itself gets a manifest row (and a pin in config.R).
in_zip <- grepl("::", files$store, fixed = TRUE)
for (i in which(in_zip)) {
  actual <- sha256_file(files$source[i])
  if (!is.na(files$sha256[i]) && !identical(actual, files$sha256[i])) stop("Source changed since pinned: ", files$source[i])
  files$sha256[i] <- actual
}
zips <- unique(sub("::.*$", "", files$store[in_zip]))
zip_rows <- lapply(zips, function(z) {
  rows <- which(in_zip & startsWith(files$store, paste0(z, "::")))
  member <- sub("^.*::", "", files$store[rows])
  zip_dest <- file.path(data_root, z)
  if (!file.exists(zip_dest)) {
    stage <- tempfile("zip_stage")
    for (k in seq_along(rows)) {
      to <- file.path(stage, member[k])
      dir.create(dirname(to), recursive = TRUE, showWarnings = FALSE)
      file.copy(files$source[rows[k]], to, copy.date = TRUE)
    }
    tmp_zip <- tempfile(fileext = ".zip")
    zip::zip(tmp_zip, files = member, root = stage)
    dir.create(dirname(zip_dest), recursive = TRUE, showWarnings = FALSE)
    if (!file.copy(tmp_zip, zip_dest)) stop("Copy failed: ", zip_dest)
    log_message("Wrote ", z, " (", length(member), " members)")
  }
  check_dir <- tempfile("zip_check")
  utils::unzip(zip_dest, exdir = check_dir)
  bad <- member[vapply(file.path(check_dir, member), sha256_file, "") != files$sha256[rows]]
  if (length(bad)) stop(z, ": members differ from their sources: ", paste(bad, collapse = ", "))
  is_adria <- identical(z, adria_zip)
  tibble::tibble(input = if (is_adria) "adria_snapshot" else paste0(files$input[rows[1]], "_zip"), group = "pinned",
                 source = if (is_adria) "ADRIA.jl (members below)" else "zip of the member below",
                 store = z, sha256 = sha256_file(zip_dest), bytes = file.size(zip_dest))
})
files <- dplyr::bind_rows(files, zip_rows)
in_zip <- grepl("::", files$store, fixed = TRUE)

for (i in which(!in_zip & !files$store %in% zips)) {
  src <- files$source[i]
  dest <- file.path(data_root, files$store[i])
  expected <- files$sha256[i]
  if (is.na(expected)) expected <- sha256_file(src)
  if (!file.exists(dest) || file.size(dest) != files$bytes[i] || !identical(sha256_file(dest), expected)) {
    dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
    if (!file.copy(src, dest, overwrite = TRUE, copy.date = TRUE)) stop("Copy failed: ", src)
    actual <- sha256_file(dest)
    if (!identical(actual, expected)) stop("Checksum mismatch after copying ", src, "\n  expected ", expected, "\n  actual   ", actual)
    log_message("Copied ", files$store[i])
  }
  files$sha256[i] <- expected
}

manifest <- files |>
  dplyr::mutate(source = vapply(source, function(p) {
    for (root in c(cfg$repo_root, adria_root, cotsmod_root)) if (startsWith(p, root)) return(file.path(basename(root), rel_to(p, root)))
    p
  }, "")) |>
  dplyr::select(store, input, group, bytes, sha256, source) |>
  dplyr::arrange(store)
utils::write.csv(manifest, file.path(data_root, "DATA_MANIFEST.csv"), row.names = FALSE)
utils::write.csv(manifest, file.path(cfg$final_dir, "DATA_MANIFEST.csv"), row.names = FALSE)
log_message("Staged and verified ", nrow(manifest), " files (", round(sum(manifest$bytes) / 2^30, 2), " GB).")
