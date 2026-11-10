# Layer: manta-trained habitat model, within reef-year, 30 m.
#
# The model 17 found best for ranking manta tows (held-out within reef-year
# AUC 0.61 vs 0.55 for the cull-trained layer; no better for cull-site
# selection). Response: manta tow with strong COTS evidence (direct COTS, >= 2
# counted scars, legacy "common" scars) vs strict negative. Predictors: the RG
# stack (mean along each tow), no survey year. Each reef-year's strong-tow
# rate is a fixed offset (xgboost base_margin), so the trees learn only where
# COTS sit within an active reef-year. Fitted here on every labelled tow.
#
#   cots_manta_within_30m.tif
#     band 1  relative log-odds of strong COTS evidence, against the reef-year
#             average (0 = typical for the reef-year; comparable within a reef,
#             not an absolute probability)
#     band 2  percentile within the GBR reef polygon
# NA where the RG layers are missing (same coverage as RG + year).
# Manta tows run along reef perimeters, so values away from tow paths (e.g.
# back reef, lagoon) are extrapolations.

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "..", "R", "common.R"))
cfg <- initialise_final()
strong_classes <- c("direct_cots", "multiple_counted_scars", "legacy_common_scars")
params <- list(objective = "binary:logistic", eta = 0.00137382379588326, max_depth = 5,
               min_child_weight = 20, gamma = 1e-10, nthread = 8, tree_method = "hist")
rounds <- 1500L

stack <- terra::rast(checked_input(cfg, "predictor_stack_rg"))
rg <- cfg$legacy$rg_layer_names
names(stack) <- rg

# --- Labelled tows with mean predictors along each tow (cache from 14) --------
tow_cache <- file.path(cfg$dirs$cache, "manta_tow_predictors.rds")
if (!file.exists(tow_cache)) stop("Run 14_manta_aggregation_cv.R first (needs ", tow_cache, ")")
tows <- readRDS(checked_input(cfg, "manta_tows")) |>
  dplyr::select("visit_id", "year", "reef_name", "manta_evidence_class")
geom <- readRDS(checked_input(cfg, "footprint_links"))$manta |>
  dplyr::group_by(.data$visit_id) |>
  dplyr::summarise(sx = dplyr::first(.data$start_x), ex = dplyr::first(.data$end_x), .groups = "drop")
tows <- dplyr::inner_join(tows, geom, by = "visit_id") |>
  dplyr::filter(is.finite(.data$sx), is.finite(.data$ex)) |>
  dplyr::mutate(label = dplyr::case_when(.data$manta_evidence_class %in% strong_classes ~ 1,
                                         .data$manta_evidence_class == "strict_negative" ~ 0),
                reef_year = paste(.data$reef_name, .data$year)) |>
  dplyr::group_by(.data$reef_year) |>
  dplyr::mutate(ry_margin = stats::qlogis((sum(.data$label %in% 1) + 0.5) / (sum(!is.na(.data$label)) + 1))) |>
  dplyr::ungroup()
tow_points <- readRDS(tow_cache)
v <- as.matrix(tow_points$values[, rg])
ok <- !is.na(v)
v[!ok] <- 0
sums <- rowsum(v, tow_points$tow)
counts <- rowsum(ok * 1, tow_points$tow)
tow_x <- matrix(NA_real_, nrow(tows), length(rg), dimnames = list(NULL, rg))
tow_x[as.integer(rownames(sums)), ] <- ifelse(counts > 0, sums / counts, NA)
train <- !is.na(tows$label)

# --- Fit on every labelled tow --------------------------------------------------------
model_path <- file.path(cfg$dirs$models, "manta_within_reef_year_xgb.json")
model <- xgboost::xgb.train(params = params, nrounds = rounds, verbose = 0,
  data = xgboost::xgb.DMatrix(tow_x[train, ], label = tows$label[train], base_margin = tows$ry_margin[train], missing = NA))
xgboost::xgb.save(model, model_path)
log_message("Fitted on ", format(sum(train), big.mark = ","), " tows (", format(sum(tows$label %in% 1), big.mark = ","),
            " strong evidence) in ", dplyr::n_distinct(tows$reef_year[train]), " reef-years")
write_table(as.data.frame(xgboost::xgb.importance(model = model)), cfg, "manta_within_importance.csv")

# --- Predict the 30 m layer ----------------------------------------------------------
score_path <- file.path(cfg$dirs$cache, "manta_within_score_30m.tif")
if (!file.exists(score_path) || identical(script_args()$force, "1")) {
  log_message("Predicting at 30 m over the RG stack")
  pred_fun <- function(m, d) stats::predict(m, xgboost::xgb.DMatrix(as.matrix(d[, rg]), missing = NA), outputmargin = TRUE)
  terra::predict(stack, model, fun = pred_fun, na.rm = TRUE, filename = score_path, overwrite = TRUE,
                 wopt = list(names = "relative_log_odds", datatype = "FLT4S", gdal = layer_gdal_options))
}
score <- terra::rast(score_path)

# Check the raster against direct predictions.
set.seed(cfg$hmsc$seed)
pts <- terra::spatSample(score, size = 2000, method = "random", na.rm = TRUE, as.points = TRUE, values = FALSE, exhaust = TRUE)
env <- terra::extract(stack, pts, ID = FALSE)
direct <- stats::predict(model, xgboost::xgb.DMatrix(as.matrix(env[, rg]), missing = NA), outputmargin = TRUE)
stored <- terra::extract(score, pts, ID = FALSE)[[1]]
check <- tibble::tibble(layer_id = "manta_within", sampled_cells = length(stored),
                        max_absolute_error = max(abs(direct - stored)), correlation = stats::cor(direct, stored))
write_table(check, cfg, "manta_within_reproduction_check.csv")
if (check$max_absolute_error > 1e-4) stop("Manta layer does not match its model.")

# --- Within-reef percentile and the delivered two-band layer ---------------------
reef_id <- terra::rast(file.path(cfg$dirs$cache, "reef_id_30m.tif"))
if (!terra::compareGeom(score, reef_id, stopOnError = FALSE)) stop("Score and reef-id grids differ.")
percentile_path <- file.path(cfg$dirs$cache, "manta_within_percentile_30m.tif")
log_message("Writing within-reef percentile")
cells <- write_reef_percentile(score, reef_id, percentile_path)
layer_path <- file.path(cfg$dirs$layers, "cots_manta_within_30m.tif")
layer <- c(score, terra::rast(percentile_path))
names(layer) <- c("relative_log_odds", "within_reef_percentile")
write_layer(layer, layer_path)

register_layer(
  cfg, "manta_within", layer_path,
  description = paste("Manta-trained habitat model, within reef-year. Band 1: relative log-odds of strong COTS evidence on a manta tow",
                      "vs the reef-year average; band 2: percentile within the GBR reef polygon. NA outside RG coverage"),
  model = "XGBoost on RG stack, strong-evidence vs strict-negative tows, reef-year offset (scripts 17-18)",
  units = "relative log-odds / percentile",
  extra = list(reefs = dplyr::n_distinct(cells$reef))
)
