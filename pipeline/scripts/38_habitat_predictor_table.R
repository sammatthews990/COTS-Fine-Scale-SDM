# Question: which variables enter each habitat model, where do they come from,
# and how much does each model use them? (Supplementary Table S1 of the
# Frontiers manuscript; Appendix B of the report.)
# Method: a fixed description of every predictor, taken from the scripts that
# built the stacks (upstream/predictors) and the Hmsc v4 fits
# (upstream/cots_sdm_v2/iterations/hmsc_corrected_v4), joined to each model's
# gain share (tables from 01-02 and 18), Hmsc component membership, value
# ranges sampled from the pinned 30 m stack, and training sizes and tuned
# settings read from the frozen XGBoost models.
# Validation: nothing is refitted; inputs are the pinned stack and models.

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "..", "R", "common.R"))
cfg <- initialise_final()
suppressPackageStartupMessages(library(tidymodels))

# Sources are shared within a group (Table S1 prints each once, on the group's
# header row); processing holds only what is specific to a variable. Elevation
# and turbidity sources confirmed by the author on 12 November 2026; the
# turbidity layer stores FNU x 10, as the Allen Coral Atlas distributes it.
group_source <- c(
  `Terrain (30 m)` = paste(
    "Allen Coral Atlas satellite-derived depth at 30 m (Li et al. 2021); in the Torres Strait, Allen Coral Atlas 15 m depth",
    "averaged to 30 m and preferred over the 30 m DeepReef model (Beaman 2017). Derivatives as in Guillaume et al. (2025);",
    "bilinear to the 30 m grid"),
  `Coral habitat` = paste("MaxEnt habitat suitability (Guillaume et al. 2025); newer predictions merged over the earlier 15 m",
                          "predictions, bilinear to 30 m"),
  ReefGuidance = paste("ReefGuidance processed layers (19 November 2025), accessed through ReefGuide; regional tiles",
                       "reprojected (bilinear) to the 30 m grid and merged"),
  Time = "Cull records",
  Geomorphology = "GBR10 geomorphic map (GBRMPA 2021; Roelfsema et al. 2021)",
  `Terrain (multiscale)` = paste(
    "Computed from the 30 m depth model: Gaussian smoothing with sigma 0, 30, 60, 120 and 180 m (nominal scales 30-360 m),",
    "each derivative recomputed at every scale, then the first principal component across the five scales"))
predictors <- tibble::tribble(
  ~code, ~variable, ~group, ~units, ~processing,
  "DEM", "Depth", "Terrain (30 m)", "m", NA,
  "SLO", "Slope", "Terrain (30 m)", "radians", NA,
  "EAST", "Eastness (sine of aspect)", "Terrain (30 m)", "-1 to 1", NA,
  "NORTH", "Northness (cosine of aspect)", "Terrain (30 m)", "-1 to 1", NA,
  "BPI", "Bathymetric position index", "Terrain (30 m)", "m", NA,
  "VRM", "Vector ruggedness measure", "Terrain (30 m)", "0 to 1", NA,
  "SVF", "Sky-view factor", "Terrain (30 m)", "0 to 1", NA,
  "HCU", "Horizontal (plan) curvature", "Terrain (30 m)", "1/m", NA,
  "VCU", "Vertical (profile) curvature", "Terrain (30 m)", "1/m", NA,
  "AhyaD", "Acropora hyacinthus suitability", "Coral habitat", "0 to 1", NA,
  "Aspat", "Acropora spathulata suitability", "Coral habitat", "0 to 1", NA,
  "Aten", "Acropora tenuis suitability", "Coral habitat", "0 to 1", NA,
  "RG_bathy", "Depth (ReefGuidance)", "ReefGuidance", "m", NA,
  "RG_slope", "Slope (ReefGuidance)", "ReefGuidance", "degrees", NA,
  "RG_turbid", "Turbidity", "ReefGuidance", "FNU x 10",
    "Allen Coral Atlas Sentinel-2 turbidity (Li et al. 2022); a value of 55 is 5.5 FNU",
  "RG_waves_Hs", "Significant wave height (Hs)", "ReefGuidance", "m", NA,
  "RG_waves_Tp", "Peak wave period (Tp)", "ReefGuidance", "s", NA,
  "Year", "Survey year", "Time", "year", "Set to 2026 for mapping",
  "geomorphic_class", "Geomorphic zone", "Geomorphology", "class", "11 levels; class at the 120 m cell centre",
  "dem_ms", "Depth (multiscale)", "Terrain (multiscale)", "PC score", NA,
  "slope_ms", "Slope (multiscale)", "Terrain (multiscale)", "PC score", NA,
  "eastness_ms", "Eastness (multiscale)", "Terrain (multiscale)", "PC score", NA,
  "northness_ms", "Northness (multiscale)", "Terrain (multiscale)", "PC score", NA,
  "hcu_ms", "Plan curvature (multiscale)", "Terrain (multiscale)", "PC score", NA,
  "vcu_ms", "Profile curvature (multiscale)", "Terrain (multiscale)", "PC score", NA,
  "bpi_ms", "Bathymetric position (multiscale)", "Terrain (multiscale)", "PC score",
    "Depth minus the mean within max(60 m, 1.5 x scale)",
  "sd_ms", "Local relief (multiscale)", "Terrain (multiscale)", "PC score",
    "Standard deviation of depth within max(60 m, 1.5 x scale)"
)
predictors$source <- unname(group_source[predictors$group])

# Model membership. XGBoost layers use the 30 m values; Hmsc uses 4 x 4 block
# means of the 30 m layers at 120 m (hmsc_corrected_v4 scripts 04, 04b, 06b).
clean12 <- cfg$legacy$rg_layer_names[1:12]
rg17 <- cfg$legacy$rg_layer_names
multiscale <- predictors$code[endsWith(predictors$code, "_ms")]
members <- list(
  clean_year = c(clean12, "Year"),
  rg_year = c(rg17, "Year"),
  manta_within = rg17,
  hmsc_enhanced = c("geomorphic_class", "AhyaD", "Aspat", "Aten", "RG_turbid", "RG_waves_Hs", "RG_waves_Tp", multiscale),
  hmsc_clean = c("AhyaD", "Aspat", "Aten", multiscale),
  hmsc_clean_native = c(setdiff(clean12, "BPI"))
)
for (m in names(members)) predictors[[paste0("in_", m)]] <- predictors$code %in% members[[m]]

# Gain shares from the fitted XGBoost models.
gain <- function(table, layer = NULL, feature = "feature", value = "gain_share") {
  d <- read.csv(file.path(cfg$dirs$tables, table), stringsAsFactors = FALSE)
  if (!is.null(layer)) d <- d[d$layer_id == layer, ]
  stats::setNames(d[[value]] / sum(d[[value]]), d[[feature]])
}
g_clean <- gain("legacy_variable_importance.csv", "legacy_clean_year2026")
g_rg <- gain("legacy_variable_importance.csv", "legacy_rg_year2026")
g_manta <- gain("manta_within_importance.csv", feature = "Feature", value = "Gain")
predictors$gain_clean_year <- unname(g_clean[predictors$code])
predictors$gain_rg_year <- unname(g_rg[predictors$code])
predictors$gain_manta_within <- unname(g_manta[predictors$code])

# Values over reef cells (5th, 50th and 95th percentiles of a regular sample).
stack <- terra::rast(checked_input(cfg, "predictor_stack_rg"))
names(stack) <- rg17
set.seed(cfg$hmsc$seed)
sample_values <- terra::spatSample(stack, size = 4000000, method = "regular", na.rm = FALSE)
q <- t(vapply(rg17, function(v) {
  x <- sample_values[[v]]
  x <- x[is.finite(x)]
  c(cells = length(x), stats::quantile(x, c(0.05, 0.5, 0.95), names = FALSE))
}, numeric(4)))
predictors$sampled_cells <- unname(q[match(predictors$code, rownames(q)), 1])
predictors$p05 <- unname(q[match(predictors$code, rownames(q)), 2])
predictors$median <- unname(q[match(predictors$code, rownames(q)), 3])
predictors$p95 <- unname(q[match(predictors$code, rownames(q)), 4])
write_table(predictors, cfg, "habitat_predictor_table.csv")

# Training sizes and the tuned XGBoost settings of the frozen cull-trained models.
xgb_spec <- function(key, label) {
  m <- readRDS(checked_input(cfg, key))
  y <- m$fit$pre$mold$outcomes[[1]]
  best <- tune::select_best(m$cv_res, metric = "roc_auc")
  tibble::tibble(model = label, training_rows = length(y), positives = sum(y == "1"),
                 trees = 1500L, tree_depth = best$tree_depth, learn_rate = best$learn_rate,
                 loss_reduction = best$loss_reduction, min_n = best$min_n,
                 tuning = "30-point space-filling grid, five-fold spatial-block CV (50 km blocks), ROC AUC")
}
hmsc <- dplyr::bind_rows(lapply(c("enhanced", "clean", "clean_native"), function(k)
  read.csv(file.path(cfg$dirs$tables, paste0("hmsc_", k, "_joint_production_summary.csv")), stringsAsFactors = FALSE)))
spec <- dplyr::bind_rows(
  xgb_spec("legacy_clean_model", "clean_year"),
  xgb_spec("legacy_rg_model", "rg_year"),
  tibble::tibble(model = paste0("hmsc_", sub("_joint$", "", hmsc$component)), training_rows = hmsc$rows,
                 positives = hmsc$target_positives, target_rows = hmsc$target_rows, covariates = hmsc$covariates,
                 samples = hmsc$samples, thin = hmsc$thin, transient = hmsc$transient, chains = hmsc$chains))
write_table(spec, cfg, "habitat_model_specification.csv")
print(as.data.frame(predictors[, c("code", "sampled_cells", "p05", "median", "p95", "gain_clean_year", "gain_rg_year", "gain_manta_within")]))
print(as.data.frame(spec))
