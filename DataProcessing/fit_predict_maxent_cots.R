############################################################
# fit_predict_maxent_cots.R
# Fit a MaxEnt (presence-only) SDM for COTS using maxnet,
# with ENMeval hyperparameter tuning.
# Presence data: Manta Tow midpoints + ReefScan detections
############################################################

# --- 0. Setup and Packages ---
library(dplyr)
library(readr)
library(readxl)
library(sf)
library(terra)
library(maxnet)
library(ENMeval)
library(ggplot2)

# Set up paths
base_dir <- "C:/Users/smatthew/Documents/GitKraken/COTS Fine Scale SDM"
dp_dir   <- file.path(base_dir, "DataProcessing")

predict_stack_full <- file.path(dp_dir, "data/predictors_terra_30m_full_extended.tif")
manta_csv          <- file.path(dp_dir, "data/COTS Program  Manta Tow Data-2026-02-04.csv")
reefscan_dir       <- file.path(base_dir, "data/ReefScan")

out_dir <- file.path(dp_dir, "outputs")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

seed <- 123
n_bg <- 50000  # Background points

# Output file paths
output_tif       <- file.path(out_dir, "COTS_maxent_suitability.tif")
output_model     <- file.path(out_dir, "maxent_model_cots.rds")
output_vip       <- file.path(out_dir, "maxent_vip_cots.png")
output_response  <- file.path(out_dir, "maxent_response_curves.png")
output_enmeval   <- file.path(out_dir, "maxent_enmeval_results.rds")

# --- 1. Load Predictors (Clean 12-layer stack for full domain coverage into Torres Strait) ---
print("Loading predictor stack...")
predictors_all <- terra::rast(predict_stack_full)

# Standardize layer names
names(predictors_all)[1:12] <- c("BPI", "EAST", "HCU", "NORTH", "SLO", "SVF",
                                  "VCU", "VRM", "DEM", "AhyaD", "Aspat", "Aten")

# Use only the 12 clean predictors (excluding ReefGuide layers) to ensure 100% coverage into Torres Strait
predictors <- predictors_all[[1:12]]

pred_names <- names(predictors)
print(paste("Clean predictor layers (12 layers):", paste(pred_names, collapse = ", ")))

# --- 2. Assemble Presence Points ---
print("Assembling presence points...")

# --- 2a. Manta Tow Presences ---
print("  Loading manta tow data...")
manta_raw <- read_csv(manta_csv, show_col_types = FALSE)

manta_pres <- manta_raw %>%
  filter(
    !is.na(StartLongitude), !is.na(EndLongitude),
    !is.na(StartLatitude),  !is.na(EndLatitude),
    StartLongitude != 0, StartLatitude != 0,
    EndLongitude != 0,   EndLatitude != 0
  ) %>%
  filter(
    CrownOfThornsStarfishCount > 0 |
    FeedingScarCountRangeCode == "c"
  ) %>%
  mutate(
    longitude = (StartLongitude + EndLongitude) / 2,
    latitude  = (StartLatitude + EndLatitude) / 2,
    source    = "manta_tow"
  ) %>%
  select(longitude, latitude, source)

print(paste("  Manta tow presences:", nrow(manta_pres)))

# --- 2b. ReefScan Presences ---
print("  Loading ReefScan monthly CSVs...")
rs_files <- list.files(reefscan_dir, pattern = "\\.csv$", full.names = TRUE)
print(paste("  Found", length(rs_files), "ReefScan CSV files"))

rs_list <- lapply(rs_files, function(f) {
  df <- read_csv(f, show_col_types = FALSE)
  df %>%
    select(latitude, longitude) %>%
    mutate(source = "reefscan")
})

reefscan_pres <- bind_rows(rs_list)
print(paste("  ReefScan presences:", nrow(reefscan_pres)))

# --- 2c. Combine and De-duplicate ---
print("  Combining and de-duplicating...")
all_pres <- bind_rows(manta_pres, reefscan_pres)
print(paste("  Total raw presences:", nrow(all_pres)))

# Convert to sf and project to raster CRS
pres_sf <- st_as_sf(all_pres, coords = c("longitude", "latitude"), crs = 4326) %>%
  st_transform(crs(predictors))

# De-duplicate to unique raster cells
pres_vect <- terra::vect(pres_sf)
cell_ids  <- terra::cellFromXY(predictors, terra::crds(pres_vect))

# Keep one point per unique cell
unique_mask <- !duplicated(cell_ids) & !is.na(cell_ids)
pres_vect   <- pres_vect[unique_mask]
pres_cells  <- cell_ids[unique_mask]

print(paste("  Unique cell-level presences:", length(pres_cells)))

# Get presence coordinates as a data.frame for ENMeval
pres_coords <- as.data.frame(terra::crds(pres_vect))
names(pres_coords) <- c("x", "y")

# --- 3. Generate Background Points ---
print(paste("Generating", n_bg, "background points..."))
set.seed(seed)

# Sample from non-NA cells (spatSample), avoiding presence cells
bg_vect <- terra::spatSample(predictors[[1]], size = n_bg * 2,
                              method = "random", na.rm = TRUE,
                              as.points = TRUE)
bg_cells <- terra::cellFromXY(predictors, terra::crds(bg_vect))

# Remove any that overlap with presence cells
bg_keep <- !(bg_cells %in% pres_cells)
bg_vect <- bg_vect[bg_keep]

# Trim to exactly n_bg
if (nrow(bg_vect) > n_bg) {
  bg_vect <- bg_vect[1:n_bg]
}

bg_coords <- as.data.frame(terra::crds(bg_vect))
names(bg_coords) <- c("x", "y")
print(paste("  Background points generated:", nrow(bg_coords)))

# --- 4. Extract Environmental Values ---
print("Extracting environmental values...")

pres_env <- terra::extract(predictors, pres_coords, ID = FALSE)
bg_env   <- terra::extract(predictors, bg_coords,   ID = FALSE)

# Remove rows with NAs
pres_complete <- complete.cases(pres_env)
bg_complete   <- complete.cases(bg_env)

pres_coords <- pres_coords[pres_complete, ]
pres_env    <- pres_env[pres_complete, ]

bg_coords <- bg_coords[bg_complete, ]
bg_env    <- bg_env[bg_complete, ]

print(paste("  Presences after NA removal:", nrow(pres_env)))
print(paste("  Background after NA removal:", nrow(bg_env)))

# Combine coords and env for ENMeval tabular mode
occs_df <- cbind(pres_coords, pres_env)
bg_df   <- cbind(bg_coords, bg_env)

# --- 5. ENMeval Hyperparameter Tuning ---
print("Running ENMeval tuning on tabular environmental data...")
print("  Feature classes: L, LQ, LQH, H")
print("  Regularization multipliers: 0.5, 1, 2, 3, 5")

set.seed(seed)
enmeval_res <- ENMevaluate(
  occs       = occs_df,
  bg         = bg_df,
  envs       = NULL,
  algorithm  = "maxnet",
  partitions = "block",
  tune.args  = list(
    fc = c("L", "LQ", "LQH", "H"),
    rm = c(0.5, 1, 2, 3, 5)
  ),
  doClamp = TRUE
)

# Print results summary
print("ENMeval Results:")
print(enmeval_res@results %>%
        select(fc, rm, auc.train, auc.val.avg, auc.val.sd,
               auc.diff.avg, or.mtp.avg, or.10p.avg, AICc) %>%
        arrange(AICc))

# Save full ENMeval results
saveRDS(enmeval_res, output_enmeval)
print(paste("Saved ENMeval results to:", output_enmeval))

# Select best model by AICc (or delta.AICc == 0)
best_row <- enmeval_res@results %>%
  filter(!is.na(AICc)) %>%
  slice_min(AICc, n = 1)

print(paste("Best model: fc =", best_row$fc, ", rm =", best_row$rm,
            ", AICc =", round(best_row$AICc, 2),
            ", AUC.val =", round(best_row$auc.val.avg, 3)))

# Extract the best model object
best_idx   <- which(enmeval_res@results$fc == best_row$fc &
                    enmeval_res@results$rm == best_row$rm)
best_model <- enmeval_res@models[[best_idx]]

# Save the best model
saveRDS(list(model = best_model,
             best_params = best_row,
             pred_names = pred_names),
        output_model)
print(paste("Saved best MaxEnt model to:", output_model))

# --- 6. Variable Importance & Response Curves ---
print("Generating variable importance plot...")

# Extract variable contributions from the maxnet model coefficients
# Group coefficient magnitudes by variable
coef_vals <- best_model$betas
if (is.null(coef_vals)) coef_vals <- coef(best_model)

# Parse variable names from coefficient names (maxnet uses e.g. "BPI", "BPI^2", "hinge(BPI)")
var_contrib <- data.frame(
  coef_name = names(coef_vals),
  abs_coef  = abs(as.numeric(coef_vals)),
  stringsAsFactors = FALSE
)

# Map each coefficient to its base variable
var_contrib$base_var <- sapply(var_contrib$coef_name, function(cn) {
  # Remove hinge(), threshold(), categorical() wrappers and exponents
  cn <- gsub("^(hinge|thresholds|categorical)\\(", "", cn)
  cn <- gsub("\\)$", "", cn)
  cn <- gsub("\\^[0-9]+$", "", cn)
  # For interaction terms (var1:var2), take the first variable
  cn <- strsplit(cn, ":")[[1]][1]
  cn
})

vip_df <- var_contrib %>%
  group_by(base_var) %>%
  summarise(importance = sum(abs_coef), .groups = "drop") %>%
  arrange(desc(importance)) %>%
  mutate(base_var = factor(base_var, levels = rev(base_var)))

vip_plot <- ggplot(vip_df, aes(x = base_var, y = importance)) +
  geom_col(fill = "#2c7fb8") +
  coord_flip() +
  theme_minimal(base_size = 14) +
  labs(
    title = "MaxEnt Variable Importance (Summed |Coefficients|)",
    x = NULL, y = "Summed Absolute Coefficient"
  )

ggsave(output_vip, plot = vip_plot, width = 8, height = 6, dpi = 300)
print(paste("Saved VIP plot to:", output_vip))

# Response curves
print("Generating response curves...")
png(output_response, width = 1200, height = 900, res = 150)
par(mfrow = c(4, 5), mar = c(3, 3, 2, 1))
plot(best_model, type = "cloglog")
dev.off()
print(paste("Saved response curves to:", output_response))

# --- 7. Predict to Spatial Grid ---
print("Predicting MaxEnt suitability to spatial grid...")
print("  (This may take a while for a GBR-wide raster)")

# Predict using terra::predict with the maxnet model
# maxnet::predict returns cloglog by default
pred_fun <- function(model, newdata) {
  newdata <- as.data.frame(newdata)
  names(newdata) <- pred_names
  p <- predict(model, newdata, type = "cloglog")
  as.numeric(p)
}

suitability <- terra::predict(
  predictors,
  best_model,
  fun   = pred_fun,
  na.rm = TRUE,
  cores = 1
)

names(suitability) <- "suitability"

output_tif_clean <- file.path(out_dir, "COTS_maxent_suitability_clean.tif")

if (file.exists(output_tif)) file.remove(output_tif)
writeRaster(
  suitability,
  filename  = output_tif,
  overwrite = TRUE,
  gdal      = c("COMPRESS=LZW", "TILED=YES", "BIGTIFF=YES")
)

if (file.exists(output_tif_clean)) file.remove(output_tif_clean)
writeRaster(
  suitability,
  filename  = output_tif_clean,
  overwrite = TRUE,
  gdal      = c("COMPRESS=LZW", "TILED=YES", "BIGTIFF=YES")
)

# Summary stats
suit_vals <- terra::global(suitability, fun = c("mean", "min", "max", "sd"), na.rm = TRUE)
n_cells   <- terra::global(suitability, fun = "notNA")
cat("\n=== MaxEnt Prediction Summary ===\n")
cat("Non-NA cells:", n_cells[[1]], "\n")
cat("Mean suitability:", round(suit_vals$mean, 4), "\n")
cat("Min suitability:", round(suit_vals$min, 4), "\n")
cat("Max suitability:", round(suit_vals$max, 4), "\n")
cat("SD suitability:", round(suit_vals$sd, 4), "\n")

print("Done! MaxEnt model fitting and prediction complete.")
