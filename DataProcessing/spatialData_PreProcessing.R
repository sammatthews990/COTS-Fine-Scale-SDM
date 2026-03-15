#######################################################
# spatialData_PreProcessing.R
# Pre-processes spatial predictors and SDM predictions
# by merging old and new spatial datasets and aligning
# them to a unified 30m template grid.
#######################################################
library(terra)
library(dplyr)
library(stringr)

library(sf)

# --- Configuration Paths ---
predict_old_dir <- "C:/Users/smatthew/Documents/GitKraken/COTS Fine Scale SDM/DataProcessing/data"
predict_new_dir <- "C:/Users/smatthew/Documents/GitKraken/COTS Fine Scale SDM/DataProcessing/data/SDMModel_Guillame/maxent_prediction_raster_new"

vars_old_dir <- "C:/Users/smatthew/Documents/GitKraken/COTS Fine Scale SDM/DataProcessing/data/GBR_bathymetryDEM_ACA_UTM55_030_crop"

out_layers_dir <- "C:/Users/smatthew/Documents/GitKraken/COTS Fine Scale SDM/DataProcessing/data/processed_layers"
if (!dir.exists(out_layers_dir)) dir.create(out_layers_dir, recursive = TRUE)

# Base 15m Data
dem_path_15 <- "C:/Users/smatthew/Documents/GitKraken/COTS Fine Scale SDM/DataProcessing/data/SDMModel_Guillame/processed_dem/TS_bathymetryDEM_15_UTM55_030UTM55_clean_030.tif"
var_dir_15 <- "C:/Users/smatthew/Documents/GitKraken/COTS Fine Scale SDM/DataProcessing/data/SDMModel_Guillame/variables/TS_bathymetryDEM_15_UTM55_030030"

# Base 30m Data (DeepReef30)
dem_path_30 <- "C:/Users/smatthew/Documents/GitKraken/COTS Fine Scale SDM/DataProcessing/data/SDMModel_Guillame/processed_dem/TS_bathymetryDEM_DeepReef30_UTM55_030.tif"
var_dir_30 <- "C:/Users/smatthew/Documents/GitKraken/COTS Fine Scale SDM/DataProcessing/data/SDMModel_Guillame/variables/TS_bathymetryDEM_DeepReef30_UTM55_030030"

shape_path <- "C:/Users/smatthew/Documents/GitKraken/COTS Fine Scale SDM/DataProcessing/data/GBR_TS_Reef_Features.gpkg"

output_stack_path <- "C:/Users/smatthew/Documents/GitKraken/COTS Fine Scale SDM/DataProcessing/data/predictors_terra_30m_combined.tif"
output_stack_masked_path <- "C:/Users/smatthew/Documents/GitKraken/COTS Fine Scale SDM/DataProcessing/data/predictors_terra_30m_combined_masked.tif"

# Define the species codes used in your files
species_codes <- c("AhyaD", "Aspat", "Aten")

# --- Load Reef Shapefile and compute 1km Buffer ---
print("Loading reef shapefile and applying 1km buffer...")
shape_ext <- st_read(shape_path)
shape_ext_buf <- sf::st_buffer(shape_ext, dist = 1000)

# --- 1. Define 30m Master Template From Full Extents ---
print("Loading old full GBR predictor stack...")
old_stack_path <- "C:/Users/smatthew/Documents/GitKraken/COTS Fine Scale SDM/DataProcessing/data/predictors_terra_30m.tif"
if (!file.exists(old_stack_path)) stop(paste("Old predictor stack missing:", old_stack_path))
env30_old <- terra::rast(old_stack_path)

print("Building 30m Master Template Extent...")
# Load the base 15m and 30m TS DEMs to get the maximum northern spatial extent
dem15 <- terra::rast(dem_path_15)
dem30 <- terra::rast(dem_path_30)

# Extend the old dataset extent to include the northern regions (Torres Strait)
template30 <- terra::extend(env30_old[[1]], ext(dem15))
template30 <- terra::extend(template30, ext(dem30))
template30 <- terra::init(template30, NA)
names(template30) <- "template"

print(paste("Template Extent:", paste(as.vector(ext(template30)), collapse = ", ")))
print(paste("Template Resolution:", paste(res(template30), collapse = "x")))


# --- 2. Merge & Align SDM Predictions ---
print("Processing model predictions...")
sdm_layers <- list()

for (sp in species_codes) {
    print(paste(" -> Harmonizing predictions for", sp))
    layer_name <- paste0("SDM_", sp)
    out_file <- file.path(out_layers_dir, paste0(layer_name, ".tif"))

    if (file.exists(out_file)) {
        print(paste("      Loading existing processed layer:", layer_name))
        sdm_layers[[sp]] <- terra::rast(out_file)
        next
    }

    # Load Old Prediction
    old_pred_file <- file.path(predict_old_dir, paste0("maxent_predrast_GBR_", sp, "_015_lq2.tif"))
    if (!file.exists(old_pred_file)) stop(paste("Old prediction missing:", old_pred_file))
    old_pred <- terra::rast(old_pred_file)
    old_pred_aligned <- terra::resample(old_pred, template30, method = "bilinear")

    # Load New Prediction
    new_pred_file <- file.path(predict_new_dir, paste0("predict_", sp, ".tif"))
    new_pred_exists <- file.exists(new_pred_file)

    if (new_pred_exists) {
        new_pred <- terra::rast(new_pred_file)
        print("      Merging old and new predictions...")
        new_pred_aligned <- terra::resample(new_pred, template30, method = "bilinear")

        # 1. Merge the spatial extents (Prioritize the new high-res predictions over the old ones)
        aligned_pred <- terra::merge(new_pred_aligned, old_pred_aligned)
    } else {
        print("      No new predictions found. Using only old prediction.")
        aligned_pred <- old_pred_aligned
    }

    names(aligned_pred) <- layer_name
    print(paste("      Saving processed predictor:", layer_name))
    terra::writeRaster(aligned_pred, filename = out_file, overwrite = TRUE, gdal = c("COMPRESS=LZW"))
    sdm_layers[[sp]] <- aligned_pred
}

coral_sdm_30 <- terra::rast(sdm_layers)


# --- 3. Merge & Align Environmental Predictor Variables ---
print("Processing environmental predictors...")

# Get the layer names from the old stack (don't extend the whole thing at once!)
old_layer_names <- names(env30_old)

# We will create a list of final harmonized environmental layers
final_env_layers <- list()
loaded_base_15 <- list()
loaded_base_30 <- list()

# Map the base core_vars we have new TS data for
core_vars_map <- list(
    "BPI" = "BPI",
    "EAST" = "EAST",
    "HCU" = "HCU",
    "NORTH" = "NORTH",
    "SLO" = "SLO",
    "SVF" = "SVF",
    "VCU" = "VCU",
    "VRM" = "VRM",
    "ext" = "DEM" # The model uses 'ext' for the DEM layer naming
)

for (layer_name in old_layer_names) {
    out_file <- file.path(out_layers_dir, paste0(layer_name, ".tif"))
    if (file.exists(out_file)) {
        print(paste(" -> Loading existing processed layer:", layer_name))
        final_env_layers[[layer_name]] <- terra::rast(out_file)
        next
    }

    # Extract and resample SINGLE layer to full template (avoids bulk extend)
    print(paste(" -> Resampling old layer to full template:", layer_name))
    old_rast <- terra::resample(env30_old[[layer_name]], template30, method = "bilinear")

    # Identify if we have TS data to merge for this layer
    core_var <- NA
    for (key in names(core_vars_map)) {
        if (grepl(paste0("_", key, "$"), layer_name)) {
            core_var <- core_vars_map[[key]]
            break
        }
    }

    # If it is a ReefGuidance or Species layer without new TS data
    if (is.na(core_var)) {
        print(paste(" -> Keeping old data only for:", layer_name))
        aligned_env <- old_rast
    } else {
        print(paste(" -> Merging new Torres Strait data for predictor:", layer_name, "(Core Var:", core_var, ")"))

        new_ts_merged <- NULL

        # 1. Load TS 15m (ACA)
        aca_file <- if (core_var == "DEM") dem_path_15 else file.path(var_dir_15, paste0("TS_bathymetryDEM_15_UTM55_030030_", core_var, ".sdat"))
        if (file.exists(aca_file)) {
            if (is.null(loaded_base_15[[core_var]])) loaded_base_15[[core_var]] <- terra::rast(aca_file)
            r_aca <- loaded_base_15[[core_var]]

            # Aggregate 15m to 30m exactly mathematically
            r_aca_30 <- terra::aggregate(r_aca, fact = 2, fun = "mean", na.rm = TRUE)
            new_ts_merged <- r_aca_30
        }

        # 2. Load TS 30m (DeepReef)
        deepreef_file <- if (core_var == "DEM") dem_path_30 else file.path(var_dir_30, paste0("TS_bathymetryDEM_DeepReef30_UTM55_030030_", core_var, ".sdat"))
        if (file.exists(deepreef_file)) {
            if (is.null(loaded_base_30[[core_var]])) loaded_base_30[[core_var]] <- terra::rast(deepreef_file)
            r_deep <- loaded_base_30[[core_var]]

            if (is.null(new_ts_merged)) {
                new_ts_merged <- r_deep
            } else {
                # Both need to be aligned to mosaic them safely. Resample DeepReef to match ACA's 30m grid before merge
                r_deep_aligned <- terra::resample(r_deep, new_ts_merged, method = "bilinear", threads = TRUE)
                # Merge ACA and DeepReef (prioritize ACA for shallow fidelity)
                new_ts_merged <- terra::merge(new_ts_merged, r_deep_aligned)
            }
        }

        # 3. Align everything perfectly to the final template and merge with the old dataset
        if (!is.null(new_ts_merged)) {
            # Align the merged TS data to the master 30m template
            new_ts_aligned <- terra::resample(new_ts_merged, template30, method = "bilinear", threads = TRUE)

            # Merge TS data with the old GBR data (prioritize new TS data over overlapping old regions)
            aligned_env <- terra::merge(new_ts_aligned, old_rast)
        } else {
            aligned_env <- old_rast
        }
    }

    names(aligned_env) <- layer_name
    print(paste(" -> Saving processed predictor:", layer_name))
    terra::writeRaster(aligned_env, filename = out_file, overwrite = TRUE, gdal = c("COMPRESS=LZW"))
    final_env_layers[[layer_name]] <- aligned_env
}

combined_env_30 <- terra::rast(final_env_layers)


# --- 4. Create and Save Final Master Stack ---
print("Assembling master predictor stack...")
predictors_terra_combined <- c(combined_env_30, coral_sdm_30)

print(predictors_terra_combined)

print(paste("Saving unified stack to:", output_stack_path, "(This may take a while...)"))
terra::writeRaster(
    predictors_terra_combined,
    filename = output_stack_path,
    overwrite = TRUE,
    gdal = c("COMPRESS=LZW") # Smaller file size
)

print(paste("Applying 1km Reef Buffer Mask..."))
# Ensure shapefile CRS matches the master template (should map exactly)
shape_ext_buf <- sf::st_transform(shape_ext_buf, crs = terra::crs(predictors_terra_combined))
predictors_masked <- terra::mask(predictors_terra_combined, terra::vect(shape_ext_buf))

print(paste("Saving masked unified stack to:", output_stack_masked_path))
terra::writeRaster(
    predictors_masked,
    filename = output_stack_masked_path,
    overwrite = TRUE,
    gdal = c("COMPRESS=LZW")
)

print("Preprocessing complete!")
