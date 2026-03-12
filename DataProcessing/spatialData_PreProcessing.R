#######################################################
# spatialData_PreProcessing.R
# Pre-processes spatial predictors and SDM predictions
# by merging old and new spatial datasets and aligning
# them to a unified 30m template grid.
#######################################################
library(terra)
library(dplyr)
library(stringr)

# --- Configuration Paths ---
predict_old_dir <- "C:/Users/smatthew/Documents/GitKraken/COTS Fine Scale SDM/DataProcessing/data"
predict_new_dir <- "C:/Users/smatthew/Documents/GitKraken/COTS Fine Scale SDM/DataProcessing/data/SDMModel_Guillame/maxent_prediction_raster_new"

vars_old_dir <- "C:/Users/smatthew/Documents/GitKraken/COTS Fine Scale SDM/DataProcessing/data/GBR_bathymetryDEM_ACA_UTM55_030_crop"
vars_new_dir <- "C:/Users/smatthew/Documents/GitKraken/COTS Fine Scale SDM/DataProcessing/data/SDMModel_Guillame/variables/TS_bathymetryDEM_15_UTM55_030030"
dem_new_path <- "C:/Users/smatthew/Documents/GitKraken/COTS Fine Scale SDM/DataProcessing/data/SDMModel_Guillame/processed_dem/TS_bathymetryDEM_15_UTM55_030UTM55_clean_030.tif"

output_stack_path <- "C:/Users/smatthew/Documents/GitKraken/COTS Fine Scale SDM/DataProcessing/data/predictors_terra_30m_combined.tif"

# Define the species codes used in your files
species_codes <- c("AhyaD", "Aspat", "Aten")


# --- 1. Define 30m Master Template ---
print("Loading 30m Master Template...")
# Load all old 30m environmental rasters to form the base environment stack
old_env_files <- list.files(vars_old_dir, pattern = "\\.tif$", full.names = TRUE)
if (length(old_env_files) == 0) stop("No old environmental predictors found!")
env30_old <- terra::rast(old_env_files)

# Use the first layer as the absolute template for extent and resolution
template30 <- env30_old[[1]]
print(paste("Template Resolution:", paste(res(template30), collapse = "x")))


# --- 2. Merge & Align SDM Predictions ---
print("Processing model predictions...")
sdm_layers <- list()

for (sp in species_codes) {
    print(paste(" -> Harmonizing predictions for", sp))

    # Load Old Prediction
    old_pred_file <- file.path(predict_old_dir, paste0("maxent_predrast_GBR_", sp, "_015_lq2.tif"))
    if (!file.exists(old_pred_file)) stop(paste("Old prediction missing:", old_pred_file))
    old_pred <- terra::rast(old_pred_file)

    # Load New Prediction
    new_pred_file <- file.path(predict_new_dir, paste0("predict_", sp, ".tif"))
    new_pred_exists <- file.exists(new_pred_file)

    if (new_pred_exists) {
        new_pred <- terra::rast(new_pred_file)
        print("      Merging old and new predictions...")

        # 1. Merge the spatial extents (this handles non-overlapping and overlapping regions seamlessly)
        merged_pred <- terra::merge(old_pred, new_pred)
    } else {
        print("      No new predictions found. Using only old prediction.")
        merged_pred <- old_pred
    }

    # 2. Resample the merged file to exactly match the 30m template grid
    print("      Resampling to 30m template...")
    aligned_pred <- terra::resample(merged_pred, template30, method = "bilinear")
    names(aligned_pred) <- paste0("SDM_", sp)

    sdm_layers[[sp]] <- aligned_pred
}

coral_sdm_30 <- terra::rast(sdm_layers)


# --- 3. Merge & Align Environmental Predictor Variables ---
print("Processing environmental predictors...")

# The old variables are already loaded in `env30_old`. We need to align them perfectly.
env30_old_aligned <- terra::resample(env30_old, template30, method = "bilinear")

# Load the new variables from .sdat files
new_var_files <- list.files(vars_new_dir, pattern = "\\.sdat$", full.names = TRUE)
print(paste("Found", length(new_var_files), "new variable files to merge."))

# In the previous script, variables like EAST, BPI, etc. were matched. We will iterate and merge overlapping types.
# First, extract the core names from the old variables to know what to look for
old_var_names <- names(env30_old)

# We will create a list of final harmonized environmental layers
final_env_layers <- list()

for (old_name in old_var_names) {
    # Extract the core variable type (e.g. "BPI", "EAST") from the old name (e.g. "ACA15_030_BPI")
    core_var <- stringr::word(old_name, 3, sep = "_")

    # Check if this variable exists in the NEW dataset
    if (core_var == "DEM") {
        new_file <- dem_new_path
    } else {
        new_file <- file.path(vars_new_dir, paste0("TS_bathymetryDEM_15_UTM55_015030_", core_var, ".sdat"))
    }

    old_rast <- env30_old_aligned[[old_name]]

    if (file.exists(new_file)) {
        print(paste(" -> Merging new data for predictor:", core_var))
        new_rast <- terra::rast(new_file)

        # Merge old and new extents together
        merged_env <- terra::merge(old_rast, new_rast)

        # Resample back to unified template
        aligned_env <- terra::resample(merged_env, template30, method = "bilinear")
    } else {
        print(paste(" -> No new data found for predictor:", core_var, "- keeping old only."))
        aligned_env <- old_rast # Already aligned
    }

    names(aligned_env) <- old_name
    final_env_layers[[old_name]] <- aligned_env
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

print("Preprocessing complete!")
