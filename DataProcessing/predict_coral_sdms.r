#######################################################
# predict_coral_sdms.R
# Standalone prediction script for existing MaxEnt Models
#######################################################
library(sf)
library(terra)
library(dplyr)
library(stringr)
library(maxnet)
library(enmSdmX)

# --- Configuration ---
# Provide paths to original models
model_dir <- "C:/Users/smatthew/Documents/GitKraken/COTS Fine Scale SDM/DataProcessing/data/SDMModel_Guillame"
models <- c(
    "AhyaD" = file.path(model_dir, "maxent_model_FINAL_maxnet_ext-targetreefs_ALLallres_AhyaD_lq2.RData"),
    "Aspat" = file.path(model_dir, "maxent_model_FINAL_maxnet_ext-targetreefs_ALLallres_Aspat_lq2.RData"),
    "Aten"  = file.path(model_dir, "maxent_model_FINAL_maxnet_ext-targetreefs_ALLallres_Aten_lq2.RData")
)

# New data paths
# Base 15m Data
dem_path_15 <- "C:/Users/smatthew/Documents/GitKraken/COTS Fine Scale SDM/DataProcessing/data/SDMModel_Guillame/processed_dem/TS_bathymetryDEM_15_UTM55_030UTM55_clean_030.tif"
var_dir_15 <- "C:/Users/smatthew/Documents/GitKraken/COTS Fine Scale SDM/DataProcessing/data/SDMModel_Guillame/variables/TS_bathymetryDEM_15_UTM55_030030"

# Base 30m Data (DeepReef30)
dem_path_30 <- "C:/Users/smatthew/Documents/GitKraken/COTS Fine Scale SDM/DataProcessing/data/SDMModel_Guillame/processed_dem/TS_bathymetryDEM_DeepReef30_UTM55_030.tif"
var_dir_30 <- "C:/Users/smatthew/Documents/GitKraken/COTS Fine Scale SDM/DataProcessing/data/SDMModel_Guillame/variables/TS_bathymetryDEM_DeepReef30_UTM55_030030"

shape_path <- "C:/Users/smatthew/Documents/GitKraken/COTS Fine Scale SDM/DataProcessing/data/GBR_TS_Reef_Features.gpkg"

output_dir <- file.path(model_dir, "maxent_prediction_raster_new")
dir.create(output_dir, showWarnings = FALSE)

# --- 1. Load Extent Shapefile ---
print("Loading reef shapefile constraints...")
shape_ext <- st_read(shape_path)

# --- 2. Create Empty Reference Raster ---
print("Loading reference DEM and creating reference raster (15m base)...")
ref_dem <- terra::rast(dem_path_15)

# Ensure matches raster CRS
shape_ext <- sf::st_transform(shape_ext, crs = terra::crs(ref_dem))

# Create a 1km buffered version used for tile extents (ensures predictions
# extend slightly beyond reef edges to avoid clipping artefacts)
shape_ext_buf <- sf::st_buffer(shape_ext, dist = 1000) # 1000m = 1 km (UTM units)

# Crop DEM to the BUFFERED shapefile extent to save time/memory
print("Cropping and masking DEM to buffered reef shapefile (1km buffer)...")
ref_dem_cropped <- terra::crop(ref_dem, shape_ext_buf)
ref_dem_masked <- terra::mask(ref_dem_cropped, shape_ext_buf)

# Create an empty template from the masked reference
empty_rast <- ref_dem_masked / ref_dem_masked # 1s where data exists, NA elsewhere
names(empty_rast) <- "reference"

# --- Visual Check ---
print("Plotting the extent of the new data and the reefs shapefile to check alignment (this might take a moment to render)...")
# Plot the full extent of the reference raster (as just a bounding box for speed)
plot(terra::ext(ref_dem), main = "Raster Extent vs Reef Shapefile", col = "lightgrey", border = NA, axes = TRUE)
# Overlay the actual new data extent (masked)
plot(empty_rast, add = TRUE, col = "lightblue", legend = FALSE)
# Overlay the shapefile
plot(sf::st_geometry(shape_ext), border = "red", add = TRUE)
print("Please review the Plot window to ensure the extents match before predictions begin.")
Sys.sleep(3) # Give the user a moment to see it flash up


# Enable Tiling for large predictions
tile_dir <- file.path(output_dir, "tiles")
dir.create(tile_dir, showWarnings = FALSE)
print("Creating tiles...")

# Calculate reasonably sized chunks (e.g., 20x20 grid of tiles)
n_tiles_x <- 20
n_tiles_y <- 20
tile_size_x <- max(1, floor(ncol(empty_rast) / n_tiles_x))
tile_size_y <- max(1, floor(nrow(empty_rast) / n_tiles_y))

# Save the base empty raster so makeTiles can process it
empty_ref_file <- file.path(tile_dir, "empty_reference.tif")
terra::writeRaster(empty_rast, empty_ref_file, overwrite = TRUE)
empty_rast_file_read <- terra::rast(empty_ref_file)

terra::makeTiles(
    x = empty_rast_file_read,
    y = c(tile_size_y, tile_size_x),
    filename = file.path(tile_dir, "empty_tile_.tif"),
    extend = TRUE,
    overwrite = TRUE
)

tiles <- list.files(tile_dir, pattern = "empty_tile_.*\\.tif$", full.names = TRUE)
print(paste("Created", length(tiles), "tiles."))

# --- 3. Loop over models to Predict ---
for (sp in names(models)) {
    # for (sp in c("Aspat", "Aten")) {
    print(paste("-----------------------------------------"))
    print(paste("Working on species:", sp))

    # Load the MaxEnt model -> loads as 'me'
    me <- NULL
    if (file.exists(models[[sp]])) {
        load(models[[sp]])
    } else {
        warning(paste("Model file not found:", models[[sp]]))
        next
    }

    if (is.null(me)) {
        warning(paste("Model 'me' failed to load from", models[[sp]]))
        next
    }

    # Identify variables expected by the model
    var.list <- names(me$levels)
    print(paste("Model expects the following variables:", paste(var.list, collapse = ", ")))

    # Process each tile
    for (t_idx in seq_along(tiles)) {
        # print(paste("  Processing tile", t_idx, "of", length(tiles)))
        empty.rast.T <- terra::rast(tiles[t_idx])

        # Check if the tile actually contains valid cells (is not fully NA)
        # Using terra::global to get the max value safely
        max_val <- suppressWarnings(terra::global(empty.rast.T, "max", na.rm = TRUE)[1, 1])

        # If the max isn't finite (i.e., NO valid cells), skip
        if (!is.finite(max_val)) {
            next
        }
        print(paste("    -> Valid cells found in tile", t_idx, "of", length(tiles)))

        # Load and map variables for this tile
        rastpred <- list()
        missing_vars <- FALSE

        # Cache memory for 15m and 30m base variables for this specific tile to prevent 72x disk read penalty
        loaded_base_15 <- list()
        loaded_base_30 <- list()

        for (VAR in var.list) {
            # Parse the model variable (e.g. ACA15_120_DEM or gbr030_060_EAST or gbr100_100_DEM)
            var_parts <- strsplit(VAR, "_")[[1]]
            dataset_prefix <- var_parts[1] # "ACA15", "gbr030", "gbr100"
            target_res_str <- var_parts[2] # "015", "030", "060", "120", "100"
            var_type <- var_parts[3] # "DEM", "BPI", "EAST", etc

            target_res <- as.numeric(target_res_str)

            # Determine which physical dataset to pull from
            is_aca <- grepl("^ACA", dataset_prefix)

            if (is_aca) {
                base_res <- 15
                base_dem <- dem_path_15
                base_dir <- var_dir_15
                file_pref <- "TS_bathymetryDEM_15_UTM55_030030_"
                base_cache <- loaded_base_15
            } else {
                # gbr030 or gbr100 use the DeepReef30 source
                base_res <- 30
                base_dem <- dem_path_30
                base_dir <- var_dir_30
                file_pref <- "TS_bathymetryDEM_DeepReef30_UTM55_030030_"
                base_cache <- loaded_base_30
            }

            # Resolve file path
            if (var_type == "DEM") {
                temp.fn <- base_dem
            } else {
                temp.fn <- file.path(base_dir, paste0(file_pref, var_type, ".sdat"))
            }

            if (!file.exists(temp.fn)) {
                stop(paste("File missing for required variable type:", VAR, "->", temp.fn))
                missing_vars <- TRUE
                break
            }

            # Load and crop base layer if not cached down to the tile yet
            if (is.null(base_cache[[var_type]])) {
                temp.r <- terra::rast(temp.fn)
                # Crop identically from geographic map to our local tile boundary
                temp.r <- terra::crop(temp.r, empty.rast.T)
                base_cache[[var_type]] <- temp.r
            }

            var_layer <- base_cache[[var_type]]

            # MULTISCALE AGGREGATION
            # If the model structurally demands a wider focal window resolution (e.g. 120m)
            # than what the base data provides (e.g. 15m), we mathematically aggregate (blur) it up in RAM.
            # Just print the variable we are working on to see the target:
            print(paste("Processing:", VAR, "| Base Res:", base_res, "| Target Res:", target_res))
            if (target_res > base_res) {
                agg_factor <- floor(target_res / base_res)
                if (agg_factor >= 2) {
                    var_layer <- terra::aggregate(var_layer, fact = agg_factor, fun = "mean", na.rm = TRUE)
                    print(paste("   -> Aggregated from", base_res, "to", res(var_layer)[1]))
                }
            }

            # STACK ALIGNMENT
            # predictMaxNet math demands all 72 outputs must share the IDENTICAL spatial coordinate
            # bounds and cellular resolutions to align matrices. We project our (aggregated or raw)
            # layer flawlessly back onto our master 15m target tile grid.
            var_layer <- terra::resample(var_layer, empty.rast.T, method = "bilinear", threads = TRUE)
            print(paste("   -> Resampled back to tile grid:", res(var_layer)[1]))


            # Store cache back to top-level lists
            if (is_aca) {
                loaded_base_15 <- base_cache
            } else {
                loaded_base_30 <- base_cache
            }

            # OVERRIDE the raster name to match the model's expected equation parameters EXACTLY
            names(var_layer) <- VAR
            rastpred <- c(rastpred, var_layer)
        }

        if (missing_vars) next

        # Stack the mathematically complete 72 dimension list
        rastpred <- terra::rast(rastpred)

        # Predict MaxNet
        pred <- enmSdmX::predictMaxNet(me, rastpred, clamp = TRUE, type = "cloglog")

        # Save output
        out_tile_name <- file.path(tile_dir, paste0("pred_", sp, "_tile_", t_idx, ".tif"))
        terra::writeRaster(pred, out_tile_name, overwrite = TRUE)

        # Clean up memory
        rm(rastpred, pred)
        gc()
    }

    # Merge all prediction tiles for this species
    print(paste("Merging predicted tiles for", sp, "..."))
    pred_files <- list.files(tile_dir, pattern = paste0("pred_", sp, "_tile_.*\\.tif$"), full.names = TRUE)
    if (length(pred_files) > 0) {
        list.temp.r <- lapply(pred_files, terra::rast)
        rsrc <- terra::sprc(list.temp.r)

        # 1. Save Full Prediction
        final_out <- file.path(output_dir, paste0("predict_", sp, ".tif"))
        merged_pred <- terra::merge(rsrc, filename = final_out, overwrite = TRUE)
        print(paste("Saved final prediction:", final_out))

        # 2. Save Masked Prediction (clipped to exact reef polygons, no buffer)
        print(paste("Generating reef-masked prediction for", sp, "..."))
        masked_pred <- terra::mask(merged_pred, terra::vect(shape_ext_buf))
        masked_out <- file.path(output_dir, paste0("predict_", sp, "_masked.tif"))
        terra::writeRaster(masked_pred, filename = masked_out, overwrite = TRUE)
        print(paste("Saved masked prediction:", masked_out))
    } else {
        print(paste("No prediction tiles generated for", sp, "- maybe all variables were masked?"))
    }
}

# Cleanup empty model tiles
print("Cleaning up temporary tiles...")
unlink(tile_dir, recursive = TRUE)

print("Done. All models predicted.")
