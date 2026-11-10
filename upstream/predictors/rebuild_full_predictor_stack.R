############################################################
# rebuild_full_predictor_stack.R
# Remosaics ReefGuide layers onto the full-extent template
# and builds a unified 17-layer predictor stack.
# CORRECTED: Uses 1km buffer for masking and GBR_TS reefs.
############################################################
library(terra)
library(sf)
library(dplyr)

# --- 1. Paths ---
clean_stack_path <- "data/predictors_clean_12layer.tif"
rg_dir <- "data/GBR-ReefGuidance-processed_2025-11-19"
reef_poly_path <- "data/GBR_TS_Reef_Features.gpkg" # Correct TS-inclusive features
output_stack_path <- "data/predictors_terra_30m_full_extended.tif"

# Variables we want from ReefGuide
rg_vars <- c("bathy", "slope", "turbid", "waves_Hs", "waves_Tp")

# --- 2. Load Template and Base Layers ---
print("Loading full-extent template from clean stack...")
base_stack <- terra::rast(clean_stack_path)
template30 <- base_stack[[1]] # Use first layer as extent/res template
template30 <- terra::init(template30, NA)

print(paste("Template Extent:", paste(as.vector(ext(template30)), collapse = ", ")))

# --- 3. Remosaic ReefGuide Layers ---
rg_layers_list <- list()

for (v in rg_vars) {
  print(paste("Processing ReefGuide variable:", v))

  files_v <- list.files(
    rg_dir,
    pattern     = paste0(v, "\\.tif$"),
    full.names  = TRUE,
    ignore.case = TRUE
  )

  if (length(files_v) == 0) {
    warning(paste("No files found for ReefGuide variable:", v))
    next
  }

  # Project/Align each regional tile to the full 30m template
  proj_tiles <- list()
  for (i in seq_along(files_v)) {
    f <- files_v[i]
    r <- terra::rast(f)
    
    # Project tile to template (handles different extents/CRS/alignment)
    proj_tiles[[i]] <- terra::project(
      r,
      y      = template30,
      method = "bilinear"
    )
  }

  # Merge projected tiles into one full layer
  r_full <- do.call(terra::merge, proj_tiles)
  names(r_full) <- paste0("RG_", v)
  
  rg_layers_list[[v]] <- r_full
  print(paste("  Done remosaicking RG_", v))
}

# Combine all new RG layers
reefguide_extended <- terra::rast(rg_layers_list)

# --- 4. Assemble Final 17-layer Stack ---
print("Assembling final extended stack (12 base + 5 RG layers)...")
final_stack <- c(base_stack, reefguide_extended)

print("Final Stack Details:")
print(final_stack)

# --- 5. Mask to Reef Features (with 1km Buffer) ---
print("Loading reef polygons and applying 1km buffer...")
# Using sf for the buffer as it is generally more robust for this
reef_sf <- sf::st_read(reef_poly_path)
print("  Buffering 1000m...")
reef_sf_buf <- sf::st_buffer(reef_sf, dist = 1000)

# Convert to terra SpatVector and project to match raster CRS
print("  Projecting mask to raster CRS...")
reef_vect_buf <- terra::vect(reef_sf_buf)
reef_vect_buf <- terra::project(reef_vect_buf, terra::crs(final_stack))

# Mask the stack (removes values outside reef buffer)
print("Applying mask...")
final_stack_masked <- terra::mask(final_stack, reef_vect_buf)

# --- 6. Save ---
print(paste("Saving extended stack to:", output_stack_path))
terra::writeRaster(
  final_stack_masked,
  filename  = output_stack_path,
  overwrite = TRUE,
  gdal      = c("COMPRESS=LZW")
)

print("Process Complete!")
