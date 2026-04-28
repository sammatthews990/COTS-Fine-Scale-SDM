############################################################
# build_clean_predictor_stack.R
# Creates a clean predictor .tif with ONLY:
#   - env layers (GBR_bathymetryDEM_*)
#   - new full-extent coral SDMs (AhyaD, Aspat, Aten)
# Excludes SDM_* (old/partial) and RG_* (restricted extent)
############################################################
library(terra)

input_stack <- "C:/Users/smatthew/Documents/GitKraken/COTS Fine Scale SDM/DataProcessing/data/predictors_terra_30m_combined_masked.tif"
output_stack <- "C:/Users/smatthew/Documents/GitKraken/COTS Fine Scale SDM/DataProcessing/data/predictors_clean_12layer.tif"

cat("Loading full stack...\n")
r <- terra::rast(input_stack)
nms <- names(r)
cat("Input layers:", nlyr(r), "\n")
cat(paste(nms, collapse = "\n"), "\n\n")

# Select layers to KEEP
keep_env   <- grep("^GBR_", nms, value = TRUE)
keep_coral <- c("AhyaD", "Aspat", "Aten")
keep <- c(keep_env, keep_coral)

# Verify all expected layers exist
missing <- setdiff(keep, nms)
if (length(missing) > 0) stop("Missing expected layers: ", paste(missing, collapse = ", "))

cat("KEEPING", length(keep), "layers:\n")
cat(paste(keep, collapse = "\n"), "\n\n")

dropped <- setdiff(nms, keep)
cat("DROPPING", length(dropped), "layers:\n")
cat(paste(dropped, collapse = "\n"), "\n\n")

# Subset
r_clean <- r[[keep]]

# Rename env layers to short names for cleaner model output
short_names <- c(
  "GBR_bathymetryDEM_ACA15_UTM55_030_BPI"   = "BPI",
  "GBR_bathymetryDEM_ACA15_UTM55_030_EAST"  = "EAST",
  "GBR_bathymetryDEM_ACA15_UTM55_030_HCU"   = "HCU",
  "GBR_bathymetryDEM_ACA15_UTM55_030_NORTH" = "NORTH",
  "GBR_bathymetryDEM_ACA15_UTM55_030_SLO"   = "SLO",
  "GBR_bathymetryDEM_ACA15_UTM55_030_SVF"   = "SVF",
  "GBR_bathymetryDEM_ACA15_UTM55_030_VCU"   = "VCU",
  "GBR_bathymetryDEM_ACA15_UTM55_030_VRM"   = "VRM",
  "GBR_bathymetryDEM_ACA_UTM55_030_ext"     = "DEM"
)

new_names <- names(r_clean)
for (i in seq_along(new_names)) {
  if (new_names[i] %in% names(short_names)) {
    new_names[i] <- short_names[new_names[i]]
  }
}
names(r_clean) <- new_names

cat("Final layer names:\n")
cat(paste(names(r_clean), collapse = "\n"), "\n\n")

# Write clean stack
cat("Writing clean stack to:", output_stack, "\n")
cat("(This may take several minutes for a large raster...)\n")
terra::writeRaster(
  r_clean,
  filename  = output_stack,
  overwrite = TRUE,
  gdal      = c("COMPRESS=LZW", "TILED=YES", "BIGTIFF=YES")
)

cat("Done! Clean 12-layer stack saved.\n")
cat("File:", output_stack, "\n")
