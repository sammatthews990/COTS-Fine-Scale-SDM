library(sf)
library(terra)

base_tif_path <- "DataProcessing/outputs/COTS_prob_0.02_cpue_year2025_clean.tif"
r_base <- terra::rast(base_tif_path)

# Target WGS84 expanded extent: Lat -15.75 to -14.50, Lon 144.90 to 146.10
ext_wgs84 <- ext(144.90, 146.10, -15.75, -14.50)

# Project extent to raster CRS
ext_proj <- terra::project(ext_wgs84, from = "EPSG:4326", to = crs(r_base))

# Crop in native CRS
crop_base_native <- terra::crop(r_base, ext_proj)
crop_base_wgs <- terra::project(crop_base_native, "EPSG:4326")

vals <- values(crop_base_wgs, mat = FALSE)
vals <- vals[!is.na(vals) & !is.nan(vals)]

cat("Summary of COTS probability values in expanded area:\n")
print(summary(vals))

cat("\nQuantiles:\n")
print(quantile(vals, probs = c(0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.85, 0.9, 0.95, 0.99)))
