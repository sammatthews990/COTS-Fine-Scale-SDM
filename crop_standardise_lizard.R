library(terra)

tif_path <- "DataProcessing/outputs/COTS_prob_0.02_cpue_year2025_clean.tif"
gpkg_path <- "DataProcessing/data/lizard_sites_june.gpkg"
out_path <- "DataProcessing/outputs/COTS_prob_0.02_cpue_year2025_Lizard_Island.tif"

cat("Loading data...\n")
r <- rast(tif_path)
v <- vect(gpkg_path)

cat("Cropping and masking...\n")
if(crs(r) != crs(v)) {
  v <- project(v, crs(r))
}

r_crop <- crop(r, v)
r_mask <- mask(r_crop, v)

cat("Standardising values (0-1)...\n")
min_val <- as.numeric(global(r_mask, "min", na.rm=TRUE))
max_val <- as.numeric(global(r_mask, "max", na.rm=TRUE))

if (is.na(min_val) || is.na(max_val)) {
    cat("Warning: All NA values in masked raster.\n")
    r_std <- r_mask
} else if(min_val == max_val) {
  cat("Warning: Min and Max are the same, setting all valid values to 1.\n")
  r_std <- r_mask
  r_std[!is.na(r_std)] <- 1
} else {
  r_std <- (r_mask - min_val) / (max_val - min_val)
}

cat("Writing output to", out_path, "...\n")
writeRaster(r_std, out_path, overwrite=TRUE)
cat("Done!\n")
