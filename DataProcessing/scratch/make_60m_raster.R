library(terra)

# Input path
in_path <- "outputs/COTS_prob_0.02_cpue_year2025.tif"
out_path <- "outputs/COTS_prob_0.02_cpue_year2025_60m.tif"

message("Loading 15m raster...")
r <- rast(in_path)

message("Aggregating to 60m (fact=4) using mean...")
r_60m <- aggregate(r, fact = 4, fun = "mean", na.rm = TRUE)

message("Saving 60m raster...")
writeRaster(r_60m, out_path, overwrite = TRUE)

message("Done!")
