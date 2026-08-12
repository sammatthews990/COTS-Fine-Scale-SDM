library(sf)
library(terra)
library(dplyr)
library(jsonlite)
library(png)
library(viridis)

cat("=== Preparing Expanded Region Leaflet Data (1-Degree South, Inshore to Offshore) ===\n")

# Expanded bounding box: Lat -15.75 to -14.50, Lon 144.80 to 146.10
ext_wgs84 <- ext(144.80, 146.10, -15.75, -14.50)

base_tif_path   <- "DataProcessing/outputs/COTS_prob_0.02_cpue_year2025_clean.tif"
eco_tif_path    <- "DataProcessing/outputs/COTS_prob_0.02_cpue_year2025_clean_ecoCentroid.tif"
pseudo_tif_path <- "DataProcessing/outputs/COTS_prob_0.02_cpue_year2025_clean_pseudorepsN5.tif"

r_base   <- terra::rast(base_tif_path)
r_eco    <- terra::rast(eco_tif_path)
r_pseudo <- terra::rast(pseudo_tif_path)

# Crop in native CRS first for maximum speed
ext_proj <- terra::project(ext_wgs84, from = "EPSG:4326", to = crs(r_base))

crop_base   <- terra::crop(r_base, ext_proj)
crop_eco    <- terra::crop(r_eco, ext_proj)
crop_pseudo <- terra::crop(r_pseudo, ext_proj)

r_base_wgs   <- terra::project(crop_base, "EPSG:4326")
r_eco_wgs    <- terra::project(crop_eco, "EPSG:4326")
r_pseudo_wgs <- terra::project(crop_pseudo, "EPSG:4326")

# Align grid rasters exactly
r_base_wgs   <- terra::resample(r_base_wgs, r_base_wgs)
r_eco_wgs    <- terra::resample(r_eco_wgs, r_base_wgs)
r_pseudo_wgs <- terra::resample(r_pseudo_wgs, r_base_wgs)

b <- ext(r_base_wgs)
bounds_json <- list(
  south = b$ymin,
  north = b$ymax,
  west  = b$xmin,
  east  = b$xmax
)

output_dir <- "DataProcessing/outputs/lizard_leaflet"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

write_json(bounds_json, file.path(output_dir, "bounds.json"), auto_unbox = TRUE)

# Color scale mapping function with low-end noise suppression
pal_vir <- viridis::viridis(256)

raster_to_rgba_png <- function(r, filename, min_val = 0.40, max_val = 0.58, noise_cutoff = 0.38, pal_cols = pal_vir, alpha_val = 0.90) {
  vals <- values(r, mat = FALSE)
  valid_idx <- !is.na(vals) & !is.nan(vals) & (vals >= noise_cutoff)
  
  nr <- nrow(r)
  nc <- ncol(r)
  
  img_array <- array(0, dim = c(nr, nc, 4))
  
  norm_vals <- pmin(pmax((vals[valid_idx] - min_val) / (max_val - min_val), 0), 1)
  col_indices <- round(norm_vals * 255) + 1
  
  pal_rgba <- col2rgb(pal_cols, alpha = TRUE) / 255
  
  r_chan <- rep(0, length(vals))
  g_chan <- rep(0, length(vals))
  b_chan <- rep(0, length(vals))
  a_chan <- rep(0, length(vals))
  
  r_chan[valid_idx] <- pal_rgba[1, col_indices]
  g_chan[valid_idx] <- pal_rgba[2, col_indices]
  b_chan[valid_idx] <- pal_rgba[3, col_indices]
  a_chan[valid_idx] <- alpha_val
  
  img_array[,,1] <- matrix(r_chan, nrow = nr, ncol = nc, byrow = TRUE)
  img_array[,,2] <- matrix(g_chan, nrow = nr, ncol = nc, byrow = TRUE)
  img_array[,,3] <- matrix(b_chan, nrow = nr, ncol = nc, byrow = TRUE)
  img_array[,,4] <- matrix(a_chan, nrow = nr, ncol = nc, byrow = TRUE)
  
  png::writePNG(img_array, file.path(output_dir, filename))
  cat("Wrote High-Signal RGBA PNG:", filename, "\n")
}

# Export primary SDM rasters with noise cutoff at 0.38 and high signal range 0.40-0.58
raster_to_rgba_png(r_base_wgs, "sdm_baseline.png", min_val = 0.40, max_val = 0.58, noise_cutoff = 0.38)
raster_to_rgba_png(r_eco_wgs, "sdm_ecoCentroid.png", min_val = 0.40, max_val = 0.58, noise_cutoff = 0.38)
raster_to_rgba_png(r_pseudo_wgs, "sdm_pseudorepsN5.png", min_val = 0.40, max_val = 0.58, noise_cutoff = 0.38)

# Difference map: Pseudo - Baseline
# Highlight differences in [-0.06, +0.06] range, suppress neutral values [-0.005, 0.005]
r_diff <- r_pseudo_wgs - r_base_wgs
pal_diff <- colorRampPalette(c("#d7191c", "#fdae61", "#ffffbf", "#abd9e9", "#2c7bb6"))(256)

raster_diff_to_rgba_png <- function(r, filename, diff_limit = 0.06, deadzone = 0.005) {
  vals <- values(r, mat = FALSE)
  valid_idx <- !is.na(vals) & !is.nan(vals) & (abs(vals) >= deadzone)
  
  nr <- nrow(r)
  nc <- ncol(r)
  
  img_array <- array(0, dim = c(nr, nc, 4))
  
  norm_vals <- pmin(pmax((vals[valid_idx] - (-diff_limit)) / (2 * diff_limit), 0), 1)
  col_indices <- round(norm_vals * 255) + 1
  
  pal_rgba <- col2rgb(pal_diff, alpha = TRUE) / 255
  
  r_chan <- rep(0, length(vals))
  g_chan <- rep(0, length(vals))
  b_chan <- rep(0, length(vals))
  a_chan <- rep(0, length(vals))
  
  r_chan[valid_idx] <- pal_rgba[1, col_indices]
  g_chan[valid_idx] <- pal_rgba[2, col_indices]
  b_chan[valid_idx] <- pal_rgba[3, col_indices]
  a_chan[valid_idx] <- 0.90
  
  img_array[,,1] <- matrix(r_chan, nrow = nr, ncol = nc, byrow = TRUE)
  img_array[,,2] <- matrix(g_chan, nrow = nr, ncol = nc, byrow = TRUE)
  img_array[,,3] <- matrix(b_chan, nrow = nr, ncol = nc, byrow = TRUE)
  img_array[,,4] <- matrix(a_chan, nrow = nr, ncol = nc, byrow = TRUE)
  
  png::writePNG(img_array, file.path(output_dir, filename))
  cat("Wrote Difference Map RGBA PNG:", filename, "\n")
}

raster_diff_to_rgba_png(r_diff, "sdm_diff_pseudo_vs_base.png", diff_limit = 0.06)

# Vector points & polygons in expanded area
bbox_sf <- st_as_sfc(st_bbox(c(xmin = 144.80, ymin = -15.75, xmax = 146.10, ymax = -14.50), crs = st_crs(4326)))

cull_polys <- st_read("DataProcessing/data/cull_site_polygons.gpkg", quiet = TRUE) %>% st_transform(4326) %>% st_filter(bbox_sf)
eco_pts    <- st_read("DataProcessing/data/cull_ecological_centroids.gpkg", quiet = TRUE) %>% st_transform(4326) %>% st_filter(bbox_sf)
pseudo_pts <- st_read("DataProcessing/data/cull_pseudo_replicates_N5.gpkg", quiet = TRUE) %>% st_transform(4326) %>% st_filter(bbox_sf)

base_pts   <- st_centroid(cull_polys) %>% mutate(PointType = "Baseline Centroid")

cat(sprintf("Expanded Area Layers: %d polys, %d baseline pts, %d eco pts, %d pseudo pts\n",
            nrow(cull_polys), nrow(base_pts), nrow(eco_pts), nrow(pseudo_pts)))

st_write(cull_polys, file.path(output_dir, "lizard_polygons.geojson"), delete_dsn = TRUE, quiet = TRUE)
st_write(base_pts, file.path(output_dir, "lizard_base_pts.geojson"), delete_dsn = TRUE, quiet = TRUE)
st_write(eco_pts, file.path(output_dir, "lizard_eco_pts.geojson"), delete_dsn = TRUE, quiet = TRUE)
st_write(pseudo_pts, file.path(output_dir, "lizard_pseudo_pts.geojson"), delete_dsn = TRUE, quiet = TRUE)

cat("Successfully generated expanded area RGBA PNG rasters and GeoJSON point layers!\n")
