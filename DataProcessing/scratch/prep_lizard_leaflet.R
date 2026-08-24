library(sf)
library(terra)
library(dplyr)
library(jsonlite)
library(png)
library(viridis)

cat("=== Preparing Multi-Region Leaflet Data ===\n")

# --- Configuration ---
# Set these before sourcing to use timestamped rasters
if (!exists("base_tif_path"))   base_tif_path   <- NULL
if (!exists("refitted_tif_path")) refitted_tif_path <- NULL  # New refitted raster to compare

# Auto-detect: find the most recent timestamped clean raster
output_dir_root <- "DataProcessing/outputs"

if (is.null(base_tif_path)) {
  base_tif_path <- file.path(output_dir_root, "COTS_prob_0.02_cpue_year2025_clean.tif")
  cat("Using baseline (legacy) raster:", base_tif_path, "\n")
}

if (is.null(refitted_tif_path)) {
  candidates <- list.files(output_dir_root, pattern = "^COTS_prob_0\\.02_cpue_year2026_clean_\\d{8}_\\d{4}\\.tif$", full.names = TRUE)
  if (length(candidates) > 0) {
    refitted_tif_path <- candidates[length(candidates)]
    cat("Auto-detected refitted raster:", refitted_tif_path, "\n")
  } else {
    cat("No timestamped refitted raster found — skipping refitted layer.\n")
  }
}

# Also load Agnostic models
old_agnostic_path <- file.path(output_dir_root, "COTS_prob_0.02_cpue_agnostic_clean.tif")

new_agnostic_path <- NULL
agnostic_candidates <- list.files(output_dir_root, pattern = "^COTS_prob_0\\.02_cpue_agnostic_clean_\\d{8}_\\d{4}\\.tif$", full.names = TRUE)
if (length(agnostic_candidates) > 0) {
  new_agnostic_path <- agnostic_candidates[length(agnostic_candidates)]
  cat("Auto-detected new agnostic raster:", new_agnostic_path, "\n")
} else {
  cat("No new agnostic raster found — skipping new agnostic layer.\n")
}

# --- Define Regions ---
# Each region is defined by a name and bounding box (WGS84)
regions <- list(
  list(
    name = "lizard",
    label = "Lizard Island to Cooktown",
    south = -15.75, north = -14.50, west = 144.80, east = 146.10
  ),
  list(
    name = "cairns",
    label = "Cairns Reef Cluster",
    south = -17.30, north = -16.30, west = 145.60, east = 146.60
  ),
  list(
    name = "townsville",
    label = "Townsville Reef Cluster",
    south = -19.30, north = -18.20, west = 146.40, east = 147.70
  ),
  list(
    name = "whitsunday",
    label = "Whitsunday Reef Cluster",
    south = -20.50, north = -19.50, west = 148.50, east = 150.00
  )
)

# Allow selecting specific regions via parameter
if (!exists("selected_regions")) selected_regions <- NULL  # NULL = all

# --- Load Rasters Once ---
cat("Loading rasters...\n")
r_base <- terra::rast(base_tif_path)

r_refitted   <- if (!is.null(refitted_tif_path) && file.exists(refitted_tif_path)) terra::rast(refitted_tif_path) else NULL
r_old_agnos  <- if (file.exists(old_agnostic_path)) terra::rast(old_agnostic_path) else NULL
r_new_agnos  <- if (!is.null(new_agnostic_path) && file.exists(new_agnostic_path)) terra::rast(new_agnostic_path) else NULL

raster_to_gray_png <- function(r, filename, out_dir) {
  vals <- values(r, mat = FALSE)
  valid_idx <- !is.na(vals) & !is.nan(vals)
  
  nr <- nrow(r)
  nc <- ncol(r)
  
  img_array <- array(0, dim = c(nr, nc, 2))
  
  vals_gray <- rep(0, length(vals))
  vals_gray[valid_idx] <- pmin(pmax(vals[valid_idx], 0), 1)
  img_array[,,1] <- matrix(vals_gray, nrow = nr, ncol = nc, byrow = TRUE)
  
  vals_alpha <- rep(0, length(vals))
  vals_alpha[valid_idx] <- 1
  img_array[,,2] <- matrix(vals_alpha, nrow = nr, ncol = nc, byrow = TRUE)
  
  png::writePNG(img_array, file.path(out_dir, filename))
  cat("  Wrote Grayscale PNG:", filename, "\n")
}

# Difference map generator (mapped to 0-1 where 0.5 is no difference)
raster_diff_to_gray_png <- function(r, filename, diff_limit = 0.06, out_dir) {
  vals <- values(r, mat = FALSE)
  valid_idx <- !is.na(vals) & !is.nan(vals)
  
  nr <- nrow(r)
  nc <- ncol(r)
  
  img_array <- array(0, dim = c(nr, nc, 2))
  
  vals_gray <- rep(0, length(vals))
  vals_gray[valid_idx] <- pmin(pmax((vals[valid_idx] - (-diff_limit)) / (2 * diff_limit), 0), 1)
  img_array[,,1] <- matrix(vals_gray, nrow = nr, ncol = nc, byrow = TRUE)
  
  vals_alpha <- rep(0, length(vals))
  vals_alpha[valid_idx] <- 1
  img_array[,,2] <- matrix(vals_alpha, nrow = nr, ncol = nc, byrow = TRUE)
  
  png::writePNG(img_array, file.path(out_dir, filename))
  cat("  Wrote Difference Grayscale PNG:", filename, "\n")
}

# --- Process each region ---
for (reg in regions) {
  if (!is.null(selected_regions) && !(reg$name %in% selected_regions)) next
  
  cat(sprintf("\n--- Processing region: %s (%s) ---\n", reg$label, reg$name))
  
  output_dir <- file.path("DataProcessing/outputs", paste0(reg$name, "_leaflet"))
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  ext_wgs84 <- ext(reg$west, reg$east, reg$south, reg$north)
  ext_proj <- terra::project(ext_wgs84, from = "EPSG:4326", to = crs(r_base))
  
  # Crop and reproject
  crop_base <- tryCatch(terra::crop(r_base, ext_proj), error = function(e) NULL)
  if (is.null(crop_base)) {
    cat("  WARNING: Baseline raster doesn't cover this region. Skipping.\n")
    next
  }
  
  r_base_wgs <- terra::project(crop_base, "EPSG:4326")
  
  # Write bounds
  b <- ext(r_base_wgs)
  bounds_json <- list(south = b$ymin, north = b$ymax, west = b$xmin, east = b$xmax)
  write_json(bounds_json, file.path(output_dir, "bounds.json"), auto_unbox = TRUE)
  
  # Export baseline
  raster_to_gray_png(r_base_wgs, "sdm_baseline.png", out_dir = output_dir)
  
  # Export refitted (if available)
  if (!is.null(r_refitted)) {
    crop_refit <- tryCatch(terra::crop(r_refitted, ext_proj), error = function(e) NULL)
    if (!is.null(crop_refit)) {
      r_refit_wgs <- terra::project(crop_refit, "EPSG:4326")
      r_refit_wgs <- terra::resample(r_refit_wgs, r_base_wgs)
      raster_to_gray_png(r_refit_wgs, "sdm_refitted.png", out_dir = output_dir)
      
      # Difference: refitted - baseline
      r_diff_refit <- r_refit_wgs - r_base_wgs
      raster_diff_to_gray_png(r_diff_refit, "sdm_diff_refitted_vs_base.png", out_dir = output_dir)
    }
  }
  
  # Export old agnostic
  if (!is.null(r_old_agnos)) {
    crop_old_agnos <- tryCatch(terra::crop(r_old_agnos, ext_proj), error = function(e) NULL)
    if (!is.null(crop_old_agnos)) {
      r_old_agnos_wgs <- terra::project(crop_old_agnos, "EPSG:4326")
      r_old_agnos_wgs <- terra::resample(r_old_agnos_wgs, r_base_wgs)
      raster_to_gray_png(r_old_agnos_wgs, "sdm_agnostic_old.png", out_dir = output_dir)
    }
  }
  
  # Export new agnostic
  if (!is.null(r_new_agnos)) {
    crop_new_agnos <- tryCatch(terra::crop(r_new_agnos, ext_proj), error = function(e) NULL)
    if (!is.null(crop_new_agnos)) {
      r_new_agnos_wgs <- terra::project(crop_new_agnos, "EPSG:4326")
      r_new_agnos_wgs <- terra::resample(r_new_agnos_wgs, r_base_wgs)
      raster_to_gray_png(r_new_agnos_wgs, "sdm_agnostic_new.png", out_dir = output_dir)
    }
  }
  
  # Export vector layers (cull polygons + points) filtered to this region
  bbox_sf <- st_as_sfc(st_bbox(c(xmin = reg$west, ymin = reg$south, xmax = reg$east, ymax = reg$north), crs = st_crs(4326)))
  
  cull_polys <- tryCatch(
    st_read("DataProcessing/data/cull_site_polygons.gpkg", quiet = TRUE) %>% st_transform(4326) %>% st_filter(bbox_sf),
    error = function(e) { cat("  No cull polygons found.\n"); NULL }
  )
  
  if (!is.null(cull_polys) && nrow(cull_polys) > 0) {
    base_pts <- st_centroid(cull_polys) %>% mutate(PointType = "Baseline Centroid")
    
    cat(sprintf("  Region layers: %d polys, %d baseline pts\n", nrow(cull_polys), nrow(base_pts)))
    
    st_write(cull_polys, file.path(output_dir, paste0(reg$name, "_polygons.geojson")), delete_dsn = TRUE, quiet = TRUE)
    st_write(base_pts, file.path(output_dir, paste0(reg$name, "_base_pts.geojson")), delete_dsn = TRUE, quiet = TRUE)
  }
  
  cat(sprintf("  Region %s complete!\n", reg$name))
}

# Write region index JSON for the HTML map builder
region_index <- lapply(regions, function(r) {
  if (!is.null(selected_regions) && !(r$name %in% selected_regions)) return(NULL)
  list(name = r$name, label = r$label, south = r$south, north = r$north, west = r$west, east = r$east)
})
region_index <- Filter(Negate(is.null), region_index)
write_json(region_index, "DataProcessing/outputs/region_index.json", auto_unbox = TRUE, pretty = TRUE)

cat("\n=== Multi-region leaflet data prep complete! ===\n")
cat(sprintf("Regions processed: %s\n", paste(sapply(region_index, `[[`, "name"), collapse = ", ")))
