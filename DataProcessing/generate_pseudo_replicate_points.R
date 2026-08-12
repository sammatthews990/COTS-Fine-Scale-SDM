############################################################
# generate_pseudo_replicate_points.R
# Extracts ecologically weighted points inside cull site polygons
# (DataProcessing/data/cull_site_polygons.gpkg).
#
# Rules:
# - Downweight shallow & flat pixels (DEM > -3m AND SLO < 2 deg)
# - Retain deeper flat pixels (DEM <= -3m) and sloped pixels
# - Output 1: Ecological Centroids (1 point/polygon, max W)
# - Output 2: Pseudo-Replicate Points (N=5 points/polygon, weighted random sample)
############################################################
library(sf)
library(terra)
library(dplyr)
library(readxl)
library(purrr)

set.seed(123)

# Paths
poly_path   <- "C:/Users/smatthew/Documents/GitKraken/COTS Fine Scale SDM/DataProcessing/data/cull_site_polygons.gpkg"
stack_path  <- "C:/Users/smatthew/Documents/GitKraken/COTS Fine Scale SDM/DataProcessing/data/predictors_clean_12layer.tif"

eco_gpkg_out <- "C:/Users/smatthew/Documents/GitKraken/COTS Fine Scale SDM/DataProcessing/data/cull_ecological_centroids.gpkg"
pseudo_gpkg_out <- "C:/Users/smatthew/Documents/GitKraken/COTS Fine Scale SDM/DataProcessing/data/cull_pseudo_replicates_N5.gpkg"

cat("Loading predictor stack...\n")
predictors <- terra::rast(stack_path)
dem_slo <- predictors[[c("DEM", "SLO")]]

cat("Loading cull site polygons...\n")
polys_sf <- st_read(poly_path, quiet = TRUE)
cat("Loaded", nrow(polys_sf), "polygons.\n")

# Transform polygons to raster CRS if needed
if (st_crs(polys_sf) != crs(predictors)) {
  cat("Transforming polygon CRS to raster CRS...\n")
  polys_sf <- st_transform(polys_sf, crs(predictors))
}

# Pre-extract centroids for fallback
centroids_sf <- st_centroid(polys_sf)
cent_coords <- st_coordinates(centroids_sf)

poly_df <- st_drop_geometry(polys_sf)
poly_df$cent_x <- cent_coords[, 1]
poly_df$cent_y <- cent_coords[, 2]

# Define Weighting Function
calc_weights <- function(dem, slo) {
  w <- rep(1.0, length(dem))
  is_na <- is.na(dem) | is.na(slo)
  w[is_na] <- 0
  valid <- !is_na

  # 1. Shallow & Flat penalty (DEM > -3m AND SLO < 2 deg) -> heavily downweight
  shallow_flat <- valid & (dem > -3.0) & (slo < 2.0)
  w[shallow_flat] <- 0.001
  
  # 2. Shallow & Sloped (Reef Crest) (DEM > -3m AND SLO >= 2 deg) -> high priority
  shallow_slope <- valid & (dem > -3.0) & (slo >= 2.0)
  w[shallow_slope] <- 1.0 + (2.0 * slo[shallow_slope])
  
  # 3. Deeper zone (-30m <= DEM <= -3m) -> flat is ok (1.0), sloped gets extra boost
  deep_zone <- valid & (dem <= -3.0) & (dem >= -30.0)
  w[deep_zone] <- 1.0 + (1.0 * slo[deep_zone])
  
  # 4. Too deep (DEM < -30m)
  too_deep <- valid & (dem < -30.0)
  w[too_deep] <- 0.001
  
  return(w)
}

cat("Extracting raster values across all polygons (chunked processing)...\n")

n_polys <- nrow(poly_df)
chunk_size <- 1000
chunks <- split(seq_len(n_polys), ceiling(seq_len(n_polys) / chunk_size))

eco_rows_list    <- vector("list", n_polys)
pseudo_rows_list <- vector("list", n_polys)

t0 <- Sys.time()
for (c_idx in seq_along(chunks)) {
  idx <- chunks[[c_idx]]
  sub_vect <- terra::vect(polys_sf[idx, ])
  
  # Extract DEM and SLO with cell coordinates
  ext_df <- terra::extract(dem_slo, sub_vect, xy = TRUE, cells = TRUE, ID = TRUE)
  ext_df$orig_idx <- idx[ext_df$ID]
  ext_df$weight <- calc_weights(ext_df$DEM, ext_df$SLO)
  
  grouped <- split(ext_df, ext_df$orig_idx)
  
  for (p_str in names(grouped)) {
    p_i <- as.integer(p_str)
    cell_data <- grouped[[p_str]]
    meta <- poly_df[p_i, ]
    
    cell_data <- cell_data %>% filter(!is.na(weight), weight > 0)
    
    if (nrow(cell_data) == 0) {
      eco_rows_list[[p_i]] <- data.frame(
        CullSiteName = meta$CullSiteName, ReefName = meta$ReefName,
        x = meta$cent_x, y = meta$cent_y, DEM = NA_real_, SLO = NA_real_,
        weight = 1.0, point_type = "Centroid_Fallback", stringsAsFactors = FALSE
      )
      
      pseudo_rows_list[[p_i]] <- data.frame(
        CullSiteName = rep(meta$CullSiteName, 5), ReefName = rep(meta$ReefName, 5),
        x = rep(meta$cent_x, 5), y = rep(meta$cent_y, 5), DEM = rep(NA_real_, 5), SLO = rep(NA_real_, 5),
        weight = rep(0.2, 5), point_id = 1:5, point_type = "Centroid_Fallback", stringsAsFactors = FALSE
      )
      next
    }
    
    # 1. Ecological Centroid: Max weight pixel
    best_row <- cell_data[which.max(cell_data$weight), ]
    eco_rows_list[[p_i]] <- data.frame(
      CullSiteName = meta$CullSiteName, ReefName = meta$ReefName,
      x = best_row$x, y = best_row$y, DEM = best_row$DEM, SLO = best_row$SLO,
      weight = best_row$weight, point_type = "Ecological_Centroid", stringsAsFactors = FALSE
    )
    
    # 2. Pseudo-Replicates: N=5 weighted random sample
    probs <- cell_data$weight / sum(cell_data$weight)
    n_sample <- 5
    replace_flag <- nrow(cell_data) < n_sample
    
    s_idx <- sample.int(nrow(cell_data), size = n_sample, replace = replace_flag, prob = probs)
    s_rows <- cell_data[s_idx, ]
    
    pseudo_rows_list[[p_i]] <- data.frame(
      CullSiteName = rep(meta$CullSiteName, n_sample), ReefName = rep(meta$ReefName, n_sample),
      x = s_rows$x, y = s_rows$y, DEM = s_rows$DEM, SLO = s_rows$SLO,
      weight = rep(1.0 / n_sample, n_sample), point_id = 1:n_sample, point_type = "Pseudo_Replicate",
      stringsAsFactors = FALSE
    )
  }
  
  if (c_idx %% 2 == 0 || c_idx == length(chunks)) {
    cat(sprintf("  Processed %d / %d polygons (%.1f%%)\n", min(c_idx * chunk_size, n_polys), n_polys, 100 * min(c_idx * chunk_size, n_polys) / n_polys))
  }
}

cat(sprintf("Extraction completed in %.2f seconds.\n", as.numeric(difftime(Sys.time(), t0, units = "secs"))))

cat("Binding rows into clean dataframes...\n")
eco_df <- dplyr::bind_rows(eco_rows_list)
pseudo_df <- dplyr::bind_rows(pseudo_rows_list)

cat("Converting to sf objects...\n")
eco_sf <- st_as_sf(eco_df, coords = c("x", "y"), crs = crs(predictors))
pseudo_sf <- st_as_sf(pseudo_df, coords = c("x", "y"), crs = crs(predictors))

cat("Saving Ecological Centroids to:", eco_gpkg_out, "\n")
st_write(eco_sf, eco_gpkg_out, delete_dsn = TRUE, quiet = TRUE)

cat("Saving Pseudo-Replicates to:", pseudo_gpkg_out, "\n")
st_write(pseudo_sf, pseudo_gpkg_out, delete_dsn = TRUE, quiet = TRUE)

# Print Summary Diagnostics
cat("\n=== Extraction Summary & Diagnostics ===\n")
cat("Ecological Centroids count:", nrow(eco_sf), "\n")
cat("Pseudo-Replicates count:", nrow(pseudo_sf), "\n")

cat("\nDEM Distribution Comparison:\n")
eco_dem <- eco_sf$DEM[!is.na(eco_sf$DEM)]
cat(sprintf("  Eco Centroids DEM: Min = %.2fm, Median = %.2fm, Mean = %.2fm, Max = %.2fm\n",
            min(eco_dem), median(eco_dem), mean(eco_dem), max(eco_dem)))

eco_slo <- eco_sf$SLO[!is.na(eco_sf$SLO)]
cat(sprintf("  Eco Centroids SLO: Min = %.2f deg, Median = %.2f deg, Mean = %.2f deg, Max = %.2f deg\n",
            min(eco_slo), median(eco_slo), mean(eco_slo), max(eco_slo)))

pseudo_dem <- pseudo_sf$DEM[!is.na(pseudo_sf$DEM)]
cat(sprintf("  Pseudo Reps DEM: Min = %.2fm, Median = %.2fm, Mean = %.2fm, Max = %.2fm\n",
            min(pseudo_dem), median(pseudo_dem), mean(pseudo_dem), max(pseudo_dem)))

cat("\nDone! Point generation completed successfully.\n")
