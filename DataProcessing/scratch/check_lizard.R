library(sf)
library(terra)

liz_gpkg <- "DataProcessing/data/lizard_sites_june.gpkg"
if (file.exists(liz_gpkg)) {
  sf_liz <- st_read(liz_gpkg, quiet = TRUE)
  cat("Lizard sites GPKG found. CRS:", st_crs(sf_liz)$epsg, "BBOX:\n")
  print(st_bbox(sf_liz))
} else {
  cat("No lizard_sites_june.gpkg found.\n")
}

base_tif <- "DataProcessing/outputs/COTS_prob_0.02_cpue_year2025_clean.tif"
eco_tif <- "DataProcessing/outputs/COTS_prob_0.02_cpue_year2025_clean_ecoCentroid.tif"
pseudo_tif <- "DataProcessing/outputs/COTS_prob_0.02_cpue_year2025_clean_pseudorepsN5.tif"

cat("Base TIF exists:", file.exists(base_tif), "\n")
cat("Eco TIF exists:", file.exists(eco_tif), "\n")
cat("Pseudo TIF exists:", file.exists(pseudo_tif), "\n")

# Let's check cull polygons geopackage
cull_polys <- "DataProcessing/data/cull_site_polygons.gpkg"
eco_pts <- "DataProcessing/data/cull_ecological_centroids.gpkg"
pseudo_pts <- "DataProcessing/data/cull_pseudo_replicates_N5.gpkg"

cat("Cull Polys exists:", file.exists(cull_polys), "\n")
cat("Eco Points exists:", file.exists(eco_pts), "\n")
cat("Pseudo Points exists:", file.exists(pseudo_pts), "\n")
