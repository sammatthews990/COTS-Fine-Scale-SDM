library(sf)
library(dplyr)

cull_polys <- st_read("DataProcessing/data/cull_site_polygons.gpkg", quiet = TRUE) %>% st_transform(4326)
eco_pts    <- st_read("DataProcessing/data/cull_ecological_centroids.gpkg", quiet = TRUE) %>% st_transform(4326)
pseudo_pts <- st_read("DataProcessing/data/cull_pseudo_replicates_N5.gpkg", quiet = TRUE) %>% st_transform(4326)

# Expanded bounding box
# Lat: -15.75 to -14.50, Lon: 144.80 to 146.10
bbox_sf <- st_as_sfc(st_bbox(c(xmin = 144.80, ymin = -15.75, xmax = 146.10, ymax = -14.50), crs = st_crs(4326)))

exp_polys  <- st_filter(cull_polys, bbox_sf)
exp_eco    <- st_filter(eco_pts, bbox_sf)
exp_pseudo <- st_filter(pseudo_pts, bbox_sf)

cat("Expanded Area Counts:\n")
cat("Polygons:", nrow(exp_polys), "\n")
cat("Eco Points:", nrow(exp_eco), "\n")
cat("Pseudo Points (N=5):", nrow(exp_pseudo), "\n")
