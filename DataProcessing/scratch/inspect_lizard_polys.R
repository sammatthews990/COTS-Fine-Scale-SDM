library(sf)
library(dplyr)

cull_polys <- st_read("DataProcessing/data/cull_site_polygons.gpkg", quiet = TRUE)
cat("Cull Polys total rows:", nrow(cull_polys), "\n")
cat("ReefNames containing Lizard:\n")
lizard_polys <- cull_polys %>% filter(grepl("Lizard", ReefName, ignore.case = TRUE))
print(table(lizard_polys$ReefName))

cat("BBOX of Lizard Reef cull site polygons (EPSG:4326):\n")
lizard_polys_wgs84 <- st_transform(lizard_polys, 4326)
print(st_bbox(lizard_polys_wgs84))
