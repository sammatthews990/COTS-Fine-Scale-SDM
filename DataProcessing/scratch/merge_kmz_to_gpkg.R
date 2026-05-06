library(sf)

kmz_path <- paste0("/vsizip/", "data/Eotr_CotsCullSites_2025_11_19.kmz")
gpkg_out <- "data/cull_site_polygons.gpkg"

layers <- st_layers(kmz_path)$name
cat("Reading", length(layers), "layers from KMZ...\n")

all_polys <- list()
for (i in seq_along(layers)) {
  lyr <- layers[i]
  tryCatch({
    d <- st_read(kmz_path, layer = lyr, quiet = TRUE)
    d <- st_zm(d, drop = TRUE, what = "ZM")  # Drop Z/M dimensions
    # Extract reef name from layer name (strip "Reef: " prefix)
    reef_name <- sub("^Reef: ", "", lyr)
    d$ReefName <- reef_name
    # Rename Name to CullSiteName for join key
    names(d)[names(d) == "Name"] <- "CullSiteName"
    # Keep only essential columns
    d <- d[, c("CullSiteName", "ReefName", "geometry")]
    all_polys[[i]] <- d
  }, error = function(e) {
    message("Skipping layer: ", lyr, " - ", conditionMessage(e))
  })
  if (i %% 50 == 0) cat("  Processed", i, "/", length(layers), "\n")
}

merged <- do.call(rbind, all_polys)
cat("Total polygons:", nrow(merged), "\n")
cat("Geometry types:", paste(unique(st_geometry_type(merged)), collapse = ", "), "\n")

# Save as GeoPackage
st_write(merged, gpkg_out, delete_dsn = TRUE, quiet = TRUE)
cat("Saved to:", gpkg_out, "\n")
