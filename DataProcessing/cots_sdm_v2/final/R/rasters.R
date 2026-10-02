# Raster helpers shared by the layer and figure scripts.

layer_gdal_options <- c("COMPRESS=DEFLATE", "PREDICTOR=3", "TILED=YES", "BIGTIFF=IF_SAFER")

write_layer <- function(x, path, datatype = "FLT4S") {
  options <- if (datatype %in% c("INT1U", "INT2S", "INT2U")) {
    c("COMPRESS=DEFLATE", "TILED=YES", "BIGTIFF=IF_SAFER")
  } else layer_gdal_options
  terra::writeRaster(x, path, overwrite = TRUE, datatype = datatype, gdal = options)
  log_message("Wrote ", path)
  invisible(path)
}

# Predict a frozen legacy tidymodels/XGBoost workflow for a fixed survey year.
predict_legacy_workflow <- function(stack, workflow, year, path) {
  pred_fun <- function(model, newdata) {
    newdata <- as.data.frame(newdata)
    newdata$Year <- year
    as.numeric(predict(model, new_data = newdata, type = "prob")$.pred_1)
  }
  terra::predict(
    stack, workflow, fun = pred_fun, na.rm = TRUE,
    filename = path, overwrite = TRUE,
    wopt = list(names = "cots_prob", datatype = "FLT4S", gdal = layer_gdal_options)
  )
}

# Paired sample of two aligned rasters over cells where both are finite.
paired_sample <- function(a, b, size = 200000L, seed = 20261001L) {
  set.seed(seed)
  pts <- terra::spatSample(a, size = size, method = "random", na.rm = TRUE,
                           as.points = TRUE, values = FALSE, exhaust = TRUE)
  va <- terra::extract(a, pts, ID = FALSE)[[1]]
  vb <- terra::extract(b, pts, ID = FALSE)[[1]]
  keep <- is.finite(va) & is.finite(vb)
  data.frame(a = va[keep], b = vb[keep])
}

# Rank each cell within its reef polygon: 0 = lowest, 1 = highest.
within_reef_percentile <- function(x, reef_ids) {
  v <- terra::values(x, mat = FALSE)
  id <- terra::values(reef_ids, mat = FALSE)
  out <- rep(NA_real_, length(v))
  ok <- is.finite(v) & !is.na(id)
  out[ok] <- stats::ave(v[ok], id[ok], FUN = function(z) {
    if (length(z) < 2L) return(rep(NA_real_, length(z)))
    (rank(z, ties.method = "average") - 1) / (length(z) - 1)
  })
  terra::setValues(x, out)
}

# Block-mean a 30 m layer onto the 0-anchored 120 m support grid (the same
# 16-point quadrature used for the Hmsc predictors). Cached at `out`.
aggregate_to_support_120 <- function(path, out, template) {
  if (!file.exists(out)) {
    r <- terra::rast(path)
    alignment <- quadrature_alignment(r)
    cropped <- terra::crop(r, alignment$crop_ext, snap = "near")
    aggregated <- terra::aggregate(cropped, fact = 4, fun = "mean", na.rm = TRUE)
    aggregated <- terra::shift(aggregated, dx = alignment$dx, dy = alignment$dy)
    terra::writeRaster(aggregated, out, overwrite = TRUE, datatype = "FLT4S", gdal = layer_gdal_options)
  }
  x <- terra::rast(out)
  if (!terra::compareGeom(x, template, stopOnError = FALSE)) stop(out, " is not on the 120 m support grid.")
  x
}
