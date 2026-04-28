############################################################
# diag_layer_extents.R
# Crop to a small far-northern reef area, compare layer
# coverage using the layer with MOST data as benchmark.
############################################################
library(terra)

stack_path <- "C:/Users/smatthew/Documents/GitKraken/COTS Fine Scale SDM/DataProcessing/data/predictors_terra_30m_combined_masked.tif"
r <- terra::rast(stack_path)
nms <- names(r)
cat("Full stack:", nlyr(r), "layers |", ncell(r), "cells\n\n")

# Crop to far northern reef area (Cooktown / Lizard Island region)
# UTM55 coords — covers a few reefs
crop_ext <- ext(300000, 350000, 8320000, 8400000)  # ~50x80km box
cat("Cropping to far-north window:", as.vector(crop_ext), "\n")
rc <- terra::crop(r, crop_ext)
cat("Cropped dims:", nrow(rc), "x", ncol(rc), "=", ncell(rc), "cells\n\n")

# Count non-NA per layer in cropped area
v <- terra::values(rc, mat = TRUE)
counts <- colSums(!is.na(v))
names(counts) <- nms

# Find benchmark (layer with most non-NA)
benchmark_name <- names(which.max(counts))
benchmark_val <- max(counts)
cat(sprintf("BENCHMARK: '%s' with %d non-NA cells\n\n", benchmark_name, benchmark_val))

# Compare all layers to benchmark
cat("=== Per-layer coverage relative to benchmark ===\n")
for (i in seq_along(counts)) {
  pct <- round(100 * counts[i] / benchmark_val, 1)
  gap <- benchmark_val - counts[i]
  flag <- if (pct >= 99) "" else if (pct >= 90) " < partial" else " <<< RESTRICTED"
  cat(sprintf("  [%2d] %-45s  %8d cells  (%5.1f%% of benchmark)%s\n",
              i, nms[i], counts[i], pct, flag))
}

# Group summary
cat("\n=== Summary by group ===\n")
groups <- ifelse(grepl("^GBR_", nms), "env (GBR_*)",
           ifelse(grepl("^SDM_", nms), "coral_old (SDM_*)",
            ifelse(grepl("^RG_", nms), "ReefGuide (RG_*)",
              "coral_new (bare)")))

for (g in unique(groups)) {
  idx <- which(groups == g)
  vals <- counts[idx]
  cat(sprintf("  %-25s  layers: %d  |  min: %8d  max: %8d  (%5.1f - %5.1f%% of benchmark)\n",
              g, length(idx), min(vals), max(vals),
              round(100 * min(vals)/benchmark_val, 1),
              round(100 * max(vals)/benchmark_val, 1)))
}

cat("\nDone.\n")
