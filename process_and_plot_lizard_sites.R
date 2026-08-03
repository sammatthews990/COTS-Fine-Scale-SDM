library(terra)

tif_path <- "DataProcessing/outputs/COTS_prob_0.02_cpue_year2025_clean.tif"
gpkg_path <- "DataProcessing/data/lizard_sites_june.gpkg"
out_tif <- "DataProcessing/outputs/COTS_prob_0.02_cpue_year2025_Lizard_Island_masked.tif"
out_png <- "DataProcessing/outputs/COTS_prob_Lizard_Island_categorised_masked.png"

cat("Loading data...\n")
r <- rast(tif_path)
v <- vect(gpkg_path, layer="sites_master") # SPECIFIC LAYER

cat("Cropping and masking...\n")
if(crs(r) != crs(v)) {
  v <- project(v, crs(r))
}

r_crop <- crop(r, v)
r_mask <- mask(r_crop, v)

cat("Standardising values (0-1)...\n")
min_val <- as.numeric(global(r_mask, "min", na.rm=TRUE))
max_val <- as.numeric(global(r_mask, "max", na.rm=TRUE))

if (is.na(min_val) || is.na(max_val)) {
    cat("Warning: All NA values in masked raster.\n")
    r_std <- r_mask
} else if(min_val == max_val) {
  cat("Warning: Min and Max are the same, setting all valid values to 1.\n")
  r_std <- r_mask
  r_std[!is.na(r_std)] <- 1
} else {
  r_std <- (r_mask - min_val) / (max_val - min_val)
}

cat("Writing new TIF...\n")
writeRaster(r_std, out_tif, overwrite=TRUE)

cat("Categorising and plotting...\n")
m <- matrix(c(-0.01, 0.2, 1,
              0.2, 0.4, 2,
              0.4, 0.6, 3,
              0.6, 0.8, 4,
              0.8, 1.01, 5), ncol=3, byrow=TRUE)

r_class <- classify(r_std, m, right=TRUE)
levels(r_class) <- data.frame(id=1:5, Category=c("Very Low (0-0.2)", 
                                                 "Low (0.2-0.4)", 
                                                 "Moderate (0.4-0.6)", 
                                                 "High (0.6-0.8)", 
                                                 "Very High (0.8-1.0)"))

cols <- c("navyblue", "dodgerblue", "khaki1", "darkorange", "firebrick")

png(out_png, width=1200, height=1000, res=150)
plot(r_class, col=cols, main="Standardised COTS Probability - Lizard Island Sites", type="classes", plg=list(x="bottomright"))
plot(v, add=TRUE, border="black", lwd=0.5)
dev.off()

cat("Done!\n")
