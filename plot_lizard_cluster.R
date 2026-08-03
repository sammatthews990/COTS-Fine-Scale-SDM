library(terra)

raster_path <- "DataProcessing/outputs/COTS_prob_0.02_cpue_year2025_Lizard_Island.tif"
gpkg_path <- "DataProcessing/data/lizard_sites_june.gpkg"
out_png <- "DataProcessing/outputs/COTS_prob_Lizard_Island_categorised.png"

cat("Loading data...\n")
r <- rast(raster_path)
v <- vect(gpkg_path)

cat("Categorising raster...\n")
# Create classification matrix
m <- matrix(c(-0.01, 0.2, 1,
              0.2, 0.4, 2,
              0.4, 0.6, 3,
              0.6, 0.8, 4,
              0.8, 1.01, 5), ncol=3, byrow=TRUE)

r_class <- classify(r, m, right=TRUE)

# Set categorical levels
levels(r_class) <- data.frame(id=1:5, Category=c("Very Low (0-0.2)", 
                                                 "Low (0.2-0.4)", 
                                                 "Moderate (0.4-0.6)", 
                                                 "High (0.6-0.8)", 
                                                 "Very High (0.8-1.0)"))

# Define colors
cols <- c("navyblue", "dodgerblue", "khaki1", "darkorange", "firebrick")

cat("Plotting and saving to PNG...\n")
png(out_png, width=1200, height=1000, res=150)
plot(r_class, col=cols, main="Standardised COTS Probability - Lizard Island", type="classes", plg=list(x="bottomright"))
plot(v, add=TRUE, border="black", lwd=0.5)
dev.off()

cat("Done! Plot saved to", out_png, "\n")
