library(terra)
library(dplyr)
library(ggplot2)
library(maxnet)

out_dir <- "c:/Users/smatthew/Documents/GitKraken/COTS Fine Scale SDM/DataProcessing/outputs"
dp_dir  <- "c:/Users/smatthew/Documents/GitKraken/COTS Fine Scale SDM/DataProcessing"

predict_stack_full <- file.path(dp_dir, "data/predictors_terra_30m_full_extended.tif")

m <- readRDS(file.path(out_dir, "maxent_model_cots_10k.rds"))
best_model <- m$model
pred_names <- m$pred_names

output_vip        <- file.path(out_dir, "COTS_maxent_vip_10k.png")
output_response   <- file.path(out_dir, "COTS_maxent_response_10k.png")
output_tif        <- file.path(out_dir, "COTS_maxent_suitability_10k.tif")
output_tif_clean  <- file.path(out_dir, "COTS_maxent_suitability_clean.tif")

print("Generating variable importance plot...")
coef_vec <- best_model[["betas"]]

var_contrib <- data.frame(
  coef_name = names(coef_vec),
  abs_coef  = abs(as.numeric(coef_vec)),
  stringsAsFactors = FALSE
)

get_base_var <- function(cn) {
  if (is.null(cn) || length(cn) == 0 || is.na(cn) || cn == "") return("Unknown")
  cn <- gsub("^(hinge|thresholds|categorical|I)\\((.*)\\)$", "\\2", cn)
  cn <- gsub("\\^[0-9]+$", "", cn)
  sp <- strsplit(cn, "[: ,]")
  if (length(sp) == 0 || length(sp[[1]]) == 0) return("Unknown")
  p <- sp[[1]][sp[[1]] != ""]
  if (length(p) == 0) return("Unknown")
  trimws(p[1])
}

var_contrib$base_var <- vapply(var_contrib$coef_name, get_base_var, FUN.VALUE = character(1))

vip_df <- var_contrib %>%
  group_by(base_var) %>%
  summarise(importance = sum(abs_coef), .groups = "drop") %>%
  arrange(desc(importance)) %>%
  mutate(base_var = factor(base_var, levels = rev(base_var)))

vip_plot <- ggplot(vip_df, aes(x = base_var, y = importance)) +
  geom_col(fill = "#2c7fb8") +
  coord_flip() +
  theme_minimal(base_size = 14) +
  labs(
    title = "MaxEnt Variable Importance (Summed |Coefficients|)",
    x = NULL, y = "Summed Absolute Coefficient"
  )

ggsave(output_vip, plot = vip_plot, width = 8, height = 6, dpi = 300)
print(paste("Saved VIP plot to:", output_vip))

print("Generating response curves...")
png(output_response, width = 1200, height = 900, res = 150)
par(mfrow = c(4, 5), mar = c(3, 3, 2, 1))
plot(best_model, type = "cloglog")
dev.off()
print(paste("Saved response curves to:", output_response))

print("Loading predictor stack...")
predictors_all <- terra::rast(predict_stack_full)
names(predictors_all)[1:12] <- c("BPI", "EAST", "HCU", "NORTH", "SLO", "SVF",
                                  "VCU", "VRM", "DEM", "AhyaD", "Aspat", "Aten")
predictors <- predictors_all[[1:12]]

print("Predicting MaxEnt suitability raster...")
pred_fun <- function(model, newdata) {
  newdata <- as.data.frame(newdata)
  names(newdata) <- pred_names
  p <- predict(model, newdata, type = "cloglog")
  as.numeric(p)
}

suitability <- terra::predict(
  predictors,
  best_model,
  fun   = pred_fun,
  na.rm = TRUE,
  cores = 1
)

names(suitability) <- "suitability"

if (file.exists(output_tif)) file.remove(output_tif)
writeRaster(
  suitability,
  filename  = output_tif,
  overwrite = TRUE,
  gdal      = c("COMPRESS=LZW", "TILED=YES", "BIGTIFF=YES")
)

if (file.exists(output_tif_clean)) file.remove(output_tif_clean)
writeRaster(
  suitability,
  filename  = output_tif_clean,
  overwrite = TRUE,
  gdal      = c("COMPRESS=LZW", "TILED=YES", "BIGTIFF=YES")
)

cat("Done generating MaxEnt Clean suitability rasters!\n")
