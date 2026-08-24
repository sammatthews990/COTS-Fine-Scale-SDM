############################################################
# refit_agnostic_updated_data.R
# Refit the Year-Agnostic model with the updated 260529
# dataset. This uses CPUE_max and aggregates across all
# years (no Year predictor).
#
# Also refits with ReefGuide stack for comparison.
############################################################

cat("=== Refitting Agnostic Model with Updated Data ===\n")

# --- Configuration ---
use_year      <- FALSE   # Year-agnostic: no Year predictor
predict_year  <- 2026    # Not used in agnostic mode but set for reference
use_reefguide <- FALSE   # Clean 12-predictor stack
cpue_metric   <- NULL    # Uses CPUE_max for agnostic mode

# No year filtering for agnostic — use all available data
# train_year_min <- NULL
# train_year_max <- NULL

# --- Run ---
source("DataProcessing/fit_predict_cots_sdm.R")

cat("\n=== Refit Complete (Agnostic Clean) ===\n")
cat(sprintf("Output TIF: %s\n", output_tif))
cat(sprintf("Output VIP: %s\n", output_vip_plot))
cat(sprintf("Model file: %s\n", model_file))
cat(sprintf("Metrics file: %s\n", metrics_file))

# --- Also refit Agnostic with ReefGuide ---
# cat("\n=== Now refitting Agnostic with ReefGuide stack ===\n")

# # Clean environment for the next run
# rm(list = setdiff(ls(), c("cull_data_file")))

# use_year      <- FALSE
# predict_year  <- 2026
# use_reefguide <- TRUE    # Full stack
# cpue_metric   <- NULL

# source("DataProcessing/fit_predict_cots_sdm.R")

# cat("\n=== Refit Complete (Agnostic Reef Guide) ===\n")
# cat(sprintf("Output TIF: %s\n", output_tif))
# cat(sprintf("Output VIP: %s\n", output_vip_plot))
# cat(sprintf("Model file: %s\n", model_file))
# cat(sprintf("Metrics file: %s\n", metrics_file))
