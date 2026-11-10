############################################################
# refit_fy2526_clean.R
# Refit the Clean (12-pred) model for Financial Year 2025/26.
# - Predicts for Year = 2026 (current operational year)
# - Option A: Train on ALL data (default — more training signal)
# - Option B: Train on FY window only (2-year rolling, uncomment below)
############################################################

cat("=== Refitting Clean Model for FY 25/26 (predict_year = 2026) ===\n")

# --- Configuration ---
use_year      <- TRUE
predict_year  <- 2026    # Predict for the current year
use_reefguide <- FALSE   # Clean 12-predictor stack
cpue_metric   <- NULL    # Uses CPUE_mean for year-specific mode

# Option B: Restrict training window to FY 25/26 (last 2 calendar years)
# Uncomment the lines below to use only 2025-2026 training data:
# train_year_min <- 2025
# train_year_max <- 2026

# --- Run ---
source("DataProcessing/fit_predict_cots_sdm.R")

cat("\n=== Refit Complete (FY 25/26 Clean, predict_year = 2026) ===\n")
cat(sprintf("Output TIF: %s\n", output_tif))
cat(sprintf("Output VIP: %s\n", output_vip_plot))
cat(sprintf("Model file: %s\n", model_file))
cat(sprintf("Metrics file: %s\n", metrics_file))
