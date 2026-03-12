# test_mapping.R
# Script to verify model loadings and mappings without running terra heavy geoprocessing

library(maxnet)
library(stringr)

model_dir <- "C:/Users/smatthew/Documents/GitKraken/COTS Fine Scale SDM/DataProcessing/data/SDMModel_Guillame"
models <- c(
    "AhyaD" = file.path(model_dir, "maxent_model_FINAL_maxnet_ext-targetreefs_ALLallres_AhyaD_lq2.RData"),
    "Aspat" = file.path(model_dir, "maxent_model_FINAL_maxnet_ext-targetreefs_ALLallres_Aspat_lq2.RData"),
    "Aten"  = file.path(model_dir, "maxent_model_FINAL_maxnet_ext-targetreefs_ALLallres_Aten_lq2.RData")
)

# New data paths
# Base 15m Data
dem_path_15 <- "C:/Users/smatthew/Documents/GitKraken/COTS Fine Scale SDM/DataProcessing/data/SDMModel_Guillame/processed_dem/TS_bathymetryDEM_15_UTM55_030UTM55_clean_030.tif"
var_dir_15 <- "C:/Users/smatthew/Documents/GitKraken/COTS Fine Scale SDM/DataProcessing/data/SDMModel_Guillame/variables/TS_bathymetryDEM_15_UTM55_030030"

# Base 30m Data (DeepReef30)
dem_path_30 <- "C:/Users/smatthew/Documents/GitKraken/COTS Fine Scale SDM/DataProcessing/data/SDMModel_Guillame/processed_dem/TS_bathymetryDEM_DeepReef30_UTM55_030.tif"
var_dir_30 <- "C:/Users/smatthew/Documents/GitKraken/COTS Fine Scale SDM/DataProcessing/data/SDMModel_Guillame/variables/TS_bathymetryDEM_DeepReef30_UTM55_030030"

print("Starting Model Verification...")

for (sp in names(models)) {
    print(paste("-----------------------------------------"))
    print(paste("Checking species:", sp))

    if (!file.exists(models[[sp]])) {
        print("FAIL: Model file not found")
        next
    }

    load(models[[sp]])
    var.list <- names(me$levels)
    print(paste("SUCCESS: Model Loaded. Expected Variables Found:", length(var.list)))
    print("Required Variables:")
    print(var.list)

    # Check mappings
    missing_vars <- c()
    found_vars <- c()

    for (VAR in var.list) {
        # Parse the model variable (e.g. ACA15_120_DEM or gbr030_060_EAST or gbr100_100_DEM)
        var_parts <- strsplit(VAR, "_")[[1]]
        dataset_prefix <- var_parts[1]
        target_res_str <- var_parts[2]
        var_type <- var_parts[3]

        target_res <- as.numeric(target_res_str)

        # Determine which physical dataset to pull from
        is_aca <- grepl("^ACA", dataset_prefix)

        if (is_aca) {
            base_res <- 15
            base_dem <- dem_path_15
            base_dir <- var_dir_15
            file_pref <- "TS_bathymetryDEM_15_UTM55_030030_"
        } else {
            base_res <- 30
            base_dem <- dem_path_30
            base_dir <- var_dir_30
            file_pref <- "TS_bathymetryDEM_DeepReef30_UTM55_030030_"
        }

        # Resolve file path
        if (var_type == "DEM") {
            temp.fn <- base_dem
        } else {
            temp.fn <- file.path(base_dir, paste0(file_pref, var_type, ".sdat"))
        }

        if (!file.exists(temp.fn)) {
            missing_vars <- c(missing_vars, paste(VAR, "->", temp.fn))
        } else {
            found_vars <- c(found_vars, temp.fn)
        }
    }

    if (length(missing_vars) > 0) {
        print(paste("FAIL: Missing following dataset variables:", paste(missing_vars, collapse = ", ")))
    } else {
        print("SUCCESS: All required variables successfully mapped to the new dataset paths!")
    }
}
print("Verification Complete.")
