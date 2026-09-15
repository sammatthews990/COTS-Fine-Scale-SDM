x <- readRDS('DataProcessing/cots_sdm_v2/models/operational_2026_models.rds')
for (n in names(x)) cat(n, ':', paste(x[[n]]$feature_names, collapse='|'), '\n')
