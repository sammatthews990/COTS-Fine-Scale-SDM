manifest <- readRDS('DataProcessing/cots_sdm_v2/artifacts/legacy_clean_year_manifest.rds')
f <- manifest$files$path[manifest$files$role == 'model']
x <- readRDS(f)
cat('FILE', f, '\nCLASS', paste(class(x), collapse='|'), '\nNAMES', paste(names(x), collapse='|'), '\n')
if (inherits(x, 'workflow')) {
  fit <- workflows::extract_fit_engine(x)
  cat('ENGINE', paste(class(fit), collapse='|'), '\nFEATURES', paste(fit$feature_names, collapse='|'), '\n')
}
