x <- readRDS('DataProcessing/cots_sdm_v2/models/operational_2026_models.rds')
for (n in names(x)) {
  z <- x[[n]]
  cat('\nNAME', n, '\n')
  cat('CLASS', paste(class(z), collapse='|'), '\n')
  cat('NAMES', paste(names(z), collapse='|'), '\n')
  if (inherits(z, 'Learner')) {
    cat('FEATURES', paste(z$state$train_task$feature_names, collapse='|'), '\n')
  }
  if (!is.null(z$model)) cat('MODEL_CLASS', paste(class(z$model), collapse='|'), '\n')
}
