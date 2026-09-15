x <- readRDS('DataProcessing/outputs/fitted_model_clean_20260820_1041.rds')
for(n in names(x)) cat(n, paste(class(x[[n]]), collapse='|'), paste(names(x[[n]]), collapse='|'), '\n')
if (inherits(x$fit, 'workflow')) {
  e <- workflows::extract_fit_engine(x$fit)
  cat('FIT_ENGINE', paste(class(e), collapse='|'), '\nFEATURES', paste(e$feature_names, collapse='|'), '\n')
}
if (inherits(x$wf, 'workflow')) cat('WF workflow\n')
