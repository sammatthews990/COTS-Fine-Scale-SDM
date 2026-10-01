# Fit production-length MCMC for one Hmsc component of the joint-model layer.
#
#   Rscript DataProcessing/cots_sdm_v2/final/scripts/03_fit_hmsc_production.R --component=enhanced_joint
#
# Components: enhanced_joint, clean_joint, clean_native_joint. The training
# table, responses, traits and terrain PCA basis are the frozen v4 artifacts;
# only the chain length changes.

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "..", "R", "common.R"))
cfg <- initialise_final()
args <- script_args()
component <- args$component
if (is.null(component) || !component %in% cfg$hmsc$components) {
  stop("--component must be one of: ", paste(cfg$hmsc$components, collapse = ", "))
}

input_name <- paste0("hmsc_", component)
artifact <- readRDS(checked_input(cfg, input_name))
model <- build_hmsc_model(artifact)
chain <- cfg$hmsc$production
smoke <- identical(args$smoke, "1")
if (smoke) {
  # Pipeline test only: tiny chains, outputs tagged so they never feed a layer.
  chain <- list(samples = 20L, thin = 1L, transient = 20L, chains = 2L)
  component_tag <- paste0(component, "_SMOKE")
} else {
  component_tag <- component
}
log_message(
  component, ": ", nrow(artifact$data), " rows, ", length(artifact$responses),
  " responses, ", ncol(model$X), " covariates; samples=", chain$samples,
  " thin=", chain$thin, " transient=", chain$transient, " chains=", chain$chains
)

set.seed(cfg$hmsc$seed)
started <- Sys.time()
model <- Hmsc::sampleMcmc(
  model,
  samples = chain$samples, thin = chain$thin, transient = chain$transient,
  nChains = chain$chains, nParallel = chain$chains, verbose = 0
)
elapsed <- as.numeric(difftime(Sys.time(), started, units = "mins"))
log_message(component, ": sampling finished in ", round(elapsed, 1), " min")

model_path <- file.path(cfg$dirs$models, paste0("hmsc_", component_tag, "_production.rds"))
saveRDS(model, model_path)

convergence <- hmsc_convergence(model, component)
write_table(convergence, cfg, paste0("hmsc_", component_tag, "_convergence.csv"))

# Apparent fit for the target, and a check that the manual linear predictor
# used for mapping reproduces Hmsc's own predictions (catches X-scaling errors).
target <- match(cfg$hmsc$target_response, artifact$responses)
target_rows <- which(!is.na(artifact$data[[cfg$hmsc$target_response]]))
design <- component_design(artifact$data[target_rows, , drop = FALSE], artifact, model)
if (length(design$rows) != length(target_rows)) stop("Design matrix dropped training rows.")
manual <- summarise_target_posterior(design$X, target_beta_draws(model, target))[, "mean"]

# The manual path re-applies the PCA basis from raw terrain columns; Hmsc's
# predict() uses the stored XData. Compare on a 200-draw subset.
pooled <- Hmsc::poolMcmcChains(model$postList)
set.seed(cfg$hmsc$seed)
check_draws <- sort(sample.int(length(pooled), min(200L, length(pooled))))
hmsc_check <- predict(
  model, post = pooled[check_draws],
  XData = as.data.frame(artifact$XData[target_rows, , drop = FALSE]),
  expected = TRUE
)
hmsc_target <- rowMeans(vapply(hmsc_check, function(p) p[, target], numeric(length(target_rows))))
manual_check <- rowMeans(stats::pnorm(design$X %*% target_beta_draws(model, target)[, check_draws]))
reproduction_error <- max(abs(manual_check - hmsc_target))
if (reproduction_error > 1e-6) {
  stop("Manual posterior prediction differs from Hmsc by ", reproduction_error)
}

truth <- artifact$data[[cfg$hmsc$target_response]][target_rows]
hmsc_target <- manual
summary <- tibble::tibble(
  component = component,
  rows = nrow(artifact$data),
  target_rows = length(target_rows),
  target_positives = sum(truth),
  covariates = ncol(model$X),
  samples = chain$samples, thin = chain$thin, transient = chain$transient, chains = chain$chains,
  elapsed_minutes = elapsed,
  min_ess = min(convergence$effective_sample_size, na.rm = TRUE),
  max_psrf_upper = max(convergence$psrf_upper, na.rm = TRUE),
  convergence_pass = min(convergence$effective_sample_size, na.rm = TRUE) >= cfg$hmsc$gates$min_ess &
    max(convergence$psrf_upper, na.rm = TRUE) <= cfg$hmsc$gates$max_psrf_upper,
  apparent_auc = as.numeric(pROC::auc(pROC::roc(truth, hmsc_target, quiet = TRUE))),
  manual_prediction_max_error = reproduction_error,
  model_sha256 = sha256_file(model_path)
)
write_table(summary, cfg, paste0("hmsc_", component_tag, "_production_summary.csv"))
print(summary)
