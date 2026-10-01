# Forward habitat refit for one cutoff year: train the clean_joint Hmsc
# component on rows with year < cutoff only.
#
#   Rscript .../09a_fit_forward_habitat.R --cutoff=2024 --variant=base|context
#
# `base` uses the delivered clean_joint terms; `context` adds the 0.5-1 km
# reef-context terms from 08. Shorter chains than production (validation only).

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "..", "R", "common.R"))
cfg <- initialise_final()
source(file.path(cfg$final_dir, "R", "forward.R"))
args <- script_args()
cutoff <- as.integer(args$cutoff)
variant <- args$variant
stopifnot(is.finite(cutoff), variant %in% c("base", "context"))

artifact <- readRDS(checked_input(cfg, "hmsc_clean_joint"))
keep <- artifact$data$year < cutoff
artifact$data <- artifact$data[keep, , drop = FALSE]
artifact$XData <- artifact$XData[keep, , drop = FALSE]

if (variant == "context") {
  context <- terra::rast(file.path(cfg$dirs$cache, "reef_context_120m.tif"))
  raw <- terra::extract(context, cbind(artifact$data$projected_x, artifact$data$projected_y))
  terms <- context_terms(raw)
  artifact$context_impute <- attr(terms, "impute")
  artifact$XData <- dplyr::bind_cols(artifact$XData, terms)
  artifact$XFormula <- stats::as.formula(paste("~", paste(names(artifact$XData), collapse = " + ")))
}
# Drop responses with no observations before the cutoff (e.g. SALAD pre-2021).
observed <- vapply(artifact$responses, function(r) sum(!is.na(artifact$data[[r]])) >= 30, logical(1))
artifact$responses <- artifact$responses[observed]
artifact$distribution <- artifact$distribution[observed]

model <- build_hmsc_model(artifact)
log_message("cutoff ", cutoff, " ", variant, ": ", nrow(artifact$data), " rows, ",
            length(artifact$responses), " responses, ", ncol(model$X), " covariates")
set.seed(cfg$hmsc$seed + cutoff)
model <- Hmsc::sampleMcmc(model, samples = 500, thin = 20, transient = 5000,
                          nChains = 2, nParallel = 2, verbose = 0)
convergence <- hmsc_convergence(model, paste(variant, cutoff))
dir.create(file.path(cfg$dirs$models, "forward"), showWarnings = FALSE)
saveRDS(list(model = model, artifact = artifact, cutoff = cutoff, variant = variant,
             min_ess = min(convergence$effective_sample_size, na.rm = TRUE),
             max_psrf_upper = max(convergence$psrf_upper, na.rm = TRUE)),
        file.path(cfg$dirs$models, "forward", sprintf("clean_joint_%s_before%d.rds", variant, cutoff)))
log_message("done: min ESS ", round(min(convergence$effective_sample_size, na.rm = TRUE)),
            ", max PSRF ", round(max(convergence$psrf_upper, na.rm = TRUE), 3))
