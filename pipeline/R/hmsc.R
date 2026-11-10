# Hmsc helpers for the final joint-model layer.
#
# The model specification (responses, traits, terrain PCA basis, conditioning
# on active reef-years) is carried over unchanged from Hmsc corrected v4
# (iterations/hmsc_corrected_v4/scripts/04, 04b and 06b). Only the MCMC length
# differs.

terrain_families <- function(names_in) {
  families <- c("dem", "slope", "eastness", "northness", "hcu", "vcu", "bpi", "sd")
  stats::setNames(lapply(families, function(x) {
    grep(paste0("^", x, "_g(30|60|120|240|360)m$"), names_in, value = TRUE)
  }), families)
}

# Project raw multiscale terrain onto the per-family PCA fitted on training data.
apply_family_pca <- function(data, fits) {
  pieces <- lapply(names(fits), function(family) {
    fit <- fits[[family]]
    z <- stats::predict(fit, newdata = as.data.frame(data[, fit$columns, drop = FALSE]))
    z <- z[, seq_len(fit$n_components), drop = FALSE]
    colnames(z) <- paste0(family, "_pc", seq_len(ncol(z)))
    tibble::as_tibble(z)
  })
  dplyr::bind_cols(pieces)
}

trait_table <- function(responses) {
  data.frame(
    method = factor(dplyr::case_when(
      grepl("^cull", responses) ~ "cull",
      grepl("^salad", responses) ~ "salad",
      grepl("^rhis", responses) ~ "rhis",
      TRUE ~ "other"
    )),
    response_block = factor(ifelse(
      responses %in% c("rhis_branching_live_asin", "rhis_table_live_asin"),
      "coral_resource", "cots_evidence"
    )),
    target_response = as.integer(responses == "cull_high_density"),
    row.names = responses
  )
}

# Rebuild the unfitted Hmsc object exactly as v4 specified it.
build_hmsc_model <- function(artifact) {
  responses <- artifact$responses
  Y <- as.matrix(artifact$data[, responses, drop = FALSE])
  storage.mode(Y) <- "double"
  Hmsc::Hmsc(
    Y = Y, XData = as.data.frame(artifact$XData), XFormula = artifact$XFormula,
    TrData = trait_table(responses),
    TrFormula = ~ method + response_block + target_response,
    distr = artifact$distribution, XScale = TRUE, YScale = TRUE
  )
}

hmsc_convergence <- function(model, component) {
  coda_beta <- Hmsc::convertToCodaObject(model)$Beta
  ess <- coda::effectiveSize(coda_beta)
  psrf <- coda::gelman.diag(coda_beta, multivariate = FALSE)$psrf
  tibble::tibble(
    component = component, parameter = names(ess),
    effective_sample_size = as.numeric(ess),
    psrf_point = as.numeric(psrf[, 1]), psrf_upper = as.numeric(psrf[, 2])
  )
}

# Columns a component needs before its design matrix can be built.
component_required_columns <- function(artifact) {
  pc_names <- grep("_pc[0-9]+$", names(artifact$XData), value = TRUE)
  base_names <- setdiff(names(artifact$XData), pc_names)
  raw_names <- unlist(lapply(artifact$pca, `[[`, "columns"), use.names = FALSE)
  list(base = base_names, required = unique(c(base_names, raw_names)))
}

# Design matrix for the rows of `data` that the component can predict.
# Returns the row positions and X aligned to the model's covariate order.
component_design <- function(data, artifact, model) {
  cols <- component_required_columns(artifact)
  complete <- stats::complete.cases(data[, cols$required, drop = FALSE])
  if (!any(complete)) return(list(rows = integer(0), X = NULL))
  z <- data[complete, , drop = FALSE]
  XData <- dplyr::select(tibble::as_tibble(z), dplyr::all_of(cols$base))
  if (length(artifact$pca)) XData <- dplyr::bind_cols(XData, apply_family_pca(z, artifact$pca))
  if ("geomorphic_class" %in% names(XData)) {
    XData$geomorphic_class <- factor(
      as.character(XData$geomorphic_class),
      levels = levels(artifact$XData$geomorphic_class)
    )
  }
  # Unseen geomorphic classes become NA here and drop to the next fallback.
  valid <- stats::complete.cases(XData)
  X <- stats::model.matrix(artifact$XFormula, XData[valid, , drop = FALSE])
  X <- X[, model$covNames, drop = FALSE]
  list(rows = which(complete)[valid], X = X)
}

# Posterior draws of the target response's coefficients (covariates x draws).
target_beta_draws <- function(model, target, n_draws = NULL, seed = 1L) {
  post <- Hmsc::poolMcmcChains(model$postList)
  draws <- vapply(post, function(s) s$Beta[, target], numeric(model$nc))
  if (!is.null(n_draws) && n_draws < ncol(draws)) {
    set.seed(seed)
    draws <- draws[, sort(sample.int(ncol(draws), n_draws)), drop = FALSE]
  }
  rownames(draws) <- model$covNames
  draws
}

# Posterior summaries of P(target) for each row of X, in chunks to bound memory.
# `transform` is applied to every draw first (used for fallback calibration).
summarise_target_posterior <- function(X, beta_draws, chunk_size = 50000L, transform = identity) {
  n <- nrow(X)
  out <- matrix(NA_real_, n, 4, dimnames = list(NULL, c("mean", "sd", "q05", "q95")))
  starts <- seq.int(1L, n, by = chunk_size)
  for (s in starts) {
    idx <- s:min(n, s + chunk_size - 1L)
    p <- transform(stats::pnorm(X[idx, , drop = FALSE] %*% beta_draws))
    out[idx, "mean"] <- rowMeans(p)
    out[idx, "sd"] <- sqrt(pmax(rowMeans(p^2) - out[idx, "mean"]^2, 0))
    q <- matrixStats::rowQuantiles(p, probs = c(0.05, 0.95))
    out[idx, "q05"] <- q[, 1]
    out[idx, "q95"] <- q[, 2]
  }
  out
}

clip_probability <- function(x) pmin(pmax(x, 1e-5), 1 - 1e-5)

# Probit-linear seam calibration of a fallback onto the primary component's
# probability scale over cells both can predict. Monotone, so it preserves the
# fallback's own ranking.
fit_probit_calibration <- function(primary, fallback) {
  keep <- is.finite(primary) & is.finite(fallback)
  fit <- stats::lm(y ~ x, data = data.frame(
    y = stats::qnorm(clip_probability(primary[keep])),
    x = stats::qnorm(clip_probability(fallback[keep]))
  ))
  list(
    intercept = unname(stats::coef(fit)[1]), slope = unname(stats::coef(fit)[2]),
    overlap_rows = sum(keep),
    spearman = stats::cor(primary[keep], fallback[keep], method = "spearman")
  )
}

apply_probit_calibration <- function(probability, calibration) {
  stats::pnorm(calibration$intercept + calibration$slope * stats::qnorm(clip_probability(probability)))
}
