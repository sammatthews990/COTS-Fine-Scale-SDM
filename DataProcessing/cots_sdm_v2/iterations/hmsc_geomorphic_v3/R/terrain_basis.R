terrain_families <- function(names_in) {
  families <- c("dem", "slope", "eastness", "northness", "hcu", "vcu", "bpi", "sd")
  stats::setNames(lapply(families, function(x) {
    grep(paste0("^", x, "_g(30|60|120|240|360)m$"), names_in, value = TRUE)
  }), families)
}

fit_family_pca <- function(data, variance = 0.80, max_components = 2L) {
  groups <- terrain_families(names(data))
  fits <- lapply(groups, function(cols) {
    x <- as.data.frame(data[, cols, drop = FALSE])
    complete <- stats::complete.cases(x)
    if (sum(complete) < 20L) return(NULL)
    fit <- stats::prcomp(x[complete, , drop = FALSE], center = TRUE, scale. = TRUE)
    cumulative <- cumsum(fit$sdev^2) / sum(fit$sdev^2)
    fit$n_components <- min(which(cumulative >= variance)[[1]], max_components)
    fit$columns <- cols
    fit
  })
  fits[!vapply(fits, is.null, logical(1))]
}

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

assert_deployable_habitat_terms <- function(terms, forbidden) {
  bad_pattern <- grepl("legacy|outbreak|effort", terms, ignore.case = TRUE) |
    grepl("^(reef(_|$)|year$|method$)", terms, ignore.case = TRUE)
  bad <- unique(terms[terms %in% forbidden | bad_pattern])
  if (length(bad)) stop("Forbidden deployable habitat terms: ", paste(bad, collapse = ", "))
  invisible(TRUE)
}
