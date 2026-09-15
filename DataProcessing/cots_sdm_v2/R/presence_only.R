fit_maxnet_lq <- function(
    presence, background, features, regmult = 2, classes = "l") {
  x <- dplyr::bind_rows(
    dplyr::select(presence, dplyr::all_of(features)),
    dplyr::select(background, dplyr::all_of(features))
  )
  y <- c(rep(1, nrow(presence)), rep(0, nrow(background)))
  keep <- vapply(x, function(z) {
    z <- z[is.finite(z)]
    length(z) > 1L && stats::sd(z) > 0
  }, logical(1))
  x <- as.data.frame(x[, keep, drop = FALSE])
  formula <- maxnet::maxnet.formula(y, x, classes = classes)
  model <- maxnet::maxnet(y, x, f = formula, regmult = regmult)
  structure(
    list(model = model, features = names(x), classes = classes),
    class = "cots_maxnet"
  )
}

predict_cots_maxnet <- function(object, new_data) {
  as.numeric(stats::predict(
    object$model,
    as.data.frame(new_data[, object$features, drop = FALSE]),
    type = "cloglog"
  ))
}

deduplicate_presence_cells <- function(data) {
  data |>
    dplyr::filter(is.finite(.data$cell_s30m)) |>
    dplyr::arrange(.data$source) |>
    dplyr::distinct(.data$cell_s30m, .keep_all = TRUE)
}

maxent_feature_sets <- function(base_names) {
  topographic <- intersect(
    c("BPI", "EAST", "HCU", "NORTH", "SLO", "SVF", "VCU", "VRM", "DEM"),
    base_names
  )
  biological <- setdiff(base_names, topographic)
  list(
    native_30m = paste0(base_names, "_s30m"),
    multiscale_topography_selected = c(
      "BPI_s30m", "EAST_s120m", "HCU_s120m", "NORTH_s120m",
      "SLO_s60m", "SVF_s120m", "VCU_s30m", "VRM_s60m", "DEM_s30m",
      paste0(biological, "_s30m")
    )
  )
}
