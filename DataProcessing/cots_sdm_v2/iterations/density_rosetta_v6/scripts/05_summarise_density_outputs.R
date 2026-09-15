args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
v6_dir <- normalizePath(file.path(dirname(script_file), ".."), winslash = "/")
v2_dir <- normalizePath(file.path(v6_dir, "..", ".."), winslash = "/")
v4_dir <- normalizePath(file.path(v6_dir, "..", "hmsc_corrected_v4"), winslash = "/")
source(file.path(v6_dir, "R", "bootstrap.R"))
source(file.path(v2_dir, "R", "mapping.R"))
ctx <- initialise_density_v6(script_file)
cfg <- ctx$cfg

cor_raster <- function(a, b) {
  av <- terra::values(a, mat = FALSE); bv <- terra::values(b, mat = FALSE)
  keep <- is.finite(av) & is.finite(bv)
  if (sum(keep) < 3L) return(NA_real_)
  stats::cor(av[keep], bv[keep], method = "spearman")
}
reef_names <- c(
  "Eyrie Reef (14-118)", "Batt Reef (16-029)",
  "(Big) Broadhurst Reef (No 1) (18-100a)", "Bait Reef (19-137)"
)
slugs <- c("eyrie", "batt", "broadhurst", "bait")
map_dir <- file.path(ctx$paths$outputs, "four_reef_density_maps_60m")
manifest <- lapply(seq_along(slugs), function(i) {
  slug <- slugs[[i]]
  cull <- terra::rast(file.path(map_dir, paste0(slug, "_cull_native_60m.tif")))
  salad <- terra::rast(file.path(map_dir, paste0(slug, "_salad_native_60m.tif")))
  balanced <- terra::rast(file.path(map_dir, paste0(slug, "_joint_balanced_60m.tif")))
  legacy <- align_raster_to_template(
    terra::rast(cfg$paths$legacy_2026_raster), salad
  )
  v2 <- align_raster_to_template(
    terra::rast(file.path(
      cfg$paths$v2_outputs_dir, "hmsc_v2_reef_maps",
      sprintf("%02d_%s_hmsc_v2_static_habitat.tif", i, slug)
    )), salad
  )
  v4 <- align_raster_to_template(
    terra::rast(file.path(
      v4_dir, "outputs", "four_reef_maps_120m",
      paste0(slug, "_v4_final_hybrid_120m.tif")
    )), salad
  )
  sv <- terra::values(salad, mat = FALSE)
  tibble::tibble(
    reef_name = reef_names[[i]], support_m = 60L,
    raster_cells = terra::ncell(salad),
    enhanced_cells = sum(is.finite(sv)),
    enhanced_coverage = mean(is.finite(sv)),
    salad_density_median = stats::median(sv, na.rm = TRUE),
    salad_density_p95 = stats::quantile(sv, 0.95, na.rm = TRUE),
    salad_vs_legacy_spearman = cor_raster(salad, legacy),
    salad_vs_v2_spearman = cor_raster(salad, v2),
    salad_vs_v4_spearman = cor_raster(salad, v4),
    cull_vs_salad_spearman = cor_raster(cull, salad),
    balanced_vs_salad_spearman = cor_raster(balanced, salad),
    expert_review_status = "pending",
    figure = file.path("figures", paste0("density_model_comparison_", slug, ".png"))
  )
}) |> dplyr::bind_rows()
utils::write.csv(
  manifest, file.path(ctx$paths$outputs, "four_reef_density_map_manifest.csv"),
  row.names = FALSE
)

coefficients <- list()
for (support_m in c(60L, 120L)) {
  artifact <- readRDS(file.path(
    ctx$paths$artifacts, sprintf("density_data_%dm.rds", support_m)
  ))
  for (model_name in c("cull_native", "salad_native", "joint_rosetta", "joint_rosetta_balanced")) {
    model <- readRDS(file.path(
      ctx$paths$models, sprintf("%s_%dm.rds", model_name, support_m)
    ))
    beta_values <- as.numeric(model$beta)
    coefficients[[length(coefficients) + 1L]] <- tibble::tibble(
      support_m, model = model_name, term = artifact$design$columns,
      coefficient = beta_values
    )
  }
}
utils::write.csv(
  dplyr::bind_rows(coefficients),
  file.path(ctx$paths$outputs, "density_model_coefficients.csv"), row.names = FALSE
)
print(manifest)
