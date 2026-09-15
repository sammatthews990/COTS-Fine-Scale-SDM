args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
v2_dir <- normalizePath(file.path(dirname(script_file), ".."), winslash = "/")
source(file.path(v2_dir, "R", "bootstrap.R"))
ctx <- initialise_v2(script_file)
cfg <- ctx$config

required_packages <- c(
  "terra", "xgboost", "png", "jsonlite", "base64enc", "viridisLite"
)
missing_packages <- required_packages[!vapply(
  required_packages, requireNamespace, logical(1), quietly = TRUE
)]
if (length(missing_packages) > 0L) {
  stop("Missing visual-comparison packages: ", paste(missing_packages, collapse = ", "))
}

required_files <- c(
  cfg$paths$predictor_stack,
  cfg$paths$predictor_stack_rg,
  cfg$paths$outbreak_forecast_file,
  cfg$paths$legacy_2026_raster,
  cfg$paths$legacy_agnostic_raster,
  file.path(cfg$paths$models_dir, "operational_2026_models.rds"),
  file.path(cfg$paths$models_dir, "additive_operational_bridge.rds"),
  file.path(cfg$paths$artifacts_dir, "site_year_first_visit.rds"),
  file.path(cfg$paths$artifacts_dir, "site_habitat_all_sites.rds"),
  file.path(cfg$paths$outputs_dir, "manta_strong_negative_reefs.csv")
)
if (any(!file.exists(required_files))) {
  stop("Missing visual-comparison inputs: ", paste(
    required_files[!file.exists(required_files)], collapse = "; "
  ))
}

visual_dir <- file.path(cfg$paths$outputs_dir, "visual_comparison")
dir.create(visual_dir, recursive = TRUE, showWarnings = FALSE)

cat("Loading visual-comparison inputs...\n")
predictor_stack <- terra::rast(cfg$paths$predictor_stack)
rg_predictor_stack <- terra::rast(cfg$paths$predictor_stack_rg)
names(rg_predictor_stack)[13:17] <- c(
  "RG_bathy", "RG_slope", "RG_turbid",
  "RG_waves_Hs", "RG_waves_Tp"
)
legacy_2026 <- terra::rast(cfg$paths$legacy_2026_raster)
legacy_agnostic <- terra::rast(cfg$paths$legacy_agnostic_raster)
forecast <- read_outbreak_probabilities(
  cfg$paths$outbreak_forecast_file, "adjusted_forecast"
)
models <- readRDS(file.path(
  cfg$paths$models_dir, "operational_2026_models.rds"
))
static_model <- models[["operational_static__train_through_2025"]]
operational_model <- models[[
  "operational_static_outbreak__train_through_2025"
]]
additive_bridge <- readRDS(file.path(
  cfg$paths$models_dir, "additive_operational_bridge.rds"
))

all_habitat_sites <- readRDS(file.path(
  cfg$paths$artifacts_dir, "site_habitat_all_sites.rds"
))
habitat_sites <- all_habitat_sites |>
  dplyr::filter(
    .data$clean_complete == 1L,
    .data$rg_complete == 1L
  )
cull_clean_habitat_sites <- all_habitat_sites |>
  dplyr::filter(
    .data$evidence_source == "cull_site_history",
    .data$clean_complete == 1L
  )
cull_habitat_sites <- habitat_sites |>
  dplyr::filter(.data$evidence_source == "cull_site_history")
manta_habitat_sites <- habitat_sites |>
  dplyr::filter(.data$evidence_source == "manta_repeated_site_zero")
clean_features <- names(predictor_stack)
rg_features <- names(rg_predictor_stack)
cull_habitat_sites$observation_weight <- 1
cull_clean_habitat_sites$observation_weight <- 1
manta_habitat_sites$observation_weight <- min(
  1,
  0.10 * sum(cull_habitat_sites$high_density == 0L) /
    max(1, nrow(manta_habitat_sites))
)

cat("Fitting full-data site-habitat models for visual QA...\n")
site_habitat_models <- list(
  clean_cull_only = fit_xgb_binary(
    cull_habitat_sites, clean_features, cfg$model,
    seed_offset = 30001L, weight_col = "observation_weight"
  ),
  clean_full_domain = fit_xgb_binary(
    cull_clean_habitat_sites, clean_features, cfg$model,
    seed_offset = 30004L, weight_col = "observation_weight"
  ),
  rg_cull_only = fit_xgb_binary(
    cull_habitat_sites, rg_features, cfg$model,
    seed_offset = 30002L, weight_col = "observation_weight"
  ),
  rg_light_manta = fit_xgb_binary(
    dplyr::bind_rows(cull_habitat_sites, manta_habitat_sites),
    rg_features, cfg$model,
    seed_offset = 30003L, weight_col = "observation_weight"
  )
)
saveRDS(
  site_habitat_models,
  file.path(cfg$paths$models_dir, "site_habitat_full_models.rds")
)

site_year <- readRDS(file.path(
  cfg$paths$artifacts_dir, "site_year_first_visit.rds"
)) |>
  dplyr::filter(.data$year == cfg$outbreak$forecast_year) |>
  dplyr::mutate(
    observed_class = dplyr::if_else(
      .data$high_density == 1L, "High density", "Below threshold"
    )
  )
manta_reefs <- utils::read.csv(
  file.path(cfg$paths$outputs_dir, "manta_strong_negative_reefs.csv"),
  stringsAsFactors = FALSE
)

regions <- v2_visual_regions()
map_data <- list()
layer_summaries <- list()
correlation_summaries <- list()
observed_summaries <- list()

base64_png <- function(path) {
  paste0("data:image/png;base64,", base64enc::base64encode(path))
}

for (region in regions) {
  cat("\nPreparing region: ", region$label, "...\n", sep = "")
  region_dir <- file.path(visual_dir, region$name)
  dir.create(region_dir, recursive = TRUE, showWarnings = FALSE)

  extent_wgs84 <- terra::ext(
    region$west, region$east, region$south, region$north
  )
  extent_native <- terra::project(
    extent_wgs84,
    from = "EPSG:4326",
    to = terra::crs(predictor_stack)
  )
  predictors <- terra::crop(predictor_stack, extent_native)
  rg_predictors <- terra::crop(rg_predictor_stack, extent_native)
  aggregation_factor <- cfg$visual$raster_aggregation_factor
  if (!is.null(aggregation_factor) && aggregation_factor > 1L) {
    cat(
      "  Aggregating the 30 m predictors by a factor of ",
      aggregation_factor, " for visual QA only...\n", sep = ""
    )
    predictors <- terra::aggregate(
      predictors,
      fact = aggregation_factor,
      fun = "mean",
      na.rm = TRUE
    )
    rg_predictors <- terra::aggregate(
      rg_predictors,
      fact = aggregation_factor,
      fun = "mean",
      na.rm = TRUE
    )
  }
  native_template <- predictors[[1]]

  static_path <- file.path(region_dir, "v2_static.tif")
  operational_path <- file.path(region_dir, "v2_operational.tif")
  additive_path <- file.path(region_dir, "v2_additive_bridge.tif")
  habitat_clean_path <- file.path(
    region_dir, "site_habitat_clean_cull_only.tif"
  )
  habitat_clean_full_path <- file.path(
    region_dir, "site_habitat_clean_full_domain.tif"
  )
  habitat_rg_path <- file.path(
    region_dir, "site_habitat_rg_cull_only.tif"
  )
  habitat_rg_manta_path <- file.path(
    region_dir, "site_habitat_rg_light_manta.tif"
  )
  habitat_rg_fallback_path <- file.path(
    region_dir, "site_habitat_rg_with_clean_fallback.tif"
  )
  v2_static <- predict_xgb_raster(predictors, static_model, static_path)
  outbreak <- make_forecast_voronoi_raster(
    forecast,
    native_template,
    region,
    cfg$visual$forecast_voronoi_buffer_degrees
  )
  operational_predictors <- c(predictors, outbreak)
  v2_operational <- predict_xgb_raster(
    operational_predictors,
    operational_model,
    operational_path
  )
  v2_additive_bridge <- predict_additive_bridge_raster(
    v2_static,
    outbreak,
    additive_bridge,
    additive_path
  )
  habitat_clean <- predict_xgb_raster(
    predictors,
    site_habitat_models$clean_cull_only,
    habitat_clean_path
  )
  habitat_clean_full <- predict_xgb_raster(
    predictors,
    site_habitat_models$clean_full_domain,
    habitat_clean_full_path
  )
  habitat_rg <- predict_xgb_raster(
    rg_predictors,
    site_habitat_models$rg_cull_only,
    habitat_rg_path
  )
  habitat_rg_manta <- predict_xgb_raster(
    rg_predictors,
    site_habitat_models$rg_light_manta,
    habitat_rg_manta_path
  )
  habitat_clean_fallback <- terra::resample(
    habitat_clean_full,
    habitat_rg,
    method = "bilinear"
  )
  habitat_rg_fallback <- terra::cover(
    habitat_rg,
    habitat_clean_fallback
  )
  names(habitat_rg_fallback) <- "prediction"
  terra::writeRaster(
    habitat_rg_fallback,
    habitat_rg_fallback_path,
    overwrite = TRUE,
    wopt = list(datatype = "FLT4S", gdal = c("COMPRESS=DEFLATE"))
  )

  old_2026 <- align_raster_to_template(legacy_2026, native_template)
  old_agnostic <- align_raster_to_template(
    legacy_agnostic, native_template
  )
  difference_operational_static <- v2_operational - v2_static
  difference_additive_static <- v2_additive_bridge - v2_static
  difference_operational_legacy <- v2_operational - old_2026

  native_layers <- list(
    legacy_2026 = old_2026,
    legacy_agnostic = old_agnostic,
    site_habitat_clean = habitat_clean,
    site_habitat_rg = habitat_rg,
    site_habitat_rg_with_clean_fallback = habitat_rg_fallback,
    site_habitat_rg_light_manta = habitat_rg_manta,
    v2_static = v2_static,
    outbreak_2026 = outbreak,
    v2_operational = v2_operational,
    v2_additive_bridge = v2_additive_bridge,
    diff_operational_static = difference_operational_static,
    diff_additive_static = difference_additive_static,
    diff_operational_legacy = difference_operational_legacy
  )
  wgs_layers <- project_layers_to_wgs84(native_layers, native_template)
  names(wgs_layers) <- names(native_layers)

  image_paths <- list()
  for (layer_name in names(wgs_layers)) {
    image_path <- file.path(region_dir, paste0(layer_name, ".png"))
    if (startsWith(layer_name, "diff_")) {
      raster_difference_to_gray_png(
        wgs_layers[[layer_name]],
        image_path,
        cfg$visual$difference_limit
      )
    } else {
      raster_to_gray_png(wgs_layers[[layer_name]], image_path)
    }
    image_paths[[layer_name]] <- image_path
    layer_summaries[[length(layer_summaries) + 1L]] <-
      raster_layer_summary(
        native_layers[[layer_name]], region$name, layer_name
      )
  }

  correlation_layers <- c(
    old_2026, old_agnostic, habitat_clean, habitat_rg,
    habitat_rg_fallback,
    habitat_rg_manta, v2_static, outbreak,
    v2_operational, v2_additive_bridge
  )
  names(correlation_layers) <- c(
    "legacy_2026", "legacy_agnostic",
    "site_habitat_clean", "site_habitat_rg",
    "site_habitat_rg_with_clean_fallback",
    "site_habitat_rg_light_manta", "v2_static",
    "outbreak_2026", "v2_operational", "v2_additive_bridge"
  )
  sampled <- terra::spatSample(
    correlation_layers,
    size = 50000,
    method = "regular",
    na.rm = TRUE,
    as.df = TRUE
  )
  correlation <- stats::cor(sampled, use = "pairwise.complete.obs")
  correlation_summaries[[region$name]] <- as.data.frame(as.table(correlation)) |>
    dplyr::transmute(
      region = region$name,
      layer_1 = as.character(.data$Var1),
      layer_2 = as.character(.data$Var2),
      correlation = as.numeric(.data$Freq)
    ) |>
    dplyr::filter(.data$layer_1 < .data$layer_2)

  observations <- site_year |>
    dplyr::filter(
      .data$longitude >= region$west,
      .data$longitude <= region$east,
      .data$latitude >= region$south,
      .data$latitude <= region$north
    )
  if (nrow(observations) > 0L) {
    extracted <- terra::extract(
      correlation_layers,
      terra::vect(
        observations,
        geom = c("longitude", "latitude"),
        crs = "EPSG:4326"
      ) |>
        terra::project(terra::crs(native_template))
    )
    observations <- dplyr::bind_cols(
      observations,
      extracted[, setdiff(names(extracted), "ID"), drop = FALSE]
    )
    observed_summaries[[region$name]] <- observations |>
      dplyr::group_by(.data$observed_class) |>
      dplyr::summarise(
        region = region$name,
        site_years = dplyr::n(),
        mean_legacy_2026 = mean(.data$legacy_2026, na.rm = TRUE),
        mean_legacy_agnostic = mean(.data$legacy_agnostic, na.rm = TRUE),
        mean_site_habitat_clean = mean(
          .data$site_habitat_clean, na.rm = TRUE
        ),
        mean_site_habitat_rg = mean(
          .data$site_habitat_rg, na.rm = TRUE
        ),
        mean_site_habitat_rg_with_clean_fallback = mean(
          .data$site_habitat_rg_with_clean_fallback, na.rm = TRUE
        ),
        mean_site_habitat_rg_light_manta = mean(
          .data$site_habitat_rg_light_manta, na.rm = TRUE
        ),
        mean_v2_static = mean(.data$v2_static, na.rm = TRUE),
        mean_outbreak_2026 = mean(.data$outbreak_2026, na.rm = TRUE),
        mean_v2_operational = mean(.data$v2_operational, na.rm = TRUE),
        mean_v2_additive_bridge = mean(
          .data$v2_additive_bridge, na.rm = TRUE
        ),
        .groups = "drop"
      )
  }

  negative_reefs <- manta_reefs |>
    dplyr::filter(
      .data$reef_longitude >= region$west,
      .data$reef_longitude <= region$east,
      .data$reef_latitude >= region$south,
      .data$reef_latitude <= region$north
    )

  wgs_template <- wgs_layers[[1]]
  bounds <- terra::ext(wgs_template)
  quicklook_path <- file.path(region_dir, paste0(region$name, "_quicklook.png"))
  grDevices::png(quicklook_path, width = 2400, height = 1200, res = 170)
  graphics::par(mfrow = c(3, 4), mar = c(1.5, 1.5, 2.7, 3.2))
  probability_layers <- c(
    "legacy_2026", "legacy_agnostic",
    "site_habitat_clean", "site_habitat_rg",
    "site_habitat_rg_with_clean_fallback",
    "site_habitat_rg_light_manta", "v2_static",
    "outbreak_2026", "v2_operational", "v2_additive_bridge"
  )
  probability_titles <- c(
    "Legacy year 2026", "Legacy agnostic",
    "Site habitat: clean", "Site habitat: RG",
    "Site habitat: RG + clean fallback",
    "Site habitat: RG + light manta", "v2 static",
    "Adjusted reef forecast", "v2 operational (unconstrained)",
    "v2 additive bridge (rank-preserving)"
  )
  for (i in seq_along(probability_layers)) {
    terra::plot(
      wgs_layers[[probability_layers[i]]],
      col = viridisLite::viridis(100),
      range = c(cfg$visual$probability_min, cfg$visual$probability_max),
      main = probability_titles[i],
      axes = FALSE
    )
    if (nrow(observations) > 0L) {
      graphics::points(
        observations$longitude,
        observations$latitude,
        pch = ifelse(observations$high_density == 1L, 19, 1),
        col = ifelse(observations$high_density == 1L, "#ef4444", "#38bdf8"),
        cex = 0.55
      )
    }
  }
  difference_layers <- c(
    "diff_operational_static", "diff_additive_static"
  )
  difference_titles <- c(
    "Unconstrained - static", "Additive bridge - static"
  )
  for (i in seq_along(difference_layers)) {
    terra::plot(
      wgs_layers[[difference_layers[i]]],
      col = grDevices::colorRampPalette(c(
        "#b2182b", "#f7f7f7", "#2166ac"
      ))(101),
      range = c(
        -cfg$visual$difference_limit,
        cfg$visual$difference_limit
      ),
      main = difference_titles[i],
      axes = FALSE
    )
    if (nrow(observations) > 0L) {
      graphics::points(
        observations$longitude,
        observations$latitude,
        pch = ifelse(observations$high_density == 1L, 19, 1),
        col = ifelse(
          observations$high_density == 1L, "#ef4444", "#38bdf8"
        ),
        cex = 0.55
      )
    }
  }
  grDevices::dev.off()

  map_data[[region$name]] <- list(
    name = region$name,
    label = region$label,
    bounds = list(
      south = bounds$ymin, north = bounds$ymax,
      west = bounds$xmin, east = bounds$xmax
    ),
    images = lapply(image_paths, base64_png),
    observations = point_feature_collection(
      observations,
      properties = intersect(c(
        "reef_name", "site_name", "observed_class", "cpue_mean",
        "legacy_2026", "legacy_agnostic",
        "site_habitat_clean", "site_habitat_rg",
        "site_habitat_rg_with_clean_fallback",
        "site_habitat_rg_light_manta", "v2_static",
        "outbreak_2026", "v2_operational", "v2_additive_bridge"
      ), names(observations))
    ),
    manta_negatives = point_feature_collection(
      negative_reefs,
      longitude = "reef_longitude",
      latitude = "reef_latitude",
      properties = intersect(c(
        "reef_name", "manta_tows", "survey_years",
        "cross_method_status", "reef_outbreak_probability_2026"
      ), names(negative_reefs))
    ),
    quicklook = basename(quicklook_path)
  )
}

layer_summary <- dplyr::bind_rows(layer_summaries)
correlation_summary <- dplyr::bind_rows(correlation_summaries)
observed_summary <- dplyr::bind_rows(observed_summaries)
utils::write.csv(
  layer_summary,
  file.path(cfg$paths$outputs_dir, "visual_layer_summary.csv"),
  row.names = FALSE
)
utils::write.csv(
  correlation_summary,
  file.path(cfg$paths$outputs_dir, "visual_layer_correlations.csv"),
  row.names = FALSE
)
utils::write.csv(
  observed_summary,
  file.path(cfg$paths$outputs_dir, "visual_observed_checks.csv"),
  row.names = FALSE
)

map_json <- jsonlite::toJSON(map_data, auto_unbox = TRUE, null = "null")
html <- '<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>COTS SDM v2 Visual Comparison</title>
<link rel="stylesheet" href="https://unpkg.com/leaflet@1.9.4/dist/leaflet.css">
<script src="https://unpkg.com/leaflet@1.9.4/dist/leaflet.js"></script>
<style>
*{box-sizing:border-box}body{margin:0;background:#07111f;color:#e2e8f0;font-family:Inter,Segoe UI,sans-serif}#map{height:100vh;width:100vw}.panel{position:absolute;z-index:1000;top:14px;right:14px;width:390px;max-height:calc(100vh - 28px);overflow:auto;background:rgba(8,20,38,.93);border:1px solid #274060;border-radius:14px;padding:17px;box-shadow:0 16px 42px #0009}.title{font-weight:750;color:#67e8f9;font-size:17px;margin-bottom:4px}.sub{font-size:11px;color:#94a3b8;line-height:1.4;margin-bottom:12px}select{width:100%;background:#17243a;color:#fff;border:1px solid #416180;border-radius:8px;padding:8px;margin-bottom:13px}.section{font-size:11px;color:#a5b4fc;text-transform:uppercase;letter-spacing:.08em;margin:12px 0 7px}.option{display:flex;gap:9px;align-items:center;padding:7px 9px;margin:4px 0;border-radius:7px;background:#17243a99;font-size:12px;cursor:pointer}.option:hover{background:#233655}.row{display:flex;align-items:center;gap:9px;margin:8px 0;font-size:11px}.row input[type=range]{flex:1}.value{color:#67e8f9;width:36px;text-align:right}.legend{height:11px;border-radius:6px;background:linear-gradient(90deg,#440154,#3b528b,#21908d,#5dc963,#fde725);margin-top:8px}.legend-labels{display:flex;justify-content:space-between;font-size:10px;color:#94a3b8;margin-top:3px}.note{font-size:10px;line-height:1.45;color:#94a3b8;margin-top:10px;border-top:1px solid #263b55;padding-top:9px}.leaflet-popup-content-wrapper,.leaflet-popup-tip{background:#0f1d31;color:#e2e8f0}.popup-title{font-weight:700;color:#67e8f9;margin-bottom:5px}.popup-row{font-size:11px;margin:2px 0}.dot{width:10px;height:10px;border-radius:50%;display:inline-block}.high{background:#ef4444}.low{border:2px solid #38bdf8}.manta{background:#facc15}
</style></head><body><div id="map"></div><div class="panel"><div class="title">COTS SDM v2 Inspector</div><div class="sub">Matched-scale comparison of legacy maps, the v2 static surface, adjusted 2026 reef forecast, and the provisional operational model.</div><select id="region"></select><div class="section">Prediction layer</div><div id="layers"></div><div class="section">Display only</div><div class="row"><span>Min</span><input id="min" type="range" min="0" max="1" step=".01" value="0.20"><span class="value" id="minv">0.20</span></div><div class="row"><span>Max</span><input id="max" type="range" min="0" max="1" step=".01" value="0.80"><span class="value" id="maxv">0.80</span></div><div class="row"><span>Opacity</span><input id="opacity" type="range" min="0" max="1" step=".05" value="0.85"><span class="value" id="opacityv">85%</span></div><label class="option"><input id="observations" type="checkbox" checked><span><i class="dot high"></i> High density / <i class="dot low"></i> below threshold</span></label><label class="option"><input id="manta" type="checkbox" checked><span><i class="dot manta"></i> Repeated manta-negative reefs</span></label><div class="legend" id="legend"></div><div class="legend-labels" id="legendLabels"><span>0.20</span><span>0.50</span><span>0.80</span></div><div class="note"><b>Display limits are not model-assessment criteria.</b> Rank tests use exact unbinned prediction values. Red filled observations have first-visit pooled CPUE >= 0.02. Blue rings are below threshold. Manta points retain cross-method conflicts. Difference maps use a fixed +/-0.30 scale.</div></div>
<script>
const data=MAP_DATA;const labels={legacy_2026:"Legacy year model (2026)",legacy_agnostic:"Legacy agnostic",site_habitat_clean:"Site habitat: clean, cull only",site_habitat_rg:"Site habitat: RG, cull only",site_habitat_rg_with_clean_fallback:"Site habitat: RG + clean fallback",site_habitat_rg_light_manta:"Site habitat: RG + light manta",v2_static:"v2 static only",outbreak_2026:"Adjusted reef forecast (2026)",v2_operational:"v2 operational (unconstrained)",v2_additive_bridge:"v2 additive bridge (rank-preserving)",diff_operational_static:"Difference: unconstrained - static",diff_additive_static:"Difference: additive - static",diff_operational_legacy:"Difference: unconstrained - legacy"};
const map=L.map("map",{preferCanvas:true,zoomControl:false});L.control.zoom({position:"bottomleft"}).addTo(map);L.tileLayer("https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}",{attribution:"Esri World Imagery"}).addTo(map);
L.GrayOverlay=L.ImageOverlay.extend({initialize:function(url,bounds,opts){L.ImageOverlay.prototype.initialize.call(this,url,bounds,opts);this.min=opts.min;this.max=opts.max;this.diff=opts.diff;this.opacityValue=opts.opacity;this.vir=this.ramp([[68,1,84],[59,82,139],[33,144,141],[93,201,99],[253,231,37]]);this.div=this.ramp([[178,24,43],[247,247,247],[33,102,172]])},ramp:function(cs){const c=document.createElement("canvas");c.width=256;c.height=1;const x=c.getContext("2d"),g=x.createLinearGradient(0,0,256,0);cs.forEach((v,i)=>g.addColorStop(i/(cs.length-1),`rgb(${v.join(",")})`));x.fillStyle=g;x.fillRect(0,0,256,1);return x.getImageData(0,0,256,1).data},onAdd:function(m){this._canvas=L.DomUtil.create("canvas","leaflet-image-layer leaflet-zoom-animated");this._ctx=this._canvas.getContext("2d");this._img=new Image();this._img.src=this._url;this._img.onload=()=>{this._canvas.width=this._img.width;this._canvas.height=this._img.height;this.redraw()};this._image=this._canvas;L.ImageOverlay.prototype.onAdd.call(this,m)},setRange:function(a,b){this.min=a;this.max=b;this.redraw()},setOpacity:function(v){this.opacityValue=v;this.redraw()},redraw:function(){if(!this._img.complete||!this._img.naturalWidth)return;this._ctx.clearRect(0,0,this._canvas.width,this._canvas.height);this._ctx.drawImage(this._img,0,0);const z=this._ctx.getImageData(0,0,this._canvas.width,this._canvas.height),d=z.data,r=this.diff?this.div:this.vir;for(let i=0;i<d.length;i+=4){if(!d[i+3])continue;let k;if(this.diff){k=d[i]}else{const v=d[i]/255;k=Math.round(255*Math.max(0,Math.min(1,(v-this.min)/(this.max-this.min))))}const j=Math.max(0,Math.min(255,k))*4;d[i]=r[j];d[i+1]=r[j+1];d[i+2]=r[j+2];d[i+3]=Math.round(255*this.opacityValue)}this._ctx.putImageData(z,0,0)}});
let region=Object.keys(data)[0],layer="site_habitat_rg_with_clean_fallback",raster=null,obs=null,manta=null;const rs=document.getElementById("region"),ls=document.getElementById("layers");Object.values(data).forEach(d=>{const o=document.createElement("option");o.value=d.name;o.textContent=d.label;rs.appendChild(o)});Object.entries(labels).forEach(([k,v])=>{const l=document.createElement("label");l.className="option";l.innerHTML=`<input type="radio" name="layer" value="${k}" ${k===layer?"checked":""}><span>${v}</span>`;ls.appendChild(l)});
const popup=(p,title)=>`<div class="popup-title">${title}</div>`+Object.entries(p).map(([k,v])=>`<div class="popup-row"><b>${k.replaceAll("_"," ")}:</b> ${typeof v==="number"?v.toFixed(3):v}</div>`).join("");function load(fit=true){const d=data[region],b=[[d.bounds.south,d.bounds.west],[d.bounds.north,d.bounds.east]];if(raster)map.removeLayer(raster);if(obs)map.removeLayer(obs);if(manta)map.removeLayer(manta);raster=new L.GrayOverlay(d.images[layer],b,{opacity:+document.getElementById("opacity").value,min:+document.getElementById("min").value,max:+document.getElementById("max").value,diff:layer.startsWith("diff_")}).addTo(map);obs=L.geoJSON(d.observations,{pointToLayer:(f,ll)=>L.circleMarker(ll,f.properties.observed_class==="High density"?{radius:4,color:"#fff",weight:1,fillColor:"#ef4444",fillOpacity:.95}:{radius:4,color:"#38bdf8",weight:2,fillOpacity:0}),onEachFeature:(f,l)=>l.bindPopup(popup(f.properties,"2026 site-year"))});manta=L.geoJSON(d.manta_negatives,{pointToLayer:(f,ll)=>L.circleMarker(ll,{radius:5,color:"#111827",weight:1,fillColor:"#facc15",fillOpacity:.95}),onEachFeature:(f,l)=>l.bindPopup(popup(f.properties,"Repeated manta-negative reef"))});if(document.getElementById("observations").checked)obs.addTo(map);if(document.getElementById("manta").checked)manta.addTo(map);if(fit)map.fitBounds(b);const diff=layer.startsWith("diff_");document.getElementById("legend").style.background=diff?"linear-gradient(90deg,#b2182b,#f7f7f7,#2166ac)":"linear-gradient(90deg,#440154,#3b528b,#21908d,#5dc963,#fde725)";document.getElementById("legendLabels").innerHTML=diff?"<span>-0.30</span><span>0</span><span>+0.30</span>":`<span>${(+document.getElementById("min").value).toFixed(2)}</span><span>${((+document.getElementById("min").value + +document.getElementById("max").value)/2).toFixed(2)}</span><span>${(+document.getElementById("max").value).toFixed(2)}</span>`}
rs.onchange=e=>{region=e.target.value;load(true)};document.querySelectorAll("input[name=layer]").forEach(x=>x.onchange=e=>{layer=e.target.value;load(false)});["min","max"].forEach(id=>document.getElementById(id).oninput=e=>{document.getElementById(id+"v").textContent=(+e.target.value).toFixed(2);if(raster&&!layer.startsWith("diff_"))raster.setRange(+document.getElementById("min").value,+document.getElementById("max").value);load(false)});document.getElementById("opacity").oninput=e=>{document.getElementById("opacityv").textContent=Math.round(100*+e.target.value)+"%";raster.setOpacity(+e.target.value)};document.getElementById("observations").onchange=e=>e.target.checked?obs.addTo(map):map.removeLayer(obs);document.getElementById("manta").onchange=e=>e.target.checked?manta.addTo(map):map.removeLayer(manta);load(true);
</script></body></html>'
html <- sub("MAP_DATA", map_json, html, fixed = TRUE)
html_path <- file.path(cfg$paths$outputs_dir, "v2_visual_comparison.html")
writeLines(html, html_path, useBytes = TRUE)

cat("\nVisual comparison complete.\n")
cat("Interactive map:", html_path, "\n")
cat("Quicklooks:", visual_dir, "\n")
print(observed_summary)
