library(jsonlite)
library(base64enc)
library(dplyr)

cat("=== Building Expanded Region Leaflet Map HTML ===\n")

input_dir   <- "DataProcessing/outputs/lizard_leaflet"
output_html <- "DataProcessing/outputs/lizard_island_sdm_map.html"

# Read bounds
bounds <- read_json(file.path(input_dir, "bounds.json"))

# Base64 encode PNG images
base64_png <- function(file_path) {
  raw_bytes <- readBin(file_path, "raw", file.info(file_path)$size)
  paste0("data:image/png;base64,", base64encode(raw_bytes))
}

img_base   <- base64_png(file.path(input_dir, "sdm_baseline.png"))
img_eco    <- base64_png(file.path(input_dir, "sdm_ecoCentroid.png"))
img_pseudo <- base64_png(file.path(input_dir, "sdm_pseudorepsN5.png"))
img_diff   <- base64_png(file.path(input_dir, "sdm_diff_pseudo_vs_base.png"))

# Read GeoJSON strings
geojson_polys  <- readLines(file.path(input_dir, "lizard_polygons.geojson"), warn = FALSE) %>% paste(collapse = "\n")
geojson_base   <- readLines(file.path(input_dir, "lizard_base_pts.geojson"), warn = FALSE) %>% paste(collapse = "\n")
geojson_eco    <- readLines(file.path(input_dir, "lizard_eco_pts.geojson"), warn = FALSE) %>% paste(collapse = "\n")
geojson_pseudo <- readLines(file.path(input_dir, "lizard_pseudo_pts.geojson"), warn = FALSE) %>% paste(collapse = "\n")

template <- '<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>COTS SDM Expanded Area: Lizard Island to Cooktown</title>
  
  <link rel="stylesheet" href="https://unpkg.com/leaflet@1.9.4/dist/leaflet.css" />
  <script src="https://unpkg.com/leaflet@1.9.4/dist/leaflet.js"></script>
  
  <link rel="preconnect" href="https://fonts.googleapis.com">
  <link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
  <link href="https://fonts.googleapis.com/css2?family=Inter:wght@300;400;500;600;700&display=swap" rel="stylesheet">

  <style>
    * { box-sizing: border-box; margin: 0; padding: 0; }
    body { font-family: "Inter", -apple-system, BlinkMacSystemFont, sans-serif; background: #0f172a; color: #f8fafc; overflow: hidden; }
    #map { width: 100vw; height: 100vh; background: #0b0f19; }

    .control-panel {
      position: absolute;
      top: 16px;
      right: 16px;
      z-index: 1000;
      width: 370px;
      background: rgba(15, 23, 42, 0.88);
      backdrop-filter: blur(12px);
      -webkit-backdrop-filter: blur(12px);
      border: 1px solid rgba(255, 255, 255, 0.12);
      border-radius: 14px;
      padding: 20px;
      box-shadow: 0 20px 40px rgba(0, 0, 0, 0.5);
      max-height: calc(100vh - 32px);
      overflow-y: auto;
    }

    .panel-header {
      border-bottom: 1px solid rgba(255, 255, 255, 0.1);
      padding-bottom: 12px;
      margin-bottom: 16px;
    }
    .panel-title { font-size: 17px; font-weight: 700; color: #38bdf8; display: flex; align-items: center; gap: 8px; }
    .panel-subtitle { font-size: 12px; color: #94a3b8; margin-top: 4px; }

    .metrics-card {
      background: rgba(30, 41, 59, 0.7);
      border-radius: 10px;
      padding: 12px;
      margin-bottom: 16px;
      border: 1px solid rgba(255, 255, 255, 0.06);
    }
    .metrics-title { font-size: 11px; text-transform: uppercase; letter-spacing: 0.05em; color: #94a3b8; font-weight: 600; margin-bottom: 8px; }
    .metrics-grid { display: grid; grid-template-columns: repeat(3, 1fr); gap: 6px; text-align: center; }
    .metric-box { background: rgba(15, 23, 42, 0.6); padding: 8px 4px; border-radius: 6px; }
    .metric-val { font-size: 15px; font-weight: 700; }
    .metric-val.baseline { color: #94a3b8; }
    .metric-val.eco { color: #38bdf8; }
    .metric-val.pseudo { color: #4ade80; }
    .metric-lbl { font-size: 10px; color: #64748b; margin-top: 2px; }

    .section-title { font-size: 12px; font-weight: 600; color: #cbd5e1; margin-bottom: 10px; text-transform: uppercase; letter-spacing: 0.05em; }
    .control-group { margin-bottom: 18px; }

    .option-btn {
      display: flex;
      align-items: center;
      gap: 10px;
      padding: 9px 12px;
      border-radius: 8px;
      background: rgba(30, 41, 59, 0.4);
      border: 1px solid rgba(255, 255, 255, 0.05);
      margin-bottom: 6px;
      cursor: pointer;
      transition: all 0.2s ease;
      font-size: 13px;
    }
    .option-btn:hover { background: rgba(51, 65, 85, 0.6); border-color: rgba(56, 189, 248, 0.3); }
    .option-btn input { accent-color: #38bdf8; cursor: pointer; }

    .slider-container { display: flex; align-items: center; gap: 12px; margin-top: 8px; }
    .slider-container input[type=range] { flex: 1; accent-color: #38bdf8; }
    .slider-val { font-size: 12px; font-weight: 600; width: 36px; text-align: right; color: #38bdf8; }

    .legend-box {
      margin-top: 14px;
      padding-top: 14px;
      border-top: 1px solid rgba(255, 255, 255, 0.1);
    }
    .legend-ramp {
      height: 12px;
      border-radius: 6px;
      background: linear-gradient(to right, #440154, #3b528b, #21908d, #5dc963, #fde725);
      margin: 8px 0 4px 0;
    }
    .legend-labels { display: flex; justify-content: space-between; font-size: 11px; color: #94a3b8; }
    .noise-badge { font-size: 10px; background: rgba(56, 189, 248, 0.15); color: #38bdf8; padding: 2px 6px; border-radius: 4px; margin-top: 4px; display: inline-block; }

    .leaflet-popup-content-wrapper {
      background: #0f172a !important;
      color: #f8fafc !important;
      border-radius: 10px !important;
      border: 1px solid rgba(56, 189, 248, 0.3);
      box-shadow: 0 10px 25px rgba(0, 0, 0, 0.6) !important;
    }
    .leaflet-popup-tip { background: #0f172a !important; }
    .popup-title { font-weight: 700; color: #38bdf8; margin-bottom: 6px; font-size: 13px; }
    .popup-row { font-size: 12px; margin-bottom: 3px; color: #cbd5e1; }
    .popup-row span { color: #94a3b8; }
  </style>
</head>
<body>

  <div id="map"></div>

  <div class="control-panel">
    <div class="panel-header">
      <div class="panel-title">
        <svg width="20" height="20" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"><circle cx="12" cy="12" r="10"/><path d="M12 2a14.5 14.5 0 0 0 0 20 14.5 14.5 0 0 0 0-20"/><path d="M2 12h20"/></svg>
        COTS SDM Expanded Area
      </div>
      <div class="panel-subtitle">Lizard Island to Cooktown (-14.5&deg; to -15.75&deg;)</div>
    </div>

    <div class="metrics-card">
      <div class="metrics-title">Out-of-Fold Spatial Cross-Validation AUC</div>
      <div class="metrics-grid">
        <div class="metric-box">
          <div class="metric-val baseline">0.610</div>
          <div class="metric-lbl">Baseline</div>
        </div>
        <div class="metric-box">
          <div class="metric-val eco">0.624</div>
          <div class="metric-lbl">EcoCentroid</div>
        </div>
        <div class="metric-box">
          <div class="metric-val pseudo">0.641</div>
          <div class="metric-lbl">Pseudo N=5</div>
        </div>
      </div>
    </div>

    <div class="control-group">
      <div class="section-title">Spatial Prediction Rasters</div>
      <label class="option-btn">
        <input type="radio" name="rasterLayer" value="pseudo" checked>
        <span><strong>PseudoReps N=5</strong> (AUC: 0.641)</span>
      </label>
      <label class="option-btn">
        <input type="radio" name="rasterLayer" value="eco">
        <span><strong>EcoCentroid</strong> (AUC: 0.624)</span>
      </label>
      <label class="option-btn">
        <input type="radio" name="rasterLayer" value="base">
        <span><strong>Baseline Centroid</strong> (AUC: 0.610)</span>
      </label>
      <label class="option-btn">
        <input type="radio" name="rasterLayer" value="diff">
        <span><strong>Difference Map</strong> (Pseudo - Base)</span>
      </label>

      <div class="slider-container">
        <span style="font-size: 12px; color: #94a3b8;">Opacity:</span>
        <input type="range" id="opacitySlider" min="0" max="1" step="0.05" value="0.90">
        <span class="slider-val" id="opacityVal">90%</span>
      </div>
    </div>

    <div class="control-group">
      <div class="section-title">Sample Points & Boundaries</div>
      <label class="option-btn">
        <input type="checkbox" id="chkPseudoPts" checked>
        <span style="display:flex; align-items:center; gap:6px;">
          <span style="width:10px; height:10px; border-radius:50%; background:#4ade80; display:inline-block;"></span>
          Pseudo-Replicates N=5 (9,344 pts)
        </span>
      </label>
      <label class="option-btn">
        <input type="checkbox" id="chkEcoPts" checked>
        <span style="display:flex; align-items:center; gap:6px;">
          <span style="width:10px; height:10px; border-radius:50%; background:#38bdf8; display:inline-block;"></span>
          Ecological Centroids (1,866 pts)
        </span>
      </label>
      <label class="option-btn">
        <input type="checkbox" id="chkBasePts" checked>
        <span style="display:flex; align-items:center; gap:6px;">
          <span style="width:10px; height:10px; border-radius:50%; background:#f43f5e; display:inline-block;"></span>
          Geometric Centroids (1,876 pts)
        </span>
      </label>
      <label class="option-btn">
        <input type="checkbox" id="chkPolys">
        <span style="display:flex; align-items:center; gap:6px;">
          <span style="width:10px; height:10px; border:1.5px dashed #fbbf24; display:inline-block;"></span>
          10 ha Cull Site Polygons (1,876)
        </span>
      </label>
    </div>

    <div class="legend-box">
      <div class="section-title" id="legendTitle">High-Signal COTS Risk Probability</div>
      <div class="legend-ramp" id="legendRamp"></div>
      <div class="legend-labels" id="legendLabels">
        <span>0.40 (Low Signal)</span>
        <span>0.49</span>
        <span>0.58+ (High Signal)</span>
      </div>
      <div class="noise-badge" id="noiseBadge">Background noise < 0.38 masked</div>
    </div>
  </div>

  <script>
    const bounds = [[__BOUNDS_SOUTH__, __BOUNDS_WEST__], [__BOUNDS_NORTH__, __BOUNDS_EAST__]];
    
    // High-performance canvas renderer for 10k+ points
    const canvasRenderer = L.canvas();

    const map = L.map("map", {
      center: [-15.125, 145.45],
      zoom: 10,
      zoomControl: false,
      preferCanvas: true
    });

    L.control.zoom({ position: "bottomright" }).addTo(map);

    const esriSat = L.tileLayer("https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}", {
      attribution: "Tiles &copy; Esri"
    }).addTo(map);

    const rasters = {
      pseudo: L.imageOverlay("__IMG_PSEUDO__", bounds, { opacity: 0.90 }).addTo(map),
      eco:    L.imageOverlay("__IMG_ECO__", bounds, { opacity: 0.90 }),
      base:   L.imageOverlay("__IMG_BASE__", bounds, { opacity: 0.90 }),
      diff:   L.imageOverlay("__IMG_DIFF__", bounds, { opacity: 0.90 })
    };

    let currentRaster = "pseudo";

    const polysData  = __GEOJSON_POLYS__;
    const basePtsData = __GEOJSON_BASE__;
    const ecoPtsData  = __GEOJSON_ECO__;
    const pseudoPtsData = __GEOJSON_PSEUDO__;

    const layerPolys = L.geoJSON(polysData, {
      style: {
        color: "#fbbf24",
        weight: 1.2,
        dashArray: "3, 3",
        fillColor: "#fbbf24",
        fillOpacity: 0.03
      },
      onEachFeature: (feature, layer) => {
        layer.bindPopup(`
          <div class="popup-title">${feature.properties.CullSiteName || "Cull Site Polygon"}</div>
          <div class="popup-row"><span>Reef:</span> ${feature.properties.ReefName}</div>
        `);
      }
    });

    const layerBasePts = L.geoJSON(basePtsData, {
      renderer: canvasRenderer,
      pointToLayer: (feature, latlng) => {
        return L.circleMarker(latlng, {
          radius: 3.5,
          fillColor: "#f43f5e",
          color: "#ffffff",
          weight: 1,
          opacity: 0.9,
          fillOpacity: 0.85
        });
      },
      onEachFeature: (feature, layer) => {
        layer.bindPopup(`
          <div class="popup-title">Baseline Geometric Centroid</div>
          <div class="popup-row"><span>Site:</span> ${feature.properties.CullSiteName}</div>
          <div class="popup-row"><span>Reef:</span> ${feature.properties.ReefName}</div>
        `);
      }
    }).addTo(map);

    const layerEcoPts = L.geoJSON(ecoPtsData, {
      renderer: canvasRenderer,
      pointToLayer: (feature, latlng) => {
        return L.circleMarker(latlng, {
          radius: 4,
          fillColor: "#38bdf8",
          color: "#ffffff",
          weight: 1,
          opacity: 0.9,
          fillOpacity: 0.85
        });
      },
      onEachFeature: (feature, layer) => {
        const props = feature.properties;
        layer.bindPopup(`
          <div class="popup-title">Ecological Centroid (Max Weight)</div>
          <div class="popup-row"><span>Site:</span> ${props.CullSiteName}</div>
          <div class="popup-row"><span>Depth (DEM):</span> ${props.DEM != null ? props.DEM.toFixed(1) + "m" : "N/A"}</div>
          <div class="popup-row"><span>Slope:</span> ${props.SLO != null ? props.SLO.toFixed(1) + "&deg;" : "N/A"}</div>
          <div class="popup-row"><span>Habitat Weight:</span> ${props.Weight != null ? props.Weight.toFixed(2) : "N/A"}</div>
        `);
      }
    }).addTo(map);

    const layerPseudoPts = L.geoJSON(pseudoPtsData, {
      renderer: canvasRenderer,
      pointToLayer: (feature, latlng) => {
        return L.circleMarker(latlng, {
          radius: 2.5,
          fillColor: "#4ade80",
          color: "#052e16",
          weight: 0.5,
          opacity: 0.8,
          fillOpacity: 0.75
        });
      },
      onEachFeature: (feature, layer) => {
        const props = feature.properties;
        layer.bindPopup(`
          <div class="popup-title">Pseudo-Replicate Point (N=5)</div>
          <div class="popup-row"><span>Site:</span> ${props.CullSiteName}</div>
          <div class="popup-row"><span>Depth (DEM):</span> ${props.DEM != null ? props.DEM.toFixed(1) + "m" : "N/A"}</div>
          <div class="popup-row"><span>Slope:</span> ${props.SLO != null ? props.SLO.toFixed(1) + "&deg;" : "N/A"}</div>
          <div class="popup-row"><span>Sample Weight:</span> w = 0.20</div>
        `);
      }
    }).addTo(map);

    document.querySelectorAll("input[name=\'rasterLayer\']").forEach(radio => {
      radio.addEventListener("change", (e) => {
        const val = e.target.value;
        map.removeLayer(rasters[currentRaster]);
        currentRaster = val;
        rasters[currentRaster].addTo(map);
        
        const ramp = document.getElementById("legendRamp");
        const title = document.getElementById("legendTitle");
        const labels = document.getElementById("legendLabels");
        const badge = document.getElementById("noiseBadge");
        
        if (val === "diff") {
          title.innerText = "Probability Shift (Pseudo - Base)";
          ramp.style.background = "linear-gradient(to right, #d7191c, #fdae61, #ffffbf, #abd9e9, #2c7bb6)";
          labels.innerHTML = "<span>-0.06 (Lower)</span><span>0.00</span><span>+0.06 (Higher)</span>";
          badge.innerText = "Neutral differences (< 0.005) masked";
        } else {
          title.innerText = "High-Signal COTS Risk Probability";
          ramp.style.background = "linear-gradient(to right, #440154, #3b528b, #21908d, #5dc963, #fde725)";
          labels.innerHTML = "<span>0.40 (Low Signal)</span><span>0.49</span><span>0.58+ (High Signal)</span>";
          badge.innerText = "Background noise < 0.38 masked";
        }
      });
    });

    const opacitySlider = document.getElementById("opacitySlider");
    const opacityVal = document.getElementById("opacityVal");
    opacitySlider.addEventListener("input", (e) => {
      const val = parseFloat(e.target.value);
      opacityVal.innerText = Math.round(val * 100) + "%";
      Object.values(rasters).forEach(r => r.setOpacity(val));
    });

    document.getElementById("chkPseudoPts").addEventListener("change", e => e.target.checked ? map.addLayer(layerPseudoPts) : map.removeLayer(layerPseudoPts));
    document.getElementById("chkEcoPts").addEventListener("change", e => e.target.checked ? map.addLayer(layerEcoPts) : map.removeLayer(layerEcoPts));
    document.getElementById("chkBasePts").addEventListener("change", e => e.target.checked ? map.addLayer(layerBasePts) : map.removeLayer(layerBasePts));
    document.getElementById("chkPolys").addEventListener("change", e => e.target.checked ? map.addLayer(layerPolys) : map.removeLayer(layerPolys));

  </script>
</body>
</html>'

out <- template %>%
  sub("__BOUNDS_SOUTH__", bounds$south, .) %>%
  sub("__BOUNDS_WEST__", bounds$west, .) %>%
  sub("__BOUNDS_NORTH__", bounds$north, .) %>%
  sub("__BOUNDS_EAST__", bounds$east, .) %>%
  sub("__IMG_PSEUDO__", img_pseudo, .) %>%
  sub("__IMG_ECO__", img_eco, .) %>%
  sub("__IMG_BASE__", img_base, .) %>%
  sub("__IMG_DIFF__", img_diff, .) %>%
  sub("__GEOJSON_POLYS__", geojson_polys, .) %>%
  sub("__GEOJSON_BASE__", geojson_base, .) %>%
  sub("__GEOJSON_ECO__", geojson_eco, .) %>%
  sub("__GEOJSON_PSEUDO__", geojson_pseudo, .)

writeLines(out, output_html)
cat("Successfully generated expanded area Leaflet HTML map:", output_html, "\n")
