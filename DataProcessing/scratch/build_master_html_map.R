library(jsonlite)
library(base64enc)
library(dplyr)

cat("=== Building Master Multi-Region Leaflet Map HTML ===\n")

region_index_path <- "DataProcessing/outputs/region_index.json"
if (!file.exists(region_index_path)) stop("Region index not found.")
regions <- read_json(region_index_path)

base64_png <- function(file_path) {
  if (!file.exists(file_path)) return(NULL)
  raw_bytes <- readBin(file_path, "raw", file.info(file_path)$size)
  paste0("data:image/png;base64,", base64encode(raw_bytes))
}

read_geojson <- function(file_path) {
  if (!file.exists(file_path)) return(NULL)
  readLines(file_path, warn = FALSE) %>% paste(collapse = "\n") %>% fromJSON(simplifyVector = FALSE)
}

map_data <- list()

for (reg in regions) {
  reg_name <- reg$name
  input_dir <- file.path("DataProcessing/outputs", paste0(reg_name, "_leaflet"))
  
  bounds_file <- file.path(input_dir, "bounds.json")
  if (!file.exists(bounds_file)) next
  
  bounds <- read_json(bounds_file)
  
  images <- list(
    base = base64_png(file.path(input_dir, "sdm_baseline.png")),
    refitted = base64_png(file.path(input_dir, "sdm_refitted.png")),
    diff_refit = base64_png(file.path(input_dir, "sdm_diff_refitted_vs_base.png")),
    agnostic_old = base64_png(file.path(input_dir, "sdm_agnostic_old.png")),
    agnostic_new = base64_png(file.path(input_dir, "sdm_agnostic_new.png"))
  )
  
  vectors <- list(
    polys = read_geojson(file.path(input_dir, paste0(reg_name, "_polygons.geojson"))),
    base_pts = read_geojson(file.path(input_dir, paste0(reg_name, "_base_pts.geojson")))
  )
  
  map_data[[reg_name]] <- list(
    name = reg_name,
    label = reg$label,
    bounds = bounds,
    images = images,
    vectors = vectors
  )
}

map_data_json <- toJSON(map_data, auto_unbox = TRUE)

html_template <- '<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>COTS SDM Master Inspector</title>
  
  <link rel="stylesheet" href="https://unpkg.com/leaflet@1.9.4/dist/leaflet.css" />
  <script src="https://unpkg.com/leaflet@1.9.4/dist/leaflet.js"></script>
  <link href="https://fonts.googleapis.com/css2?family=Inter:wght@300;400;500;600;700&display=swap" rel="stylesheet">

  <style>
    * { box-sizing: border-box; margin: 0; padding: 0; }
    body { font-family: "Inter", sans-serif; background: #0f172a; color: #f8fafc; overflow: hidden; }
    #map { width: 100vw; height: 100vh; background: #0b0f19; }

    .control-panel {
      position: absolute; top: 16px; right: 16px; z-index: 1000; width: 370px;
      background: rgba(15, 23, 42, 0.88); backdrop-filter: blur(12px); border-radius: 14px; padding: 20px;
      border: 1px solid rgba(255, 255, 255, 0.12); box-shadow: 0 20px 40px rgba(0, 0, 0, 0.5);
      max-height: calc(100vh - 32px); overflow-y: auto;
    }

    .panel-header { border-bottom: 1px solid rgba(255, 255, 255, 0.1); padding-bottom: 12px; margin-bottom: 16px; }
    .panel-title { font-size: 17px; font-weight: 700; color: #38bdf8; display: flex; align-items: center; gap: 8px; }
    
    .region-select { width: 100%; margin-top: 12px; padding: 8px 10px; border-radius: 8px; background: rgba(30, 41, 59, 0.9); border: 1px solid rgba(56, 189, 248, 0.4); color: #f8fafc; font-family: inherit; font-size: 14px; cursor: pointer; outline: none; }

    .section-title { font-size: 12px; font-weight: 600; color: #cbd5e1; margin-bottom: 10px; text-transform: uppercase; letter-spacing: 0.05em; }
    .control-group { margin-bottom: 18px; }

    .option-btn { display: flex; align-items: center; gap: 10px; padding: 9px 12px; border-radius: 8px; background: rgba(30, 41, 59, 0.4); border: 1px solid rgba(255, 255, 255, 0.05); margin-bottom: 6px; cursor: pointer; font-size: 13px; }
    .option-btn:hover { background: rgba(51, 65, 85, 0.6); border-color: rgba(56, 189, 248, 0.3); }
    .option-btn input { accent-color: #38bdf8; cursor: pointer; }

    .slider-container { display: flex; align-items: center; gap: 12px; margin-top: 8px; }
    .slider-container input[type=range] { flex: 1; accent-color: #38bdf8; cursor: pointer; }
    .slider-val { font-size: 12px; font-weight: 600; width: 36px; text-align: right; color: #38bdf8; }

    .legend-box { margin-top: 14px; padding-top: 14px; border-top: 1px solid rgba(255, 255, 255, 0.1); }
    .legend-ramp { height: 12px; border-radius: 6px; margin: 8px 0 4px 0; background: linear-gradient(to right, #440154, #3b528b, #21908d, #5dc963, #fde725); }
    .legend-labels { display: flex; justify-content: space-between; font-size: 11px; color: #94a3b8; }
    .noise-badge { font-size: 10px; background: rgba(56, 189, 248, 0.15); color: #38bdf8; padding: 2px 6px; border-radius: 4px; margin-top: 4px; display: inline-block; }

    .leaflet-popup-content-wrapper { background: #0f172a !important; color: #f8fafc !important; border-radius: 10px !important; border: 1px solid rgba(56, 189, 248, 0.3); }
    .leaflet-popup-tip { background: #0f172a !important; }
    .popup-title { font-weight: 700; color: #38bdf8; margin-bottom: 6px; font-size: 13px; }
    .popup-row { font-size: 12px; margin-bottom: 3px; color: #cbd5e1; }
  </style>
</head>
<body>

  <div id="map"></div>

  <div class="control-panel">
    <div class="panel-header">
      <div class="panel-title">
        <svg width="20" height="20" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"><circle cx="12" cy="12" r="10"/><path d="M12 2a14.5 14.5 0 0 0 0 20 14.5 14.5 0 0 0 0-20"/><path d="M2 12h20"/></svg>
        COTS SDM Master Inspector
      </div>
      <select class="region-select" id="regionSelect"></select>
    </div>

    <div class="control-group">
      <div class="section-title">Spatial Prediction Rasters</div>
      <label class="option-btn"><input type="radio" name="rasterLayer" value="base" checked><span><strong>Baseline (2025)</strong></span></label>
      <label class="option-btn"><input type="radio" name="rasterLayer" value="refitted"><span><strong>Refitted (2026)</strong></span></label>
      <label class="option-btn"><input type="radio" name="rasterLayer" value="diff_refit"><span><strong>Difference</strong> (Refitted - Baseline)</span></label>
      <label class="option-btn"><input type="radio" name="rasterLayer" value="agnostic_old"><span><strong>Agnostic (2025)</strong></span></label>
      <label class="option-btn" id="lblAgnosticNew"><input type="radio" name="rasterLayer" value="agnostic_new"><span><strong>Agnostic (2026)</strong></span></label>

      <div class="slider-container" style="margin-top: 14px;">
        <span style="font-size: 12px; color: #94a3b8;">Ramp Min:</span>
        <input type="range" id="rampMinSlider" min="0" max="1" step="0.01" value="0.45">
        <span class="slider-val" id="rampMinVal">0.45</span>
      </div>
      <div class="slider-container">
        <span style="font-size: 12px; color: #94a3b8;">Ramp Max:</span>
        <input type="range" id="rampMaxSlider" min="0" max="1" step="0.01" value="1.00">
        <span class="slider-val" id="rampMaxVal">1.00</span>
      </div>
      <div class="slider-container">
        <span style="font-size: 12px; color: #94a3b8;">Opacity:</span>
        <input type="range" id="opacitySlider" min="0" max="1" step="0.05" value="0.90">
        <span class="slider-val" id="opacityVal">90%</span>
      </div>
    </div>

    <div class="control-group">
      <div class="section-title">Sample Points & Boundaries</div>
      <label class="option-btn">
        <input type="checkbox" id="chkBasePts" checked>
        <span style="display:flex; align-items:center; gap:6px;">
          <span style="width:10px; height:10px; border-radius:50%; background:#f43f5e; display:inline-block;"></span>
          Geometric Centroids
        </span>
      </label>
      <label class="option-btn">
        <input type="checkbox" id="chkPolys">
        <span style="display:flex; align-items:center; gap:6px;">
          <span style="width:10px; height:10px; border:1.5px dashed #fbbf24; display:inline-block;"></span>
          Cull Site Polygons
        </span>
      </label>
    </div>

    <div class="legend-box">
      <div class="section-title" id="legendTitle">COTS Risk Probability</div>
      <div class="legend-ramp" id="legendRamp"></div>
      <div class="legend-labels" id="legendLabels">
        <span>0.40 (Low)</span><span>0.49</span><span>0.58+ (High)</span>
      </div>
      <div class="noise-badge" id="noiseBadge">Color ramp stretches between min & max</div>
    </div>
  </div>

  <script>
    const mapData = MAP_DATA_PLACEHOLDER;
    
    const map = L.map("map", { zoomControl: false, preferCanvas: true });
    L.control.zoom({ position: "bottomright" }).addTo(map);
    L.tileLayer("https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}", {
      attribution: "Tiles &copy; Esri"
    }).addTo(map);

    // Custom Canvas Overlay
    L.CanvasImageOverlay = L.ImageOverlay.extend({
      initialize: function (url, bounds, options) {
        L.ImageOverlay.prototype.initialize.call(this, url, bounds, options);
        this.rampMin = options.rampMin || 0.45;
        this.rampMax = (options.rampMax !== undefined) ? options.rampMax : 1.0;
        this.isDiff = options.isDiff || false;
        
        this.rampViridis = this._buildRamp([[68,1,84], [59,82,139], [33,144,141], [93,201,99], [253,231,37]]);
        this.rampDiff = this._buildRamp([[215,25,28], [253,174,97], [255,255,191], [171,217,233], [44,123,182]]);
      },
      _buildRamp: function(colors) {
        const canvas = document.createElement("canvas");
        canvas.width = 256; canvas.height = 1;
        const ctx = canvas.getContext("2d");
        const grad = ctx.createLinearGradient(0,0,256,0);
        for(let i=0; i<colors.length; i++) {
           grad.addColorStop(i/(colors.length-1), `rgb(${colors[i].join(",")})`);
        }
        ctx.fillStyle = grad; ctx.fillRect(0,0,256,1);
        return ctx.getImageData(0,0,256,1).data;
      },
      onAdd: function (map) {
        this._canvas = L.DomUtil.create("canvas", "leaflet-image-layer leaflet-zoom-animated");
        this._ctx = this._canvas.getContext("2d");
        
        this._img = new Image();
        this._img.src = this._url;
        this._img.onload = () => {
          this._canvas.width = this._img.width;
          this._canvas.height = this._img.height;
          this._redrawCanvas();
        };
        
        this._image = this._canvas;
        L.ImageOverlay.prototype.onAdd.call(this, map);
      },
      setRampRange: function(min, max) {
        this.rampMin = min;
        this.rampMax = max;
        this._redrawCanvas();
      },
      setOpacity: function(val) {
        this.options.opacity = val;
        this._redrawCanvas();
      },
      _redrawCanvas: function () {
        if (!this._img.complete || this._img.naturalWidth === 0) return;
        this._ctx.clearRect(0, 0, this._canvas.width, this._canvas.height);
        this._ctx.drawImage(this._img, 0, 0);
        let imgData = this._ctx.getImageData(0, 0, this._canvas.width, this._canvas.height);
        let data = imgData.data;
        let rampMin = this.rampMin * 255;
        let rampMax = this.rampMax * 255;
        if (rampMin >= rampMax) rampMax = rampMin + 1; // prevent div zero
        
        let ramp = this.isDiff ? this.rampDiff : this.rampViridis;
        let alphaMod = (this.options.opacity !== undefined ? this.options.opacity : 0.9);
        
        for (let i = 0; i < data.length; i += 4) {
           let origAlpha = data[i+3];
           if (origAlpha === 0) continue;
           
           let val = data[i]; 
           let rIdx = 0;
           
           if (!this.isDiff) {
             let v = val <= rampMin ? 0 : val >= rampMax ? 255 : ((val - rampMin) / (rampMax - rampMin)) * 255;
             rIdx = Math.floor(v) * 4;
           } else {
             rIdx = val * 4;
           }
           
           rIdx = Math.max(0, Math.min(1020, rIdx));
           
           data[i] = ramp[rIdx];
           data[i+1] = ramp[rIdx+1];
           data[i+2] = ramp[rIdx+2];
           data[i+3] = Math.floor(255 * alphaMod);
        }
        this._ctx.putImageData(imgData, 0, 0);
      }
    });

    let currentRegion = Object.keys(mapData)[0];
    let currentRasterType = "base";
    let activeRaster = null;
    let activePolys = null;
    let activeBasePts = null;
    
    // UI Elements
    const regionSelect = document.getElementById("regionSelect");
    const rampMinSlider = document.getElementById("rampMinSlider");
    const rampMinVal = document.getElementById("rampMinVal");
    const rampMaxSlider = document.getElementById("rampMaxSlider");
    const rampMaxVal = document.getElementById("rampMaxVal");
    const opacitySlider = document.getElementById("opacitySlider");
    const opacityVal = document.getElementById("opacityVal");
    const chkBasePts = document.getElementById("chkBasePts");
    const chkPolys = document.getElementById("chkPolys");
    
    // Populate dropdown
    Object.values(mapData).forEach(r => {
      let opt = document.createElement("option");
      opt.value = r.name; opt.innerText = r.label;
      regionSelect.appendChild(opt);
    });
    
    function loadRegion(regName) {
      currentRegion = regName;
      const d = mapData[regName];
      const bounds = [[d.bounds.south, d.bounds.west], [d.bounds.north, d.bounds.east]];
      
      if (!map._loaded) {
        map.fitBounds(bounds);
      } else {
        map.flyToBounds(bounds, { duration: 1.5 });
      }
      
      if (activeRaster) map.removeLayer(activeRaster);
      if (activePolys) map.removeLayer(activePolys);
      if (activeBasePts) map.removeLayer(activeBasePts);
      
      const imgSrc = d.images[currentRasterType];
      if (imgSrc) {
        const isDiff = currentRasterType.startsWith("diff");
        activeRaster = new L.CanvasImageOverlay(imgSrc, bounds, { 
          opacity: parseFloat(opacitySlider.value),
          rampMin: parseFloat(rampMinSlider.value),
          rampMax: parseFloat(rampMaxSlider.value),
          isDiff: isDiff
        }).addTo(map);
      }
      
      if (d.vectors.polys) {
        activePolys = L.geoJSON(d.vectors.polys, {
          style: { color: "#fbbf24", weight: 1.2, dashArray: "3, 3", fillColor: "#fbbf24", fillOpacity: 0.03 },
          onEachFeature: (f, l) => l.bindPopup(`<div class="popup-title">${f.properties.CullSiteName || "Cull Site"}</div><div class="popup-row"><span>Reef:</span> ${f.properties.ReefName || ""}</div>`)
        });
        if (chkPolys.checked) activePolys.addTo(map);
      }
      
      if (d.vectors.base_pts) {
        activeBasePts = L.geoJSON(d.vectors.base_pts, {
          renderer: L.canvas(),
          pointToLayer: (f, ll) => L.circleMarker(ll, { radius: 3.5, fillColor: "#f43f5e", color: "#fff", weight: 1, opacity: 0.9, fillOpacity: 0.85 }),
          onEachFeature: (f, l) => l.bindPopup(`<div class="popup-title">Baseline Centroid</div><div class="popup-row"><span>Site:</span> ${f.properties.CullSiteName || ""}</div><div class="popup-row"><span>Reef:</span> ${f.properties.ReefName || ""}</div>`)
        });
        if (chkBasePts.checked) activeBasePts.addTo(map);
      }
      
      document.getElementById("lblAgnosticNew").style.display = d.images.agnostic_new ? "flex" : "none";
    }
    
    regionSelect.addEventListener("change", e => loadRegion(e.target.value));
    
    document.querySelectorAll("input[name=\'rasterLayer\']").forEach(radio => {
      radio.addEventListener("change", e => {
        currentRasterType = e.target.value;
        loadRegion(currentRegion);
        
        const ramp = document.getElementById("legendRamp");
        const title = document.getElementById("legendTitle");
        const labels = document.getElementById("legendLabels");
        
        if (currentRasterType.startsWith("diff")) {
          title.innerText = "Probability Difference";
          ramp.style.background = "linear-gradient(to right, #d7191c, #fdae61, #ffffbf, #abd9e9, #2c7bb6)";
          labels.innerHTML = "<span>-0.06 (Lower)</span><span>0.00</span><span>+0.06 (Higher)</span>";
        } else {
          title.innerText = "COTS Risk Probability";
          ramp.style.background = "linear-gradient(to right, #440154, #3b528b, #21908d, #5dc963, #fde725)";
          labels.innerHTML = "<span>0.40 (Low)</span><span>0.49</span><span>0.58+ (High)</span>";
        }
      });
    });

    function updateRamps() {
      const minVal = parseFloat(rampMinSlider.value);
      const maxVal = parseFloat(rampMaxSlider.value);
      rampMinVal.innerText = minVal.toFixed(2);
      rampMaxVal.innerText = maxVal.toFixed(2);
      if (activeRaster) activeRaster.setRampRange(minVal, maxVal);
    }
    rampMinSlider.addEventListener("input", updateRamps);
    rampMaxSlider.addEventListener("input", updateRamps);

    opacitySlider.addEventListener("input", e => {
      const val = parseFloat(e.target.value);
      opacityVal.innerText = Math.round(val * 100) + "%";
      if (activeRaster) activeRaster.setOpacity(val);
    });

    chkBasePts.addEventListener("change", e => { if (activeBasePts) e.target.checked ? activeBasePts.addTo(map) : map.removeLayer(activeBasePts); });
    chkPolys.addEventListener("change", e => { if (activePolys) e.target.checked ? activePolys.addTo(map) : map.removeLayer(activePolys); });
    
    // Initialize
    loadRegion(currentRegion);
  </script>
</body>
</html>'

html_out <- sub("MAP_DATA_PLACEHOLDER", map_data_json, html_template, fixed = TRUE)
writeLines(html_out, "DataProcessing/outputs/master_sdm_map.html")
cat("Successfully generated Master Leaflet HTML map: DataProcessing/outputs/master_sdm_map.html\n")
