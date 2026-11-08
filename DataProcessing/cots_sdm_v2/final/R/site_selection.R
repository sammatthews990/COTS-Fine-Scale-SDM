# Ranking designated cull sites and scoring the ranking against dive outcomes
# (used by 19 and 20; 15-17 carry their own copies).

# Percentile within a group, ties broken at random (0 = lowest, 1 = highest).
pct_rank <- function(v) {
  ok <- is.finite(v)
  out <- rep(NA_real_, length(v))
  out[ok] <- (rank(v[ok], ties.method = "random") - 1) / max(sum(ok) - 1, 1)
  out
}

# Model terms for the site-evidence model. `habitat` names the habitat rank
# column; `habitat2` optionally adds a second habitat rank.
site_evidence_terms <- function(f, habitat, habitat2 = NULL) {
  fill <- function(x) ifelse(is.finite(x), x, 0.5)
  out <- data.frame(habitat = fill(f[[habitat]]),
                    site_problematic_prev = f$site_problematic_prev, site_quiet_prev = f$site_quiet_prev,
                    cull_rate = f$cull_rate, log_cull_effort = log1p(f$cull_effort),
                    log_manta_strong = log1p(f$manta_strong), log_manta_effort = log1p(f$manta_effort))
  if (!is.null(habitat2)) out$habitat2 <- fill(f[[habitat2]])
  out
}
site_evidence_formula <- function(habitat = TRUE, habitat2 = FALSE) {
  rhs <- c(if (habitat) "habitat", if (habitat2) "habitat2", "site_problematic_prev", "site_quiet_prev",
           "cull_rate", "log_cull_effort", "log_manta_strong", "log_manta_effort")
  stats::as.formula(paste("problematic ~", paste(rhs, collapse = " + ")))
}

# Share of problematic dives kept, and of low-density and ghost dives avoided,
# when only sites ranked in the top fraction `a` of their reef are dived.
site_metrics <- function(d, col, a) {
  kept <- d[[col]] >= 1 - a - 1e-9
  c(problematic_dives_kept = sum(kept & d$outcome == "problematic") / sum(d$outcome == "problematic"),
    low_dives_avoided = sum(!kept & d$outcome == "low") / sum(d$outcome == "low"),
    ghost_dives_avoided = sum(!kept & d$outcome == "ghost") / sum(d$outcome == "ghost"),
    dives_kept = mean(kept))
}

# Sites searched and dives avoided when a ranking keeps a fixed share of
# problematic dives (linear interpolation along its curve over site fractions).
matched_retention <- function(d, col, targets, fractions = round(seq(0.05, 1, 0.05), 2)) {
  m <- vapply(fractions, function(a) site_metrics(d, col, a), numeric(4))
  pk <- c(0, m["problematic_dives_kept", ])
  curves <- list(sites_searched = c(0, fractions), low_dives_avoided = c(1, m["low_dives_avoided", ]),
                 ghost_dives_avoided = c(1, m["ghost_dives_avoided", ]))
  unlist(lapply(targets, function(t) {
    v <- vapply(curves, function(y) stats::approx(pk, y, xout = t, ties = mean)$y, numeric(1))
    stats::setNames(v, paste0(names(curves), "_at_", round(100 * t)))
  }))
}

# The EoTR designated cull-site KML with each polygon coloured by its priority
# class and `fields` (named list of character vectors, one value per row of
# `sites`) added to its data and pop-up. Returns the KML text and, per
# placemark, the EoTR "Cull site ID" and the reef||site key.
priority_kml <- function(cfg, sites, priority, fields, title) {
  kml_path <- list.files(file.path(cfg$dirs$cache, "cull_sites_kmz"), pattern = "\\.kml$", full.names = TRUE)[1]
  kml <- paste(readLines(kml_path, encoding = "UTF-8", warn = FALSE), collapse = "\n")
  parts <- strsplit(kml, "<Placemark>", fixed = TRUE)[[1]]
  head_part <- parts[1]
  pm <- parts[-1]
  data_value <- function(x, name) {
    r <- regexpr(paste0("<Data name=\"", name, "\">\\s*<value>[^<]*</value>"), x)
    out <- rep(NA_character_, length(x))
    out[r > 0] <- sub(".*<value>([^<]*)</value>", "\\1", regmatches(x, r))
    out
  }
  xml_unescape <- function(x) gsub("&amp;", "&", gsub("&lt;", "<", gsub("&gt;", ">", gsub("&apos;", "'", gsub("&quot;", "\"", x)))))
  pm_id <- data_value(pm, "Cull site ID")
  pm_key <- paste(xml_unescape(data_value(pm, "Reef")), clean_site_name(xml_unescape(data_value(pm, "Name"))), sep = "||")
  idx <- match(pm_key, paste(sites$reef_name, clean_site_name(sites$site_name), sep = "||"))
  log_message(sprintf("KMZ: %d placemarks, %.1f%% matched to a ranked site", length(pm), 100 * mean(!is.na(idx))))
  # Priority colours (KML aabbggrr): A red, B orange, C yellow, D pale yellow,
  # E grey (the last class is always grey).
  breaks <- cfg$site_ranking$priority_breaks
  labels <- c(names(breaks), LETTERS[length(breaks) + 1L])
  colours <- stats::setNames(c("cc4849e3", "cc3468eb", "9900a1ed", "8099e6f5", "4d818789")[c(seq_along(breaks), 5L)], labels)
  pri <- ifelse(is.na(idx), NA, priority[idx])
  normal_colour <- "<color>33669900</color>"
  pos <- regexpr(normal_colour, pm, fixed = TRUE)
  has <- pos > 0 & !is.na(pri)
  pm[has] <- paste0(substr(pm[has], 1, pos[has] - 1), "<color>", colours[pri[has]], "</color>",
                    substr(pm[has], pos[has] + nchar(normal_colour), nchar(pm[has])))
  new_data <- ifelse(is.na(idx), "", paste0("  ", do.call(paste0, lapply(names(fields), function(f) paste0(
    "<Data name=\"", f, "\">\n            <value>", fields[[f]][idx], "</value>\n          </Data>\n          "))), ""))
  pos <- regexpr("</ExtendedData>", pm, fixed = TRUE)
  has <- pos > 0 & nzchar(new_data)
  pm[has] <- paste0(substr(pm[has], 1, pos[has] - 1), new_data[has], substr(pm[has], pos[has], nchar(pm[has])))
  out_kml <- paste0(head_part, paste0("<Placemark>", pm, collapse = ""))
  balloon_rows <- paste0(vapply(names(fields), function(f)
    paste0("<tr>\n<td class='property-name'>", f, "</td>\n<td class='property-value'>$[", f, "]</td>\n</tr>\n"), character(1)), collapse = "")
  out_kml <- gsub("</tbody>", paste0(balloon_rows, "</tbody>"), out_kml, fixed = TRUE)
  out_kml <- sub("<name>Eotr_CotsCullSites_[^<]*</name>", paste0("<name>", title, "</name>"), out_kml)
  list(kml = out_kml, pm_id = pm_id, pm_key = pm_key)
}

# Stage a KML as cache/<stage>/doc.kml (21 adds exploration zones to it) and
# zip it to `kmz_path`.
write_staged_kmz <- function(cfg, kml, stage, kmz_path) {
  stage_dir <- file.path(cfg$dirs$cache, stage)
  dir.create(stage_dir, showWarnings = FALSE)
  writeLines(kml, file.path(stage_dir, "doc.kml"), useBytes = TRUE)
  if (file.exists(kmz_path)) file.remove(kmz_path)
  utils::zip(kmz_path, file.path(stage_dir, "doc.kml"), flags = "-j9Xq")
  log_message("Wrote ", kmz_path)
}
