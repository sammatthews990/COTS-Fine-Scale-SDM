# Ranked designated cull sites for the next season (2027).
#
# Scores every EoTR designated cull site with the site-evidence model: a
# logistic model of "a dive here is problematic (CPUE >= 0.02)" on the habitat
# rank(s) and the previous season's evidence (dives at and near the site,
# strong and total manta tows nearby). Model choice and training window come
# from cfg$site_ranking, set from the forward test in 19.
#
# Inputs for 2027:
#   evidence  the most recent 12 months of Cull dives and manta tows
#             (cfg$site_ranking$evidence_end back 365 days). Calendar 2026
#             holds only January-May, and the model was trained on a full
#             season of previous-season evidence.
#   habitat   deployed layers, mean over each site polygon, ranked among the
#             designated sites on the reef: cull-trained RG + year (clean
#             fallback) and the within reef-year manta model (18).
# Outputs (outputs/site_ranking/):
#   cots_cull_site_priority_2027.csv   one row per designated site
#   cots_cull_site_priority_2027.kmz   the EoTR cull-site KMZ with each polygon
#                                      coloured by priority and the ranking
#                                      added to its data and pop-up
#   cots_cull_site_priority_2027_reefs.csv   per-reef summary
#   cots_cull_site_priority_classes_2027.csv what each class retained in the
#                                      forward test (19)

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "..", "R", "common.R"))
cfg <- initialise_final()
source(file.path(cfg$final_dir, "R", "forward.R"))
source(file.path(cfg$final_dir, "R", "cull_sites.R"))
source(file.path(cfg$final_dir, "R", "site_selection.R"))
sr <- cfg$site_ranking
set.seed(cfg$hmsc$seed + 20L)
strong_classes <- c("direct_cots", "multiple_counted_scars", "legacy_common_scars")
key <- c("reef_name", "year", "site_row")
out_dir <- file.path(cfg$dirs$outputs, "site_ranking")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# --- Designated sites and their habitat ranks -----------------------------------
sites_sf <- read_designated_sites(cfg)
sites <- sf::st_drop_geometry(sites_sf)
sites$site_row <- seq_len(nrow(sites))
site_xy <- sf::st_coordinates(suppressWarnings(sf::st_point_on_surface(sf::st_geometry(sites_sf))))[, 1:2]
v <- terra::vect(sites_sf)
cull_layer <- terra::rast(file.path(cfg$dirs$layers, "cots_legacy_rg_year2026_clean_fallback_30m.tif"))
manta_layer <- terra::rast(file.path(cfg$dirs$layers, "cots_manta_within_30m.tif"))[["relative_log_odds"]]
polygon_mean <- function(r) terra::extract(r, v, fun = mean, na.rm = TRUE, touches = TRUE, ID = FALSE)[[1]]
sites$habitat_cull <- polygon_mean(cull_layer)
sites$habitat_manta <- polygon_mean(manta_layer)
sites <- sites |>
  dplyr::group_by(.data$reef_name) |>
  dplyr::mutate(rank_cull = pct_rank(.data$habitat_cull), rank_within = pct_rank(.data$habitat_manta),
                rank_ensemble = pct_rank((.data$rank_cull + .data$rank_within) / 2)) |>
  dplyr::ungroup()

# --- Previous-season evidence (most recent 12 months) ---------------------------
evidence_start <- sr$evidence_end - sr$evidence_days + 1
dives <- matched_cull_dives(cfg, sites_sf)
has_xy <- is.finite(dives$longitude) & is.finite(dives$latitude)
xy <- matrix(NA_real_, nrow(dives), 2)
xy[has_xy, ] <- sf::st_coordinates(sf::st_transform(sf::st_as_sf(dives[has_xy, ], coords = c("longitude", "latitude"), crs = 4326), sf::st_crs(sites_sf)))
fill <- !has_xy & !is.na(dives$site_row)
xy[fill, ] <- site_xy[dives$site_row[fill], ]
dives$x <- xy[, 1]
dives$y <- xy[, 2]
recent <- dives[dives$survey_date >= evidence_start & dives$survey_date <= sr$evidence_end & is.finite(dives$x), ]

tows <- readRDS(checked_input(cfg, "manta_tows")) |>
  dplyr::select("visit_id", "survey_date", "reef_name", "manta_evidence_class")
geom <- readRDS(checked_input(cfg, "footprint_links"))$manta |>
  dplyr::group_by(.data$visit_id) |>
  dplyr::summarise(x = (dplyr::first(.data$start_x) + dplyr::first(.data$end_x)) / 2,
                   y = (dplyr::first(.data$start_y) + dplyr::first(.data$end_y)) / 2, .groups = "drop")
recent_tows <- dplyr::inner_join(tows, geom, by = "visit_id") |>
  dplyr::filter(is.finite(.data$x), .data$survey_date >= evidence_start, .data$survey_date <= sr$evidence_end) |>
  dplyr::mutate(strong = .data$manta_evidence_class %in% strong_classes)
log_message(sprintf("Evidence %s to %s: %s dives, %s manta tows (%s strong)", evidence_start, sr$evidence_end,
                    format(nrow(recent), big.mark = ","), format(nrow(recent_tows), big.mark = ","), sum(recent_tows$strong)))

ck <- history_kernel(cbind(recent$x, recent$y), as.numeric(recent$outcome == "problematic"), rep(0, nrow(recent)),
                     site_xy, sr$kernel_sigma_m, Inf)
mk <- history_kernel(cbind(recent_tows$x, recent_tows$y), as.numeric(recent_tows$strong), rep(0, nrow(recent_tows)),
                     site_xy, sr$kernel_sigma_m, Inf)
at_site <- recent[!is.na(recent$site_row), ] |>
  dplyr::group_by(.data$site_row) |>
  dplyr::summarise(dives_last_season = dplyr::n(), problematic_last_season = sum(.data$outcome == "problematic"),
                   max_cpue_last_season = max(.data$cpue), .groups = "drop")
sites <- sites |>
  dplyr::left_join(at_site, by = "site_row") |>
  dplyr::mutate(dives_last_season = dplyr::coalesce(.data$dives_last_season, 0L),
                problematic_last_season = dplyr::coalesce(.data$problematic_last_season, 0L),
                site_problematic_prev = as.numeric(.data$problematic_last_season > 0),
                site_quiet_prev = as.numeric(.data$dives_last_season > 0 & .data$problematic_last_season == 0),
                cull_rate = ifelse(ck$W > 0, ck$H, 0), cull_effort = ck$W,
                manta_strong = ifelse(mk$W > 0, mk$H * mk$W, 0), manta_effort = mk$W,
                reef_dived_last_season = .data$reef_name %in% recent$reef_name,
                reef_towed_last_season = .data$reef_name %in% recent_tows$reef_name)

# --- Fit the site-evidence model on the training seasons ------------------------
habitat <- switch(sr$habitat, cull = list(h = "rank_cull"), ensemble = list(h = "rank_ensemble"),
                  both = list(h = "rank_cull", h2 = "rank_within"))
ev <- readRDS(file.path(cfg$dirs$cache, "cull_site_evidence_scores.rds")) |>
  dplyr::filter(.data$window == "last_season") |>
  dplyr::select(dplyr::all_of(key), "site_problematic_prev", "site_quiet_prev", "cull_rate", "cull_effort",
                "manta_strong", "manta_effort") |>
  dplyr::inner_join(readRDS(file.path(cfg$dirs$cache, "manta_response_site_ranks.rds")) |>
                      dplyr::select(dplyr::all_of(key), "rank_cull", "rank_within", "rank_ensemble"), by = key)
train <- dives[!is.na(dives$site_row) & dives$year >= sr$training_years[1] & dives$year <= sr$training_years[2], c(key, "outcome")] |>
  dplyr::inner_join(ev, by = key)
terms <- site_evidence_terms(train, habitat$h, habitat$h2)
terms$problematic <- as.numeric(train$outcome == "problematic")
model <- stats::glm(site_evidence_formula(habitat2 = !is.null(habitat$h2)), family = stats::binomial(), data = terms)
saveRDS(model, file.path(cfg$dirs$models, "site_ranking_2027_glm.rds"))
coef_table <- tibble::tibble(term = names(stats::coef(model)), estimate = unname(stats::coef(model)),
                             std_error = unname(sqrt(diag(stats::vcov(model)))),
                             training_dives = nrow(train), training_years = paste(sr$training_years, collapse = "-"),
                             habitat = sr$habitat)
write_table(coef_table, cfg, "site_ranking_2027_coefficients.csv")
log_message("Model: ", sr$habitat, " habitat + evidence, ", format(nrow(train), big.mark = ","), " dives ",
            paste(sr$training_years, collapse = "-"))

# --- Score and rank every designated site ----------------------------------------
sites$score <- stats::predict(model, site_evidence_terms(sites, habitat$h, habitat$h2))
breaks <- sr$priority_breaks
sites <- sites |>
  dplyr::group_by(.data$reef_name) |>
  dplyr::arrange(dplyr::desc(.data$score), dplyr::desc(.data$rank_cull), .by_group = TRUE) |>
  dplyr::mutate(sites_on_reef = dplyr::n(), rank_on_reef = dplyr::row_number(),
                top_fraction = .data$rank_on_reef / .data$sites_on_reef,
                priority = dplyr::case_when(.data$top_fraction <= breaks[[1]] ~ names(breaks)[1],
                                            .data$top_fraction <= breaks[[2]] ~ names(breaks)[2],
                                            .data$top_fraction <= breaks[[3]] ~ names(breaks)[3],
                                            TRUE ~ "D")) |>
  dplyr::ungroup() |>
  dplyr::mutate(p_problematic_if_dived = stats::plogis(.data$score),
                evidence = dplyr::case_when(
                  .data$site_problematic_prev == 1 ~ "problematic dive here last season",
                  .data$manta_strong >= 0.5 ~ "strong manta evidence nearby last season",
                  .data$site_quiet_prev == 1 ~ "dived last season, nothing problematic",
                  .data$cull_effort > 0.1 | .data$manta_effort > 0.5 ~ "searched nearby last season, little found",
                  .data$reef_dived_last_season | .data$reef_towed_last_season ~ "reef visited, nothing nearby",
                  TRUE ~ "no recent evidence on this reef (habitat only)"))

table_out <- sites |>
  dplyr::transmute(reef_name = .data$reef_name, site_name = .data$site_name, priority = .data$priority,
                   rank_on_reef = .data$rank_on_reef, sites_on_reef = .data$sites_on_reef,
                   top_fraction = round(.data$top_fraction, 3), score_log_odds = round(.data$score, 3),
                   p_problematic_if_dived = round(.data$p_problematic_if_dived, 3), evidence = .data$evidence,
                   dives_last_season = .data$dives_last_season, problematic_dives_last_season = .data$problematic_last_season,
                   max_cpue_last_season = round(.data$max_cpue_last_season, 4),
                   nearby_problematic_rate = round(.data$cull_rate, 3), nearby_dive_effort = round(.data$cull_effort, 2),
                   nearby_strong_manta = round(.data$manta_strong, 2), nearby_manta_effort = round(.data$manta_effort, 2),
                   habitat_rank_cull = round(.data$rank_cull, 3), habitat_rank_manta = round(.data$rank_within, 3),
                   reef_dived_last_season = .data$reef_dived_last_season, reef_towed_last_season = .data$reef_towed_last_season) |>
  dplyr::arrange(.data$reef_name, .data$rank_on_reef)
csv_path <- file.path(out_dir, "cots_cull_site_priority_2027.csv")
utils::write.csv(table_out, csv_path, row.names = FALSE, na = "")
log_message("Wrote ", csv_path)

reefs <- table_out |>
  dplyr::group_by(.data$reef_name) |>
  dplyr::summarise(sites = dplyr::n(), priority_a = sum(.data$priority == "A"),
                   flagged_last_season = sum(.data$evidence %in% c("problematic dive here last season", "strong manta evidence nearby last season")),
                   dived_last_season = dplyr::first(.data$reef_dived_last_season),
                   towed_last_season = dplyr::first(.data$reef_towed_last_season), .groups = "drop")
utils::write.csv(reefs, file.path(out_dir, "cots_cull_site_priority_2027_reefs.csv"), row.names = FALSE)

# What each class kept in the forward test (19), for the chosen score.
fw <- readRDS(file.path(cfg$dirs$cache, "site_ranking_forward.rds"))
col <- paste(sr$habitat, sr$forward_window, sep = "_")
cum <- c(breaks, D = 1)
classes <- tibble::tibble(
  priority = names(cum), searched_through = c(sprintf("top %d%% of the reef's sites", round(100 * breaks)), "all sites"),
  problematic_dives_kept = vapply(cum, function(a) site_metrics(fw, col, a)[["problematic_dives_kept"]], numeric(1)),
  low_dives_avoided = vapply(cum, function(a) site_metrics(fw, col, a)[["low_dives_avoided"]], numeric(1)),
  ghost_dives_avoided = vapply(cum, function(a) site_metrics(fw, col, a)[["ghost_dives_avoided"]], numeric(1)),
  forward_test = paste0("dives 2019-2026, score ", col))
utils::write.csv(classes, file.path(out_dir, "cots_cull_site_priority_classes_2027.csv"), row.names = FALSE)

# --- KMZ: the EoTR file with priority colours and fields added ------------------------
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
pm_name <- data_value(pm, "Name")
pm_reef <- data_value(pm, "Reef")
pm_id <- data_value(pm, "Cull site ID")
xml_unescape <- function(x) gsub("&amp;", "&", gsub("&lt;", "<", gsub("&gt;", ">", gsub("&apos;", "'", gsub("&quot;", "\"", x)))))
pm_key <- paste(xml_unescape(pm_reef), clean_site_name(xml_unescape(pm_name)), sep = "||")
idx <- match(pm_key, paste(sites$reef_name, clean_site_name(sites$site_name), sep = "||"))
log_message(sprintf("KMZ: %d placemarks, %.1f%% matched to a ranked site", length(pm), 100 * mean(!is.na(idx))))
# Priority colours (KML aabbggrr): A red, B orange, C yellow, D grey.
colours <- c(A = "cc4849e3", B = "cc3468eb", C = "9900a1ed", D = "4d818789")
pri <- ifelse(is.na(idx), NA, sites$priority[idx])
normal_colour <- "<color>33669900</color>"
pos <- regexpr(normal_colour, pm, fixed = TRUE)
has <- pos > 0 & !is.na(pri)
pm[has] <- paste0(substr(pm[has], 1, pos[has] - 1), "<color>", colours[pri[has]], "</color>",
                  substr(pm[has], pos[has] + nchar(normal_colour), nchar(pm[has])))
fmt <- function(x, d) ifelse(is.finite(x), formatC(x, format = "f", digits = d), "")
new_data <- ifelse(is.na(idx), "", paste0(
  "  <Data name=\"Priority 2027\">\n            <value>", pri, "</value>\n          </Data>\n",
  "          <Data name=\"Rank on reef\">\n            <value>", sites$rank_on_reef[idx], " of ", sites$sites_on_reef[idx], "</value>\n          </Data>\n",
  "          <Data name=\"P(problematic if dived)\">\n            <value>", fmt(sites$p_problematic_if_dived[idx], 2), "</value>\n          </Data>\n",
  "          <Data name=\"Last season\">\n            <value>", sites$evidence[idx], "</value>\n          </Data>\n        "))
pos <- regexpr("</ExtendedData>", pm, fixed = TRUE)
has <- pos > 0 & nzchar(new_data)
pm[has] <- paste0(substr(pm[has], 1, pos[has] - 1), new_data[has], substr(pm[has], pos[has], nchar(pm[has])))
out_kml <- paste0(head_part, paste0("<Placemark>", pm, collapse = ""))
balloon_rows <- paste0(vapply(c("Priority 2027", "Rank on reef", "P(problematic if dived)", "Last season"), function(f)
  paste0("<tr>\n<td class='property-name'>", f, "</td>\n<td class='property-value'>$[", f, "]</td>\n</tr>\n"), character(1)), collapse = "")
out_kml <- gsub("</tbody>", paste0(balloon_rows, "</tbody>"), out_kml, fixed = TRUE)
out_kml <- sub("<name>Eotr_CotsCullSites_[^<]*</name>",
               paste0("<name>COTS cull site priority ", sr$target_season, " (A red, B orange, C yellow, D grey)</name>"), out_kml)
stage <- file.path(cfg$dirs$cache, "site_ranking_kmz")
dir.create(stage, showWarnings = FALSE)
writeLines(out_kml, file.path(stage, "doc.kml"), useBytes = TRUE)
kmz_path <- file.path(out_dir, "cots_cull_site_priority_2027.kmz")
if (file.exists(kmz_path)) file.remove(kmz_path)
utils::zip(kmz_path, file.path(stage, "doc.kml"), flags = "-j9Xq")
log_message("Wrote ", kmz_path)

# Add the EoTR "Cull site ID" to the table (first placemark per site).
table_key <- paste(table_out$reef_name, clean_site_name(table_out$site_name), sep = "||")
table_out <- dplyr::bind_cols(table_out[, c("reef_name", "site_name")],
                              cull_site_id = pm_id[match(table_key, pm_key)],
                              table_out[, setdiff(names(table_out), c("reef_name", "site_name"))])
utils::write.csv(table_out, csv_path, row.names = FALSE, na = "")
log_message("Added Cull site ID to ", csv_path, " (", sum(!is.na(table_out$cull_site_id)), " of ", nrow(table_out), ")")

print(table(table_out$priority))
print(table(table_out$evidence))
print(classes)
print(coef_table)
