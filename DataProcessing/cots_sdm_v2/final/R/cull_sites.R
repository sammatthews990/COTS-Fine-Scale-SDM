# Designated cull sites and Cull dives matched to them (as built by 15).

clean_site_name <- function(x) toupper(gsub("^[[:space:](]+|[[:space:]]+$", "", x))

# Designated sites in the same row order as 15 (cache written by 15).
read_designated_sites <- function(cfg) {
  path <- file.path(cfg$dirs$cache, "cull_sites_designated.gpkg")
  if (!file.exists(path)) stop("Run 15_cull_search_area_cv.R first (needs ", path, ")")
  sites_sf <- sf::st_read(path, quiet = TRUE)
  sites_sf$site_key <- paste(sites_sf$reef_name, clean_site_name(sites_sf$site_name), sep = "||")
  sites_sf[!duplicated(sites_sf$site_key), ]
}

# Cull dives with outcome class and the designated site they were made at
# (site name, else the nearest site on the same reef within 100 m).
matched_cull_dives <- function(cfg, sites_sf) {
  threshold <- cfg$legacy$cpue_threshold
  dives <- readRDS(file.path(cfg$exploration_dir, "artifacts", "cull_observations_clean.rds")) |>
    dplyr::filter(is.finite(.data$total_cots), is.finite(.data$bottom_time), .data$bottom_time > 0) |>
    dplyr::mutate(cpue = .data$total_cots / .data$bottom_time,
                  outcome = dplyr::case_when(.data$total_cots == 0 ~ "ghost",
                                             .data$cpue < threshold ~ "low",
                                             TRUE ~ "problematic"),
                  site_key = paste(.data$reef_name, clean_site_name(.data$site_name), sep = "||"))
  dives$site_row <- match(dives$site_key, sites_sf$site_key)
  unmatched <- which(is.na(dives$site_row) & is.finite(dives$longitude) & is.finite(dives$latitude))
  if (length(unmatched)) {
    dpts <- sf::st_transform(sf::st_as_sf(dives[unmatched, ], coords = c("longitude", "latitude"), crs = 4326),
                             sf::st_crs(sites_sf))
    near <- sf::st_nearest_feature(dpts, sites_sf)
    dist <- as.numeric(sf::st_distance(dpts, sites_sf[near, ], by_element = TRUE))
    ok <- dist <= 100 & sites_sf$reef_name[near] == dives$reef_name[unmatched]
    dives$site_row[unmatched[ok]] <- near[ok]
  }
  dives
}
