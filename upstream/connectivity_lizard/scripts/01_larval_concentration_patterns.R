# 01_larval_concentration_patterns.R
#
# Question: Do modelled COTS larval concentrations, 2022/23 spawning season,
#   show (a) elevated concentration/retention around the Lizard Island reef
#   group relative to other mid-shelf reefs, and (b) relatively low
#   concentration on the outer-barrier Ribbon Reefs?
# Method: Sub-reef polygon daily mean concentrations (fngbr_reef_daily_means_
#   dayNN) for 10 releases (1 Dec 2022 - 28 Jan 2023) at larval ages 5, 9 and
#   10 days. Polygons are aggregated to sub-reef sites (indexV2, the node set of
#   the far-north connectivity matrices) by area-weighted mean, then to reefs
#   (GBRMPA_ID). Concentrations are compared as within-release percentile
#   ranks, so release-to-release differences in total larval load do not
#   dominate. Reefs are assigned to five shelf zones (see assign_zones()).
# Validation: Descriptive. Zone contrasts are tested by permuting zone labels
#   among reefs (ignores spatial autocorrelation, so p-values are optimistic)
#   and by counting releases where the contrast holds. Treat as hypothesis-
#   generating; the network analysis in ANALYSIS_PLAN.md is the formal test.
#
# Run from the repository root:
#   Rscript DataProcessing/connectivity_lizard/scripts/01_larval_concentration_patterns.R

suppressPackageStartupMessages({
  library(sf); library(dplyr); library(tidyr); library(ggplot2)
  library(patchwork); library(ggrepel)
})
sf_use_s2(FALSE)

source("DataProcessing/cots_sdm_v2/final/R/figures.R")  # palette + themes

root <- "DataProcessing/connectivity_lizard"
dir_fig <- file.path(root, "outputs", "figures")
dir_tab <- file.path(root, "outputs", "tables")
dir.create(dir_fig, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_tab, recursive = TRUE, showWarnings = FALSE)

conc_dir <- "DataProcessing/data/fngbr_LarvalConcentrations"
ages <- c(5, 9, 10)
conc_files <- file.path(conc_dir, sprintf("fngbr_reef_daily_means_day%02d.shp", ages))
features_file <- "DataProcessing/data/Great_Barrier_Reef_Features.gpkg"

lizard_xy <- c(145.455, -14.675)  # Lizard Island
lizard_radius_km <- 20
utm <- 32755

write_tab <- function(x, name) {
  path <- file.path(dir_tab, name)
  write.csv(x, path, row.names = FALSE)
  message("Wrote ", path)
}
save_fig <- function(p, name, w, h) {
  path <- file.path(dir_fig, name)
  ggsave(path, p, width = w, height = h, dpi = 300, device = ragg::agg_png, bg = ink$surface)
  message("Wrote ", path)
}

# ---- Inputs --------------------------------------------------------------

shp_parts <- function(f) sub("\\.shp$", c(".shp", ".dbf", ".shx", ".prj"), f)
manifest <- tibble(file = c(unlist(lapply(conc_files, shp_parts)), features_file)) |>
  mutate(sha256 = vapply(file, \(f) digest::digest(file = f, algo = "sha256"), ""))
write_tab(manifest, "input_manifest.csv")

polys <- st_read(conc_files[1], quiet = TRUE) |>
  select(indexV3, indexV2, GBRMPA_ID, area)

read_conc <- function(f, age) {
  st_read(f, quiet = TRUE) |>
    st_drop_geometry() |>
    select(indexV3, indexV2, GBRMPA_ID, area, starts_with("d20")) |>
    pivot_longer(starts_with("d20"), names_to = "release", values_to = "conc") |>
    mutate(release = as.Date(sub("^d", "", release), "%Y%m%d"), age = age)
}
long <- bind_rows(Map(read_conc, conc_files, ages))

features <- st_read(features_file, quiet = TRUE)
reef_names <- features |>
  st_drop_geometry() |>
  filter(!grepl("^U/N", GBR_NAME)) |>
  group_by(LABEL_ID) |>
  summarise(name = first(GBR_NAME), .groups = "drop")

bbox_ll <- st_as_sfc(st_bbox(c(xmin = 144.9, xmax = 146.1, ymin = -15.95, ymax = -14.15), crs = 4326))
land <- features |>
  filter(FEAT_NAME %in% c("Mainland", "Island")) |>
  st_transform(4326) |>
  st_make_valid() |>
  st_intersection(bbox_ll)
mainland_utm <- land |> filter(FEAT_NAME == "Mainland") |> st_transform(utm)

# ---- Sites and reefs ------------------------------------------------------

# Area-weighted mean over the polygons that make up each sub-reef site.
site_conc <- long |>
  mutate(w = pmax(area, 1)) |>
  group_by(age, release, indexV2, GBRMPA_ID) |>
  summarise(conc = sum(conc * w) / sum(w), area = sum(area), .groups = "drop") |>
  group_by(age, release) |>
  mutate(pct = percent_rank(conc)) |>
  ungroup()

reef_geom <- polys |>
  st_transform(utm) |>
  group_by(GBRMPA_ID) |>
  summarise(area = sum(area), .groups = "drop")
reef_cent <- st_centroid(st_geometry(reef_geom))
ll <- st_coordinates(st_transform(reef_cent, 4326))
lizard_utm <- st_transform(st_sfc(st_point(lizard_xy), crs = 4326), utm)

reefs <- reef_geom |>
  st_drop_geometry() |>
  mutate(
    lon = ll[, 1], lat = ll[, 2],
    dist_coast_km = as.numeric(st_distance(reef_cent, st_union(mainland_utm))) / 1000,
    dist_lizard_km = as.numeric(st_distance(reef_cent, lizard_utm)) / 1000
  ) |>
  left_join(reef_names, by = c(GBRMPA_ID = "LABEL_ID")) |>
  mutate(name = coalesce(name, sub("[a-z]$", "", GBRMPA_ID) |>
                           (\(id) reef_names$name[match(id, reef_names$LABEL_ID)])()),
         name = coalesce(name, paste0("U/N (", GBRMPA_ID, ")")))

# Shelf zones.
# - Ribbon Reefs: within 4 km (east-west) of, or seaward of, the line through
#   the named Ribbon No 1-10 / No 10 Patches reefs, south of Lizard.
# - Outer barrier N of Lizard: Yonge, No Name, Carter, Day, Hicks, Sydney
#   Schubert (outer barrier > 44 km offshore, north of -14.67, east of 145.48).
# - Lizard group: other reefs within lizard_radius_km of Lizard Island.
# - Inner shelf: < 18 km from the mainland. Everything else: other mid-shelf.
assign_zones <- function(r) {
  rib <- r |> filter(grepl("^Ribbon No|^No 10 Patches", name))
  rib_line <- approxfun(rib$lat, rib$lon, rule = 2)
  r |>
    mutate(
      ribbon_offset_km = (rib_line(lat) - lon) * 111.32 * cos(lat * pi / 180),
      zone = case_when(
        lat < -14.78 & ribbon_offset_km <= 4 ~ "Ribbon Reefs",
        lat > -14.67 & lon > 145.48 & dist_coast_km > 44 ~ "Outer barrier N of Lizard",
        dist_lizard_km <= lizard_radius_km ~ "Lizard group",
        dist_coast_km < 18 ~ "Inner shelf",
        TRUE ~ "Other mid-shelf"
      )
    )
}
zone_levels <- c("Inner shelf", "Lizard group", "Other mid-shelf",
                 "Ribbon Reefs", "Outer barrier N of Lizard")
zone_colours <- setNames(c("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4"), zone_levels)

reefs <- assign_zones(reefs) |> mutate(zone = factor(zone, zone_levels))
write_tab(reefs |> arrange(lat), "reef_zones.csv")

# Reef-level concentration per release and age (area-weighted over sites),
# ranked among reefs within each release.
reef_conc <- site_conc |>
  mutate(w = pmax(area, 1)) |>
  group_by(age, release, GBRMPA_ID) |>
  summarise(conc = sum(conc * w) / sum(w), .groups = "drop") |>
  group_by(age, release) |>
  mutate(pct = percent_rank(conc)) |>
  ungroup() |>
  left_join(reefs |> select(GBRMPA_ID, zone, lat, lon, name), by = "GBRMPA_ID")

reef_summary <- reef_conc |>
  group_by(age, GBRMPA_ID, zone, lat, lon, name) |>
  summarise(mean_pct = mean(pct), top20_share = mean(pct >= 0.8),
            mean_conc = mean(conc), .groups = "drop")
write_tab(reef_summary, "reef_concentration_summary.csv")

site_summary <- site_conc |>
  group_by(age, indexV2) |>
  summarise(mean_pct = mean(pct), top20_share = mean(pct >= 0.8), .groups = "drop")

# ---- Zone statistics -------------------------------------------------------

zone_release <- reef_conc |>
  group_by(age, release) |>
  mutate(domain_mean = mean(conc)) |>
  group_by(age, release, zone) |>
  summarise(n_reefs = n(), mean_pct = mean(pct),
            rel_conc = mean(conc) / first(domain_mean), .groups = "drop") |>
  mutate(log2_rel = log2(pmax(rel_conc, 1 / 64)))
write_tab(zone_release, "zone_by_release.csv")

perm_contrast <- function(d, a, b, var = "mean_pct", n_perm = 9999) {
  d <- d |> filter(zone %in% c(a, b))
  v <- d[[var]]
  obs <- mean(v[d$zone == a]) - mean(v[d$zone == b])
  na <- sum(d$zone == a)
  set.seed(20221201)
  null <- replicate(n_perm, {
    i <- sample(nrow(d), na)
    mean(v[i]) - mean(v[-i])
  })
  tibble(zone_a = a, zone_b = b, n_a = na, n_b = nrow(d) - na,
         measure = var, diff = obs,
         p_greater = (1 + sum(null >= obs)) / (n_perm + 1),
         p_less = (1 + sum(null <= obs)) / (n_perm + 1))
}
contrasts <- list(c("Lizard group", "Other mid-shelf"),
                  c("Ribbon Reefs", "Other mid-shelf"),
                  c("Ribbon Reefs", "Lizard group"),
                  c("Outer barrier N of Lizard", "Lizard group"))
zone_tests <- reef_summary |>
  group_by(age) |>
  group_modify(\(d, k) bind_rows(lapply(contrasts, \(z) perm_contrast(d, z[1], z[2])))) |>
  ungroup()
release_wins <- zone_release |>
  select(age, release, zone, mean_pct) |>
  pivot_wider(names_from = zone, values_from = mean_pct)
zone_tests$releases_a_gt_b <- mapply(function(ag, a, b) {
  r <- release_wins[release_wins$age == ag, ]
  sum(r[[a]] > r[[b]])
}, zone_tests$age, zone_tests$zone_a, zone_tests$zone_b)
zone_tests$n_releases <- length(unique(zone_release$release))
write_tab(zone_tests, "zone_contrasts.csv")
print(zone_tests, width = 200)

# Concentration rank declines strongly with latitude (north = low), which
# confounds zone contrasts. Remove a smooth latitude trend and re-test zones
# on the residuals, i.e. compare each reef with reefs at the same latitude.
reef_summary <- reef_summary |>
  group_by(age) |>
  mutate(lat_resid = residuals(loess(mean_pct ~ lat, span = 0.5))) |>
  ungroup()
lat_tests <- reef_summary |>
  group_by(age) |>
  group_modify(\(d, k) bind_rows(lapply(contrasts, \(z) perm_contrast(d, z[1], z[2], "lat_resid")))) |>
  ungroup()
write_tab(lat_tests, "zone_contrasts_latitude_adjusted.csv")
print(lat_tests, width = 200)

lat_trend <- reef_summary |>
  group_by(age) |>
  summarise(spearman_rank_vs_lat = cor(mean_pct, lat, method = "spearman"), .groups = "drop")
write_tab(lat_trend, "latitude_trend.csv")
print(lat_trend)

# Do big reefs simply hold more larvae (local release + retention)?
area_check <- reef_summary |>
  left_join(reefs |> select(GBRMPA_ID, area), by = "GBRMPA_ID") |>
  group_by(age) |>
  summarise(spearman_rank_vs_log_area = cor(mean_pct, log(area), method = "spearman"),
            partial_after_lat = cor(lat_resid, log(area), method = "spearman"),
            .groups = "drop")
write_tab(area_check, "area_check.csv")
print(area_check)

# ---- Figures ----------------------------------------------------------------

age_lab <- function(a) paste0("Larval age ", a, " d")
map_limits <- list(x = c(145.1, 145.9), y = c(-15.8, -14.27))
zoom <- list(x = c(145.28, 145.74), y = c(-14.95, -14.45))

site_polys <- polys |>
  group_by(indexV2) |>
  summarise(.groups = "drop")

base_map <- function(xlim, ylim) {
  list(
    geom_sf(data = land, fill = ink$land, colour = NA),
    theme_map(),
    theme(panel.background = element_rect(fill = ink$surface, colour = NA))
  )
}
map_coord <- function(xlim, ylim) coord_sf(xlim = xlim, ylim = ylim, expand = FALSE)
pct_scale <- scale_fill_gradientn(
  colours = ramp_blue, limits = c(0, 1), breaks = c(0, 0.5, 1),
  labels = c("Lowest", "Median", "Highest"),
  name = "Within-release\nconcentration rank", aesthetics = c("fill", "colour")
)
key_labels <- reefs |>
  filter(grepl("^Lizard Island Reef \\(Lagoon|^Ribbon No (1|3|5|7|10) Reef|^Yonge|^Day Reef|^Osterlund|^Williamson|^Mackay|^Jewell|^Hope Islands Reef \\(West", name)) |>
  mutate(label = sub(" Reef.*$| Reefs$| \\(.*$", "", name),
         label = sub("^Lizard Island", "Lizard Is.", label))

# Fig 1: mean rank by larval age (small multiples).
f1_data <- site_polys |>
  inner_join(site_summary, by = "indexV2") |>
  mutate(age = factor(age_lab(age), age_lab(ages)))
fig1 <- ggplot() +
  base_map(map_limits$x, map_limits$y) +
  geom_sf(data = f1_data, aes(fill = mean_pct, colour = mean_pct), linewidth = 0.35) +
  pct_scale +
  map_coord(map_limits$x, map_limits$y) +
  facet_wrap(~age, nrow = 1) +
  labs(title = "Where larvae concentrate changes with larval age",
       subtitle = "Sub-reef sites, mean within-release rank of concentration over 10 releases (Dec 2022 - Jan 2023)",
       caption = "Rank removes differences in total larval load between releases. Grey = land.")
save_fig(fig1, "lc1_concentration_rank_by_age.png", 10, 6.2)

# Fig 2: persistence of day-10 hotspots, domain + Lizard zoom.
f2_data <- site_polys |>
  inner_join(site_summary |> filter(age == 10), by = "indexV2")
persist_scale <- scale_fill_gradientn(
  colours = ramp_orange, limits = c(0, 1), labels = scales::percent,
  name = "Releases in\ntop 20 %", aesthetics = c("fill", "colour")
)
zoom_box <- st_as_sfc(st_bbox(c(xmin = zoom$x[1], xmax = zoom$x[2],
                                ymin = zoom$y[1], ymax = zoom$y[2]), crs = 4326))
f2a <- ggplot() +
  base_map(map_limits$x, map_limits$y) +
  geom_sf(data = f2_data, aes(fill = top20_share, colour = top20_share), linewidth = 0.35) +
  geom_sf(data = zoom_box, fill = NA, colour = ink$secondary, linewidth = 0.3, linetype = 2, inherit.aes = FALSE) +
  geom_text_repel(data = key_labels, aes(lon, lat, label = label), size = 2.4,
                  colour = ink$primary, family = report_font, min.segment.length = 0,
                  segment.colour = ink$muted, segment.size = 0.2, box.padding = 0.3, seed = 1) +
  persist_scale +
  map_coord(map_limits$x, map_limits$y) +
  labs(title = "Domain")
zoom_labels <- reefs |>
  filter(lon > zoom$x[1], lon < zoom$x[2], lat > zoom$y[1], lat < zoom$y[2],
         !grepl("^U/N", name)) |>
  group_by(sub("[a-z]$", "", GBRMPA_ID)) |>
  slice_max(area, n = 1) |>
  ungroup() |>
  mutate(label = sub(" \\(.*$", "", name))
f2b <- ggplot() +
  base_map(zoom$x, zoom$y) +
  geom_sf(data = f2_data, aes(fill = top20_share, colour = top20_share), linewidth = 0.25) +
  geom_text_repel(data = zoom_labels, aes(lon, lat, label = label), size = 2.4,
                  colour = ink$primary, family = report_font, min.segment.length = 0,
                  segment.colour = ink$muted, segment.size = 0.2, box.padding = 0.35, seed = 2) +
  persist_scale +
  map_coord(zoom$x, zoom$y) +
  labs(title = "Lizard Island group")
fig2 <- (f2a | f2b) + plot_layout(widths = c(1, 1.25), guides = "collect") +
  plot_annotation(
    title = "Persistent larval hotspots at competency age (day 10)",
    subtitle = "Share of the 10 releases in which a sub-reef site ranks in the domain's top 20 % for concentration",
    caption = "Model: far-north GBR sub-reef larval concentration, 2022/23 spawning season.",
    theme = theme_report()
  )
save_fig(fig2, "lc2_day10_hotspot_persistence.png", 10, 7)

# Fig 3: reef-level zone comparison.
f3_data <- reef_summary |> mutate(age = factor(age_lab(age), age_lab(ages)))
f3_med <- f3_data |> group_by(age, zone) |>
  summarise(med = median(mean_pct), n = n(), .groups = "drop")
fig3 <- ggplot(f3_data, aes(zone, mean_pct)) +
  geom_hline(yintercept = 0.5, colour = ink$axis, linewidth = 0.3) +
  geom_point(aes(colour = zone), position = position_jitter(width = 0.18, seed = 3),
             size = 1.8, alpha = 0.85) +
  geom_crossbar(data = f3_med, aes(y = med, ymin = med, ymax = med),
                width = 0.55, colour = ink$primary, linewidth = 0.4) +
  geom_text(data = f3_med, aes(y = 1.02, label = paste0("n=", n)), size = 2.4,
            colour = ink$muted, family = report_font) +
  scale_colour_manual(values = zone_colours, guide = "none") +
  scale_x_discrete(labels = \(x) sub(" N of ", "\nN of ", sub(" mid-shelf", "\nmid-shelf", x))) +
  scale_y_continuous(limits = c(0, 1.05), breaks = c(0, 0.5, 1),
                     labels = c("Lowest", "Median", "Highest")) +
  facet_wrap(~age, ncol = 1) +
  labs(x = NULL, y = "Mean within-release rank among reefs",
       title = "Concentration by shelf zone",
       subtitle = "Each point is a reef (mean of 10 releases); bar = zone median",
       caption = paste0("Zones: Ribbon Reefs within 4 km of the Ribbon No 1-10 line; Lizard group within 20 km of Lizard Island; inner shelf < 18 km from the mainland.\n",
                        "Not latitude-adjusted: the low 'other mid-shelf' reefs lie north of Lizard (Jewell, Hilder, Parke). See lc6 for the same-latitude comparison.")) +
  theme_report() +
  theme(panel.grid.major.x = element_blank())
save_fig(fig3, "lc3_zone_comparison.png", 7.5, 8)

# Fig 4: release x zone relative concentration (episodic clustering).
f4_data <- zone_release |>
  mutate(age = factor(age_lab(age), age_lab(ages)),
         zone = factor(zone, rev(zone_levels)),
         rel_lab = ifelse(rel_conc < 0.1, "<0.1", formatC(rel_conc, digits = 1, format = "f")))
fig4 <- ggplot(f4_data, aes(release, zone, fill = log2_rel)) +
  geom_tile(colour = ink$surface, linewidth = 0.8) +
  geom_text(aes(label = rel_lab), size = 2.3, family = report_font,
            colour = ifelse(abs(f4_data$log2_rel) > 1.6, "white", ink$primary)) +
  scale_fill_gradientn(colours = rev(diverging_blue_red), limits = c(-3, 3),
                       oob = scales::squish, breaks = -3:3,
                       labels = c("1/8", "1/4", "1/2", "1", "2", "4", "8"),
                       name = "Zone mean /\ndomain mean") +
  scale_x_date(breaks = sort(unique(zone_release$release)), date_labels = "%d %b") +
  facet_wrap(~age, ncol = 1) +
  labs(x = "Release date", y = NULL,
       title = "The Lizard group retains larvae in some releases, rarely to competency age",
       subtitle = "Mean reef concentration in each zone relative to the domain mean, by release",
       caption = "Values > 1 mean the zone holds more larvae per unit reef area than the domain average for that release.") +
  theme_report() +
  theme(panel.grid.major = element_blank(), axis.text.x = element_text(angle = 0, size = 7))
save_fig(fig4, "lc4_zone_by_release.png", 10, 7.5)

# Fig 5: latitudinal profile.
f5_data <- reef_summary |> mutate(age = factor(age_lab(age), age_lab(ages)))
fig5 <- ggplot(f5_data, aes(mean_pct, lat, colour = zone)) +
  annotate("rect", xmin = -Inf, xmax = Inf,
           ymin = lizard_xy[2] - lizard_radius_km / 111, ymax = lizard_xy[2] + lizard_radius_km / 111,
           fill = ink$grid, alpha = 0.5) +
  annotate("text", x = 0.02, y = lizard_xy[2] + 0.03, label = "Lizard Is.", hjust = 0,
           size = 2.5, colour = ink$secondary, family = report_font) +
  geom_point(size = 1.6, alpha = 0.9) +
  scale_colour_manual(values = zone_colours, name = NULL) +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1),
                     labels = c("Lowest", "Median", "Highest")) +
  facet_wrap(~age, nrow = 1) +
  labs(x = "Mean within-release rank among reefs", y = "Latitude (°)",
       title = "North-south profile of larval concentration",
       subtitle = "Each point is a reef; shaded band = 20 km around Lizard Island",
       caption = "Domain spans -14.30 to -15.76°. Reefs near either edge lack upstream sources beyond the domain, so edge values may be biased low.") +
  theme_report() +
  theme(legend.position = "bottom")
save_fig(fig5, "lc5_latitude_profile.png", 10, 6.5)

# Fig 6: cross-shelf profiles within latitude bands (holds latitude constant).
band_breaks <- c(-15.8, -15.35, -14.9, -14.6, -14.25)
band_labels <- c("South: Cooktown to Ribbon 4", "Central: Ribbon 5 to 9",
                 "Lizard: Ribbon 10 to Lizard", "North of Lizard")
f6_data <- reef_summary |>
  left_join(reefs |> select(GBRMPA_ID, dist_coast_km), by = "GBRMPA_ID") |>
  mutate(band = factor(cut(lat, band_breaks, labels = band_labels), rev(band_labels)),
         age = factor(age_lab(age), age_lab(ages)))
fig6 <- ggplot(f6_data, aes(dist_coast_km, mean_pct, colour = zone)) +
  geom_hline(yintercept = 0.5, colour = ink$axis, linewidth = 0.3) +
  geom_point(size = 1.7, alpha = 0.9) +
  scale_colour_manual(values = zone_colours, name = NULL) +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1),
                     labels = c("Lowest", "Median", "Highest")) +
  facet_grid(band ~ age, labeller = labeller(band = label_wrap_gen(14))) +
  labs(x = "Distance from mainland (km)", y = "Mean within-release rank among reefs",
       title = "Across the shelf at the same latitude",
       subtitle = "Each point is a reef; rows are latitude bands, north at top",
       caption = "Compares inner, mid-shelf and outer-barrier reefs without the north-south gradient.") +
  theme_report() +
  theme(legend.position = "bottom", strip.text.y = element_text(angle = 0, hjust = 0))
save_fig(fig6, "lc6_cross_shelf_by_latitude.png", 10, 9)

message("Done.")
