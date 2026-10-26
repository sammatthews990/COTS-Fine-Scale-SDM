# Question: what can the completed concentration maps say about larval dispersal?
# Method: import the 23 October outputs; prepare a portrait version of lc2.
# Validation: pinned inputs (including DBF); day-10 reef means replay the source
# table. Descriptive concentration, not a settlement/connectivity matrix test.
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "..", "R", "common.R"))
cfg <- initialise_final()
suppressPackageStartupMessages({library(sf); library(dplyr); library(tidyr)})
sf_use_s2(FALSE)
keys <- grep("^larval_report_", names(cfg$inputs), value = TRUE)
paths <- setNames(lapply(keys, function(k) checked_input(cfg, k)), keys)
tables <- c(summary = "larval_concentration_reef_summary.csv", zones = "larval_concentration_reef_zones.csv",
            releases = "larval_concentration_releases.csv", latitude = "larval_concentration_latitude.csv")
for (k in names(tables)) write_table(read.csv(paths[[paste0("larval_report_", k)]]), cfg, tables[[k]])
p <- st_read(paths$larval_report_day10_shp, quiet = TRUE)
d <- p |> st_drop_geometry() |> select(indexV2, GBRMPA_ID, area, starts_with("d20")) |>
  pivot_longer(starts_with("d20"), names_to = "release", values_to = "conc") |>
  mutate(w = pmax(area, 1)) |> group_by(release, indexV2, GBRMPA_ID) |>
  summarise(conc = sum(conc * w) / sum(w), area = sum(area), .groups = "drop")
replay <- d |> mutate(w = pmax(area, 1)) |> group_by(release, GBRMPA_ID) |>
  summarise(conc = sum(conc * w) / sum(w), .groups = "drop") |>
  group_by(release) |> mutate(pct = percent_rank(conc)) |> ungroup() |>
  group_by(GBRMPA_ID) |> summarise(mean_pct = mean(pct), .groups = "drop")
source_summary <- read.csv(paths$larval_report_summary) |> filter(age == 10)
check <- inner_join(replay, source_summary, by = "GBRMPA_ID", suffix = c("_replay", "_source"))
stopifnot(nrow(check) == nrow(source_summary), max(abs(check$mean_pct_replay - check$mean_pct_source)) < 1e-10)
log_message("Day-10 concentration ranks reproduce the frozen source reef summary.")
site <- d |> group_by(release) |> mutate(pct = percent_rank(conc)) |> ungroup() |>
  group_by(indexV2) |> summarise(top20_share = mean(pct >= .8), .groups = "drop")
geometry <- p |> select(indexV2) |> group_by(indexV2) |> summarise(.groups = "drop") |>
  st_transform(4326) |> left_join(site, by = "indexV2")
features <- st_read(checked_input(cfg, "reef_features"), quiet = TRUE)
bbox <- st_as_sfc(st_bbox(c(xmin = 144.9, xmax = 146.1, ymin = -15.95, ymax = -14.15), crs = 4326))
land <- features |> filter(FEAT_NAME %in% c("Mainland", "Island")) |> st_transform(4326) |> st_make_valid()
land <- suppressWarnings(st_intersection(land, bbox))
saveRDS(list(sites = geometry, land = land), file.path(cfg$dirs$cache, "larval_concentration_report.rds"))
releases <- as.Date(sub("^d", "", unique(d$release)), "%Y%m%d")
write_table(data.frame(sites = nrow(site), reefs = nrow(source_summary), releases = length(releases),
                       first_release = min(releases), last_release = max(releases)), cfg, "larval_concentration_scope.csv")
write_table(data.frame(input = keys, path = unlist(paths), sha256 = vapply(paths, sha256_file, "")),
            cfg, "larval_concentration_report_inputs.csv")
