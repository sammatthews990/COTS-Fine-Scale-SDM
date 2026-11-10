# Density-aware v2 priority classes for the 2027 cull-site ranking.
#
# Question: how should the 2027 product look under the layer that did best in
# 36's replay (2020-2026, matched effort)? The density-aware classes (20b)
# fixed how many sites each reef gets from a 12-month density fitted on dived
# reef-years; reefs culled down months ago (Big Broadhurst) or towed clean and
# left (John Brewer) kept the level of reefs crews chose to dive, and a dense
# reef's every site rose together.
#
# Method (fitted and calibrated in 36, applied to evidence to cfg evidence_end):
#   reefs with dives or manta tows in the year to evidence_end:
#     site + visit model, the sum of two logistic models' log odds:
#       site model    P(dive problematic) on the reef's recent two-part
#                     density, the site's logit rank on the reef (20) and its
#                     own history (dives within 300 m and tows linked to it:
#                     recent, 90-day half-life, and 1-6 years before), fitted
#                     on dives 2019-2026
#       visit model   x 0.5: the same terms plus log counts, fitted on every
#                     designated site of each 2019-2026 cull visit, outcome "a
#                     problematic dive here on this visit" (sites crews towed
#                     past and never dived count as zeros; this encodes where
#                     crews chose to dive as well as COTS)
#     untowed sites (no linked tow in the year) get the mean recent strong-
#     tow rate of the reef's towed sites, as if towed with the reef's average
#     result (10 November): scored with no tow data they had a 50 % strong
#     rate and outranked towed ground with strong scars (Tongue, John Brewer)
#     classes from cut-offs set on the 2019-2026 visits' designated sites so
#     that the pooled shares match the within-reef classes (A 20, B 50, C 66,
#     D 80 %)
#   reefs with neither (known only from the GBR outbreak forecast): the
#     recent two-part density's scaled quota, sites ordered by the site +
#     visit score. Reef density there is P(survey needs culling | forecast) x
#     P(problematic | dived); the site models would give these reefs about
#     60 % A-C against 27 % under the quota.
#   recent two-part density: P(a survey occasion needs culling) x the share
#     of problematic dives if dived, from the reef's dives and tows weighted by
#     a 90-day half-life and the outbreak forecast for the year before.
#
# Validation: 36 (strictly forward replay of 2020-2026 visits; paired
# bootstrap of reef-years against the within-reef and density-aware classes).
# The weight 0.5 was chosen on those test seasons (grid 0.5-5), so its replay
# margin is optimistic. Never-dived reefs cannot be checked by the replay.
#
# Outputs (outputs/site_ranking/):
#   cots_cull_site_priority_density_v2_2027.csv        one row per designated site
#   cots_cull_site_priority_density_v2_2027_reefs.csv  per-reef density and classes
#   cots_cull_site_priority_density_v2_2027.kmz        EoTR KMZ coloured by these
#                                                      classes (21 adds exploration)

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "..", "R", "common.R"))
cfg <- initialise_final()
source(file.path(cfg$final_dir, "R", "cull_sites.R"))
source(file.path(cfg$final_dir, "R", "site_selection.R"))
sr <- cfg$site_ranking
layer <- "density-aware v2, untowed neutral"
density_layer <- "recent, two-part"
out_dir <- file.path(cfg$dirs$outputs, "site_ranking")
labels <- c(names(sr$priority_breaks), LETTERS[length(sr$priority_breaks) + 1L])

layers <- utils::read.csv(file.path(cfg$dirs$tables, "density_layers_2027_sites.csv"), check.names = FALSE,
                          stringsAsFactors = FALSE)
if (!layer %in% names(layers)) stop("Run 36_protocol_density.R first (needs the '", layer, "' column)")
reef_dens <- utils::read.csv(file.path(cfg$dirs$tables, "density_methods_2027_reefs.csv"), check.names = FALSE,
                             stringsAsFactors = FALSE)
v1 <- utils::read.csv(file.path(out_dir, "cots_cull_site_priority_density_2027.csv"), stringsAsFactors = FALSE)
key <- function(d) paste(d$reef_name, d$site_name, sep = "||")
i <- match(key(v1), key(layers))
stopifnot(!anyNA(i))
sites <- v1 |>
  dplyr::mutate(priority_density_v2 = layers[[layer]][i],
                class_basis = ifelse(layers$reef_recent_evidence[i], "site + visit model (reef has recent dives or tows)",
                                     "reef quota (forecast only), site + visit order"),
                site_visit_rank = layers$site_visit_rank_neutral[i], site_visit_score = layers$site_visit_score_neutral[i],
                untowed_12m = layers$untowed_12m[i],
                searched_clean_12m = layers$searched_clean_12m[i],
                reef_density_v2 = round(reef_dens[[paste0(density_layer, ".density")]][match(.data$reef_name, reef_dens$reef_name)], 3))
table_out <- sites |>
  dplyr::transmute(reef_name = .data$reef_name, site_name = .data$site_name, cull_site_id = .data$cull_site_id,
                   priority_density_v2 = .data$priority_density_v2, priority_density = .data$priority_density,
                   priority_within_reef = .data$priority_within_reef, class_basis = .data$class_basis,
                   site_visit_rank = .data$site_visit_rank, rank_on_reef = .data$rank_on_reef, sites_on_reef = .data$sites_on_reef,
                   site_visit_score = .data$site_visit_score, reef_density_v2 = .data$reef_density_v2,
                   reef_density_v1 = .data$reef_expected_problematic, density_source = .data$density_source,
                   searched_clean_12m = .data$searched_clean_12m, untowed_12m = .data$untowed_12m, reef_dives_12m = .data$reef_dives_12m,
                   reef_problematic_dives_12m = .data$reef_problematic_dives_12m, reef_tows_12m = .data$reef_tows_12m,
                   reef_strong_tows_12m = .data$reef_strong_tows_12m,
                   reef_outbreak_forecast_2026 = .data$reef_outbreak_forecast_2026, evidence = .data$evidence)
utils::write.csv(table_out, file.path(out_dir, "cots_cull_site_priority_density_v2_2027.csv"), row.names = FALSE, na = "")
log_message("Wrote ", file.path(out_dir, "cots_cull_site_priority_density_v2_2027.csv"))

count_class <- function(p, cls) sum(p %in% cls)
reefs <- table_out |>
  dplyr::group_by(.data$reef_name) |>
  dplyr::summarise(sites = dplyr::n(), class_basis = dplyr::first(.data$class_basis), density_source = dplyr::first(.data$density_source),
                   reef_density_v2 = dplyr::first(.data$reef_density_v2), reef_density_v1 = dplyr::first(.data$reef_density_v1),
                   dives_12m = dplyr::first(.data$reef_dives_12m), problematic_dives_12m = dplyr::first(.data$reef_problematic_dives_12m),
                   tows_12m = dplyr::first(.data$reef_tows_12m), strong_tows_12m = dplyr::first(.data$reef_strong_tows_12m),
                   priority_a = count_class(.data$priority_density_v2, "A"), priority_b = count_class(.data$priority_density_v2, "B"),
                   priority_c = count_class(.data$priority_density_v2, "C"), priority_d = count_class(.data$priority_density_v2, "D"),
                   priority_e = count_class(.data$priority_density_v2, "E"),
                   density_v1_a_to_c = count_class(.data$priority_density, c("A", "B", "C")),
                   within_reef_a_to_c = count_class(.data$priority_within_reef, c("A", "B", "C")), .groups = "drop") |>
  dplyr::arrange(dplyr::desc(.data$reef_density_v2))
utils::write.csv(reefs, file.path(out_dir, "cots_cull_site_priority_density_v2_2027_reefs.csv"), row.names = FALSE, na = "")

# How the classes moved, by where each reef's evidence came from.
shift <- table_out |>
  dplyr::group_by(.data$density_source) |>
  dplyr::summarise(reefs = dplyr::n_distinct(.data$reef_name), sites = dplyr::n(),
                   within_a_to_c = count_class(.data$priority_within_reef, c("A", "B", "C")),
                   density_v1_a_to_c = count_class(.data$priority_density, c("A", "B", "C")),
                   density_v2_a = count_class(.data$priority_density_v2, "A"),
                   density_v2_a_to_c = count_class(.data$priority_density_v2, c("A", "B", "C")),
                   reefs_without_a_to_c_v2 = sum(tapply(.data$priority_density_v2 %in% c("A", "B", "C"), .data$reef_name, sum) == 0),
                   median_reef_density_v2 = stats::median(tapply(.data$reef_density_v2, .data$reef_name, dplyr::first)),
                   .groups = "drop")
shift <- dplyr::bind_rows(shift, shift |>
                            dplyr::summarise(density_source = "all reefs", dplyr::across(-c("density_source", "median_reef_density_v2"), sum),
                                             median_reef_density_v2 = stats::median(.env$reefs$reef_density_v2)))
write_table(shift, cfg, "site_ranking_density_v2_2027_classes.csv")

# --- KMZ: the EoTR file coloured by the v2 classes ------------------------------------
fmt <- function(x, d) ifelse(is.finite(x), formatC(x, format = "f", digits = d), "")
k <- priority_kml(cfg, table_out, table_out$priority_density_v2,
                  fields = list(`Priority 2027 (density-aware v2)` = table_out$priority_density_v2,
                                `Density-aware priority (v1)` = table_out$priority_density,
                                `Within-reef priority` = table_out$priority_within_reef,
                                `Classes from` = table_out$class_basis,
                                `Site rank on reef (site record)` = paste(table_out$site_visit_rank, "of", table_out$sites_on_reef),
                                `Reef density` = paste0(round(100 * table_out$reef_density_v2), " % of dives expected problematic (",
                                                        table_out$density_source, ", recent evidence weighted)"),
                                `Last season` = table_out$evidence),
                  title = paste0("COTS cull site priority ", sr$target_season,
                                 ", density-aware v2 (A red, B orange, C yellow, D pale yellow, E grey)"))
write_staged_kmz(cfg, k$kml, "site_ranking_density_v2_kmz", file.path(out_dir, "cots_cull_site_priority_density_v2_2027.kmz"))

print(shift, width = Inf)
print(utils::head(reefs[, c("reef_name", "class_basis", "reef_density_v2", "reef_density_v1", "sites", "priority_a", "priority_b",
                            "priority_c", "density_v1_a_to_c", "within_reef_a_to_c")], 15), width = Inf)
print(as.data.frame(reefs[grepl("18-100a|18-075|16-029|14-118", reefs$reef_name),
                          c("reef_name", "reef_density_v2", "reef_density_v1", "sites", "priority_a", "priority_b", "priority_c",
                            "density_v1_a_to_c", "within_reef_a_to_c")]))
