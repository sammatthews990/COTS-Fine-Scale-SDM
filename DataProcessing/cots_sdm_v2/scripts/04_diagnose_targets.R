args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
v2_dir <- normalizePath(file.path(dirname(script_file), ".."), winslash = "/")
source(file.path(v2_dir, "R", "bootstrap.R"))
ctx <- initialise_v2(script_file)
cfg <- ctx$config

required_artifacts <- c(
  "site_year_all.rds",
  "site_year_first_visit.rds",
  "site_history_diagnostics.rds"
)
artifact_paths <- file.path(cfg$paths$artifacts_dir, required_artifacts)
if (any(!file.exists(artifact_paths))) {
  stop("Run scripts/01_prepare_data.R before target diagnostics.")
}

site_year_all <- readRDS(artifact_paths[1])
site_year_first <- readRDS(artifact_paths[2])
site_history <- readRDS(artifact_paths[3])
threshold <- cfg$response$cpue_threshold

summary_row <- function(name, truth, cpue) {
  tibble::tibble(
    definition = name,
    rows = length(truth),
    positives = sum(truth == 1L),
    prevalence = mean(truth),
    median_cpue = stats::median(cpue, na.rm = TRUE),
    mean_cpue = mean(cpue, na.rm = TRUE)
  )
}

definition_summary <- dplyr::bind_rows(
  summary_row(
    "all_dives_pooled_cpue",
    site_year_all$high_density,
    site_year_all$cpue_mean
  ),
  summary_row(
    "first_visit_pooled_cpue",
    site_year_first$high_density,
    site_year_first$cpue_mean
  ),
  summary_row(
    "any_dive_max_cpue_diagnostic",
    as.integer(site_year_all$cpue_max_diagnostic >= threshold),
    site_year_all$cpue_max_diagnostic
  )
)

comparison <- site_year_all |>
  dplyr::select(dplyr::all_of(c(
    "site_year_id", "first_survey_date", "last_survey_date", "n_dives",
    "cpue_mean", "cpue_max_diagnostic", "high_density"
  ))) |>
  dplyr::rename(
    cpue_all_dives = "cpue_mean",
    high_all_dives = "high_density"
  ) |>
  dplyr::left_join(
    site_year_first |>
      dplyr::select(dplyr::all_of(c(
        "site_year_id", "cpue_mean", "high_density"
      ))) |>
      dplyr::rename(
        cpue_first_visit = "cpue_mean",
        high_first_visit = "high_density"
      ),
    by = "site_year_id"
  ) |>
  dplyr::mutate(
    multiple_visit_dates = .data$first_survey_date != .data$last_survey_date,
    high_any_dive_max = as.integer(
      .data$cpue_max_diagnostic >= threshold
    )
  )

disagreement_summary <- comparison |>
  dplyr::summarise(
    site_years = dplyr::n(),
    multiple_visit_date_site_years = sum(.data$multiple_visit_dates),
    first_visit_target_changes = sum(
      .data$high_all_dives != .data$high_first_visit
    ),
    pooled_to_first_positive = sum(
      .data$high_all_dives == 0L & .data$high_first_visit == 1L
    ),
    pooled_to_first_negative = sum(
      .data$high_all_dives == 1L & .data$high_first_visit == 0L
    ),
    max_vs_pooled_target_changes = sum(
      .data$high_any_dive_max != .data$high_all_dives
    ),
    max_positive_pooled_negative = sum(
      .data$high_any_dive_max == 1L & .data$high_all_dives == 0L
    ),
    cpue_correlation_all_vs_first = stats::cor(
      .data$cpue_all_dives, .data$cpue_first_visit,
      use = "complete.obs"
    )
  )

effort_breaks <- c(0, 1, 2, 5, 10, 20, 50, Inf)
effort_labels <- c("1", "2", "3-5", "6-10", "11-20", "21-50", ">50")
history_effort <- site_history |>
  dplyr::mutate(
    lifetime_dive_bin = cut(
      .data$lifetime_dives,
      breaks = effort_breaks,
      labels = effort_labels,
      include.lowest = TRUE,
      right = TRUE
    )
  ) |>
  dplyr::group_by(.data$lifetime_dive_bin) |>
  dplyr::summarise(
    sites = dplyr::n(),
    median_dives = stats::median(.data$lifetime_dives),
    ever_high_by_single_dive_max = mean(.data$ever_high_by_dive_max),
    any_high_by_pooled_site_year = mean(.data$any_high_density_year),
    median_high_year_fraction = stats::median(
      .data$high_density_year_fraction
    ),
    .groups = "drop"
  )

utils::write.csv(
  definition_summary,
  file.path(cfg$paths$outputs_dir, "target_definition_summary.csv"),
  row.names = FALSE
)
utils::write.csv(
  disagreement_summary,
  file.path(cfg$paths$outputs_dir, "target_definition_disagreements.csv"),
  row.names = FALSE
)
utils::write.csv(
  history_effort,
  file.path(cfg$paths$outputs_dir, "history_effort_bias.csv"),
  row.names = FALSE
)

cat("Target diagnostics complete.\n")
print(definition_summary)
print(disagreement_summary)
print(history_effort)
