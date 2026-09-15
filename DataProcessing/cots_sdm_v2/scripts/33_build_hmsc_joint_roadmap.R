args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
v2_dir <- normalizePath(file.path(dirname(script_file), ".."), winslash = "/")
source(file.path(v2_dir, "R", "bootstrap.R"))
ctx <- initialise_v2(script_file)
cfg <- ctx$config

input <- c(
  maps = file.path(
    cfg$paths$outputs_dir,
    "hmsc_fixed_protocol_visual_map_manifest.csv"
  ),
  settings = file.path(
    cfg$paths$outputs_dir, "hmsc_fixed_protocol_map_settings.csv"
  ),
  cv = file.path(
    cfg$paths$outputs_dir, "hmsc_joint_200m_cv_performance.csv"
  ),
  decomposition = file.path(
    cfg$paths$outputs_dir,
    "hmsc_joint_200m_high_density_decomposition.csv"
  ),
  ablation = file.path(
    cfg$paths$outputs_dir, "hmsc_effort_ablation_metrics.csv"
  ),
  variance = file.path(
    cfg$paths$outputs_dir, "hmsc_joint_200m_variance_partition.csv"
  ),
  beta = file.path(
    cfg$paths$outputs_dir, "hmsc_joint_200m_beta.csv"
  ),
  convergence = file.path(
    cfg$paths$outputs_dir, "hmsc_joint_200m_audit.csv"
  ),
  inventory = file.path(
    cfg$paths$outputs_dir, "joint_response_data_inventory.csv"
  ),
  alignment = file.path(
    cfg$paths$outputs_dir, "joint_response_alignment_audit.csv"
  )
)
if (any(!file.exists(input))) {
  stop("Run scripts 23, 26, 27 and 32 before building this report.")
}

maps <- utils::read.csv(input[["maps"]], stringsAsFactors = FALSE)
settings <- utils::read.csv(input[["settings"]])
cv <- utils::read.csv(input[["cv"]])
decomposition <- utils::read.csv(input[["decomposition"]])
ablation <- utils::read.csv(input[["ablation"]])
variance <- utils::read.csv(input[["variance"]])
beta <- utils::read.csv(input[["beta"]])
audit <- utils::read.csv(input[["convergence"]])
inventory <- utils::read.csv(
  input[["inventory"]], stringsAsFactors = FALSE
)
alignment <- utils::read.csv(
  input[["alignment"]], stringsAsFactors = FALSE
)

if (nrow(maps) < 2L || nrow(maps) > 4L) {
  stop("The Hmsc report requires two to four reef comparisons.")
}
map_files <- file.path(v2_dir, "reports", maps$figure)
if (any(!file.exists(map_files))) {
  stop("One or more Hmsc comparison images are missing.")
}

fmt <- function(x, digits = 3L) {
  ifelse(is.finite(x), formatC(x, digits = digits, format = "f"), "NA")
}
markdown_table <- function(headers, rows) {
  c(
    paste0("| ", paste(headers, collapse = " | "), " |"),
    paste0("|", paste(rep("---", length(headers)), collapse = "|"), "|"),
    apply(rows, 1, function(x) {
      paste0("| ", paste(x, collapse = " | "), " |")
    })
  )
}

high_cv <- cv[cv$response == "cull_high_density", ]
high_var <- variance[variance$response == "cull_high_density", ]
habitat_fraction <- high_var$variance_fraction[
  high_var$component == "Multiscale habitat"
]
cull_effort_fraction <- high_var$variance_fraction[
  high_var$component == "Cull effort"
]
identity_fraction <- sum(high_var$variance_fraction[
  high_var$component %in% c("Random: cell", "Random: reef")
])
static_ablation <- ablation[ablation$model == "glm_static_habitat", ]
effort_ablation <- ablation[ablation$model == "glm_effort_year", ]
high_beta <- beta[beta$response == "cull_high_density", ]
spatial_terms <- c(
  "BPI_s30m", "SLO_s60m", "VRM_s60m", "DEM_s30m", "SVF_s120m",
  "coral_suitability_mean"
)
high_beta <- high_beta[high_beta$covariate %in% spatial_terms, ]

map_rows <- maps |>
  dplyr::transmute(
    Reef = .data$reef_name,
    `Raw 5%-median-95%` = paste0(
      fmt(.data$raw_probability_q05), "-",
      fmt(.data$raw_probability_median), "-",
      fmt(.data$raw_probability_q95)
    ),
    `Rank correlation` = fmt(.data$spearman_rank_correlation),
    `Mean absolute rank difference` = fmt(
      .data$mean_absolute_percentile_difference
    )
  )
map_blocks <- unlist(lapply(seq_len(nrow(maps)), function(i) {
  c(
    paste0("### ", maps$reef_name[[i]]),
    "",
    paste0(
      "![Hmsc fixed-protocol versus clean 2026: ", maps$reef_name[[i]],
      "](", maps$figure[[i]], ")"
    ),
    ""
  )
}))
beta_rows <- data.frame(
  Predictor = high_beta$covariate,
  `Posterior mean` = fmt(high_beta$posterior_mean),
  `Support positive` = fmt(high_beta$posterior_support_positive),
  check.names = FALSE
)
inventory_rows <- inventory |>
  dplyr::transmute(
    Method = .data$sheet,
    Rows = format(.data$rows, big.mark = ","),
    Years = paste(.data$date_min, .data$date_max, sep = " to "),
    Reefs = .data$reefs,
    `Named sites` = format(.data$named_sites, big.mark = ","),
    `Rows with COTS/scar evidence` = format(
      .data$response_nonzero_rows, big.mark = ","
    )
  )
alignment_rows <- alignment |>
  dplyr::transmute(
    Comparison = .data$comparison,
    `Named units` = ifelse(
      is.na(.data$units), "-", format(.data$units, big.mark = ",")
    ),
    `Units repeated across days` = ifelse(
      is.na(.data$units_with_repeat_days), "-",
      format(.data$units_with_repeat_days, big.mark = ",")
    ),
    `Shared reef-days` = ifelse(
      is.na(.data$shared_reef_days), "-",
      format(.data$shared_reef_days, big.mark = ",")
    )
  )

report <- c(
  "# Hmsc visual assessment and joint-modelling roadmap",
  "",
  "## Decision",
  "",
  "**Do not promote the current Hmsc surface. Retain the frozen legacy",
  "clean-year model as the current best operational map.**",
  "",
  "The rejected successor was not Hmsc: it was the static multiscale XGBoost",
  "candidate. The maps below are newly derived from the fitted production Hmsc",
  "model and use exactly the same four reef comparisons. They are a diagnostic",
  "of what Hmsc has learned about transferable fine-scale habitat.",
  "",
  "## What was mapped",
  "",
  "The mapped response is Cull high density (pooled CPUE at least 0.02). The",
  "surface uses posterior-mean fixed effects for 2026, fixes Cull and manta",
  "effort at the median matched-data protocol, and excludes reef and 200 m cell",
  "random effects. It therefore asks: before a survey, and without being told",
  "reef identity, how does the fitted Hmsc model rank habitat within a reef?",
  "",
  sprintf(
    "The fixed protocol is %.0f Cull bottom-time units, %.0f Cull dives, %.0f manta tows, %.1f m manta distance and a hard-coral observation fraction of %.2f.",
    settings$cull_bottom_time, settings$cull_n_dives, settings$manta_n_tows,
    settings$manta_total_distance_m, settings$manta_hard_coral_fraction
  ),
  "",
  "The full held-out Hmsc score cannot be attached to this deployable map.",
  sprintf(
    "The full detection-conditioned fit achieved overall AUC %s and within-reef AUC %s, but %.1f%% of Cull high-density variance was assigned to multiscale habitat, %.1f%% to Cull effort, and %.1f%% to reef/cell identity effects.",
    fmt(high_cv$roc_auc), fmt(decomposition$within_reef_roc_auc),
    100 * habitat_fraction, 100 * cull_effort_fraction,
    100 * identity_fraction
  ),
  sprintf(
    "As a check, the simpler held-out static-habitat ablation had AUC %s, whereas effort plus year alone had AUC %s.",
    fmt(static_ablation$roc_auc), fmt(effort_ablation$roc_auc)
  ),
  "",
  "## Four-reef visual comparison",
  "",
  "All first and second panels use within-reef percentiles; this reveals rank",
  "structure without treating the compressed raw Hmsc scale as a failure. The",
  "third panel is Hmsc minus incumbent percentile. Crosses are historical Cull",
  "locations and are not response symbols.",
  "",
  markdown_table(names(map_rows), as.matrix(map_rows)),
  "",
  "The rank correlation is negative on every candidate reef. This does not by",
  "itself prove the incumbent is ecologically correct, but it confirms that the",
  "current Hmsc habitat component has learned a substantially different spatial",
  "ordering and requires expert scrutiny rather than metric-led promotion.",
  "",
  map_blocks,
  "## Why the current joint model is not yet delivering a useful habitat map",
  "",
  "1. The four responses are observation channels for one focal animal",
  "   (Cull direct presence/high density and manta direct/scar presence), not a",
  "   genuine biological community. Shared survey history can create residual",
  "   association without identifying a habitat mechanism.",
  "2. Reef identity, cell identity and realised effort explain almost all of",
  "   the fitted Cull high-density variation. Those terms help reconstruct the",
  "   historical sample but cannot create transferable 30-120 m habitat detail.",
  "3. The observation methods have different footprints, effort and error.",
  "   Aggregating them to 200 m cell-years does not explicitly model imperfect",
  "   detection, scar false positives, depletion after culling, or preferential",
  "   selection of sites.",
  "4. The terrain representation is sparse: each terrain family is represented",
  "   at only one selected support. It does not yet implement the paper's full",
  "   DEM-first Gaussian multiscale derivation and scale-selection workflow.",
  sprintf(
    "5. Random-effect convergence remains inadequate (maximum upper PSRF %s; minimum ESS %s), so near-perfect residual correlations must not be interpreted as established ecological coupling.",
    fmt(audit$maximum_psrf_upper), fmt(audit$minimum_ess, 1L)
  ),
  "",
  "### Current Cull high-density habitat coefficients",
  "",
  markdown_table(names(beta_rows), as.matrix(beta_rows)),
  "",
  "Posterior support is one-sided for coral suitability, but weak for most",
  "terrain terms. These coefficients should be treated as conditional",
  "associations within the current detection-conditioned specification.",
  "",
  "## Data already available locally",
  "",
  "The combined 2026 Cull-Manta-RHIS workbook is much richer than the current",
  "four-response Hmsc table. It already contains Cull size cohorts; manta",
  "direct counts, scar counts, tow distance and coral-cover classes; and RHIS",
  "adult/juvenile COTS counts, scars by coral growth form, live coral by growth",
  "form, recent mortality, visibility, depth and survey habitat.",
  "",
  markdown_table(names(inventory_rows), as.matrix(inventory_rows)),
  "",
  markdown_table(names(alignment_rows), as.matrix(alignment_rows)),
  "",
  "These are preliminary compatibility counts, not proof of matched sampling.",
  "Same reef/day is deliberately weaker than same footprint and protocol. The",
  "large number of repeated named sites and shared reef-days nevertheless shows",
  "that a visit-level alignment audit is worthwhile before seeking new data.",
  "",
  "The older ReefScan script expects a data/ReefScan CSV directory, but that",
  "directory is not present in the current repository. Its previous ingestion",
  "also discarded survey date, effort, segment/image identity, confidence, coral",
  "composition and separate direct/scar outcomes. Those fields are required if",
  "ReefScan is to strengthen detection modelling rather than add pseudo-replicates.",
  "",
  "## Data that would make a joint approach materially stronger",
  "",
  "| Priority | Compatible sampling unit and fields | Role in model |",
  "|---|---|---|",
  "| 1 | Individual Cull dive/track: date-time, observer, pre/post control, counts by size, bottom time, footprint, visibility, depth, scars, coral taxa/cover | Separates abundance from detection and removal |",
  "| 1 | Manta tow segment/track: direct count, scar count, tow time/distance, observer, visibility, coral cover/taxa and georeferenced footprint | Method-specific detection and scar error |",
  "| 1 | ReefScan image/segment effort with validated direct-COTS and scar labels, negatives, confidence and coral composition | Independent observation channel with quantified classification error |",
  "| 1 | Repeat visits to the same cells before intervention, ideally paired across methods | Identifies non-detection versus true absence |",
  "| 2 | Adult and juvenile COTS, fresh scars, susceptible coral taxa/cover, recent coral mortality and standardized predator-guild observations on compatible units | Creates a genuinely multivariate ecological response |",
  "| 2 | Culling history, time since first detection/outbreak, removals and revisit interval | Models depletion and outbreak phase |",
  "| 2 | DEM-first Gaussian terrain at 30/60/120/240/360 m plus geomorphic zone, wave/current exposure and larval connectivity | Multiscale habitat and dispersal mechanism |",
  "| 3 | SST/heat stress, cyclones, chlorophyll/productivity, turbidity and time-varying coral condition aligned before each visit | Environmental state and disturbance context |",
  "",
  "A response matrix becomes genuinely useful when the columns are compatible",
  "ecological states observed on common units. Adult COTS, juvenile COTS, fresh",
  "scars, susceptible coral groups and coral mortality are stronger candidates",
  "than treating multiple detections of the same animal as if they were species.",
  "",
  "## Recommended modelling sequence",
  "",
  "1. **H1: transferable fixed-effect benchmark.** Refit without realised effort",
  "   and reef/cell identity, using complete Gaussian multiscale terrain and",
  "   coral composition. Validate held-out cells, reefs and future years.",
  "2. **H2: observation-aware visit model.** Retain individual visits and add",
  "   method-specific detection, effort, observer, scar error and Cull-removal",
  "   processes. If Hmsc cannot express the required observation model, use an",
  "   integrated multivariate model in Stan, NIMBLE or TMB and keep Hmsc as a",
  "   benchmark rather than forcing the data into Hmsc.",
  "3. **H3: random-structure sensitivity.** Compare reef only, cell only, reef",
  "   plus spatial field, and regularised temporal structures. Require stable",
  "   chains before interpreting residual associations.",
  "4. **H4: genuine multivariate ecology.** Add compatible ReefScan, coral taxon",
  "   and COTS life-stage responses; test whether shared latent gradients improve",
  "   out-of-sample Cull high-density ordering, not merely fit to sampled reefs.",
  "5. **H5: prospective operational test.** Pre-register fixed-effort predictions",
  "   for stratified random sites within reefs, then score discrimination and",
  "   calibration without using realised survey effort.",
  "",
  "Every iteration keeps four gates: within-reef rank validation, leave-reef",
  "transfer, forward-year validation, and two-to-four reef expert map review.",
  "A higher global AUC cannot override a failed spatial-plausibility gate.",
  "",
  "## Immediate next implementation",
  "",
  "The next bounded experiment should be H1: refit the Hmsc response structure",
  "with the expanded Gaussian multiscale predictors, no reef/cell random effects",
  "in the deployable habitat score, and no realised effort. In parallel, build a",
  "visit-level data audit showing which Cull, manta and ReefScan records can be",
  "aligned on cell, date and pre-intervention status. H2 should start only after",
  "that audit identifies enough repeated compatible visits to estimate detection.",
  "",
  "## Reproducibility",
  "",
  "- scripts/23_fit_hmsc_joint_production.R fits the production Hmsc model.",
  "- scripts/26_summarise_hmsc_joint.R summarizes performance and parameters.",
  "- scripts/27_diagnose_hmsc_effort_signal.R separates habitat and effort signal.",
  "- scripts/32_build_hmsc_reef_maps.R creates these fixed-protocol maps.",
  "- scripts/34_audit_joint_response_data.R audits available joint responses.",
  "- scripts/33_build_hmsc_joint_roadmap.R regenerates this report.",
  "- Hmsc method: https://doi.org/10.1111/2041-210X.13345"
)

report_dir <- file.path(v2_dir, "reports")
dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)
writeLines(
  report,
  file.path(report_dir, "HMSC_VISUAL_JOINT_MODELLING_ROADMAP.md")
)
message("Hmsc visual and joint-modelling roadmap complete.")
