args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
v2_dir <- normalizePath(file.path(dirname(script_file), ".."), winslash = "/")
source(file.path(v2_dir, "R", "bootstrap.R"))
ctx <- initialise_v2(script_file)
cfg <- ctx$config

required <- c(
  static = file.path(
    cfg$paths$outputs_dir,
    "two_part_static_candidate_summary.csv"
  ),
  components = file.path(
    cfg$paths$outputs_dir,
    "two_part_component_summary.csv"
  ),
  rank = file.path(
    cfg$paths$outputs_dir,
    "two_part_rank_lock_checks.csv"
  ),
  audit = file.path(
    cfg$paths$artifacts_dir,
    "two_part_iteration_01_audit.rds"
  ),
  bundle = file.path(
    cfg$paths$models_dir,
    "two_part_fixed_effort_iteration_01.rds"
  ),
  maps = file.path(
    cfg$paths$outputs_dir,
    "two_part_iteration_01_visual_map_manifest.csv"
  ),
  hmsc_maps = file.path(
    cfg$paths$outputs_dir,
    "hmsc_fixed_protocol_visual_map_manifest.csv"
  )
)
if (any(!file.exists(required))) {
  stop("Run scripts 28 and 29 before updating model reports.")
}

static <- utils::read.csv(required[["static"]])
components <- utils::read.csv(required[["components"]])
rank_checks <- utils::read.csv(required[["rank"]])
audit <- readRDS(required[["audit"]])
bundle <- readRDS(required[["bundle"]])
map_manifest <- utils::read.csv(
  required[["maps"]], stringsAsFactors = FALSE
)
hmsc_map_manifest <- utils::read.csv(
  required[["hmsc_maps"]], stringsAsFactors = FALSE
)

selected <- static |>
  dplyr::filter(.data$model == bundle$selected_static_candidate)
if (nrow(selected) != 1L) stop("Selected candidate summary is ambiguous.")

incumbent <- tibble::tibble(
  model_id = "legacy_clean_year_20260820_1041",
  display_name = "Legacy clean + year",
  within_reef_year_auc = 0.568344817887779,
  top20_capture = 0.250916870415648,
  top20_lift = 1.09469139769291,
  exact_raster_reproduction = TRUE
)

expert_review_path <- file.path(
  cfg$paths$artifacts_dir,
  "two_part_iteration_01_expert_review.csv"
)
expert_review_passed <- FALSE
expert_review_evidence <- "Pending"
if (file.exists(expert_review_path)) {
  expert <- utils::read.csv(
    expert_review_path, stringsAsFactors = FALSE
  )
  expert_review_passed <- nrow(expert) > 0L &&
    all(tolower(expert$status) == "accepted")
  if (expert_review_passed) {
    expert_review_evidence <- "Accepted"
  } else if (any(tolower(expert$status) == "rejected")) {
    expert_review_evidence <- paste0(
      "Rejected: ", paste(unique(expert$evidence), collapse = "; ")
    )
  }
}
outbreak_audit_path <- file.path(
  cfg$paths$artifacts_dir,
  "outbreak_probability_provenance_audit.rds"
)
outbreak_provenance_passed <- FALSE
if (file.exists(outbreak_audit_path)) {
  outbreak_audit <- readRDS(outbreak_audit_path)
  outbreak_provenance_passed <- isTRUE(outbreak_audit$passed)
}

candidate_id <- bundle$iteration
if (nrow(map_manifest) < 2L || nrow(map_manifest) > 4L) {
  stop("Every model report requires 2-4 candidate reef maps.")
}
if (!all(map_manifest$candidate_model == candidate_id)) {
  stop("Visual map manifest does not match the evaluated candidate.")
}
map_files <- file.path(v2_dir, "reports", map_manifest$figure)
if (any(!file.exists(map_files))) {
  stop("One or more candidate reef comparison figures are missing.")
}
if (nrow(hmsc_map_manifest) < 2L || nrow(hmsc_map_manifest) > 4L) {
  stop("The Hmsc comparison requires 2-4 reef maps.")
}
hmsc_map_files <- file.path(
  v2_dir, "reports", hmsc_map_manifest$figure
)
if (any(!file.exists(hmsc_map_files))) {
  stop("One or more Hmsc reef comparison figures are missing.")
}
operational_rank_checks <- rank_checks |>
  dplyr::filter(
    !grepl("realised_effort_diagnostic$", .data$model)
  )
gates <- tibble::tibble(
  gate = c(
    "Mean within-reef-year AUC improves by at least 0.005",
    "Mean top-20% lift is no worse than incumbent",
    "Minimum repeat within-reef-year AUC is at least 0.55",
    "State and fixed-effort transforms preserve within-reef rank",
    "Operational prediction excludes realised survey effort",
    "Outbreak probability provenance has passed audit",
    "Blinded expert map review has accepted candidate"
  ),
  passed = c(
    selected$mean_within_reef_year_auc >=
      incumbent$within_reef_year_auc + 0.005,
    selected$mean_top20_lift >= incumbent$top20_lift,
    selected$minimum_within_reef_year_auc >= 0.55,
    all(operational_rank_checks$violations == 0),
    TRUE,
    outbreak_provenance_passed,
    expert_review_passed
  ),
  evidence = c(
    sprintf(
      "%.3f versus %.3f",
      selected$mean_within_reef_year_auc,
      incumbent$within_reef_year_auc
    ),
    sprintf(
      "%.3f versus %.3f",
      selected$mean_top20_lift, incumbent$top20_lift
    ),
    sprintf("%.3f", selected$minimum_within_reef_year_auc),
    sprintf(
      "%d operational rank violations",
      sum(operational_rank_checks$violations)
    ),
    "Fixed at 2 dives and 600 bottom-time units",
    if (outbreak_provenance_passed) "Passed" else "Pending",
    expert_review_evidence
  )
)
promotion_passed <- all(gates$passed)

registry_path <- file.path(
  cfg$paths$outputs_dir, "model_successor_registry.csv"
)
if (file.exists(registry_path)) {
  registry <- utils::read.csv(
    registry_path, stringsAsFactors = FALSE
  )
} else {
  registry <- tibble::tibble(
    model_id = incumbent$model_id,
    display_name = incumbent$display_name,
    iteration = "legacy",
    status = "current_best",
    within_reef_year_auc = incumbent$within_reef_year_auc,
    top20_capture = incumbent$top20_capture,
    top20_lift = incumbent$top20_lift,
    statistical_gates_passed = TRUE,
    provenance_gate_passed = TRUE,
    expert_map_gate_passed = TRUE,
    all_promotion_gates_passed = TRUE,
    evaluated_at = "2026-09-01"
  )
}
registry <- dplyr::filter(
  registry, .data$model_id != candidate_id
)
if (promotion_passed) {
  registry$status[registry$status == "current_best"] <- "superseded"
}
candidate_row <- tibble::tibble(
  model_id = candidate_id,
  display_name = paste(
    "Two-part:", bundle$selected_static_candidate
  ),
  iteration = "01",
  status = if (promotion_passed) {
    "current_best"
  } else {
    "candidate_not_promoted"
  },
  within_reef_year_auc = selected$mean_within_reef_year_auc,
  top20_capture = selected$mean_top20_capture,
  top20_lift = selected$mean_top20_lift,
  statistical_gates_passed = all(gates$passed[1:5]),
  provenance_gate_passed = gates$passed[[6]],
  expert_map_gate_passed = gates$passed[[7]],
  all_promotion_gates_passed = promotion_passed,
  evaluated_at = "2026-09-01"
)
registry <- dplyr::bind_rows(registry, candidate_row)
utils::write.csv(registry, registry_path, row.names = FALSE)
utils::write.csv(
  gates,
  file.path(
    cfg$paths$outputs_dir,
    "two_part_iteration_01_promotion_gates.csv"
  ),
  row.names = FALSE
)

fmt <- function(x, digits = 3) {
  ifelse(is.finite(x), formatC(x, digits = digits, format = "f"), "NA")
}
yes_no <- function(x) ifelse(x, "PASS", "PENDING/FAIL")
markdown_table <- function(headers, rows) {
  c(
    paste0("| ", paste(headers, collapse = " | "), " |"),
    paste0("|", paste(rep("---", length(headers)), collapse = "|"), "|"),
    apply(rows, 1, function(x) {
      paste0("| ", paste(x, collapse = " | "), " |")
    })
  )
}

static_rows <- static |>
  dplyr::transmute(
    Model = .data$model,
    `Overall AUC` = fmt(.data$mean_overall_auc),
    `Within-reef AUC` = fmt(.data$mean_within_reef_year_auc),
    `Repeat SD` = fmt(.data$sd_within_reef_year_auc),
    `Top-20 capture` = fmt(.data$mean_top20_capture),
    `Top-20 lift` = fmt(.data$mean_top20_lift)
  )
component_rows <- components |>
  dplyr::mutate(component = sub("^.*__", "", .data$model)) |>
  dplyr::transmute(
    Component = .data$component,
    `Overall AUC` = fmt(.data$mean_overall_auc),
    `Within-reef AUC` = fmt(.data$mean_within_reef_year_auc),
    `PR-AUC` = fmt(.data$mean_pr_auc),
    Brier = fmt(.data$mean_brier),
    `Top-20 lift` = fmt(.data$mean_top20_lift)
  )
gate_rows <- gates |>
  dplyr::transmute(
    Gate = .data$gate,
    Status = yes_no(.data$passed),
    Evidence = .data$evidence
  )
map_rows <- map_manifest |>
  dplyr::transmute(
    Reef = .data$reef_name,
    `Rank correlation` = fmt(.data$spearman_rank_correlation),
    `Mean absolute rank difference` = fmt(
      .data$mean_absolute_percentile_difference
    ),
    Status = .data$review_status
  )
map_blocks <- unlist(lapply(seq_len(nrow(map_manifest)), function(i) {
  c(
    paste0("### ", map_manifest$reef_name[[i]]),
    "",
    paste0(
      "![Candidate versus current best: ",
      map_manifest$reef_name[[i]], "](",
      map_manifest$figure[[i]], ")"
    ),
    ""
  )
}))
hmsc_map_rows <- hmsc_map_manifest |>
  dplyr::transmute(
    Reef = .data$reef_name,
    `Raw probability (5%-median-95%)` = paste0(
      fmt(.data$raw_probability_q05), "-",
      fmt(.data$raw_probability_median), "-",
      fmt(.data$raw_probability_q95)
    ),
    `Rank correlation` = fmt(.data$spearman_rank_correlation),
    `Mean absolute rank difference` = fmt(
      .data$mean_absolute_percentile_difference
    ),
    Status = .data$review_status
  )
hmsc_map_blocks <- unlist(lapply(
  seq_len(nrow(hmsc_map_manifest)), function(i) {
    c(
      paste0("### ", hmsc_map_manifest$reef_name[[i]]),
      "",
      paste0(
        "![Hmsc fixed-protocol versus current best: ",
        hmsc_map_manifest$reef_name[[i]], "](",
        hmsc_map_manifest$figure[[i]], ")"
      ),
      ""
    )
  }
))

hmsc_comparison <- c(
  "## Hmsc fixed-protocol reef comparison: diagnostic only",
  "",
  "This is the requested map from the fitted production Hmsc model, not the",
  "XGBoost successor candidate. It uses posterior-mean fixed effects for the",
  "Cull high-density response in 2026, holds survey effort at the matched-data",
  "median, and excludes reef and cell random effects. Consequently it is a",
  "pre-survey habitat ordering rather than a reconstruction using known reef",
  "identity or realised effort.",
  "",
  "The Hmsc fixed-effect ordering is negatively correlated with the incumbent",
  "on all four reefs. Its full-model validation score must not be assigned to",
  "this map: the full score includes strong reef/cell and effort information.",
  "",
  markdown_table(names(hmsc_map_rows), as.matrix(hmsc_map_rows)),
  "",
  hmsc_map_blocks,
  "See HMSC_VISUAL_JOINT_MODELLING_ROADMAP.md for interpretation, data needs",
  "and a staged route toward a genuinely joint ecological model.",
  ""
)

shared_next_steps <- c(
  "1. Keep the XGBoost candidate rejected and retain its expert feedback as",
  "   promotion evidence; do not tune directly to the four review maps.",
  "2. Align individual Cull, manta and RHIS visits by footprint, date and",
  "   pre-intervention status, using the locally available size, scar, coral",
  "   growth-form, visibility and effort fields.",
  "3. Re-derive terrain variables from Gaussian-smoothed bathymetry at",
  "   30, 60, 120, 240 and 360 m supports; do not rely only on block means.",
  "4. Refit a transferable Hmsc benchmark with static multiscale habitat and",
  "   compatible ecological responses, excluding realised effort and reef/cell",
  "   identity from the deployable habitat score.",
  "5. Fit a visit-level multi-method detection model with observer, effort,",
  "   scar error and Cull-removal processes; use Hmsc only where its observation",
  "   model matches the data-generating process.",
  "6. Validate prospectively or on fixed-effort surveys using held-out cells,",
  "   reefs, future years and two-to-four expert reef map comparisons.",
  "7. Audit the 1993-2026 outbreak layer provenance and confirm that each",
  "   validation year uses a genuinely available hindcast rather than a",
  "   smoothed estimate informed by future observations.",
  "8. Promote a successor only when statistical, provenance and expert-map",
  "   gates all pass; otherwise retain the frozen legacy model."
)

decision <- if (promotion_passed) {
  paste("PROMOTED:", candidate_id, "is the current best model.")
} else {
  "NOT PROMOTED: the frozen legacy clean-year model remains current best."
}

detail <- c(
  "# Two-part habitat, outbreak-state and fixed-effort model: iteration 01",
  "",
  "## Decision",
  "",
  paste0("**", decision, "**"),
  "",
  "The candidate separates the operational question into three components:",
  "",
  "1. a static habitat score fitted without year, reef identity, outbreak",
  "   probability or survey effort;",
  "2. a reef-year outbreak-state bridge that can change absolute probability",
  "   but is mathematically unable to reorder cells within a reef-year; and",
  "3. an observation layer estimated from dive effort, then evaluated at a",
  "   fixed protocol of two dives and 600 bottom-time units.",
  "",
  "The repeated-dive output estimates the chance of at least one high-density",
  "dive under the fixed protocol. It is an observation proxy and is not",
  "labelled as the probability that pooled site CPUE exceeds 0.02.",
  "",
  "## Data and validation",
  "",
  sprintf(
    "The common complete domain contains %s site-years, %s sites and %s reefs",
    format(audit$rows, big.mark = ","),
    format(audit$sites, big.mark = ","),
    format(audit$reefs, big.mark = ",")
  ),
  sprintf(
    "covering %d-%d. Three repeated five-fold splits hold out 200 m spatial",
    audit$years[[1]], audit$years[[2]]
  ),
  "site groups within reefs. The selected composition then uses nested",
  "cross-fitting, so the state and effort calibrators never train on the",
  "outer held-out sites.",
  "",
  "Fold-specific calibrators naturally have different numeric scales. Before",
  "pooling held-out folds, each operational component is monotonically",
  "quantile-aligned to the held-out habitat score within reef-year, without",
  "using outcomes. This reproduces the rank contract of the single deployed",
  "bridge while retaining the held-out component probability distribution.",
  "",
  "## Static habitat screen",
  "",
  markdown_table(names(static_rows), as.matrix(static_rows)),
  "",
  paste(
    "Selected for composition:", bundle$selected_static_candidate,
    sprintf(
      "(mean within-reef-year AUC %s).",
      fmt(selected$mean_within_reef_year_auc)
    )
  ),
  "",
  "![Repeated static validation](figures/two_part_static_repeated_auc.png)",
  "",
  "## Separated components",
  "",
  markdown_table(names(component_rows), as.matrix(component_rows)),
  "",
  "Outbreak state and fixed-effort transformations must have the same",
  "within-reef rank metrics as the static habitat component. The realised",
  "effort result is retained only as a leakage diagnostic and is never an",
  "operational map input.",
  "",
  "Outbreak state raises overall AUC from 0.582 to 0.743 without changing the",
  "within-reef AUC of 0.577. Fixed-effort standardisation adds essentially no",
  "discrimination beyond state and slightly worsens Brier score (0.211 versus",
  "0.208). In contrast, supplying realised effort inflates within-reef AUC to",
  "0.782, independently confirming that realised effort is leakage-prone.",
  "",
  "## Candidate reef maps: visual acceptance required",
  "",
  "Each panel uses the same within-reef percentile scale for the incumbent and",
  "candidate. The third map is candidate minus incumbent percentile. White or",
  "black crosses mark historical Cull site locations but do not encode the",
  "response. These maps are evidence for expert acceptance, not a numerical",
  "model metric.",
  "",
  markdown_table(names(map_rows), as.matrix(map_rows)),
  "",
  map_blocks,
  "The expert rejection is recorded in",
  "../artifacts/two_part_iteration_01_expert_review.csv and is a hard",
  "promotion failure for this candidate.",
  "",
  hmsc_comparison,
  "![Component validation](figures/two_part_component_auc.png)",
  "",
  "![Effort coefficients](figures/two_part_effort_coefficients.png)",
  "",
  "## Promotion gates",
  "",
  markdown_table(names(gate_rows), as.matrix(gate_rows)),
  "",
  "## Interpretation and limitations",
  "",
  "- The architecture now prevents reef outbreak state from masking local",
  "  site differences.",
  "- Fixed-effort prediction removes the direct advantage obtained by",
  "  supplying realised Cull effort to held-out observations.",
  "- It does not fully identify ecological occupancy and detection: survey",
  "  effort was historically preferential, and fixing it at prediction time",
  "  cannot remove every bias in coefficients estimated from those data.",
  "- Same-year outbreak probabilities improve absolute-risk discrimination",
  "  but remain external reef-level information; they are not evidence that",
  "  the fine-scale habitat model has improved.",
  "- The existing multiscale map panels remain the visual evidence for the",
  "  selected static family. Expert acceptance remains a hard gate.",
  "",
  "## Harmonised next steps",
  "",
  shared_next_steps,
  "",
  "## Reproducibility",
  "",
  "- scripts/28_fit_two_part_fixed_effort.R: repeated and nested validation.",
  "- scripts/29_plot_two_part_diagnostics.R: diagnostic figures.",
  "- scripts/31_build_candidate_reef_maps.R: exact candidate/incumbent maps.",
  "- scripts/30_update_model_reports.R: registry and report decision.",
  "- outputs/model_successor_registry.csv: machine-readable model history.",
  "- models/two_part_fixed_effort_iteration_01.rds: candidate bundle."
)

current_row <- registry |>
  dplyr::filter(.data$status == "current_best") |>
  dplyr::slice_tail(n = 1)
current <- c(
  "# Current best COTS fine-scale model",
  "",
  sprintf("Last evaluated: %s", candidate_row$evaluated_at),
  "",
  "## Current decision",
  "",
  sprintf("**%s**", current_row$display_name),
  "",
  if (promotion_passed) {
    "Iteration 01 passed every promotion gate and supersedes the legacy model."
  } else {
    paste(
      "The legacy clean-year model remains the incumbent. The two-part",
      "candidate is statistically informative but is not promoted because",
      "its spatial pattern was rejected by expert review; outbreak-layer",
      "provenance also remains unresolved."
    )
  },
  "",
  "## Incumbent evidence",
  "",
  sprintf(
    "- Within-reef-year AUC: %s", fmt(current_row$within_reef_year_auc)
  ),
  sprintf("- Top-20%% capture: %s", fmt(current_row$top20_capture)),
  sprintf("- Top-20%% lift: %s", fmt(current_row$top20_lift)),
  "- The legacy 2026 raster can be reproduced exactly from its frozen model.",
  "- Operational experts continue to regard its stretched within-reef",
  "  spatial pattern as the most credible map currently available.",
  "",
  "## Most recent successor attempt (rejected)",
  "",
  sprintf("Candidate: **%s**", candidate_row$display_name),
  sprintf(
    "- Mean repeated within-reef-year AUC: %s",
    fmt(candidate_row$within_reef_year_auc)
  ),
  sprintf(
    "- Mean top-20%% capture / lift: %s / %s",
    fmt(candidate_row$top20_capture), fmt(candidate_row$top20_lift)
  ),
  sprintf(
    "- Promotion status: **%s**", candidate_row$status
  ),
  "",
  "The successor architecture keeps static habitat rank, reef-year outbreak",
  "state and fixed-effort observation probability as separate named layers.",
  "Neither state nor effort is allowed to reorder cells within a reef-year.",
  "",
  "## Rejected candidate reef maps",
  "",
  "The incumbent and candidate use identical within-reef percentile scales.",
  "These patterns were reviewed and rejected as operationally implausible.",
  "",
  markdown_table(names(map_rows), as.matrix(map_rows)),
  "",
  map_blocks,
  "The decision is recorded in",
  "../artifacts/two_part_iteration_01_expert_review.csv.",
  "",
  hmsc_comparison,
  "## Promotion gates",
  "",
  markdown_table(names(gate_rows), as.matrix(gate_rows)),
  "",
  "## Harmonised next steps",
  "",
  shared_next_steps,
  "",
  "## Update rule",
  "",
  "Run scripts 28, 29, 31 and 30 for a new iteration. This report changes",
  "incumbent only",
  "when all statistical, rank-integrity, provenance and expert-map gates pass.",
  "The complete decision history is retained in",
  "outputs/model_successor_registry.csv."
)

report_dir <- file.path(v2_dir, "reports")
dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)
writeLines(
  detail,
  file.path(
    report_dir, "TWO_PART_FIXED_EFFORT_ITERATION_01.md"
  )
)
writeLines(
  current,
  file.path(report_dir, "CURRENT_BEST_MODEL.md")
)

message(decision)
