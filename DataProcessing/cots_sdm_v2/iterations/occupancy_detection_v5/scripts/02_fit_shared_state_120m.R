args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
v5_dir <- normalizePath(file.path(dirname(script_file), ".."), winslash = "/")
artifact_dir <- file.path(v5_dir, "artifacts")
model_dir <- file.path(v5_dir, "models")
output_dir <- file.path(v5_dir, "outputs")
source_dir <- file.path(v5_dir, "src")
dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

artifact <- readRDS(file.path(artifact_dir, "shared_state_data_120m.rds"))
cpp <- file.path(source_dir, "shared_state_track.cpp")
old <- setwd(source_dir)
on.exit(setwd(old), add = TRUE)
dll <- TMB::dynlib("shared_state_track")
if (!file.exists(dll) || file.mtime(dll) < file.mtime(cpp)) {
  TMB::compile(basename(cpp), flags = "-O2")
}
dyn.load(dll)
on.exit(try(dyn.unload(dll), silent = TRUE), add = TRUE)

parameters <- list(
  beta_state = rep(0, ncol(artifact$data$X)),
  alpha_direct = c(-1.5, -1.5, -1.5),
  effort_direct = rep(0.2, 3),
  depletion_post = -0.5,
  depletion_log_days = -0.1,
  alpha_scar = c(-1.5, -1.5),
  effort_scar = rep(0, 2),
  logit_fp_scar = c(-4, -4),
  scar_prior_effect = 0.5,
  scar_decay = -0.1,
  scar_post_cull_effect = 0,
  alpha_manta_direct = -2,
  effort_manta_direct = 0.2,
  alpha_manta_scar = -2,
  effort_manta_scar = 0.2,
  logit_fp_manta_scar = -4
)
objective <- TMB::MakeADFun(
  data = artifact$data, parameters = parameters,
  DLL = "shared_state_track", silent = TRUE
)
started <- Sys.time()
fit <- stats::nlminb(
  objective$par, objective$fn, objective$gr,
  control = list(iter.max = 1000, eval.max = 1500, rel.tol = 1e-9)
)
elapsed <- as.numeric(difftime(Sys.time(), started, units = "mins"))
report <- objective$report(fit$par)
sdreport <- tryCatch(TMB::sdreport(objective, par.fixed = fit$par), error = identity)
gradient <- objective$gr(fit$par)

parameter_table <- tibble::tibble(
  parameter = names(fit$par), estimate = as.numeric(fit$par)
)
if (!inherits(sdreport, "error")) {
  fixed <- as.data.frame(summary(sdreport, "fixed")) |>
    tibble::rownames_to_column("parameter") |>
    dplyr::rename(estimate_sdreport = "Estimate", standard_error = "Std. Error")
  parameter_table <- dplyr::left_join(parameter_table, fixed, by = "parameter")
}
state_rows <- grepl("^beta_state", parameter_table$parameter)
state_names <- artifact$design_columns
parameter_table$component <- "observation"
parameter_table$component[state_rows] <- "state"
parameter_table$term <- parameter_table$parameter
parameter_table$term[state_rows] <- state_names[seq_len(sum(state_rows))]

psi <- as.numeric(report$psi)
sites <- artifact$sites |>
  dplyr::mutate(shared_state_probability = psi)
saveRDS(
  list(
    fit = fit, report = report,
    sdreport_summary = if (inherits(sdreport, "error")) NULL else summary(sdreport),
    parameter_table = parameter_table, sites = sites,
    design = list(
      pca = artifact$pca, base_terms = artifact$state_base_terms,
      columns = artifact$design_columns, center = artifact$design_center,
      scale = artifact$design_scale
    ),
    likelihood_note = artifact$likelihood_note,
    fitted_at = Sys.time()
  ),
  file.path(model_dir, "shared_state_track_120m.rds")
)
utils::write.csv(parameter_table, file.path(output_dir, "shared_state_parameters_120m.csv"), row.names = FALSE)
saveRDS(sites, file.path(output_dir, "shared_state_site_probabilities_120m.rds"))

audit <- tibble::tibble(
  support_m = 120L,
  convergence_code = fit$convergence,
  convergence_message = fit$message,
  objective = fit$objective,
  max_absolute_gradient = max(abs(gradient)),
  positive_definite_hessian = if (inherits(sdreport, "error")) FALSE else isTRUE(sdreport$pdHess),
  elapsed_minutes = elapsed,
  psi_min = min(psi), psi_q05 = stats::quantile(psi, 0.05),
  psi_median = stats::median(psi), psi_q95 = stats::quantile(psi, 0.95),
  psi_max = max(psi), psi_sd = stats::sd(psi)
)
utils::write.csv(audit, file.path(output_dir, "shared_state_fit_audit_120m.csv"), row.names = FALSE)
print(audit)
print(dplyr::filter(parameter_table, .data$component == "observation"))
