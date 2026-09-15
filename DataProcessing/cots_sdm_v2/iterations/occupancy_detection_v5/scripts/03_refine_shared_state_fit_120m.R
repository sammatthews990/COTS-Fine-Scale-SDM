args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
v5_dir <- normalizePath(file.path(dirname(script_file), ".."), winslash = "/")
artifact <- readRDS(file.path(v5_dir, "artifacts", "shared_state_data_120m.rds"))
initial <- readRDS(file.path(v5_dir, "models", "shared_state_track_120m.rds"))
source_dir <- file.path(v5_dir, "src")
old <- setwd(source_dir); on.exit(setwd(old), add = TRUE)
dll <- TMB::dynlib("shared_state_track")
dyn.load(dll); on.exit(try(dyn.unload(dll), silent = TRUE), add = TRUE)

parameters <- list(
  beta_state = rep(0, ncol(artifact$data$X)),
  alpha_direct = c(-1.5, -1.5, -1.5), effort_direct = rep(0.2, 3),
  depletion_post = -0.5, depletion_log_days = -0.1,
  alpha_scar = c(-1.5, -1.5), effort_scar = rep(0, 2),
  logit_fp_scar = c(-4, -4), scar_prior_effect = 0.5,
  scar_decay = -0.1, scar_post_cull_effect = 0,
  alpha_manta_direct = -2, effort_manta_direct = 0.2,
  alpha_manta_scar = -2, effort_manta_scar = 0.2,
  logit_fp_manta_scar = -4
)
objective <- TMB::MakeADFun(artifact$data, parameters, DLL = "shared_state_track", silent = TRUE)
started <- Sys.time()
fit <- stats::nlminb(
  initial$fit$par, objective$fn, objective$gr,
  control = list(iter.max = 3000, eval.max = 5000, rel.tol = 1e-12, x.tol = 1e-10)
)
gradient <- objective$gr(fit$par)
sdreport <- TMB::sdreport(objective, par.fixed = fit$par)
fixed <- as.data.frame(summary(sdreport, "fixed"))
parameter_table <- tibble::tibble(
  parameter = names(fit$par), estimate = as.numeric(fit$par),
  standard_error = fixed[["Std. Error"]]
)
state_rows <- grepl("^beta_state", parameter_table$parameter)
parameter_table$component <- ifelse(state_rows, "state", "observation")
parameter_table$term <- parameter_table$parameter
parameter_table$term[state_rows] <- artifact$design_columns[seq_len(sum(state_rows))]
state_par <- fit$par[state_rows]
psi <- stats::plogis(as.numeric(artifact$data$X %*% state_par))
sites <- dplyr::mutate(artifact$sites, shared_state_probability = psi)
elapsed <- as.numeric(difftime(Sys.time(), started, units = "mins"))

model <- list(
  fit = fit, parameter_table = parameter_table, sites = sites,
  design = initial$design, likelihood_note = artifact$likelihood_note,
  max_absolute_gradient = max(abs(gradient)), pdHess = sdreport$pdHess,
  fitted_at = Sys.time()
)
saveRDS(model, file.path(v5_dir, "models", "shared_state_track_refined_120m.rds"))
saveRDS(sites, file.path(v5_dir, "outputs", "shared_state_site_probabilities_refined_120m.rds"))
utils::write.csv(parameter_table, file.path(v5_dir, "outputs", "shared_state_parameters_refined_120m.csv"), row.names = FALSE)
audit <- tibble::tibble(
  support_m = 120L, convergence_code = fit$convergence,
  convergence_message = fit$message, objective = fit$objective,
  max_absolute_gradient = max(abs(gradient)),
  positive_definite_hessian = sdreport$pdHess, elapsed_minutes = elapsed,
  psi_min = min(psi), psi_q05 = stats::quantile(psi, 0.05),
  psi_median = stats::median(psi), psi_q95 = stats::quantile(psi, 0.95),
  psi_max = max(psi), psi_sd = stats::sd(psi)
)
utils::write.csv(audit, file.path(v5_dir, "outputs", "shared_state_fit_audit_refined_120m.csv"), row.names = FALSE)
print(audit)
