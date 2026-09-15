fit_density_weighted_tmb_v6 <- function(
    artifact, dll_name, cull_ids = NULL, salad_ids = NULL, start = NULL) {
  input <- subset_tmb_data_v6(
    artifact, "joint_rosetta", cull_ids, salad_ids, calibration_draw = NULL
  )
  n_cull <- nrow(input$cull)
  n_salad <- nrow(input$salad)
  # Give each monitoring method equal total likelihood weight. This is a
  # pre-specified sensitivity analysis, not a deployable sampling correction.
  cull_weight <- if (n_cull > 0L && n_salad > 0L) n_salad / n_cull else 1
  input$data$cull_likelihood_weight <- cull_weight
  input$data$salad_likelihood_weight <- 1

  p <- ncol(artifact$design$X)
  if (is.null(start)) {
    beta <- rep(0, p)
    beta[[1]] <- log(max(
      stats::median(input$salad$count / input$salad$effort_area_ha), 1
    ))
    start <- list(beta = beta, log_phi_cull = log(5), log_phi_salad = log(5))
  }
  objective <- TMB::MakeADFun(input$data, start, DLL = dll_name, silent = TRUE)
  fit <- stats::nlminb(
    objective$par, objective$fn, objective$gr,
    control = list(iter.max = 1500, eval.max = 2500, rel.tol = 1e-10)
  )
  gradient <- objective$gr(fit$par)
  hessian <- try(stats::optimHess(fit$par, objective$fn, objective$gr), silent = TRUE)
  pd_hessian <- !inherits(hessian, "try-error") && all(is.finite(hessian)) &&
    all(eigen(hessian, symmetric = TRUE, only.values = TRUE)$values > 0)
  list(
    model_type = "joint_rosetta", candidate = "joint_rosetta_balanced",
    fit = fit, beta = fit$par[grepl("^beta", names(fit$par))],
    max_gradient = max(abs(gradient)), pd_hessian = pd_hessian,
    calibration_draw = NULL,
    likelihood_weights = c(cull = cull_weight, salad = 1),
    input_counts = c(cull = n_cull, salad = n_salad)
  )
}
