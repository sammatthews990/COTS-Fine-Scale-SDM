# RosettaCOTS 1.0.2 contains conflicting prose about the SALAD input unit in
# the Cull-from-SALAD path. The posterior coefficients and the published
# 0.04 Cull CPUE ~= 20.6 COTS-plus-scars/ha result reconcile when the predictor
# is SALAD density per hectare. Override only the calibration density offset.
.subset_tmb_data_v6_unfixed <- subset_tmb_data_v6
subset_tmb_data_v6 <- function(artifact, model_type, cull_ids = NULL,
                               salad_ids = NULL, calibration_draw = NULL) {
  result <- .subset_tmb_data_v6_unfixed(
    artifact, model_type, cull_ids, salad_ids, calibration_draw
  )
  if (model_type == "joint_rosetta") {
    result$data$calibration_log_density_offset <- 0
  }
  result
}
