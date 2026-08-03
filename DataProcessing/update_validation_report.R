# R script to cleanly integrate MaxEnt model evaluation into SDM_Validation_Report.qmd

# Restore original file first
system("git checkout DataProcessing/SDM_Validation_Report.qmd")

qmd_path <- "DataProcessing/SDM_Validation_Report.qmd"
content <- readLines(qmd_path)

# 1. Add predict_maxent_path to configuration block
idx_cfg <- grep('predict_year2025_rg_path <-', content, fixed = TRUE)
if (length(idx_cfg) > 0 && !any(grepl('predict_maxent_path', content, fixed = TRUE))) {
  content <- append(content, 'predict_maxent_path    <- "outputs/COTS_maxent_suitability.tif"', after = idx_cfg[1])
}

# 2. Add prob_maxent to load_data block
idx_load <- grep('prob_2025_rg  <-', content, fixed = TRUE)
if (length(idx_load) > 0 && !any(grepl('prob_maxent   <-', content, fixed = TRUE))) {
  content <- append(content, 'prob_maxent   <- if(file.exists(predict_maxent_path))    terra::rast(predict_maxent_path) else NULL', after = idx_load[1])
}

# 3. Add crop_maxent and thresholds
idx_thresh <- grep('median_thresh_agnostic <-', content, fixed = TRUE)
if (length(idx_thresh) > 0 && !any(grepl('crop_maxent <-', content, fixed = TRUE))) {
  add_lines <- c(
    'crop_maxent <- if(!is.null(prob_maxent)) terra::crop(prob_maxent, ext_val) else NULL',
    'median_thresh_maxent <- if(!is.null(crop_maxent)) as.numeric(terra::global(crop_maxent, fun=function(x) median(x, na.rm=TRUE))) else NA',
    'thresh_maxent_p88    <- if(!is.null(crop_maxent)) as.numeric(terra::global(crop_maxent, fun=function(x) quantile(x, probs=0.88, na.rm=TRUE))) else NA',
    'thresh_maxent_p95    <- if(!is.null(crop_maxent)) as.numeric(terra::global(crop_maxent, fun=function(x) quantile(x, probs=0.95, na.rm=TRUE))) else NA'
  )
  content <- append(content, add_lines, after = idx_thresh[1])
}

# Update cat() output for thresholds
idx_cat_ag <- grep('cat("Median Threshold (Agnostic):', content, fixed = TRUE)
if (length(idx_cat_ag) > 0 && !any(grepl('Median Threshold (MaxEnt)', content, fixed = TRUE))) {
  add_cat_lines <- c(
    'cat("Median Threshold (MaxEnt):", round(median_thresh_maxent, 3), "\\n")',
    'cat("88th Percentile (MaxEnt):", round(thresh_maxent_p88, 3), "\\n")',
    'cat("95th Percentile (MaxEnt):", round(thresh_maxent_p95, 3), "\\n")'
  )
  content <- append(content, add_cat_lines, after = idx_cat_ag[1])
}

# 4. Add Section 1.4 MaxEnt Model Fitting Performance AFTER chunk model_fitting_metrics
idx_sec1_start <- grep('```{r model_fitting_metrics', content, fixed = TRUE)
idx_sec1_end <- grep('```', content, fixed = TRUE)
idx_chunk1_end <- idx_sec1_end[idx_sec1_end > idx_sec1_start[1]][1]

if (length(idx_chunk1_end) > 0 && !any(grepl('### 1.4 MaxEnt Presence-Only Model', content, fixed = TRUE))) {
  sec14_block <- c(
    '',
    '### 1.4 MaxEnt Presence-Only Model (ENMeval Tuning & Feature Importance)',
    '',
    '```{r maxent_eval_summary, results=\'asis\'}',
    'maxent_rds <- "outputs/maxent_enmeval_results.rds"',
    'if (file.exists(maxent_rds)) {',
    '  me_res <- readRDS(maxent_rds)',
    '  res_df <- me_res@results %>%',
    '    select(fc, rm, auc.train, auc.val.avg, auc.val.sd, AICc) %>%',
    '    arrange(AICc)',
    '  ',
    '  best_row <- res_df[1, ]',
    '  cat(paste0("**Optimal MaxEnt Configuration selected by AICc**:\\n",',
    '             "- **Feature Classes (fc)**: ", best_row$fc, "\\n",',
    '             "- **Regularization Multiplier (rm)**: ", best_row$rm, "\\n",',
    '             "- **AICc**: ", round(best_row$AICc, 2), "\\n",',
    '             "- **Validation AUC (4-fold Spatial Block CV)**: ", round(best_row$auc.val.avg, 3), " ± ", round(best_row$auc.val.sd, 3), "\\n\\n"))',
    '  ',
    '  print(kable(head(res_df, 10), digits = 3, caption = "Top 10 MaxEnt Candidate Model Configurations (ENMeval)"))',
    '  cat("\\n\\n")',
    '} else {',
    '  cat("*MaxEnt ENMeval tuning results not found.*")',
    '}',
    '```',
    '',
    '```{r maxent_plots, echo=FALSE, fig.show=\'hold\', out.width=\'48%\'}',
    'vip_img <- "outputs/maxent_vip_cots.png"',
    'rc_img  <- "outputs/maxent_response_curves.png"',
    '',
    'if(file.exists(vip_img)) knitr::include_graphics(vip_img)',
    'if(file.exists(rc_img))  knitr::include_graphics(rc_img)',
    '```'
  )
  content <- append(content, sec14_block, after = idx_chunk1_end)
}

# 5. Grid-based 100m evaluation: Add MaxEnt raster
idx_val_df <- grep('prob_agnostic = terra::values(pred_60m_agn)[,1]', content, fixed = TRUE)
if (length(idx_val_df) > 0 && !any(grepl('prob_maxent = terra::values(pred_60m_max)[,1]', content, fixed = TRUE))) {
  idx_resample <- grep('pred_60m_agn <- terra::resample', content, fixed = TRUE)
  if (length(idx_resample) > 0) {
    content <- append(content, '  pred_60m_max <- if(!is.null(crop_maxent)) terra::resample(crop_maxent, crop_2025, method = "bilinear") else NULL', after = idx_resample[1])
  }
  idx_val_df2 <- grep('prob_agnostic = terra::values(pred_60m_agn)[,1]', content, fixed = TRUE)
  content[idx_val_df2[1]] <- '    prob_agnostic = terra::values(pred_60m_agn)[,1],'
  content <- append(content, '    prob_maxent   = if(!is.null(pred_60m_max)) terra::values(pred_60m_max)[,1] else NA', after = idx_val_df2[1])
}

# Add ROC and AUC for MaxEnt in grid evaluation
idx_roc_agn <- grep('roc_agn <- roc(val_df$observed_presence, val_df$prob_agnostic, quiet=TRUE)', content, fixed = TRUE)
if (length(idx_roc_agn) > 0 && !any(grepl('roc_maxent <-', content, fixed = TRUE))) {
  content <- append(content, '  roc_maxent <- if("prob_maxent" %in% names(val_df)) roc(val_df$observed_presence, val_df$prob_maxent, quiet=TRUE) else NULL', after = idx_roc_agn[1])
  idx_cat_auc <- grep('cat("AUC (Agnostic Model):', content, fixed = TRUE)
  if (length(idx_cat_auc) > 0) {
    content <- append(content, '  if(!is.null(roc_maxent)) cat("AUC (MaxEnt Model):", round(auc(roc_maxent), 3), "\\n")', after = idx_cat_auc[1])
  }
}

# 6. Cull Dive Validation: Extract MaxEnt probability
idx_cull_ext_agn <- grep('cull_extract_agn <- terra::extract(prob_agnostic, terra::vect(cull_pts_2025))', content, fixed = TRUE)
if (length(idx_cull_ext_agn) > 0 && !any(grepl('cull_extract_max <-', content, fixed = TRUE))) {
  add_cull_max_lines <- c(
    '# Extract MaxEnt model probability',
    'cull_extract_max <- if(!is.null(prob_maxent)) terra::extract(prob_maxent, terra::vect(cull_pts_2025)) else NULL',
    'val_col_max <- if(!is.null(cull_extract_max)) setdiff(names(cull_extract_max), "ID")[1] else NULL'
  )
  content <- append(content, add_cull_max_lines, after = idx_cull_ext_agn[1])
}

# Add prob_maxent to cull_val_df
idx_cull_df_agn <- grep('prob_agnostic = cull_extract_agn[[val_col_agn]]', content, fixed = TRUE)
if (length(idx_cull_df_agn) > 0 && !any(grepl('prob_maxent = cull_extract_max', content, fixed = TRUE))) {
  content[idx_cull_df_agn[1]] <- '  prob_agnostic = cull_extract_agn[[val_col_agn]],'
  content <- append(content, '  prob_maxent   = if(!is.null(cull_extract_max)) cull_extract_max[[val_col_max]] else NA', after = idx_cull_df_agn[1])
}

# Add roc_cull_maxent
idx_roc_cull_agn <- grep('roc_cull_agn  <- roc(cull_val_df$observed_presence, cull_val_df$prob_agnostic, quiet = TRUE)', content, fixed = TRUE)
if (length(idx_roc_cull_agn) > 0 && !any(grepl('roc_cull_max  <-', content, fixed = TRUE))) {
  content <- append(content, 'roc_cull_max  <- if("prob_maxent" %in% names(cull_val_df)) roc(cull_val_df$observed_presence, cull_val_df$prob_maxent, quiet = TRUE) else NULL', after = idx_roc_cull_agn[1])
}

# Add roc_df_max to roc_df
idx_roc_df_agn_end <- grep('Model = paste0("Agnostic (AUC = ", round(auc(roc_cull_agn), 3), ")")', content, fixed = TRUE)
if (length(idx_roc_df_agn_end) > 0 && !any(grepl('roc_df_max', content, fixed = TRUE))) {
  # find closing parenthesis line after idx_roc_df_agn_end
  close_paren_idx <- idx_roc_df_agn_end[1] + 1
  roc_max_lines <- c(
    'roc_df_max <- if(!is.null(roc_cull_max)) data.frame(',
    '  FPR = 1 - roc_cull_max$specificities,',
    '  TPR = roc_cull_max$sensitivities,',
    '  Model = paste0("MaxEnt (AUC = ", round(auc(roc_cull_max), 3), ")")',
    ') else NULL'
  )
  content <- append(content, roc_max_lines, after = close_paren_idx)
  
  idx_bind_roc <- grep('roc_df <- bind_rows(roc_df_2025, roc_df_agn)', content, fixed = TRUE)
  if (length(idx_bind_roc) > 0) {
    content[idx_bind_roc[1]] <- 'roc_df <- bind_rows(roc_df_2025, roc_df_agn, roc_df_max)'
  }
}

# Add MaxEnt to cull_roc_summary and cull_threshold_summary
idx_cull_roc_sum <- grep('summarise_roc(roc_cull_agn, "Agnostic")', content, fixed = TRUE)
if (length(idx_cull_roc_sum) > 0 && !any(grepl('summarise_roc(roc_cull_max, "MaxEnt")', content, fixed = TRUE))) {
  content[idx_cull_roc_sum[1]] <- '  summarise_roc(roc_cull_agn, "Agnostic"),'
  content <- append(content, '  if(!is.null(roc_cull_max)) summarise_roc(roc_cull_max, "MaxEnt")', after = idx_cull_roc_sum[1])
}

idx_cull_thresh_sum <- grep('calc_binary_metrics(cull_val_df$prob_agnostic, cull_val_df$observed_presence, opt_agn$threshold, "Agnostic", "Youden")', content, fixed = TRUE)
if (length(idx_cull_thresh_sum) > 0 && !any(grepl('"MaxEnt"', content, fixed = TRUE))) {
  opt_max_lines <- c(
    'if(!is.null(roc_cull_max)) {',
    '  opt_max <- coords(roc_cull_max, "best", best.method = "youden", best.policy = "first", ret = c("threshold", "sensitivity", "specificity"))',
    '  cull_threshold_summary <- bind_rows(cull_threshold_summary, calc_binary_metrics(cull_val_df$prob_maxent, cull_val_df$observed_presence, opt_max$threshold, "MaxEnt", "Youden"))',
    '  if(!is.na(thresh_maxent_p88)) cull_threshold_summary <- bind_rows(cull_threshold_summary, calc_binary_metrics(cull_val_df$prob_maxent, cull_val_df$observed_presence, thresh_maxent_p88, "MaxEnt", "P88"))',
    '  if(!is.na(thresh_maxent_p95)) cull_threshold_summary <- bind_rows(cull_threshold_summary, calc_binary_metrics(cull_val_df$prob_maxent, cull_val_df$observed_presence, thresh_maxent_p95, "MaxEnt", "P95"))',
    '}'
  )
  content <- append(content, opt_max_lines, after = idx_cull_thresh_sum[1])
}

writeLines(content, qmd_path)
cat("Successfully updated SDM_Validation_Report.qmd with MaxEnt integration!\n")
