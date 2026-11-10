# Question: what supports the ranking comparison and experimental delivery?
# Method: paired habitat-summary contrasts, cohort/support counts and product
# changes, all from scripts 29-31. Validation: unchanged candidate coverage,
# constant class counts per reef and exact archive replay recorded by 29-30.
source(file.path(dirname(sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE))),'..','R','common.R'))
cfg <- initialise_final()
source(file.path(cfg$final_dir,'R','cull_sites.R'))
cache <- file.path(cfg$dirs$cache,'ranking_candidates')
a <- readRDS(file.path(cache,'assessment.rds'));obj <- readRDS(file.path(cache,'features.rds'))
read <- function(n) read.csv(file.path(cfg$dirs$tables,n))
diff <- read('ranking_candidate_differences.csv')
diff <- diff[!grepl('^h_',diff$model),]
extra <- list()
for(period in unique(a$results$period)) for(h in c('h_q50','h_q75','h_q90','h_q95')) for(metric in c('auc_high_all','auc_useful')) {
 key <- paste(h,period,'habitat_only',0,metric,sep='::')
 ref <- paste('h_mean',period,'habitat_only',0,metric,sep='::')
 v <- a$bootstrap_values[[key]]-a$bootstrap_values[[ref]]
 extra[[length(extra)+1L]] <- data.frame(model=h,period,rule='habitat_only',target=0,metric,reference='h_mean',
  value=v[1],low=unname(quantile(v[-1],.025)),high=unname(quantile(v[-1],.975)))
}
diff <- dplyr::bind_rows(diff,dplyr::bind_rows(extra));a$differences <- diff
write_table(diff,cfg,'ranking_candidate_differences.csv')
saveRDS(a,file.path(cache,'assessment.rds'))
cover <- read('ranking_candidate_habitat_coverage.csv')
points <- readRDS(checked_input(cfg,'ranking_site_sample_points'))
cell_data <- readRDS(file.path(cache,'full_polygon_cells.rds'))$cells
scope <- data.frame(metric=c('evaluated_dives','evaluated_reef_years','candidate_site_years','designated_sites',
 'sample_points','median_sample_points','polygon_cell_memberships','distinct_cells','median_full_cells',
 'half_life_days','habitat_quantile','sites_without_predictor_values','spatial_replay_max_error'),value=c(nrow(a$dives),length(a$groups),nrow(a$full),
 nrow(obj$sites),length(points$site),median(tabulate(points$site,nbins=nrow(obj$sites))),sum(cover$cells),
 length(unique(cell_data$cell[is.finite(cell_data$cell)])),median(cover$cells),cfg$ranking_candidates$half_life_days,
 cfg$ranking_candidates$habitat_quantile,sum(cover$valid_cells==0),max(read('ranking_candidate_habitat_replay.csv')$max_abs_error)))
write_table(scope,cfg,'ranking_candidate_scope.csv')
manifest <- read('ranking_candidate_product_manifest.csv')
new <- read.csv(file.path(cfg$dirs$outputs,'site_ranking_candidates',manifest$file[grepl('\\.csv$',manifest$file)]))
old <- read.csv(file.path(cfg$dirs$outputs,'site_ranking','cots_cull_site_priority_2027.csv'))
old$site_key <- paste(old$reef_name,clean_site_name(old$site_name),sep='||')
new$site_key <- paste(new$reef_name,clean_site_name(new$site_name),sep='||')
z <- dplyr::inner_join(new,old[,c('site_key','priority','rank_on_reef')],by='site_key',suffix=c('_candidate','_current'))
stopifnot(nrow(z)==nrow(obj$sites),!anyDuplicated(z$site_key))
changes <- z |> dplyr::group_by(.data$reef_name) |>
 dplyr::summarise(sites=dplyr::n(),priority_changed=sum(.data$priority_candidate!=.data$priority_current),
  crosses_D=sum((.data$priority_candidate=='D')!=(.data$priority_current=='D')),
  D_current=sum(.data$priority_current=='D'),D_candidate=sum(.data$priority_candidate=='D'),.groups='drop')
stopifnot(all(changes$D_current==changes$D_candidate),max(scope$value[scope$metric=='spatial_replay_max_error'])<1e-5)
write_table(changes,cfg,'ranking_candidate_product_changes.csv')
inputs <- c('assessment_dives','assessment_sites','assessment_evidence','assessment_manta_sites',
 'assessment_forward','assessment_forward_coefficients','ranking_train_site_years','ranking_site_sample_points',
 'predictor_stack_rg','ranking_deployed_habitat','manta_tows','footprint_links')
write_table(data.frame(input=inputs,sha256=vapply(inputs,function(n) cfg$inputs[[n]]$sha256,'')),cfg,'ranking_candidate_input_checksums.csv')
log_message('Candidate summary and product checks passed; ',sum(changes$crosses_D),' sites cross the deferred-class boundary.')
