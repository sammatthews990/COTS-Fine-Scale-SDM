# Question: is a candidate worth a separate field-comparison layer?
# Method: apply the pre-specified strong gate; select the best positive recent
# and overall comparison as experimental only, then fit on the deployment years.
# Validation: preserve incumbent hashes, retain all designated sites, verify
# temporal bounds, probabilities and ranks; reread the exported GeoPackage.
source(file.path(dirname(sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE))),'..','R','common.R'))
cfg <- initialise_final()
source(file.path(cfg$final_dir,'R','cull_sites.R'))
source(file.path(cfg$final_dir,'R','site_selection.R'))
source(file.path(cfg$final_dir,'R','ranking_candidates.R'))
source(file.path(cfg$final_dir,'R','forward.R'))
cache <- file.path(cfg$dirs$cache,'ranking_candidates')
assessment <- readRDS(file.path(cache,'assessment.rds'))
obj <- readRDS(file.path(cache,'features.rds'))
diff <- assessment$differences
delta <- function(model,period,metric,reference='archived',rule='matched_high',target=.9,col='value') {
 r <- diff[diff$model==model & diff$period==period & diff$metric==metric &
  diff$reference==reference & diff$rule==rule & abs(diff$target-target)<1e-9,]
 stopifnot(nrow(r)==1);r[[col]]
}
decisions <- lapply(setdiff(candidate_model_names,'mean_control'),function(model) {
 primary <- vapply(c('archived','mean_control'),function(ref) delta(model,'2019-2026','dives_avoided',ref),numeric(1))
 lower <- vapply(c('archived','mean_control'),function(ref) delta(model,'2019-2026','dives_avoided',ref,col='low'),numeric(1))
 recent <- vapply(c('archived','mean_control'),function(ref) delta(model,'2024-2026','dives_avoided',ref),numeric(1))
 useful <- vapply(c('archived','mean_control'),function(ref) delta(model,'2019-2026','useful_retained',ref,'A-C',.8),numeric(1))
 data.frame(model,min_primary_gain=min(primary),min_primary_lower=min(lower),
  min_recent_gain=min(recent),min_useful_retention_gain=min(useful),
  strong_gate_met=all(primary>=.01 & lower>0 & recent>=0) && all(useful>=0),
  eligible_experiment=all(primary>0 & recent>=0 & useful>=0))
})
decisions <- dplyr::bind_rows(decisions)
eligible <- decisions[decisions$eligible_experiment,]
decisions$exported_experiment <- FALSE
if(!nrow(eligible)) {
 write_table(decisions,cfg,'ranking_candidate_decisions.csv')
 log_message('No candidate supports even a positive experimental comparison; no layer exported.');quit(save='no')
}
model <- eligible$model[which.max(eligible$min_primary_gain)]
decisions$exported_experiment <- decisions$model==model
write_table(decisions,cfg,'ranking_candidate_decisions.csv')
log_message('Experimental export: ',model,'. The incumbent remains the default.')
incumbent <- file.path(cfg$dirs$outputs,'site_ranking',c('cots_cull_site_priority_2027.csv','cots_cull_site_priority_2027.kmz'))
before <- vapply(incumbent,sha256_file,'')
sites <- read_designated_sites(cfg);sites$site_row <- seq_len(nrow(sites))
d <- obj$dives;sr <- cfg$site_ranking
start <- sr$evidence_end-sr$evidence_days+1;end <- sr$evidence_end+1
recent <- d[d$survey_date>=start & d$survey_date<end,]
site_xy <- sf::st_coordinates(suppressWarnings(sf::st_point_on_surface(sites)))[,1:2]
has <- is.finite(d$longitude) & is.finite(d$latitude)
xy <- matrix(NA_real_,nrow(d),2)
xy[has,] <- sf::st_coordinates(sf::st_transform(sf::st_as_sf(d[has,],coords=c('longitude','latitude'),crs=4326),sf::st_crs(sites)))
fill <- !has & !is.na(d$site_row);xy[fill,] <- site_xy[d$site_row[fill],]
ri <- which(d$survey_date>=start & d$survey_date<end & is.finite(xy[,1]))
ck <- history_kernel(xy[ri,,drop=FALSE],as.numeric(d$outcome[ri]=='problematic'),rep(0,length(ri)),site_xy,sr$kernel_sigma_m,Inf)
tows <- readRDS(checked_input(cfg,'manta_tows'))[,c('visit_id','survey_date','manta_evidence_class')]
geom <- readRDS(checked_input(cfg,'footprint_links'))$manta |>
 dplyr::group_by(.data$visit_id) |>
 dplyr::summarise(x=(dplyr::first(.data$start_x)+dplyr::first(.data$end_x))/2,
  y=(dplyr::first(.data$start_y)+dplyr::first(.data$end_y))/2,.groups='drop')
tows <- dplyr::inner_join(tows,geom,by='visit_id') |>
 dplyr::filter(.data$survey_date>=start,.data$survey_date<end,is.finite(.data$x),is.finite(.data$y))
strong <- tows$manta_evidence_class %in% c('direct_cots','multiple_counted_scars','legacy_common_scars')
mk <- history_kernel(cbind(tows$x,tows$y),as.numeric(strong),rep(0,nrow(tows)),site_xy,sr$kernel_sigma_m,Inf)
at_site <- recent[!is.na(recent$site_row),] |>
 dplyr::mutate(w=2^(-as.numeric(end-.data$survey_date)/cfg$ranking_candidates$half_life_days)) |>
 dplyr::group_by(.data$site_row) |>
 dplyr::summarise(site_problematic_prev=as.numeric(any(.data$cpue>=.02)),
  site_high_prev=as.numeric(any(.data$cpue>.04)),site_quiet_prev=as.numeric(!any(.data$cpue>=.02)),
  recent_cpue=sum(.data$w*.data$total_cots)/sum(.data$w*.data$bottom_time),
  recent_weight=sum(.data$w),recent_dives=dplyr::n(),.groups='drop')
f <- data.frame(site_row=sites$site_row,reef_name=sites$reef_name,site_name=sites$site_name,
 cull_rate=ifelse(ck$W>0,ck$H,0),cull_effort=ck$W,
 manta_strong=ifelse(mk$W>0,mk$H*mk$W,0),manta_effort=mk$W) |>
 dplyr::left_join(at_site,by='site_row') |>
 dplyr::mutate(dplyr::across(dplyr::all_of(c('site_problematic_prev','site_high_prev','site_quiet_prev','recent_cpue','recent_weight','recent_dives')),
  ~dplyr::coalesce(.x,0)))
cells <- readRDS(file.path(cache,'full_polygon_cells.rds'))$cells
grid <- terra::rast(checked_input(cfg,'predictor_stack_rg'))[[1]]
layer <- terra::rast(checked_input(cfg,'ranking_deployed_habitat'))[[1]]
ids <- sort(unique(cells$cell))
if(terra::compareGeom(grid,layer,stopOnError=FALSE)) {
 v <- terra::extract(layer,ids,raw=TRUE)[,1]
} else {
 stopifnot(sf::st_crs(terra::crs(grid))==sf::st_crs(terra::crs(layer)))
 v <- terra::extract(layer,terra::xyFromCell(grid,ids))[,1]
}
habitat <- split(v[match(cells$cell,ids)],cells$ID)
mean_h <- vapply(habitat,function(z) if(any(is.finite(z))) mean(z,na.rm=TRUE) else NA_real_,numeric(1))
tail_h <- vapply(habitat,function(z) if(any(is.finite(z))) unname(quantile(z,.9,na.rm=TRUE)) else NA_real_,numeric(1))
f$mean <- mean_h[match(f$site_row,as.integer(names(mean_h)))];f$q90 <- tail_h[match(f$site_row,as.integer(names(tail_h)))]
rank_h <- function(x) {ok <- is.finite(x);z <- rep(NA_real_,length(x));z[ok] <- (rank(x[ok],ties.method='average')-1)/max(sum(ok)-1,1);z}
f <- f |> dplyr::group_by(.data$reef_name) |>
 dplyr::mutate(h_mean=rank_h(.data$mean),h_q90=rank_h(.data$q90)) |> dplyr::ungroup()
train <- obj$dives[obj$dives$year>=sr$training_years[1] & obj$dives$year<=sr$training_years[2] & !is.na(obj$dives$site_row),] |>
 dplyr::inner_join(obj$features,by=c('reef_name','year','site_row'),relationship='many-to-one',suffix=c('','_feature'))
x <- candidate_terms(train,model);x$problematic <- as.numeric(train$outcome=='problematic')
fit <- glm(candidate_formula(model),family=binomial(),data=x);stopifnot(fit$converged)
f$score <- as.numeric(predict(fit,candidate_terms(f,model),type='link'))
f$p_useful_if_dived <- plogis(f$score)
f <- f |> dplyr::group_by(.data$reef_name) |>
 dplyr::arrange(dplyr::desc(.data$score),dplyr::desc(.data$h_mean),.data$site_row,.by_group=TRUE) |>
 dplyr::mutate(rank_on_reef=dplyr::row_number(),sites_on_reef=dplyr::n(),
  top_fraction=.data$rank_on_reef/.data$sites_on_reef,
  priority=dplyr::case_when(.data$top_fraction<=.2~'A',.data$top_fraction<=.5~'B',.data$top_fraction<=.8~'C',TRUE~'D')) |>
 dplyr::ungroup()
f$model <- model;f$status <- 'EXPERIMENTAL; incumbent remains default'
f$evidence_start <- as.character(start);f$evidence_end <- as.character(sr$evidence_end)
stopifnot(nrow(f)==nrow(sites),all(is.finite(f$score)),all(f$p_useful_if_dived>=0 & f$p_useful_if_dived<=1))
out <- file.path(cfg$dirs$outputs,'site_ranking_candidates');dir.create(out,showWarnings=FALSE)
stem <- paste0('cots_cull_site_priority_2027_',model,'_experimental')
write.csv(f,file.path(out,paste0(stem,'.csv')),row.names=FALSE,na='')
spatial <- dplyr::left_join(sites[,c('site_row')],f,by='site_row')
gpkg <- file.path(out,paste0(stem,'.gpkg'))
sf::st_write(spatial,gpkg,layer='experimental_priority',delete_layer=TRUE,quiet=TRUE)
check <- sf::st_read(gpkg,quiet=TRUE);stopifnot(nrow(check)==nrow(sites))
coef_table <- data.frame(model,term=names(coef(fit)),estimate=unname(coef(fit)),training_dives=nrow(train),
 training_years=paste(sr$training_years,collapse='-'),status='experimental')
coef_table$input_definition <- ifelse(coef_table$term=='site_problematic_prev',
 'Internal alias for site_high_prev: any prior dive with CPUE >0.04',coef_table$term)
write_table(coef_table,cfg,'ranking_candidate_2027_coefficients.csv')
saveRDS(fit,file.path(cache,paste0('2027_',model,'_glm.rds')))

# Preserve the source KMZ geometry and pop-ups, adding experimental fields.
source_dir <- file.path(cache,'kmz_source');dir.create(source_dir,showWarnings=FALSE)
unzip(checked_input(cfg,'cull_sites_kmz'),exdir=source_dir)
kml_file <- list.files(source_dir,pattern='\\.kml$',recursive=TRUE,full.names=TRUE)[1]
kml <- paste(readLines(kml_file,encoding='UTF-8',warn=FALSE),collapse='\n')
parts <- strsplit(kml,'<Placemark>',fixed=TRUE)[[1]];pm <- parts[-1]
field <- function(x,name) {
 p <- regexpr(paste0('<Data name="',name,'">\\s*<value>[^<]*</value>'),x)
 z <- rep(NA_character_,length(x));z[p>0] <- sub('.*<value>([^<]*)</value>','\\1',regmatches(x,p));z
}
unescape <- function(x) gsub('&amp;','&',gsub('&lt;','<',gsub('&gt;','>',gsub('&apos;',"'",gsub('&quot;','"',x)))))
key <- paste(unescape(field(pm,'Reef')),clean_site_name(unescape(field(pm,'Name'))),sep='||')
index <- match(key,paste(f$reef_name,clean_site_name(f$site_name),sep='||'))
colours <- c(A='cc4849e3',B='cc3468eb',C='9900a1ed',D='4d818789')
for(i in which(!is.na(index))) {
 j <- index[i]
 pm[i] <- sub('<color>33669900</color>',paste0('<color>',colours[f$priority[j]],'</color>'),pm[i],fixed=TRUE)
 data <- paste0('<Data name="Experimental priority"><value>',f$priority[j],'</value></Data>',
  '<Data name="Candidate rank"><value>',f$rank_on_reef[j],' of ',f$sites_on_reef[j],'</value></Data>',
  '<Data name="Candidate model"><value>',model,'; NOT THE DEFAULT</value></Data>',
  '<Data name="Recency weighted CPUE"><value>',formatC(f$recent_cpue[j],digits=4,format='f'),'</value></Data>')
 pm[i] <- sub('</ExtendedData>',paste0(data,'</ExtendedData>'),pm[i],fixed=TRUE)
}
new <- paste0(parts[1],paste0('<Placemark>',pm,collapse=''))
new <- sub('<name>Eotr_CotsCullSites_[^<]*</name>',paste0('<name>EXPERIMENTAL 2027: ',model,'</name>'),new)
balloon <- paste0(vapply(c('Experimental priority','Candidate rank','Candidate model','Recency weighted CPUE'),
 function(n) paste0('<tr><td>',n,'</td><td>$[',n,']</td></tr>'),''),collapse='')
new <- gsub('</tbody>',paste0(balloon,'</tbody>'),new,fixed=TRUE)
stage <- file.path(cache,'candidate_kmz');dir.create(stage,showWarnings=FALSE)
writeLines(new,file.path(stage,'doc.kml'),useBytes=TRUE)
kmz <- file.path(out,paste0(stem,'.kmz'))
if(file.exists(kmz)) file.remove(kmz)
zip_status <- system2('python',c(shQuote(file.path(cfg$final_dir,'scripts','package_site_kmz.py')),
 shQuote(file.path(stage,'doc.kml')),shQuote(kmz)))
stopifnot(zip_status==0)
stopifnot(file.exists(kmz),all(vapply(incumbent,sha256_file,'')==before))
manifest <- data.frame(model,status='experimental, not default',file=list.files(out,pattern=paste0('^',stem),full.names=FALSE),
 evidence_start=as.character(start),evidence_end=as.character(sr$evidence_end),half_life_days=cfg$ranking_candidates$half_life_days,
 strong_gate_met=decisions$strong_gate_met[decisions$model==model],kmz_match_fraction=mean(!is.na(index)))
manifest$sha256 <- vapply(file.path(out,manifest$file),sha256_file,'')
write_table(manifest,cfg,'ranking_candidate_product_manifest.csv')
gate_text <- if(any(decisions$strong_gate_met)) 'A candidate passed the statistical gate; field validation is still required.' else
 'No candidate passed the pre-specified strong replacement gate. Use this layer for comparison and field trials.'
writeLines(c('Experimental cull-site ranking; the incumbent remains the default.',
 'Selection was exploratory after comparing three fixed alternatives. Paired intervals do not adjust for model selection.',
 gate_text,
 paste('Model:',model),paste('Evidence:',start,'to',sr$evidence_end),
 'The probabilities concern CPUE >=0.02 conditional on diving, not CPUE >0.04 or occurrence at unvisited sites.',
 'Spatial habitat validation and forward evidence testing retain the limitations stated in the analysis plan.',
 'Annual checks, exploration effort, tow detection and vessel-time savings are not established by this layer.'),file.path(out,'README.txt'))
log_message('Experimental CSV, GeoPackage and KMZ exported; incumbent checksums unchanged.')
