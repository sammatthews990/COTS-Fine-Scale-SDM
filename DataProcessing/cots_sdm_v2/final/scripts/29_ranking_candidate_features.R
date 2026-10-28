# Question: do recent CPUE and habitat hotspots add information to site ranks?
# Method: derive prior-season evidence and full-cell habitat quantiles; retain
# the original five reef folds, hyperparameters and all designated candidates.
# Validation: replay archived 16-point predictions before scoring new summaries;
# verify frozen inputs, dates, candidate keys and point-to-site coverage.
source(file.path(dirname(sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE))),'..','R','common.R'))
cfg <- initialise_final()
source(file.path(cfg$final_dir,'R','cull_sites.R'))
rc <- cfg$ranking_candidates
cache <- file.path(cfg$dirs$cache,'ranking_candidates')
dir.create(cache,showWarnings=FALSE)
keys <- c('reef_name','year','site_row')
input_names <- c('assessment_dives','assessment_sites','assessment_evidence',
 'assessment_manta_sites','assessment_forward','assessment_forward_coefficients',
 'ranking_train_site_years','ranking_site_sample_points','predictor_stack_rg',
 'manta_tows','footprint_links')
paths <- setNames(lapply(input_names,function(n) checked_input(cfg,n)),input_names)
signature <- list(version=2L,inputs=lapply(input_names,function(n) cfg$inputs[[n]]$sha256),
 half_life_days=rc$half_life_days)
final_cache <- file.path(cache,'features.rds')
if(file.exists(final_cache)) {
 previous <- readRDS(final_cache)
 if(identical(previous$signature,signature)) {
  log_message('Candidate features already match frozen inputs and settings.');quit(save='no')
 }
}
sites <- read_designated_sites(cfg)
sites$site_row <- seq_len(nrow(sites))
dives <- matched_cull_dives(cfg,sites)
old <- readRDS(paths$assessment_manta_sites)
sy <- old[,keys]
stopifnot(!anyDuplicated(sy[keys]),'cull_cv' %in% names(old))
stack <- terra::rast(paths$predictor_stack_rg)
names(stack) <- cfg$legacy$rg_layer_names
rg <- cfg$legacy$rg_layer_names;features <- c(rg,'Year')
points <- readRDS(paths$ranking_site_sample_points)
pts_by_site <- split(seq_along(points$site),points$site)

# Same spatial folds as script 17. Tow ordering is retained by the join.
tows <- readRDS(paths$manta_tows)[,c('visit_id','reef_name')]
geom <- readRDS(paths$footprint_links)$manta |>
 dplyr::group_by(.data$visit_id) |>
 dplyr::summarise(sx=dplyr::first(.data$start_x),ex=dplyr::first(.data$end_x),.groups='drop')
tows <- dplyr::inner_join(tows,geom,by='visit_id') |>
 dplyr::filter(is.finite(.data$sx),is.finite(.data$ex))
train <- readRDS(paths$ranking_train_site_years) |>
 dplyr::filter(!is.na(.data$high_density),is.finite(.data$longitude),is.finite(.data$latitude))
xy <- sf::st_coordinates(sf::st_transform(sf::st_as_sf(train,coords=c('longitude','latitude'),crs=4326),terra::crs(stack)))
tx <- terra::extract(stack,xy)[,rg];tx$Year <- train$year
reefs <- unique(c(train$reef_name,tows$reef_name,sy$reef_name))
set.seed(cfg$hmsc$seed+17L)
fold <- setNames(sample(rep(1:5,length.out=length(reefs))),reefs)
params <- list(objective='binary:logistic',eta=0.00137382379588326,max_depth=5,
 min_child_weight=20,gamma=1e-10,nthread=8,tree_method='hist')

cell_cache <- file.path(cache,'full_polygon_cells.rds')
if(file.exists(cell_cache) && identical(readRDS(cell_cache)$signature,signature)) {
 cells <- readRDS(cell_cache)$cells
} else {
 log_message('Extracting all intersecting 30 m predictor cells for ',nrow(sites),' polygons.')
 # Find polygon membership on an empty grid. This avoids rereading 17 raster
 # bands for every overlapping polygon. Then read distinct cells in row order.
 grid <- terra::rast(nrows=nrow(stack),ncols=ncol(stack),extent=terra::ext(stack),crs=terra::crs(stack))
 owners <- as.data.frame(terra::cells(grid,terra::vect(sites),touches=TRUE))
 names(owners)[1:2] <- c('ID','cell')
 owners <- owners[!duplicated(owners[c('ID','cell')]) & is.finite(owners$cell),c('ID','cell')]
 absent <- setdiff(seq_len(nrow(sites)),owners$ID)
 if(length(absent)) {
  xy_fallback <- sf::st_coordinates(suppressWarnings(sf::st_point_on_surface(sites[absent,])))
  owners <- rbind(owners,data.frame(ID=absent,cell=terra::cellFromXY(grid,xy_fallback)))
 }
 ids <- sort(unique(owners$cell));stopifnot(all(is.finite(ids)))
 log_message(length(ids),' distinct cells identified in ',nrow(owners),' polygon-cell memberships.')
 blocks <- split(seq_along(ids),ceiling(seq_along(ids)/100000L))
 values <- vector('list',length(blocks))
 for(j in seq_along(blocks)) {
  file <- file.path(cache,sprintf('cell_values_%03d.rds',j))
  ii <- blocks[[j]]
  previous <- if(file.exists(file)) readRDS(file) else NULL
  if(!is.null(previous) && identical(previous$signature,signature) && identical(previous$ids,ids[ii])) {
   values[[j]] <- previous$values
  } else {
   values[[j]] <- as.data.frame(terra::extract(stack,ids[ii],raw=TRUE));names(values[[j]]) <- rg
   saveRDS(list(signature=signature,ids=ids[ii],values=values[[j]]),file)
  }
  log_message('Raster value block ',j,'/',length(blocks),' complete.')
 }
 values <- dplyr::bind_rows(values)
 cells <- cbind(owners,values[match(owners$cell,ids),,drop=FALSE]);rownames(cells) <- NULL
 # Check the faster geometry/value route against direct polygon extraction.
 sample_ids <- unique(round(seq(1,nrow(sites),length.out=5)))
 probe <- terra::extract(stack,terra::vect(sites[sample_ids,]),cells=TRUE,touches=TRUE)
 probe$ID <- sample_ids[probe$ID]
 key <- function(x) paste(x$ID,x$cell,sep=':')
 ci <- match(key(probe),key(cells));stopifnot(!anyNA(ci))
 a <- as.matrix(probe[,rg]);b <- as.matrix(cells[ci,rg])
 stopifnot(identical(as.vector(is.na(a)),as.vector(is.na(b))),max(abs(a-b),na.rm=TRUE)<1e-12)
 log_message('Full-cell extraction verified against direct extraction at five distributed polygons.')
 saveRDS(list(signature=signature,cells=cells),cell_cache)
}
cell_counts <- tabulate(cells$ID,nbins=nrow(sites))
stopifnot(all(cell_counts>0),all(cells$ID %in% sites$site_row))
log_message(format(nrow(cells),big.mark=','),' site-cell records; median ',median(cell_counts),' per site.')
write_table(data.frame(site_row=sites$site_row,reef_name=sites$reef_name,
 cells=cell_counts,valid_cells=tabulate(cells$ID[rowSums(is.finite(as.matrix(cells[,rg])))>0],nbins=nrow(sites))),
 cfg,'ranking_candidate_habitat_coverage.csv')
summaries <- list();replay <- list()
for(k in 1:5) {
 tr <- fold[train$reef_name]!=k
 model <- xgboost::xgb.train(params=params,nrounds=1500,verbose=0,
  data=xgboost::xgb.DMatrix(as.matrix(tx[tr,features]),label=train$high_density[tr],missing=NA))
 xgboost::xgb.save(model,file.path(cache,paste0('habitat_fold_',k,'.json')))
 held <- which(fold[sy$reef_name]==k)
 pi <- pts_by_site[as.character(sy$site_row[held])]
 xx <- points$values[unlist(pi,use.names=FALSE),rg];xx$Year <- rep(sy$year[held],lengths(pi))
 pp <- predict(model,xgboost::xgb.DMatrix(as.matrix(xx[,features]),missing=NA))
 means <- as.numeric(tapply(pp,rep(seq_along(held),lengths(pi)),mean))
 error <- max(abs(means-old$cull_cv[held]),na.rm=TRUE)
 replay[[k]] <- data.frame(fold=k,site_years=length(held),max_abs_error=error)
 log_message('Fold ',k,' archived mean replay maximum error: ',signif(error,5))
 if(error>1e-5) stop('Archived habitat predictions did not replay; inspect folds/library before proceeding.')
 for(year in sort(unique(sy$year[held]))) {
  rows <- held[sy$year[held]==year];ids <- sy$site_row[rows]
  ci <- which(cells$ID %in% ids)
  x <- cells[ci,rg];x$Year <- year
  p <- predict(model,xgboost::xgb.DMatrix(as.matrix(x[,features]),missing=NA))
  z <- lapply(split(p,cells$ID[ci]),function(v) c(mean=mean(v),q50=unname(quantile(v,.5)),
   q75=unname(quantile(v,.75)),q90=unname(quantile(v,.9)),q95=unname(quantile(v,.95))))
  zz <- as.data.frame(do.call(rbind,z));zz$site_row <- as.integer(rownames(zz));zz$year <- year
  summaries[[length(summaries)+1L]] <- zz
 }
 log_message('Fold ',k,' full-cell summaries complete.')
}
write_table(dplyr::bind_rows(replay),cfg,'ranking_candidate_habitat_replay.csv')
hab <- dplyr::bind_rows(summaries) |>
 dplyr::left_join(sf::st_drop_geometry(sites)[,c('site_row','reef_name')],by='site_row') |>
 dplyr::group_by(.data$reef_name,.data$year) |>
 dplyr::mutate(dplyr::across(dplyr::all_of(c('mean','q50','q75','q90','q95')),
  ~ (rank(.x,ties.method='average')-1)/max(dplyr::n()-1,1),.names='h_{.col}')) |>
 dplyr::ungroup()
stopifnot(nrow(hab)==nrow(sy),!anyDuplicated(hab[keys]))
ev <- readRDS(paths$assessment_evidence) |>
 dplyr::filter(.data$window=='last_season') |>
 dplyr::inner_join(old[,c(keys,'rank_cull')],by=keys) |>
 dplyr::inner_join(hab,by=keys)
direct <- list();date_audit <- list()
for(year in sort(unique(ev$year))) {
 start <- as.Date(sprintf('%d-01-01',year-1));end <- as.Date(sprintf('%d-01-01',year))
 prior <- dives[dives$year==year-1 & !is.na(dives$site_row),]
 high <- prior |> dplyr::group_by(.data$site_row) |>
  dplyr::summarise(site_high_prev=as.numeric(any(.data$cpue>.04)),.groups='drop')
 dated <- prior[!is.na(prior$survey_date) & prior$survey_date>=start & prior$survey_date<end,]
 dated$age_days <- as.numeric(end-dated$survey_date)
 stopifnot(all(dated$age_days>0))
 dated$w <- 2^(-dated$age_days/rc$half_life_days)
 continuous <- dated |> dplyr::group_by(.data$site_row) |>
  dplyr::summarise(recent_cpue=sum(.data$w*.data$total_cots)/sum(.data$w*.data$bottom_time),
   recent_weight=sum(.data$w),latest_prior_date=if(dplyr::n()) max(.data$survey_date) else as.Date(NA),.groups='drop')
 direct[[as.character(year)]] <- dplyr::full_join(high,continuous,by='site_row') |>
  dplyr::mutate(year=year)
 date_audit[[as.character(year)]] <- data.frame(year,prior_records=nrow(prior),
  records_with_valid_dates=nrow(dated),excluded_dates=nrow(prior)-nrow(dated))
}
ev <- dplyr::left_join(ev,dplyr::bind_rows(direct),by=c('site_row','year')) |>
 dplyr::mutate(dplyr::across(dplyr::all_of(c('site_high_prev','recent_cpue','recent_weight')),
  ~dplyr::coalesce(.x,0)))
write_table(dplyr::bind_rows(date_audit),cfg,'ranking_candidate_date_audit.csv')
stopifnot(all(is.finite(ev$recent_cpue)),all(ev$recent_cpue>=0),
 all(is.na(ev$latest_prior_date) | ev$latest_prior_date<as.Date(sprintf('%d-01-01',ev$year))))
saveRDS(list(signature=signature,sites=sf::st_drop_geometry(sites),features=ev,
 dives=dives,habitat=hab,folds=fold,cell_counts=cell_counts),final_cache)
log_message('Candidate features saved: ',nrow(ev),' site-years. No deployed files changed.')
