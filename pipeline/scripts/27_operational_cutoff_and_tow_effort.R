# Question: what does deferring class D mean for cull dives and manta effort?
# Method: replay the frozen forward rankings without refitting; evaluate the
# exact A-C product cutoff and overlay recorded straight tow segments on sites.
# Validation: reproduce every saved forward rank; paired reef-year bootstrap;
# exact random-site expectations at the same integer site count. Tow savings
# are conditional spatial proxies, not a test of detection after skipping tows.
source(file.path(dirname(sub('^--file=', '', grep('^--file=', commandArgs(FALSE), value=TRUE))), '..', 'R', 'common.R'))
cfg <- initialise_final()
source(file.path(cfg$final_dir,'R/cull_sites.R'))
source(file.path(cfg$final_dir,'R/site_selection.R'))
keys <- c('reef_name','year','site_row')
paths <- lapply(c('assessment_dives','assessment_sites','assessment_evidence',
 'assessment_manta_sites','assessment_forward','assessment_forward_coefficients',
 'manta_tows','footprint_links'),function(n) checked_input(cfg,n))
names(paths) <- c('dives','sites','evidence','habitat','forward','coefficients','tows','links')

# Reproduce script 19, including its random tie order across ALL ten scores.
# Replaying only the selected model would consume a different random stream.
ev <- readRDS(paths$evidence) |>
 dplyr::filter(.data$window=='last_season') |>
 dplyr::select(dplyr::all_of(keys),'site_problematic_prev','site_quiet_prev',
 'cull_rate','cull_effort','manta_strong','manta_effort','reef_cull_prev','reef_manta_prev') |>
 dplyr::inner_join(readRDS(paths$habitat) |>
 dplyr::select(dplyr::all_of(keys),'rank_cull','rank_within','rank_ensemble'),by=keys)
coefs <- read.csv(paths$coefficients)
habitats <- list(cull=list(h='rank_cull'),ensemble=list(h='rank_ensemble'),
 both=list(h='rank_cull',h2='rank_within'))
set.seed(cfg$hmsc$seed+19L)
full <- list()
for(T in 2019:2026) {
 test <- ev[ev$year==T,]
 scores <- data.frame(test[,keys],habitat_only=test$rank_cull)
 for(hn in names(habitats)) for(wn in c('all','last5','last3')) {
  hab <- habitats[[hn]]
  x <- model.matrix(delete.response(terms(site_evidence_formula(habitat2=!is.null(hab$h2)))),
   site_evidence_terms(test,hab$h,hab$h2))
  b <- coefs[coefs$test_year==T & coefs$habitat==hn & coefs$window==wn,]
  scores[[paste(hn,wn,sep='_')]] <- as.vector(x[,b$term,drop=FALSE] %*% b$estimate)
 }
 score_names <- setdiff(names(scores),keys)
 full[[as.character(T)]] <- scores |> dplyr::group_by(.data$reef_name,.data$year) |>
  dplyr::mutate(dplyr::across(dplyr::all_of(score_names),pct_rank)) |> dplyr::ungroup()
}
full <- dplyr::bind_rows(full)
old <- readRDS(paths$forward)
chk <- dplyr::inner_join(old,full,by=keys,suffix=c('_old','_new'),relationship='many-to-one')
stopifnot(nrow(chk)==nrow(old))
for(s in score_names) stopifnot(identical(chk[[paste0(s,'_old')]],chk[[paste0(s,'_new')]]))
log_message('All saved forward ranks reproduced exactly; no model refit.')
full <- full |> dplyr::semi_join(unique(old[c('reef_name','year')]),by=c('reef_name','year')) |>
 dplyr::select(dplyr::all_of(keys),habitat='habitat_only',combined='cull_all') |>
 dplyr::group_by(.data$reef_name,.data$year) |> dplyr::mutate(N=dplyr::n()) |> dplyr::ungroup()
full$reef_year <- paste(full$reef_name,full$year,sep='||')
groups <- sort(unique(full$reef_year)); ng <- length(groups)
full$g <- match(full$reef_year,groups)
ns <- as.numeric(table(factor(full$reef_year,levels=groups)))
stopifnot(!anyDuplicated(full[keys]),all(full$site_row>0))
sites <- read_designated_sites(cfg)
dives <- dplyr::inner_join(matched_cull_dives(cfg,sites),full,by=keys,relationship='many-to-one')
stopifnot(nrow(dives)==nrow(old))
bands <- c('high','low','ghost')
dives$band <- ifelse(dives$cpue>.04,'high',ifelse(dives$total_cots==0,'ghost',
 ifelse(dives$cpue<.02,'low','intermediate')))
counts <- function(keep) sapply(bands,function(b) tabulate(dives$g[keep & dives$band==b],nbins=ng))
totals <- counts(rep(TRUE,nrow(dives)))

# Manta route proxy: retain every part outside designated polygons, and every
# part shared with ANY selected site. Only deferred-only portions can be saved.
# Endpoints define a straight segment, not a measured vessel track. Exclude
# source geometry fallbacks rather than treating point records as zero effort.
tows <- readRDS(paths$tows)[,c('visit_id','reef_name','year','survey_date')]
stopifnot(!anyDuplicated(tows$visit_id))
ln <- readRDS(paths$links)$manta |>
 dplyr::group_by(.data$visit_id) |>
 dplyr::summarise(dplyr::across(dplyr::all_of(c('start_longitude','start_latitude',
 'end_longitude','end_latitude','geometry_fallback')),dplyr::first),.groups='drop')
tows <- tows |> dplyr::semi_join(unique(full[c('reef_name','year')]),by=c('reef_name','year')) |>
 dplyr::left_join(ln,by='visit_id',relationship='one-to-one')
n_all_tows <- nrow(tows)
ok <- !is.na(tows$geometry_fallback) & !tows$geometry_fallback &
 apply(as.matrix(tows[,c('start_longitude','start_latitude','end_longitude','end_latitude')]),1,
 function(x) all(is.finite(x)))
tows <- tows[ok,]
sx <- sf::sf_project('EPSG:4326',sf::st_crs(sites)$wkt,as.matrix(tows[,c('start_longitude','start_latitude')]))
ex <- sf::sf_project('EPSG:4326',sf::st_crs(sites)$wkt,as.matrix(tows[,c('end_longitude','end_latitude')]))
tows$length <- sqrt(rowSums((ex-sx)^2))
ok <- is.finite(tows$length) & tows$length>0
tows <- tows[ok,]; sx <- sx[ok,,drop=FALSE]; ex <- ex[ok,,drop=FALSE]
tows$tow <- seq_len(nrow(tows)); tows$g <- match(paste(tows$reef_name,tows$year,sep='||'),groups)
stopifnot(nrow(tows)>0,!anyNA(tows$g))

# Split each straight line at all site boundaries. Owners can overlap; each
# resulting piece is counted once. Cache is keyed to the frozen geometry
# inputs, evaluation reef-years and geometry algorithm version. Increment the
# version below if changing the geometry calculation.
signature <- digest::digest(list(geometry_version=2L,vapply(paths[c('sites','tows','links')],sha256_file,character(1)),unique(full[c('reef_name','year')])) )
cache <- file.path(cfg$dirs$cache,'operational_tow_pieces.rds')
saved <- if(file.exists(cache)) readRDS(cache) else NULL
if(!is.null(saved) && identical(saved$signature,signature)) pieces <- saved$pieces else {
 pieces <- list()
 for(reef in unique(tows$reef_name)) {
  ti <- which(tows$reef_name==reef); si <- which(sites$reef_name==reef)
  lines <- sf::st_sf(tow=ti,geometry=sf::st_sfc(lapply(ti,function(i)
   sf::st_linestring(rbind(sx[i,],ex[i,]))),crs=sf::st_crs(sites)))
  polygons <- sf::st_sf(site_row=si,geometry=sf::st_make_valid(sf::st_geometry(sites[si,])))
  hit <- suppressWarnings(sf::st_intersection(lines,polygons))
  hit <- hit[as.numeric(sf::st_length(hit))>1e-7,]
  if(!nrow(hit)) next
  hit <- suppressWarnings(sf::st_cast(hit,'LINESTRING'))
  ends <- lapply(sf::st_geometry(hit),function(x) rbind(x[1,1:2],x[nrow(x),1:2]))
  hit <- sf::st_drop_geometry(hit)
  fs <- Map(function(x,i) {
   v <- ex[i,]-sx[i,]
   f <- as.vector((x-matrix(sx[i,],2,2,byrow=TRUE)) %*% v)/sum(v^2)
   c(max(0,min(f)),min(1,max(f)))
  },ends,hit$tow)
  hit$lo <- vapply(fs,`[`,numeric(1),1L)
  hit$hi <- vapply(fs,`[`,numeric(1),2L)
  by_tow <- split(seq_len(nrow(hit)),hit$tow)
  for(ix in by_tow) {
   i <- hit$tow[ix[1]]; h <- hit[ix,]
   cuts <- sort(unique(c(0,1,h$lo,h$hi)))
   mids <- (head(cuts,-1)+tail(cuts,-1))/2
   owners <- lapply(mids,function(m) sort(unique(h$site_row[h$lo<m & h$hi>m])))
   # Outside pieces do not contribute to the avoidable length; their presence
   # nevertheless prevents a whole tow being classified as avoidable below.
   inside <- lengths(owners)>0
   if(any(inside)) pieces[[length(pieces)+1L]] <- tibble::tibble(tow=i,
    fraction=diff(cuts)[inside],owners=owners[inside])
  }
  log_message('Tow geometry: ',reef,' (',length(ti),' segments)')
 }
 pieces <- dplyr::bind_rows(pieces)
 saveRDS(list(signature=signature,pieces=pieces),cache)
}
pieces$g <- tows$g[pieces$tow]
pieces$length <- pieces$fraction*tows$length[pieces$tow]
owner_keys <- paste(rep(pieces$g,lengths(pieces$owners)),unlist(pieces$owners,use.names=FALSE),sep=':')
owner_index <- match(owner_keys,paste(full$g,full$site_row,sep=':'))
stopifnot(!anyNA(owner_index))
owners <- unname(split(owner_index,rep(seq_len(nrow(pieces)),lengths(pieces$owners))))
stopifnot(all(lengths(owners)==lengths(pieces$owners)))
best <- lapply(c('habitat','combined'),function(sc) vapply(owners,function(i) max(full[[sc]][i]),numeric(1)))
names(best) <- c('habitat','combined')
m <- lengths(owners)
# Weighted totals with explicit zeroes for groups without observations.
wsum <- function(v,g,n) { z <- numeric(n); q <- tapply(v,g,sum); z[as.integer(names(q))] <- q; z }
inside <- wsum(pieces$fraction,pieces$tow,nrow(tows))
stopifnot(all(inside<=1+1e-7))
whole <- which(inside>=1-1e-7)
union_owners <- lapply(split(pieces$owners,pieces$tow),function(x) unique(unlist(x)))
whole_m <- integer(nrow(tows)); whole_m[as.integer(names(union_owners))] <- lengths(union_owners)
tow_den <- tabulate(tows$g,nbins=ng); length_den <- wsum(tows$length,tows$g,ng)

set.seed(cfg$hmsc$seed+27L)
weights <- cbind(rep(1,ng),replicate(500L,tabulate(sample.int(ng,ng,replace=TRUE),nbins=ng)))
aggregate_ratio <- function(num,den) as.vector(crossprod(weights,num))/as.vector(crossprod(weights,den))
interval <- function(x) c(value=x[1],low=unname(quantile(x[-1],.025)),high=unname(quantile(x[-1],.975)))
random_none <- function(N,K,m) ifelse(N-K<m,0,exp(lchoose(N-m,K)-lchoose(N,K)))
fixed <- matched <- policy <- list(); curves <- list()
fractions <- seq(.05,1,.05)
evaluate <- function(sc,a,operational=FALSE) {
 K <- if(operational) floor(a*ns+1e-9) else vapply(ns,function(n) sum(seq(0,1,length.out=n)>=1-a-1e-9),numeric(1))
 # Minimum percentile rank retained by the exact integer top-K rule.
 cut <- ifelse(K>0,(ns-K)/pmax(ns-1,1),Inf)
 if(sc=='random') {
  nums <- totals*(K/ns)
  prob <- random_none(ns[pieces$g],K[pieces$g],m)
  removed_length <- wsum(prob*pieces$length,pieces$g,ng)
  prob_whole <- random_none(ns[tows$g[whole]],K[tows$g[whole]],whole_m[whole])
  removed_tows <- wsum(prob_whole,tows$g[whole],ng)
 } else {
  keep <- dives[[sc]]>=cut[dives$g]-1e-9
  nums <- counts(keep)
  removed <- best[[sc]]<cut[pieces$g]-1e-9
  removed_length <- wsum(pieces$length*removed,pieces$g,ng)
  removed_fraction <- wsum(pieces$fraction*removed,pieces$tow,nrow(tows))
  removed_tows <- tabulate(tows$g[removed_fraction>=1-1e-7],nbins=ng)
 }
 list(sites_selected=aggregate_ratio(K,ns),high_retained=aggregate_ratio(nums[,1],totals[,1]),
  low_avoided=1-aggregate_ratio(nums[,2],totals[,2]),ghost_avoided=1-aggregate_ratio(nums[,3],totals[,3]),
  manta_tows_avoided=aggregate_ratio(removed_tows,tow_den),
  manta_distance_avoided=aggregate_ratio(removed_length,length_den))
}
row_result <- function(values,sc,rule,a) dplyr::bind_rows(lapply(names(values),function(metric) {
 v <- interval(values[[metric]])
 data.frame(score=sc,rule,cutoff=a,metric,value=v[1],low=v[2],high=v[3])
}))
for(sc in c('habitat','combined','random')) {
 cs <- lapply(fractions,function(a) evaluate(sc,a))
 curves[[sc]] <- lapply(names(cs[[1]]),function(metric) {
  origin <- switch(metric,sites_selected=rep(0,501),high_retained=rep(0,501),
   manta_tows_avoided=aggregate_ratio(tabulate(tows$g[whole],nbins=ng),tow_den),
   manta_distance_avoided=aggregate_ratio(wsum(pieces$length,pieces$g,ng),length_den),rep(1,501))
  cbind(origin,do.call(cbind,lapply(cs,`[[`,metric)))
 })
 names(curves[[sc]]) <- names(cs[[1]])
 for(j in seq_along(cs)) fixed[[length(fixed)+1L]] <- row_result(cs[[j]],sc,'percentile',fractions[j])
 for(target in c(.8,.9,.95,.99)) {
  cv <- curves[[sc]]
  vals <- lapply(cv,function(y) vapply(seq_len(501),function(b) {
   x <- cv$high_retained[b,]; j <- which(x>=target)[1]
   if(x[j]==target) y[b,j] else y[b,j-1]+(y[b,j]-y[b,j-1])*(target-x[j-1])/(x[j]-x[j-1])
  },numeric(1)))
  matched[[length(matched)+1L]] <- row_result(vals,sc,'matched_retention',target)
 }
 policy[[length(policy)+1L]] <- row_result(evaluate(sc,.8,TRUE),sc,'A-C',.8)
}
fixed <- dplyr::bind_rows(fixed); matched <- dplyr::bind_rows(matched); policy <- dplyr::bind_rows(policy)
stopifnot(all(is.finite(unlist(fixed[c('value','low','high')]))),all(fixed$value>=-1e-9 & fixed$value<=1+1e-9))
for(sc in c('habitat','combined','random')) {
 z <- fixed[fixed$score==sc & abs(fixed$cutoff-1)<1e-8,]
 stopifnot(all(abs(z$value[z$metric %in% c('manta_tows_avoided','manta_distance_avoided','low_avoided','ghost_avoided')])<1e-8))
}
# Old matched-retention point estimates must be unchanged. Intervals are newly
# paired across cull and tow metrics, and intentionally use a separate seed.
old_match <- read.csv(file.path(cfg$dirs$tables,'cull_high_density_matched.csv'))
for(sc in c('habitat','combined','random')) for(t in c(.8,.9,.95,.99)) {
 a <- matched$value[matched$score==sc & matched$cutoff==t & matched$metric=='sites_selected']
 b <- old_match$value[old_match$design=='forward' & old_match$period=='2019-2026' & old_match$score==sc & old_match$target==t & old_match$metric=='sites_searched']
 stopifnot(length(a)==1,length(b)==1,abs(a-b)<1e-10)
}
write_table(fixed,cfg,'operational_tow_effort_fixed.csv')
write_table(matched,cfg,'operational_tow_effort_matched.csv')
write_table(policy,cfg,'operational_priority_cutoff.csv')
write_table(data.frame(reef_years=ng,dives=nrow(dives),site_years=nrow(full),
 manta_records=n_all_tows,manta_valid_segments=nrow(tows),manta_excluded=n_all_tows-nrow(tows),
 manta_reef_years=sum(tow_den>0),manta_straight_distance_km=sum(tows$length)/1000,
 distance_inside_sites=sum(pieces$length)/sum(tows$length),tows_entirely_inside_sites=length(whole)/nrow(tows)),
 cfg,'operational_tow_effort_coverage.csv')
print(policy[policy$score=='combined',],row.names=FALSE)
print(matched[matched$score=='combined' & grepl('manta',matched$metric),],row.names=FALSE)
